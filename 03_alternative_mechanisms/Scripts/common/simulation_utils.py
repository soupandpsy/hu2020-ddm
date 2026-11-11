"""Native HDDM forward and prior-predictive simulation; no custom solver."""
import random
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from common.baseline_v import trial_parameters, build_model, natural_values, cell_key
from common.data_utils import encode_conditions


def make_design(subjects=2, trials_per_cell=8):
    if subjects < 2 or trials_per_cell < 1:
        raise ValueError("At least two subjects and one trial per cell required")
    rows = [(subj, val, ident, stim) for subj in range(subjects)
            for val in ("Good", "Bad") for ident in ("Self", "Other")
            for stim in (0, 1) for _ in range(trials_per_cell)]
    data = pd.DataFrame(rows, columns=["subj_idx", "val", "id", "stim"])
    data["response"] = data.stim
    data["rt"] = 0.6
    data["event_id"] = ["sim_%08d" % i for i in range(len(data))]
    return encode_conditions(data)


def example_truth(data, name, seed=20261106):
    """Transparent recovery scenario, NOT a prior distribution."""
    rng = np.random.default_rng(seed)
    truth = {}
    for subj in sorted(data.subj_idx.unique()):
        suffix=str(int(subj))
        truth['a_subj.'+suffix]=1.25+rng.normal(0,.04)
        truth['z_subj.'+suffix]=.55+rng.normal(0,.015)
        for ident in ['Self','Other']:
            for match in ['match','mismatch']:
                for val in ['Good','Bad']:
                    # Preserve separate match/mismatch effects: .5 only on matching Self.
                    effect=.5 if (ident=='Self' and match=='match' and val=='Good') else 0.
                    truth[cell_key('v',ident,match,val,subj)]=1.5+effect+rng.normal(0,.08)
                    truth[cell_key('t',ident,match,val,subj)]=.25+rng.normal(0,.015)
        if name=='M1':truth['beta_z.'+suffix]=.8+rng.normal(0,.1)
    return truth


def simulate_from_parameters(data, parameters, seed=20261106,
                             decision_window=8.0, grid_step=0.001):
    """Native CDF simulation in unit-boundary coordinates, then restore seconds.

    For dX=mu*dt+dW, Y=X/a and u=t/a**2 give drift mu*a and boundary 1.
    This exact change of units resolves tiny-a grid failure without clipping a.
    The numerical window/grid are in u, not total RT or experiment deadline.
    A finite mass check rejects inadequate grids instead of silently normalizing
    missing probability. Extreme RT addition rounding is explicitly recorded.
    """
    import hddm
    if (not np.isfinite(decision_window) or decision_window <= 0
            or not np.isfinite(grid_step) or not 0 < grid_step < decision_window):
        raise ValueError("Positive normalized decision window and smaller grid step required")
    if 2*decision_window/grid_step > 2_000_000:
        raise ValueError("Numerical grid exceeds the explicit two-million-point limit")
    frame = encode_conditions(data)
    pars = parameters.reset_index(drop=True)
    if len(pars) != len(frame) or not np.isfinite(pars[["a","t","z","mu"]]).all().all():
        raise ValueError("Finite aligned trial parameters required")
    if not (pars.a.gt(0).all() and pars.t.ge(0).all() and pars.z.between(0,1,inclusive="neither").all()):
        raise ValueError("Invalid native DDM parameter domain")
    random.seed(seed); np.random.seed(seed)
    result = frame.copy()
    groups = pars.groupby(["a", "t", "z", "mu"], sort=True).groups
    seeds = np.random.SeedSequence(seed).spawn(len(groups))
    grid = np.arange(-decision_window, decision_window, grid_step)
    mass_checks = []
    for ((a,t,z,mu), indices), child in zip(groups.items(),seeds):
        scale = float(a)*float(a)
        drift = float(mu)*float(a)
        if not np.isfinite(scale) or scale <= 0 or not np.isfinite(drift):
            raise ValueError("Boundary rescaling exceeds floating-point range; no parameter clipping")
        density = hddm.wfpt.pdf_array(grid, drift, 0., 1., float(z), 0., 0., 0., err=1e-8)
        mass = float(density[1:].sum()*grid_step)
        if (not np.isfinite(density).all() or (density < 0).any()
                or abs(mass-1.) > 1e-4):
            raise ValueError("Normalized native CDF grid has inadequate mass: %r; refine settings explicitly" % mass)
        mass_checks.append(mass)
        np.random.seed(int(child.generate_state(1)[0]))
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", RuntimeWarning)
                generated = hddm.generate.gen_rts(size=len(indices), range_=(-decision_window,decision_window),
                    dt=grid_step, a=1., t=0., z=float(z), v=drift, sv=0, sz=0, st=0)
        except RuntimeWarning as error:
            raise ValueError("Native CDF simulation failed for a=%r, t=%r, z=%r, mu=%r; "
                             "do not treat this draw as a valid prior prediction" % (a,t,z,mu)) from error
        if (not np.isfinite(generated.rt).all() or (generated.rt <= 0).any()
                or not generated.response.isin([0, 1]).all()):
            raise ValueError("Native simulator returned invalid RT/choice values")
        decision_time = generated.rt.to_numpy()*scale
        rt = float(t)+decision_time
        if not np.isfinite(rt).all() or (decision_time <= 0).any():
            raise ValueError("Restored RT exceeds floating-point range")
        result.loc[indices,"rt"] = rt
        result.loc[indices,"decision_time"] = decision_time
        result.loc[indices,"nondecision_time"] = float(t)
        result.loc[indices,"rt_rounding_lost_decision_time"] = rt == float(t)
        result.loc[indices,"response"] = generated.response.to_numpy().astype(int)
    result["accuracy"] = (result.response == result.stim).astype(int)
    # Original response keys do not describe generated choices.
    result = result.drop(columns=[c for c in ("response_key","upper_key") if c in result])
    result = encode_conditions(result)
    result.attrs['simulation_numerics'] = dict(method='native_cdf_unit_boundary',
        decision_window=decision_window,grid_step=grid_step,min_grid_mass=min(mass_checks),
        max_grid_mass=max(mass_checks),rt_addition_rounding_count=int(result.rt_rounding_lost_decision_time.sum()))
    return result


def simulate_dataset(data, name, truth=None, seed=20261106):
    truth = example_truth(data,name,seed) if truth is None else dict(truth)
    return simulate_from_parameters(data,trial_parameters(truth,data,name),seed), truth


def prior_predictive(data, name, draws=10, seed=20261106, priors=None, simulation_options=None):
    """Draw from the actual native model prior, without fitting observations.

    Observed nodes are excluded from the forward prior traversal; data supply
    only design and subjects. Returns (replicated_data, stochastic_prior_draws).
    """
    if draws < 1:
        raise ValueError("draws must be positive")
    import pymc as pm
    model = build_model(data,name,priors)
    graph = pm.MCMC(list(model.get_stochastics()["node"]))
    if not graph.stochastics:
        raise RuntimeError("Native prior graph contains no stochastic nodes")
    random.seed(seed); np.random.seed(seed)
    replicas, priors_out = [], []
    for draw_id in range(draws):
        for generation in graph.generations:
            for node in sorted(generation,key=lambda n:n.__name__):
                if not node.observed:
                    node.random()
        values = natural_values(model)
        generated = simulate_from_parameters(data,trial_parameters(values,data,name),seed+draw_id+1,
                                             **(simulation_options or {}))
        generated["replicate"] = draw_id
        replicas.append(generated);priors_out.append(values)
    return pd.concat(replicas,ignore_index=True),pd.DataFrame(priors_out)


def save_simulation(data, truth, output_dir):
    from common.data_utils import write_json, dataset_hash
    out=Path(output_dir);out.mkdir(parents=True,exist_ok=False)
    data.to_csv(out/"data.csv.gz",index=False)
    write_json(out/"truth.json",truth)
    write_json(out/"manifest.json",{"data_hash":dataset_hash(data),"scope":"uncensored synthetic DDM"})


def simulate_dynamic(data, name, truth, seed=20261107, settings=None):
    """Forward sample from PyDDM's joint distribution, then add cell t0.

    Undecided draws are an explicit numerical-window failure, never discarded.
    """
    from common.baseline_v import dynamic_solution, dynamic_trial_parameters
    frame=encode_conditions(data).reset_index(drop=True)
    pars=dynamic_trial_parameters(truth,frame,name)
    result=frame.copy(); masses=[]
    groups=pars.groupby(['mu','a','z','t','effect','contrast'],sort=True).groups
    seeds=np.random.SeedSequence(seed).spawn(len(groups))
    for ((mu,a,z,t0,effect,contrast),ix),child in zip(groups.items(),seeds):
        sol,scale,mass=dynamic_solution(name,mu,a,z,effect,contrast,settings)
        draw_seed=int(child.generate_state(1)[0])
        sample=sol.sample(k=len(ix),seed=draw_seed)
        if sample.undecided:
            raise FloatingPointError('Finite-window undecided draw: extend numerical duration; no filtering')
        decisions=np.r_[sample.choice_upper,sample.choice_lower]*scale
        choice=np.r_[np.ones(len(sample.choice_upper)),np.zeros(len(sample.choice_lower))].astype(int)
        order=np.random.default_rng(draw_seed).permutation(len(ix))
        result.loc[ix,'rt']=t0+decisions[order]
        result.loc[ix,'decision_time']=decisions[order]
        result.loc[ix,'response']=choice[order]
        masses.append(mass)
    if (result.decision_time<=0).any(): raise FloatingPointError('Nonpositive generated decision time')
    result=result.drop(columns=[c for c in ('response_key','upper_key') if c in result])
    result['accuracy']=(result.response==result.stim).astype(int)
    result.attrs['simulation_numerics']=dict(method='PyDDM implicit sample',min_mass=min(masses),max_mass=max(masses))
    return encode_conditions(result)


def dynamic_prior_draw(data, name, seed, settings=None):
    """Return actual hierarchical prior draw before simulation (failures retain it)."""
    import pymc as pm
    from common.baseline_v import build_dynamic_model
    model=build_dynamic_model(data,name,settings)
    graph=pm.MCMC(list(model.get_stochastics()['node']))
    random.seed(seed);np.random.seed(seed)
    for generation in graph.generations:
        for node in sorted(generation,key=lambda n:n.__name__):
            if not node.observed: node.random()
    return natural_values(model)


def prior_parameters(data, name, draws=2000, seed=20261106, priors=None, *, dynamic=None, noncentered=False):
    """Joint parameter-only prior draws, parents before children on every draw.

    Uses native PyMC nodes. Does not generate trial replicas, condition on
    responses, truncate failed draws, or call MCMC.sample().
    """
    if draws < 1:
        raise ValueError("draws must be positive")
    import pymc as pm
    random.seed(seed); np.random.seed(seed)
    if dynamic is None:
        model = build_model(data, name, priors, noncentered=noncentered)
    else:
        from common.baseline_v import build_dynamic_model
        model = build_dynamic_model(data, name, dynamic)
    graph = pm.MCMC(list(model.get_stochastics()["node"]))
    rows = []
    for _ in range(draws):
        for generation in graph.generations:
            for node in sorted(generation, key=lambda n: n.__name__):
                if not node.observed:
                    node.random()
        values = natural_values(model)
        rows.append(values)
    return pd.DataFrame(rows)
