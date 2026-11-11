"""Part 1 vtz-derived M0 and its conditional-starting-point extension.
Native HDDMStimCoding; eight independent v/t cells, shared a/z, no dc.
This is a version-explicit reanalysis baseline, not an identical 0.8 posterior.
"""
from types import MappingProxyType
DESIGN_COLUMNS=("intercept","identity","self_valence","other_valence")
CONDITION_DESIGN=MappingProxyType({
    ("Good","Self"):(1.,0.,.5,0.),("Bad","Self"):(1.,0.,-.5,0.),
    ("Good","Other"):(1.,1.,0.,.5),("Bad","Other"):(1.,1.,0.,-.5)})
MODELS=("M0","M1")

def canonical_model(name):
    """Bv is a read-compatible historical label for the unchanged M0."""
    return "M0" if name in ("Bv", "bv") else name


def parameter_label(name):
    """Research labels at the read/display boundary; native HDDM nodes stay intact."""
    if name == 'bz': return 'mu_beta_z'
    if name == 'bz_std': return 'sigma_beta_z'
    if name.startswith('bz_eta_subj.'): return 'eta_beta_z.' + name[len('bz_eta_subj.'):]
    if name.startswith('bz_subj.'): return 'beta_z.' + name[len('bz_subj.'):]
    for raw in ('beta_dc', 'beta_v', 'log_kappa_a'):
        if name == raw: return 'mu_' + raw
        if name == raw + '_std': return 'sigma_' + raw
        if name.startswith(raw + '_subj.'): return raw + '.' + name.split('_subj.', 1)[1]
    return name


def canonical_parameters(values):
    """Rename a parameter mapping without dropping aliases or changing values."""
    labels = [parameter_label(str(key)) for key in values]
    if len(set(labels)) != len(labels):
        raise ValueError('Duplicate native/canonical parameter aliases')
    return dict(zip(labels, values.values()))


DEFAULT_PRIORS={"native_informative":True,"bz_group_sd":0.5,"bz_subject_sd":0.25}


def cell_key(parameter, identity, match, valence, subject):
    return f"{parameter}_subj({identity}.{match}.{valence}).{int(subject)}"


def build_model(data,name,priors=None,noncentered=False):
    """Native model; optional M1 noncentering preserves its joint scientific prior."""
    import numpy as np
    import pymc as pm
    import hddm
    from kabuki.hierarchical import Knode
    from common.data_utils import encode_conditions,validate_data
    name=canonical_model(name)
    if name not in MODELS:raise ValueError("Expected M0 or M1; historical H3 is not this baseline")
    if not isinstance(noncentered,bool):raise ValueError('noncentered must be boolean')
    if noncentered and name!='M1':raise ValueError('Noncentering applies only to the M1 starting effect')
    settings=dict(DEFAULT_PRIORS if priors is None else priors)
    if set(settings)!=set(DEFAULT_PRIORS) or settings['native_informative'] is not True:
        raise ValueError("This baseline uses the locked native informative prior family")
    if any(not np.isfinite(settings[k]) or settings[k]<=0 for k in ['bz_group_sd','bz_subject_sd']):
        raise ValueError("Positive finite starting-point prior scales required")
    frame=encode_conditions(data);validate_data(frame)
    class ConditionalStart(Knode):
        def create_node(self,node_name,kwargs,data):
            row=data.iloc[0]
            kwargs['parents']['contrast']=float(CONDITION_DESIGN[(row['val'],row['id'])][2])
            return self.pymc_node(name=node_name,**kwargs)
    class StartingPointModel(hddm.HDDMStimCoding):
        def _create_family_normal_normal_hnormal(self, parameter, value=0., g_mu=None,
                g_tau=15**-2, std_std=2., std_value=.1):
            if parameter!='bz' or not noncentered:
                return super()._create_family_normal_normal_hnormal(
                    parameter,value,g_mu,g_tau,std_std,std_value)
            from collections import OrderedDict
            group=Knode(pm.Normal,'bz',mu=g_mu,tau=g_tau,value=value,depends=self.depends[parameter])
            scale=Knode(pm.HalfNormal,'bz_std',tau=std_std**-2,value=std_value,
                depends=self.depends[parameter] if self.std_depends else ())
            eta=Knode(pm.Normal,'bz_eta_subj',mu=0.,tau=1.,value=0.,depends=('subj_idx',),subj=True)
            beta=Knode(pm.Deterministic,'bz_subj',eval=lambda mu,sigma,eta:mu+sigma*eta,
                mu=group,sigma=scale,eta=eta,depends=('subj_idx',),subj=True,trace=False,plot=False)
            return OrderedDict(bz=group,bz_std=scale,bz_eta=eta,bz_bottom=beta)

        def get_traces(self):
            traces=super().get_traces()
            if name=='M1' and noncentered:
                for subject in self.data.subj_idx.unique():
                    suffix='.'+str(int(subject))
                    traces['bz_subj'+suffix]=traces.bz+traces.bz_std*traces['bz_eta_subj'+suffix]
            return traces

        def _create_stochastic_knodes(self,include):
            nodes=super()._create_stochastic_knodes(include)
            if name=='M1':
                nodes.update(self._create_family_normal_normal_hnormal('bz',value=0.,g_mu=0.,
                    g_tau=settings['bz_group_sd']**-2,std_std=settings['bz_subject_sd']))
                def conditional_z(z,b,contrast):
                    return 1/(1+np.exp(-(np.log(z/(1-z))+b*contrast)))
                nodes['z_condition_bottom']=ConditionalStart(pm.Deterministic,'z_condition',
                    eval=conditional_z,z=nodes['z_bottom'],b=nodes['bz_bottom'],
                    depends=['val','id'],trace=False,plot=False)
            return nodes
        def _create_wfpt_parents_dict(self,nodes):
            parents=super()._create_wfpt_parents_dict(nodes)
            if name=='M1':parents['z']=nodes['z_condition_bottom']
            return parents
    # 1.0.1RC requires explicit v/a/t inclusion; historical 0.8 include='z' is not portable.
    model=StartingPointModel(frame,include=['v','a','t','z'],stim_col='stim',split_param='v',
        depends_on={'v':['match','val','id'],'t':['match','val','id']},
        drift_criterion=False,p_outlier=0,informative=True)
    model.research_name=name;model.research_priors=settings
    model.starting_effect_parameterization='noncentered' if noncentered else 'centered'
    return model


def natural_values(model):
    """Snapshot actual stochastic values plus the native transformed z value."""
    values={str(k):float(r['node'].value) for k,r in model.iter_stochastics()}
    for subject in model.data.subj_idx.unique():
        key='z_subj.'+str(int(subject));values[key]=float(model.nodes_db.loc[key,'node'].value)
        beta='bz_subj.'+str(int(subject))
        if beta in model.nodes_db.index:values[beta]=float(model.nodes_db.loc[beta,'node'].value)
    return canonical_parameters(values)


def trial_parameters(draw,data,name):
    import numpy as np
    import pandas as pd
    from common.data_utils import encode_conditions
    name=canonical_model(name)
    if name not in MODELS:raise ValueError(name)
    draw=canonical_parameters(dict(draw))
    frame=encode_conditions(data);result=pd.DataFrame(index=frame.index)
    for (subj,ident,match,val),rows in frame.groupby(['subj_idx','id','match','val']).groups.items():
        suffix=str(int(subj))
        result.loc[rows,'v']=float(draw[cell_key('v',ident,match,val,subj)])
        result.loc[rows,'t']=float(draw[cell_key('t',ident,match,val,subj)])
        result.loc[rows,'a']=float(draw['a_subj.'+suffix])
        z=float(draw['z_subj.'+suffix])
        if name=='M1':z=1/(1+np.exp(-(np.log(z/(1-z))+float(draw['beta_z.'+suffix])*CONDITION_DESIGN[(val,ident)][2])))
        result.loc[rows,'z']=z
    result['dc']=0.
    result['mu']=frame.sign*result.v
    if not np.isfinite(result.to_numpy()).all():raise ValueError('Nonfinite parameters')
    return result


def load_part1_reference(archive=None, load_posterior=True, candidate_data=None):
    """Read the original four VTZ chains, never refit or relabel old outputs.

    ``archive`` is the original Matching archive (or HU2020_MATCHING_ARCHIVE).
    Only load trusted local PyMC databases recorded in the reference manifest.
    Posterior availability is distinct from compatibility with a new analysis.
    """
    import os, json, pickle, hashlib
    from pathlib import Path
    import numpy as np
    import pandas as pd
    from common.data_utils import PROJECT, REPOSITORY, load_data, dataset_hash
    manifest = json.loads((PROJECT/'Output/tables/m0_part1_reference.json').read_text(encoding='utf-8'))
    root = Path(archive or os.environ.get('HU2020_MATCHING_ARCHIVE',
        str(REPOSITORY.parent.parent/'history/legacy_history/02_reproduction_workspace/nnu_hu_repro/matching')))
    original = REPOSITORY/'01_reproduction_2020/Data/matching_input.csv'
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    if sha(original) != manifest['input_sha256']:
        raise ValueError('Part 1 input differs from the recorded four-chain reference')
    for record in manifest['source_records']:
        if sha(REPOSITORY/record['path']) != record['sha256']:
            raise ValueError('Part 1 configuration/diagnostic record changed: '+record['path'])
    data = load_data(original, rt_min=None)  # Exact historical reference, not current empirical input.
    candidate = load_data() if candidate_data is None else candidate_data
    columns = ['subj_idx','match','val','id','stim','response','rt']
    same_input = dataset_hash(data[columns]) == dataset_hash(candidate[columns])
    info = dict(manifest, archive=str(root), input_matches_candidate=same_input,
        candidate_rows=len(candidate), reference_rows=len(data),
        common_prior_equivalence='not_verified_HDDM_0.8_to_1.0.1RC',
        observation_equivalence='requires_explicit_empirical_alignment',
        ready_for_empirical_comparison=False, posterior_loaded=False)
    frames = []
    if load_posterior:
        for record in manifest['chains']:
            path = root/record['database']
            if sha(path) != record['sha256']:
                raise ValueError('Archived trace differs from reference: '+str(path))
            with path.open('rb') as stream:
                db = pickle.load(stream, encoding='latin1')
            traces = {k: np.asarray(v[0]) for k,v in db.items()
                      if isinstance(v,dict) and set(v)=={0} and np.asarray(v[0]).ndim==1}
            frame = pd.DataFrame(traces).sort_index(axis=1)
            if len(frame)!=record['draws'] or not np.isfinite(frame.to_numpy()).all():
                raise ValueError('Invalid complete Part 1 posterior: '+str(path))
            if frames and list(frame.columns)!=list(frames[0].columns):
                raise ValueError('Part 1 posterior columns differ across chains')
            # Exercise every subject/cell, including natural-scale z, before reuse.
            trial_parameters(frame.iloc[0], data, 'M0')
            frames.append(frame)
        info.update(posterior_loaded=True, posterior_columns=len(frames[0].columns),
                    retained_draws=[len(f) for f in frames])
    return {'info':info, 'data':data, 'chains':frames}

# Dynamic extensions: same native hierarchy, mature PyDDM density and PyMC steps.
# These are explicit validation settings, not estimates from the data.
DYNAMIC_MODELS = ('M2', 'M3', 'M4')
DYNAMIC_DEFAULTS = dict(t_star=0.5, dx=0.0025, dt=0.0004, duration=8.0,
    mass_tolerance=1e-4, beta_group_sd=0.5, beta_subject_sd=0.25,
    log_kappa_group_mean=-0.6931471805599453, log_kappa_group_sd=1.0,
    log_kappa_subject_sd=0.5)


def dynamic_settings(settings=None, *, time_blocks=False):
    import numpy as np
    values = dict(DYNAMIC_DEFAULTS)
    if settings is not None:
        if set(settings) - set(values): raise ValueError('Unknown dynamic setting')
        values.update(settings)
    if not all(np.isfinite(v) for v in values.values()): raise ValueError('Nonfinite setting')
    for key in set(values) - {'log_kappa_group_mean'}:
        if values[key] <= 0: raise ValueError('Positive setting required: ' + key)
    if (values['dx'] > .02 or values['dt'] > .01
            or (not time_blocks and values['duration']/values['dt'] > 1000000)):
        raise ValueError('Grid outside explicit numerical budget')
    return values


def dynamic_solution(name, mu, a, z, effect=0., contrast=0., settings=None, require_unit_mass=True):
    """PyDDM implicit solve in u=tau/a^2, y=(X-a/2)/a coordinates.

    No nondecision-time overlay: t0 is added/subtracted outside the solver.
    A two-bin point source interpolates location without quantizing the z parameter.
    Density and simulation use the identical returned numerical distribution.
    """
    import numpy as np
    import pyddm
    from pyddm.models import Drift, InitialCondition, NoiseConstant, BoundConstant, BoundCollapsingExponential, OverlayNone
    from pathlib import Path
    import json
    marker=Path(pyddm.__file__).with_name('hu2020_tail_patch.json')
    if not marker.exists() or json.loads(marker.read_text()).get('implicit_remaining_mass_cutoff') != 1e-300:
        raise RuntimeError('Use the documented PyDDM 0.9.0 tail-precision environment')
    # Constant boundaries allow continuation from the previous undecided density.
    # The original million-step budget applies per C call, not per legal proposal.
    s = dynamic_settings(settings, time_blocks=name in ('M2', 'M4'))
    if name not in ('M0',) + DYNAMIC_MODELS: raise ValueError(name)
    if not np.isfinite([mu,a,z,effect,contrast]).all() or a <= 0 or not 0 < z < 1:
        raise ValueError('Invalid dynamic DDM parameters')
    if name == 'M3' and effect < 0: raise ValueError('kappa_a cannot be negative')
    scale = a*a
    if not np.isfinite(scale) or scale <= 0: raise ValueError('Invalid time scaling')
    class DriftTime(Drift):
        name = 'VTZ plus one independent time mechanism'
        required_parameters = ['base', 'coefficient', 'time_scale', 'kind', 'contrast', 't_star']
        def get_drift(self, t, **kwargs):
            g = -np.expm1(-t*self.time_scale/self.t_star)
            if self.kind == 'M2':
                return self.base + np.sqrt(self.time_scale)*self.coefficient*self.contrast*g
            if self.kind == 'M4':
                return self.base*np.exp(self.coefficient*self.contrast*g)
            return self.base
    class PointSource(InitialCondition):
        name = 'Linearly interpolated point source (grid approximation)'
        required_parameters = ['z']
        def get_IC(self, x, dx, **kwargs):
            target = self.z - .5
            upper = int(np.searchsorted(x, target, side='right'))
            if not 0 < upper < len(x): raise ValueError('Point source outside initial grid')
            weight = (target-x[upper-1])/(x[upper]-x[upper-1])
            result = np.zeros(len(x)); result[upper-1] = 1-weight; result[upper] = weight
            return result
    model = pyddm.Model(drift=DriftTime(base=mu*a, coefficient=effect, time_scale=scale,
        kind=name, contrast=contrast, t_star=s['t_star']), noise=NoiseConstant(noise=1.),
        bound=(BoundCollapsingExponential(B=.5, tau=effect*scale) if name=='M3'
               else BoundConstant(B=.5)),
        IC=PointSource(z=z), overlay=OverlayNone(), dx=s['dx'], dt=s['dt'],
        T_dur=s['duration'], choice_names=('upper','lower'))
    # Vectorize inputs to the very same pinned upstream C solver. The public
    # Model wrapper builds these arrays one scalar call at a time; no PDE code
    # or integration algorithm is implemented here. Equivalence is tested.
    grid=model.t_domain()
    drift=np.asarray(model.get_dependence('drift').get_drift(t=grid))
    if drift.ndim==0: drift=np.full(len(grid),float(drift))
    bound=np.asarray(.5*np.exp(-effect*scale*grid) if name=='M3' else [.5])
    from pyddm import csolve, Solution
    if name in ('M2', 'M4') and len(grid) > 1000001:
        # Continue the same upstream implicit solver, without restarting mass,
        # resetting physical time, coarsening the grid, or changing the likelihood.
        upper_raw = np.zeros(len(grid)); lower_raw = np.zeros(len(grid))
        undecided = model.IC({})
        for start in range(0, len(grid)-1, 1000000):
            stop = min(start+1000000, len(grid)-1)
            chunk = csolve.implicit_time(drift[start:stop+1],1,np.asarray([1.]),0,bound,0,
                undecided,(stop-start)*s['dt'],s['dt'],s['dx'],stop-start+1)
            # Index zero is the initial state, already covered by the prior block.
            upper_raw[start+1:stop+1] = chunk[0][1:]
            lower_raw[start+1:stop+1] = chunk[1][1:]
            undecided = chunk[2]
        raw = (upper_raw, lower_raw, undecided)
    else:
        raw=csolve.implicit_time(drift,1,np.asarray([1.]),0,bound,1 if name=='M3' else 0,
            model.IC({}),s['duration'],s['dt'],s['dx'],len(grid))
    if any(not np.isfinite(x).all() or (x < -1e-14).any() for x in raw):
        raise FloatingPointError('Invalid raw PyDDM solver output')
    solution=Solution(raw[0]*s['dt'],raw[1]*s['dt'],model,{},pdf_undec=raw[2])
    upper, lower = solution.pdf('upper'), solution.pdf('lower')
    mass = float((upper.sum()+lower.sum())*s['dt'])
    if (not np.isfinite(upper).all() or not np.isfinite(lower).all()
            or (upper < 0).any() or (lower < 0).any()
            or (require_unit_mass and abs(mass-1)>s['mass_tolerance'])):
        raise FloatingPointError('PyDDM mass/grid failure (not filtered): mass=%r' % mass)
    return solution, scale, mass


def dynamic_cell_loglik(signed_rt, name, mu, a, z, t0, effect=0., contrast=0., settings=None,
                        solution=None):
    """Joint choice/RT log density; never conditions on having responded.

    Linear density interpolation, Jacobian 1/a^2. Values outside support are -inf;
    numerical solver failures raise and must not silently reject a prior draw.
    """
    import numpy as np
    rt = np.asarray(signed_rt, dtype=float).reshape(-1)
    if not np.isfinite(rt).all() or not np.isfinite(t0) or t0 < 0: return np.full(len(rt), -np.inf)
    decision = np.abs(rt)-t0
    if (decision <= 0).any():
        result = np.full(len(rt), -np.inf)
        valid = decision > 0
        if valid.any(): result[valid] = dynamic_cell_loglik(rt[valid],name,mu,a,z,t0,effect,contrast,settings,solution)
        return result
    options=dynamic_settings(settings, time_blocks=name in ('M2', 'M4'))
    if solution is None:
        # Density up to observed times needs no renormalization and no solution
        # beyond those times. This is a finite integration domain, not truncation.
        horizon=min(160., max(.25, .25*np.ceil((float(decision.max())/(a*a)+options['dt'])/.25)))
        options['duration']=horizon
        solution=dynamic_solution(name,mu,a,z,effect,contrast,options,require_unit_mass=False)
    sol, scale, _ = solution
    time = sol.t_domain
    if (decision/scale > time[-1]).any():
        # Extend the numerical domain, never reject a valid parameter for the
        # arbitrary simulation window. At 160 normalized seconds, verify that
        # the remaining density has reached the explicit solver precision.
        options=dynamic_settings(settings, time_blocks=name in ('M2', 'M4'))
        required=float(decision.max()/scale)+options['dt']
        options['duration']=min(160.,max(options['duration'],8.*2**np.ceil(np.log2(max(1.,required/8.)))))
        sol,scale,_=dynamic_solution(name,mu,a,z,effect,contrast,options)
        time=sol.t_domain
        if (decision/scale>time[-1]).any() and np.sum(sol.undec)>1e-300:
            raise FloatingPointError('Tail remains unresolved after explicit domain extension')
    density = np.where(rt > 0,
        np.interp(decision/scale,time,sol.pdf('upper'),right=0.),
        np.interp(decision/scale,time,sol.pdf('lower'),right=0.))/scale
    with np.errstate(divide='ignore'): return np.log(density)


def dynamic_trial_parameters(draw, data, name):
    import numpy as np
    from common.data_utils import encode_conditions
    draw = canonical_parameters(dict(draw))
    frame = encode_conditions(data)
    result = trial_parameters(draw,frame,'M0')
    result['contrast'] = [.5 if ident=='Self' and val=='Good' else
                           -.5 if ident=='Self' else 0. for val,ident in zip(frame.val,frame.id)]
    effects = []
    for subj in frame.subj_idx:
        suffix = str(int(subj))
        if name == 'M0': effect = 0.
        elif name == 'M3': effect = np.exp(draw['log_kappa_a.'+suffix])
        elif name in ('M2','M4'): effect = draw[('beta_dc.' if name=='M2' else 'beta_v.')+suffix]
        else: raise ValueError(name)
        effects.append(effect)
    result['effect'] = effects
    return result


def build_dynamic_model(data, name, settings=None):
    """Keep HDDM's v/t/a/z prior graph; replace only observed first-passage density.

    PyDDM supplies the solver, HDDM supplies common priors and PyMC SliceStep /
    conjugate Normal steps. This is hierarchical inference, not PyDDM's MLE fit.
    M0 can use the same numerical engine for a fair paired comparison.
    """
    import numpy as np
    import pymc as pm
    import hddm
    from kabuki.hierarchical import Knode
    from common.data_utils import encode_conditions, validate_data
    from functools import lru_cache
    s = dynamic_settings(settings)
    if name not in ('M0',) + DYNAMIC_MODELS: raise ValueError(name)
    frame = encode_conditions(data); validate_data(frame)
    family = {'M2':'beta_dc', 'M3':'log_kappa_a', 'M4':'beta_v'}.get(name)
    # Small per-observed-node cache: t0 proposals reuse the decision distribution.
    class DynamicNode(Knode):
        def create_node(self, node_name, kwargs, data):
            row = data.iloc[0]
            sign = 2*int(row.stim)-1
            contrast = CONDITION_DESIGN[(row['val'],row['id'])][2]
            @lru_cache(maxsize=4)
            def solve(v,a,z,effect,horizon):
                return dynamic_solution(name,sign*v,a,z,effect,contrast,dict(s,duration=horizon),require_unit_mass=False)
            def logp(value, v, a, z, t, extension):
                if a <= 0 or not 0 < z < 1 or t < 0 or (np.abs(value)<=t).any(): return -np.inf
                effect = np.exp(extension) if name=='M3' else extension
                horizon=min(160.,max(.25,.25*np.ceil(((np.abs(value).max()-t)/(a*a)+s['dt'])/.25)))
                values = dynamic_cell_loglik(value,name,sign*float(v),float(a),float(z),float(t),
                    float(effect),contrast,s,solve(float(v),float(a),float(z),float(effect),horizon))
                return float(values.sum())
            distribution = pm.stochastic_from_dist('DynamicWiener', logp=logp, dtype=float, mv=True)
            # Initial values only; the native t prior and all sampled t nodes stay intact.
            kwargs['t'].value=min(.2, .5*float(np.abs(data.rt).min()))
            return distribution(node_name, **kwargs)
    class DynamicModel(hddm.HDDMStimCoding):
        def _create_stochastic_knodes(self, include):
            nodes = super()._create_stochastic_knodes(include)
            if family:
                center = s['log_kappa_group_mean'] if name=='M3' else 0.
                group_sd = s['log_kappa_group_sd'] if name=='M3' else s['beta_group_sd']
                subject_sd = s['log_kappa_subject_sd'] if name=='M3' else s['beta_subject_sd']
                nodes.update(self._create_family_normal_normal_hnormal(family,value=center,g_mu=center,
                    g_tau=group_sd**-2,std_std=subject_sd))
            return nodes
        def _create_wfpt_knode(self, nodes):
            return DynamicNode(pm.Stochastic,'wfpt',observed=True,col_name='rt',
                depends=['stim','val','id'],v=nodes['v_bottom'],a=nodes['a_bottom'],
                z=nodes['z_bottom'],t=nodes['t_bottom'],
                extension=nodes[family+'_bottom'] if family else 0.)
    model = DynamicModel(frame,include=['v','a','t','z'],stim_col='stim',split_param='v',
        depends_on={'v':['match','val','id'],'t':['match','val','id']},
        drift_criterion=False,p_outlier=0,informative=True)
    model.research_name=name
    model.research_priors=dict(DEFAULT_PRIORS,dynamic=s)
    return model
