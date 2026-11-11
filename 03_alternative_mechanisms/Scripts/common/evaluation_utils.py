"""Bounded native fitting and evaluation for the M0/M1 pair.
Independent snapshot of the reviewed V2 utilities; no runtime Part 2 imports.
No empirical fitting is launched on import.
"""
"""Explicitly gated, bounded native-chain fitting and recovery experiments."""
import os
for _key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_key,"1")
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import random
import sys
import time
import traceback
import numpy as np
import pandas as pd
from common.data_utils import dataset_hash, validate_data, write_json
from common.baseline_v import build_model, DEFAULT_PRIORS, MODELS, parameter_label, canonical_parameters

SCOPE = "valid-response subset; untruncated descriptive DDM, not a full or deadline-conditioned observation model"


def code_fingerprint():
    scripts=Path(__file__).resolve().parents[1]
    return {f.relative_to(scripts).as_posix():hashlib.sha256(f.read_bytes()).hexdigest()
            for f in sorted(scripts.rglob("*.py"))}


def _chain_job(data, name, chain_dir, config):
    out=Path(chain_dir)
    if out.exists():
        if not (out/"config.json").exists() or json.loads((out/"config.json").read_text()) != config:
            raise ValueError("Existing chain has different configuration: "+str(out))
        if (out/"complete.json").exists():
            marker=json.loads((out/"complete.json").read_text())
            if not (out/"traces.csv.gz").exists() or marker.get("trace_sha256") != hashlib.sha256((out/"traces.csv.gz").read_bytes()).hexdigest():
                raise ValueError("Completed trace missing or changed; use a new run directory")
            expected = len(range(((config["burn"] + config["thin"] - 1) // config["thin"]) * config["thin"],
                                 config["samples"], config["thin"]))
            if marker.get("draws") != expected:
                raise ValueError("Saved chain did not finish its requested sampling budget: " + str(out))
            return {"model":name,"chain":config["chain"],"status":"reused","path":str(out)}
        # Explicit rerun of an incomplete chain. No successful results removed.
    else:
        out.mkdir(parents=True)
        write_json(out/"config.json",config)
    started=time.time()
    try:
        import hddm, arviz
        random.seed(config["seed"]);np.random.seed(config["seed"])
        if config.get("dynamic") is not None:
            from common.baseline_v import build_dynamic_model
            model=build_dynamic_model(data,name,config["dynamic"])
        else:
            if config.get('starting_effect_parameterization')=='noncentered':
                model=build_model(data,name,config["priors"],noncentered=True)
            else:
                model=build_model(data,name,config["priors"])
        prior_manifest={}
        for label,row in model.iter_stochastics():
            node=row["node"]
            prior_manifest[str(label)]={"distribution":type(node).__name__,
                "parents":{str(k):str(v) for k,v in node.parents.items()}}
        write_json(out/"priors.json",prior_manifest)
        from common.baseline_v import cell_key
        for (subj,ident,match,val),part in data.groupby(['subj_idx','id','match','val']):
            model.nodes_db.loc[cell_key('t',ident,match,val,subj),'node'].value=min(.2,float(part.rt.min())*.5)
        model.mcmc(db="pickle", dbname=str(out/"model.db"))
        if not np.isfinite(model.logp):
            raise ValueError("Nonfinite initialized native model density")
        # Bounded native MAP initialization; dynamic PDE fits keep their checked start.
        if config.get("dynamic") is None:
            model.find_starting_values(individual_subjs=False,
                minimizer_kwargs={"maxiter": 10, "maxfev": 100})
        if not np.isfinite(model.logp):
            raise ValueError("Nonfinite density after initialization")
        model.sample(config["samples"],burn=config["burn"],thin=config["thin"],
            chains=1,parallel=False,return_infdata=False,loglike=False,ppc=False,
            db="pickle", dbname=str(out/"model.db"), save_name=str(out/"model"),
            progress_bar=False)
        if not (out/"model.hddm").is_file() or not (out/"model.db").is_file():
            raise RuntimeError("Native model/database save failed; traces are not marked complete")
        traces=model.get_traces()
        for subj in data.subj_idx.unique():
            suffix='.'+str(int(subj))
            if 'z_subj'+suffix not in traces:
                traces['z_subj'+suffix]=1/(1+np.exp(-traces['z_subj_trans'+suffix]))
        if len(traces)<1 or not np.isfinite(traces.to_numpy()).all():
            raise ValueError("Empty or nonfinite posterior traces")
        expected = len(range(((config["burn"] + config["thin"] - 1) // config["thin"]) * config["thin"],
                             config["samples"], config["thin"]))
        completed_iterations = getattr(model.mc, "_current_iter", None)
        if len(traces) != expected or completed_iterations != config["samples"]:
            traces.to_csv(out/"partial_traces.csv.gz", index=False)
            write_json(out/"interrupted.json", dict(model=name, chain=config["chain"],
                requested_iterations=config["samples"], completed_iterations=completed_iterations,
                expected_draws=expected, saved_draws=len(traces), complete=False))
            raise RuntimeError("Sampling interrupted: partial output preserved, not marked complete")
        traces.to_csv(out/"traces.csv.gz",index=False)
        runtime={}
        if config.get('dynamic') is not None:
            import pyddm
            runtime={'pyddm':pyddm.__version__,'tail_precision':json.loads(Path(pyddm.__file__).with_name('hu2020_tail_patch.json').read_text())}
        write_json(out/"complete.json",{"runtime":runtime,"model":name,"chain":config["chain"],
            "seed":config["seed"],"data_hash":config["data_hash"],"draws":len(traces),
            "elapsed_seconds":time.time()-started,"scope":SCOPE,"production":config["production"],
            "trace_sha256":hashlib.sha256((out/"traces.csv.gz").read_bytes()).hexdigest(),
            "hddm":hddm.__version__,"arviz":arviz.__version__,"python":sys.version,
            "scientific_validation_passed":False})
        return {"model":name,"chain":config["chain"],"status":"complete","path":str(out)}
    except Exception:
        write_json(out/"failed.json",{"error":traceback.format_exc(),"elapsed_seconds":time.time()-started})
        raise


def fit_models(data, output_root, models=("M0","M1"), chains=4,
               samples=80, burn=40, thin=1, workers=1, seed=20261106,
               allow_sampling=False, production=False, max_models=1, priors=None, dynamic=None,
               noncentered=None, incremental=False):
    """Each process runs exactly one native chain; never nested parallelism.

    Defaults are a smoke configuration, requiring explicit allow_sampling.
    More than 200 draws requires production=True; that flag is never implicit.
    max_models bounds the number of simultaneously active model families.
    Failed chains can be retried with identical settings; successful chains reused.
    incremental=True treats chains as the target total, permits a one-chain pilot,
    and keeps each chain's configuration independent of that target total.
    """
    if not allow_sampling:
        raise PermissionError("Sampling disabled: explicitly set allow_sampling=True")
    models=tuple(models)
    if noncentered is not None and not isinstance(noncentered,bool):
        raise ValueError('noncentered must be boolean or None')
    if noncentered and (dynamic is not None or 'M1' not in models):
        raise ValueError('Noncentering is available only for native M1')
    allowed=MODELS if dynamic is None else ('M0','M2','M3','M4')
    if dynamic is not None:
        from common.baseline_v import dynamic_settings
        dynamic=dynamic_settings(dynamic)
        if priors is not None and dict(priors)!=DEFAULT_PRIORS:
            raise ValueError('Dynamic models preserve locked common priors')
    if not models or len(set(models))!=len(models) or any(m not in allowed for m in models):
        raise ValueError("Unique M0/M1 model names required")
    if any(not isinstance(v,int) or v<1 for v in (chains,samples,thin,workers,max_models)) or not 0<=burn<samples:
        raise ValueError("Invalid chain/draw/concurrency settings")
    if workers > 16:
        raise ValueError("This study permits at most 16 concurrent MCMC chains")
    if not production and (samples>200 or chains>4):
        raise PermissionError("Beyond smoke limits: explicitly enable production")
    if production and chains<4 and not incremental:
        raise ValueError("Production requires at least four independent chains")
    validate_data(data)
    data=data.reset_index(drop=True)
    priors=dict(DEFAULT_PRIORS if priors is None else priors)
    root=Path(output_root);root.mkdir(parents=True,exist_ok=True)
    digest=dataset_hash(data);fingerprint=code_fingerprint()
    manifest={"data_hash":digest,"scope":SCOPE,"rows":len(data),"priors":priors}
    if dynamic is not None: manifest["dynamic"]=dynamic
    manifest_path=root/"data_manifest.json"
    if manifest_path.exists() and json.loads(manifest_path.read_text())!=manifest:
        raise ValueError("Output root already belongs to different data/prior settings")
    write_json(manifest_path,manifest)
    data.to_csv(root/"used_data.csv.gz",index=False)
    results=[]
    for start in range(0,len(models),max_models):
        jobs=[]
        for chain in range(chains):
            for name in models[start:start+max_models]:
                chain_seed=int(np.random.SeedSequence([seed,list(allowed).index(name),chain]).generate_state(1)[0])
                config=dict(model=name,chain=chain,seed=chain_seed,seed_base=seed,
                    chains=chains,samples=samples,burn=burn,thin=thin,priors=priors,
                    data_hash=digest,scope=SCOPE,production=production,code=fingerprint)
                if incremental:
                    config.pop("chains")
                    config["chain_layout"] = "incremental-v1"
                if dynamic is not None: config['dynamic']=dynamic
                if name=='M1' and noncentered is not None:
                    config['starting_effect_parameterization']='noncentered' if noncentered else 'centered'
                jobs.append((data,name,str(root/name/("chain_%02d"%chain)),config))
        # Reject incompatible existing chains before starting any missing chain.
        for _, _, directory, config in jobs:
            saved = Path(directory)
            if saved.exists() and (not (saved/"config.json").exists() or
                    json.loads((saved/"config.json").read_text()) != config):
                raise ValueError("Existing chain has different configuration: " + str(saved))
        for name in models[start:start+max_models]:
            if (root/name).exists() and list((root/name).glob("chain_*/complete.json")):
                posterior_stamp(root/name, allow_incomplete=True)
        if workers==1:
            for args in jobs:
                try:results.append(_chain_job(*args))
                except Exception as e:results.append({"model":args[1],"chain":args[3]["chain"],"status":"failed","error":str(e)})
        else:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                pending={pool.submit(_chain_job,*args):args for args in jobs}
                for future in as_completed(pending):
                    args=pending[future]
                    try:results.append(future.result())
                    except Exception as e:results.append({"model":args[1],"chain":args[3]["chain"],"status":"failed","error":str(e)})
    write_json(root/"run_status.json",{"results":results,"scope":SCOPE})
    failures=[r for r in results if r["status"]=="failed"]
    if failures:
        raise RuntimeError("Some chains failed; other models were retained: "+str(failures))
    return results


def parameter_recovery(run_dir, truth, output=None):
    """Coverage/errors for explicitly known parameters; no success inference."""
    files=sorted(Path(run_dir).glob("chain_*/traces.csv.gz"))
    if not files or any(not (f.parent/"complete.json").exists() for f in files):
        raise ValueError("Complete chains required")
    from common.evaluation_utils import load_chains
    frames,_=load_chains(run_dir)
    posterior=pd.concat(frames,ignore_index=True)
    rows=[]
    for parameter,value in canonical_parameters(truth).items():
        if parameter not in posterior:
            raise KeyError("Truth parameter absent from fit: "+parameter)
        x=posterior[parameter].to_numpy();lo,hi=np.quantile(x,[.025,.975])
        rows.append(dict(parameter=parameter,truth=float(value),mean=x.mean(),error=x.mean()-value,
                         lower=lo,upper=hi,covered=bool(lo<=value<=hi)))
    result=pd.DataFrame(rows)
    if output is not None:result.to_csv(output,index=False)
    return result

"""Post-hoc multi-chain diagnostics and response-subset PPC; never launches MCMC."""
from pathlib import Path
import json
import numpy as np
import pandas as pd


def load_chains(run_dir):
    """Validate original chains/configs; expose canonical labels without rewriting files."""
    paths = sorted(Path(run_dir).glob("chain_*"))
    if not paths:
        raise ValueError("No chain_* directories in " + str(run_dir))
    frames, configs = [], []
    for path in paths:
        if not (path / "complete.json").is_file():
            raise ValueError("Incomplete chain: " + str(path))
        config = json.loads((path / "config.json").read_text())
        if not all(k in config for k in ("seed", "data_hash")):
            raise ValueError("Chain config requires seed and data_hash")
        frame = pd.read_csv(path / "traces.csv.gz", float_precision="round_trip")
        labels = [parameter_label(str(c)) for c in frame.columns]
        if len(set(labels)) != len(labels):
            raise ValueError("Duplicate native/canonical parameter aliases in chain")
        frame.columns = labels
        marker = json.loads((path / "complete.json").read_text())
        for key in ("model", "chain", "seed", "data_hash"):
            if marker.get(key) != config.get(key):
                raise ValueError("Completion marker mismatch: " + key)
        expected = len(range(((config["burn"] + config["thin"] - 1) // config["thin"]) * config["thin"],
                             config["samples"], config["thin"]))
        if len(frame) != expected:
            raise ValueError("Incomplete retained chain: " + str(path))
        if marker.get("draws") != len(frame):
            raise ValueError("Completion marker draw count mismatch")
        if len(frame) < 4 or frame.empty or not np.isfinite(frame.to_numpy(dtype=float)).all():
            raise ValueError("Need at least four finite retained draws")
        if frames and (list(frame.columns) != list(frames[0].columns) or len(frame) != len(frames[0])):
            raise ValueError("Chain lengths or parameter columns differ")
        frames.append(frame)
        configs.append(config)
    if len({c["seed"] for c in configs}) != len(configs):
        raise ValueError("Chains must use independent seeds")
    expected_chains = len(configs) if configs[0].get("chain_layout") == "incremental-v1" else configs[0]["chains"]
    if len(configs) != expected_chains or {c["chain"] for c in configs} != set(range(expected_chains)):
        raise ValueError("Missing or duplicate configured chains")
    for key in ("data_hash", "model", "priors", "samples", "burn", "thin", "chains", "chain_layout", "code", "production", "scope", "dynamic", "starting_effect_parameterization"):
        if len({json.dumps(c.get(key), sort_keys=True) for c in configs}) > 1:
            raise ValueError("Inconsistent chain configuration: " + key)
    # Incremental runs carry immutable checksums; historical runs remain readable.
    if configs[0].get("chain_layout") == "incremental-v1":
        posterior_stamp(run_dir)
    return frames, configs


def legacy_part1_ppc(archive):
    """Read historical aggregated PPC with corrected labels; never edit its source.

    The historical 'accuracy' statistic is mean(response), not correctness.
    Aggregates cannot recover predictive accuracy without trial-level predictions.
    """
    path=Path(archive)/'ppc/m1_p0_ppc_s500_compare.csv'
    result=pd.read_csv(path)
    result['historical_stat']=result['stat']
    result['stat']=result['stat'].replace({'accuracy':'upper_boundary_choice_proportion'})
    result=result.rename(columns={'MSE':'statistic_prediction_mse'})
    result['source']=str(path)
    return result


def inference_data(run_dir):
    import arviz as az
    frames, configs = load_chains(run_dir)
    posterior = {col: np.stack([f[col].to_numpy() for f in frames]) for col in frames[0]}
    result = az.from_dict(posterior=posterior)
    result.attrs.update(model_id=configs[0]["model"], data_hash=configs[0]["data_hash"],
                        chain_configs=json.dumps(configs, sort_keys=True))
    return result, configs


def summarize_chains(run_dir, smoke=False, output_dir=None, plots=True):
    """Rank R-hat/ESS gates: smoke runs are always scientifically unqualified."""
    import arviz as az
    idata, configs = inference_data(run_dir)
    summary = az.summary(idata, kind="all", round_to=None)
    # Explicit high-precision metrics: never decide a threshold using rounded summary.
    rhat = az.rhat(idata, method="rank")
    bulk = az.ess(idata, method="bulk")
    tail = az.ess(idata, method="tail")
    for col in summary.index:
        summary.loc[col, "r_hat"] = float(rhat[col])
        summary.loc[col, "ess_bulk"] = float(bulk[col])
        summary.loc[col, "ess_tail"] = float(tail[col])
    metrics = summary[["r_hat", "ess_bulk", "ess_tail"]]
    passed = (np.isfinite(metrics).all(axis=1) & (metrics.r_hat < 1.01)
              & (metrics.ess_bulk >= 400) & (metrics.ess_tail >= 400))
    smoke = bool(smoke or any(not c.get("production", False) or c.get("smoke", False) or c.get("mode") == "smoke" for c in configs))
    gate = dict(chains=len(configs), smoke=smoke, diagnostic_pass=bool(passed.all()),
                scientific_pass=bool(len(configs) >= 4 and passed.all() and not smoke),
                failed_parameters=summary.index[~passed].tolist(),
                scope="valid-response subset; convergence does not establish model adequacy")
    out = Path(output_dir) if output_dir else Path(run_dir) / "diagnostics"
    out.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out / "summary.csv")
    (out / "gate.json").write_text(json.dumps(gate, indent=2))
    if plots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        # Bound figure size for large participant-level models; all remain in summary.
        names = [v for v in idata.posterior.data_vars if "_subj." not in v][:12]
        if not names:
            names = list(idata.posterior.data_vars)[:12]
        az.plot_trace(idata, var_names=names, compact=True)
        plt.savefig(out / "trace.png", dpi=120, bbox_inches="tight")
        plt.close("all")
    return summary, gate


def behavioral_metrics(data, fast_threshold=0.3):
    """Joint choice/RT diagnostics, including errors and empty fast-response cells."""
    records = []
    for (val, identity, stim), part in data.groupby(["val", "id", "stim"], sort=True):
        correct = part.response.to_numpy() == stim
        rt = part.rt.to_numpy()
        metrics = {"p_match": float(part.response.mean()), "accuracy": float(correct.mean()),
                   "mean_rt": float(rt.mean()), "rt_q10": float(np.quantile(rt, .1)),
                   "rt_q50": float(np.quantile(rt, .5)), "rt_q90": float(np.quantile(rt, .9)),
                   "mean_correct_rt": float(rt[correct].mean()) if correct.any() else np.nan,
                   "mean_error_rt": float(rt[~correct].mean()) if (~correct).any() else np.nan,
                   "p_fast": float((rt <= fast_threshold).mean())}
        fast = part[part.rt <= fast_threshold]
        metrics["fast_p_match"] = float(fast.response.mean()) if len(fast) else np.nan
        for metric, value in metrics.items():
            records.append(dict(val=val, id=identity, stim=stim, metric=metric, value=value))
    # Preserve direction stratification: do not confound stimulus mix with bias.
    for stim in sorted(data.stim.unique()):
        fast = data[(data.id == "Self") & (data.stim == stim) & (data.rt <= fast_threshold)]
        probs = fast.groupby("val").response.mean()
        value = float(probs["Good"]-probs["Bad"]) if {"Good", "Bad"}.issubset(probs.index) else np.nan
        records.append(dict(val="Good-minus-Bad", id="Self", stim=stim,
                            metric="fast_match_advantage", value=value))
    return pd.DataFrame(records)


def posterior_predictive(run_dir, data, name, output_dir, draws=100, seed=9182,
                         fast_threshold=0.3, save_replicates=False):
    """Simulate one replicate at a time; retain compact metrics, not trial arrays."""
    from common.data_utils import dataset_hash
    from common.simulation_utils import simulate_from_parameters
    from common.baseline_v import trial_parameters
    if not isinstance(draws, int) or draws < 1 or fast_threshold <= 0:
        raise ValueError("Positive draws and fast_threshold required")
    frames, configs = load_chains(run_dir)
    if configs[0]["data_hash"] != dataset_hash(data):
        raise ValueError("PPC data differs from fitted data")
    if configs[0].get("model", name) != name:
        raise ValueError("Model name differs from chain configuration")
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    summaries, selections, subject_summaries = [], [], []
    subject = sorted(data.subj_idx.unique())[0]  # Fixed example, not selected for fit quality.
    for replicate in range(draws):
        chain = int(rng.integers(len(frames)))
        index = int(rng.integers(len(frames[chain])))
        sim_seed = int(rng.integers(1, 2**31-1))
        if configs[0].get('dynamic') is not None:
            from common.simulation_utils import simulate_dynamic
            simulated=simulate_dynamic(data,name,frames[chain].iloc[index],sim_seed,configs[0]['dynamic'])
        else:
            pars = trial_parameters(frames[chain].iloc[index], data, name)
            simulated = simulate_from_parameters(data, pars, seed=sim_seed)
        summary = behavioral_metrics(simulated, fast_threshold)
        summary["replicate"] = replicate
        summaries.append(summary)
        subject_summaries.append(behavioral_metrics(
            simulated.loc[simulated.subj_idx.eq(subject)], fast_threshold).assign(replicate=replicate))
        selections.append(dict(replicate=replicate, chain=chain, draw=index, seed=sim_seed))
        if save_replicates:
            simulated.to_csv(out / ("replicate_%04d.csv.gz" % replicate), index=False)
    predictions = pd.concat(summaries, ignore_index=True)
    keys = ["val", "id", "stim", "metric"]
    interval = predictions.groupby(keys).value.agg(
        predicted_mean="mean", predicted_lower=lambda x: x.quantile(.025),
        predicted_upper=lambda x: x.quantile(.975), valid_replicates="count").reset_index()
    observed = behavioral_metrics(data, fast_threshold).rename(columns={"value": "observed"})
    result = observed.merge(interval, on=keys, how="left", validate="one_to_one")
    result["outside_interval"] = ((result.observed < result.predicted_lower)
                                   | (result.observed > result.predicted_upper))
    result.loc[result[["observed", "predicted_lower", "predicted_upper"]].isna().any(axis=1), "outside_interval"] = np.nan
    result.to_csv(out / "ppc_summary.csv", index=False)
    subject_interval = pd.concat(subject_summaries).groupby(keys).value.agg(
        predicted_mean="mean", predicted_lower=lambda x: x.quantile(.025),
        predicted_upper=lambda x: x.quantile(.975), valid_replicates="count").reset_index()
    subject_observed = behavioral_metrics(data.loc[data.subj_idx.eq(subject)], fast_threshold).rename(columns={"value":"observed"})
    subject_observed.merge(subject_interval, on=keys, validate="one_to_one").assign(
        subj_idx=subject).to_csv(out/"ppc_example_subject.csv", index=False)
    predictions.to_csv(out / "ppc_metric_draws.csv.gz", index=False)
    (out / "ppc_manifest.json").write_text(json.dumps(dict(model=name, seed=seed,
        dataset_hash=dataset_hash(data), draws=draws, fast_threshold_seconds=fast_threshold,
        selections=selections, scope="untruncated valid-response subset; no timeout predictions"), indent=2))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    for ax, metric in zip(axes, ("p_match", "mean_rt")):
        part = result[result.metric == metric].reset_index(drop=True)
        x = np.arange(len(part))
        ax.plot(x, part.observed, "ko", label="Observed")
        ax.vlines(x, part.predicted_lower, part.predicted_upper, color="tab:blue")
        ax.plot(x, part.predicted_mean, "o", color="tab:blue", label="PPC mean / 95%")
        ax.set_xticks(x)
        ax.set_xticklabels([f"{r.val}/{r.id}/{r.stim}" for r in part.itertuples()], rotation=60, ha="right")
        ax.set_title(metric)
        ax.legend()
    fig.tight_layout()
    fig.savefig(out / "ppc.png", dpi=120)
    plt.close(fig)
    return result

def pointwise_loglik(draw, data, model, dynamic=None):
    """Native HDDM joint choice/RT density, evaluated without truncation."""
    import hddm
    from common.baseline_v import trial_parameters
    if dynamic is not None: return dynamic_loglik(draw,data,model,dynamic)
    validate_data(data, require_design=False)
    params = trial_parameters(draw, data, model)
    signed_rt = data.rt.to_numpy() * (2*data.response.to_numpy()-1)
    result = np.empty(len(data))
    groups = params.reset_index(drop=True).groupby(["mu", "a", "t", "z"], sort=False).groups
    for (mu, a, t, z), ix in groups.items():
        ix = np.asarray(ix, dtype=int)
        result[ix] = hddm.wfpt.pdf_array(signed_rt[ix], float(mu), 0., float(a),
                                        float(z), 0., float(t), 0., logp=True, p_outlier=0)
    if np.isnan(result).any() or np.isposinf(result).any():
        raise ValueError("Invalid native likelihood values")
    return result


def score_heldout(traces, train_data, test_data, model, max_draws=500,
                  seed=20261106, loglik_path=None, dynamic=None):
    """Streaming posterior log-mean density; optional explicitly requested disk array."""
    if set(test_data.subj_idx)-set(train_data.subj_idx):
        raise ValueError("New subjects require a different predictive target")
    if "event_id" not in train_data or "event_id" not in test_data:
        raise ValueError("Stable event IDs required")
    if set(train_data.event_id.astype(str)) & set(test_data.event_id.astype(str)):
        raise ValueError("Train/test leakage: overlapping events")
    if not len(traces) or max_draws < 1:
        raise ValueError("Posterior draws required")
    draws = np.sort(np.random.default_rng(seed).choice(len(traces), min(max_draws, len(traces)), replace=False))
    disk = None
    if loglik_path is not None:
        path = Path(loglik_path)
        if path.exists(): raise FileExistsError(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        disk = np.lib.format.open_memmap(path, mode="w+", dtype="float64", shape=(len(draws),len(test_data)))
    total = np.full(len(test_data), -np.inf)
    for j, row in enumerate(draws):
        values = pointwise_loglik(traces.iloc[row], test_data, model, dynamic)
        total = np.logaddexp(total, values)
        if disk is not None: disk[j] = values
    if disk is not None: disk.flush()
    result = test_data[["event_id", "subj_idx", "val", "id", "stim"]].copy()
    result["model"] = model
    result["log_predictive_density"] = total-np.log(len(draws))
    result["posterior_draws"] = len(draws)
    return result


def compare_scores(scores, diagnostics_pass, allow_unvalidated=False):
    """Paired ELPD differences; report both event and subject-cluster uncertainty."""
    if not diagnostics_pass and not allow_unvalidated:
        raise ValueError("Refusing scientific ranking: sampling diagnostics have not passed")
    if scores.duplicated(["model", "event_id"]).any():
        raise ValueError("Duplicate predictions")
    pivot = scores.pivot(index="event_id", columns="model", values="log_predictive_density")
    if pivot.isna().any().any() or np.isposinf(pivot.to_numpy()).any():
        raise ValueError("Incomplete or invalid common-event predictions")
    if np.isneginf(pivot.to_numpy()).any():
        if not allow_unvalidated:
            raise ValueError("Zero posterior predictive density: no scientific ranking")
        return pd.DataFrame([dict(model=model, elpd=float(pivot[model].sum()),
            reference=None, delta_elpd=np.nan, se_event=np.nan, se_subject_cluster=np.nan,
            n_events=len(pivot), n_subjects=scores.subj_idx.nunique(),
            valid_for_inference=False, status="ZERO_PREDICTIVE_SUPPORT_NO_RANKING")
            for model in pivot.columns])
    subjects = scores.groupby("event_id").subj_idx.agg(lambda x: x.nunique())
    if not subjects.eq(1).all(): raise ValueError("Event subject mismatch")
    subject_ids = scores.drop_duplicates("event_id").set_index("event_id").subj_idx.reindex(pivot.index)
    totals = pivot.sum()
    reference = totals.idxmax()
    rows = []
    for model in pivot:
        delta = pivot[model]-pivot[reference]
        cluster = delta.groupby(subject_ids).sum()
        rows.append(dict(model=model, elpd=float(totals[model]), reference=reference,
            delta_elpd=float(delta.sum()), se_event=float(np.sqrt(len(delta)*delta.var(ddof=1))) if len(delta)>1 else None,
            se_subject_cluster=float(np.sqrt(len(cluster)*cluster.var(ddof=1))) if len(cluster)>1 else None,
            n_events=len(delta), n_subjects=len(cluster),
            valid_for_inference=bool(diagnostics_pass),
            status="diagnostic_gate_passed" if diagnostics_pass else "SMOKE_ONLY_NO_SCIENTIFIC_RANKING"))
    return pd.DataFrame(rows).sort_values("elpd", ascending=False).reset_index(drop=True)


def load_run_traces(model_dir):
    """Validate completed chains against actual fit data before prediction."""
    model_dir = Path(model_dir)
    train = pd.read_csv(model_dir.parent/"used_data.csv.gz", float_precision="round_trip")
    frames, configs = load_chains(model_dir)
    if configs[0]["data_hash"] != dataset_hash(train) or configs[0]["model"] != model_dir.name:
        raise ValueError("Chain fingerprint/model mismatch")
    return pd.concat(frames,ignore_index=True), train


def posterior_stamp(run_dir, allow_incomplete=False):
    """Identify the actual completed chains, not the requested future chain count."""
    entries, configs = [], []
    for directory in sorted(Path(run_dir).glob("chain_*")):
        marker_path = directory/"complete.json"
        if not marker_path.exists():
            if allow_incomplete:
                continue
            raise ValueError("Incomplete chain: " + str(directory))
        config = json.loads((directory/"config.json").read_text())
        marker = json.loads(marker_path.read_text())
        filename = "posterior.nc" if config.get("engine") == "hssm" else "traces.csv.gz"
        artifact = directory/filename
        digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
        if digest != marker.get("trace_sha256"):
            raise ValueError("Saved posterior has changed: " + str(artifact))
        for key in ("model", "chain", "seed", "data_hash"):
            if marker.get(key) != config.get(key):
                raise ValueError("Completion marker mismatch: " + key)
        if config.get("engine") == "hssm":
            expected = config["draws"]
        else:
            expected = len(range(((config["burn"]+config["thin"]-1)//config["thin"])*config["thin"],
                                 config["samples"], config["thin"]))
        if marker.get("draws") != expected:
            raise ValueError("Incomplete sampling budget: " + str(directory))
        entries.append(dict(config=config, sha256=digest))
        configs.append(config)
    if not configs:
        raise ValueError("No completed chains: " + str(run_dir))
    if len({c["seed"] for c in configs}) != len(configs):
        raise ValueError("Chain seeds are not independent")
    if not allow_incomplete and [c["chain"] for c in configs] != list(range(len(configs))):
        raise ValueError("Missing chain indices")
    shared = [{k: v for k, v in c.items() if k not in ("chain", "seed")} for c in configs]
    if any(c != shared[0] for c in shared[1:]):
        raise ValueError("Chains have incompatible configurations")
    return dict(n_chains=len(configs), posterior_hash=hashlib.sha256(
        json.dumps(entries, sort_keys=True).encode()).hexdigest())


def fit_hssm_chains(model, run_dir, config, chains, initial, allow_sampling=False):
    """Append independent HSSM chains; an existing complete chain is never sampled again.

    Chains within this notebook run sequentially. Different model notebooks may
    run concurrently. Config excludes the requested total number of chains.
    """
    if not allow_sampling:
        raise PermissionError("Sampling disabled")
    if not isinstance(chains, int) or chains < 1 or config["draws"] < 4 or config["tune"] < 1:
        raise ValueError("Invalid HSSM chain budget")
    root = Path(run_dir)
    jobs = []
    for chain in range(chains):
        cfg = dict(config, engine="hssm", chain_layout="incremental-v1", chain=chain,
                   seed=int(np.random.SeedSequence([config["seed_base"], 3, chain]).generate_state(1)[0]))
        path = root/("chain_%02d" % chain)
        if path.exists() and (not (path/"config.json").exists() or
                json.loads((path/"config.json").read_text()) != cfg):
            raise ValueError("Existing HSSM chain has different configuration: " + str(path))
        jobs.append((path, cfg))
    if list(root.glob("chain_*/complete.json")):
        posterior_stamp(root, allow_incomplete=True)
    results = []
    for path, cfg in jobs:
        if (path/"complete.json").exists():
            results.append(dict(chain=cfg["chain"], status="reused"))
            continue
        path.mkdir(parents=True, exist_ok=True)
        write_json(path/"config.json", cfg)
        import copy
        trace = model.sample(sampler="pymc", initvals=copy.deepcopy(initial),
            random_seed=cfg["seed"], chains=1, cores=1, tune=cfg["tune"], draws=cfg["draws"],
            target_accept=cfg["target_accept"], idata_kwargs={"log_likelihood": False},
            include_response_params=False)
        p = trace["posterior"].ds
        if p.sizes["chain"] != 1 or p.sizes["draw"] != cfg["draws"] or any(
                not np.isfinite(v.values).all() for v in p.data_vars.values()):
            raise ValueError("Incomplete or nonfinite HSSM posterior; not marked complete")
        trace.attrs.update(chain_config=json.dumps(cfg, sort_keys=True))
        temporary = path/"posterior.partial.nc"
        trace.to_netcdf(temporary)
        temporary.replace(path/"posterior.nc")
        write_json(path/"complete.json", dict(model=cfg["model"], chain=cfg["chain"], seed=cfg["seed"],
            data_hash=cfg["data_hash"], draws=cfg["draws"],
            trace_sha256=hashlib.sha256((path/"posterior.nc").read_bytes()).hexdigest()))
        results.append(dict(chain=cfg["chain"], status="complete"))
    return results


def load_hssm_chains(run_dir, config):
    """Read immutable per-chain DataTrees and join only groups carrying a chain axis."""
    import xarray as xr
    posterior_stamp(run_dir)
    traces = []
    for path in sorted(Path(run_dir).glob("chain_*")):
        cfg = json.loads((path/"config.json").read_text())
        if any(cfg.get(k) != v for k, v in config.items()):
            raise ValueError("HSSM posterior does not match current configuration")
        with xr.open_datatree(path/"posterior.nc", engine="h5netcdf") as source:
            trace = source.load().copy(deep=True)
        if json.loads(trace.attrs["chain_config"]) != cfg:
            raise ValueError("HSSM file configuration mismatch")
        if trace["posterior"].sizes["draw"] != cfg["draws"]:
            raise ValueError("HSSM saved draw count mismatch")
        traces.append(trace)
    groups = {}
    if any(t.groups != traces[0].groups for t in traces):
        raise ValueError("HSSM posterior groups differ")
    for group in traces[0].groups:
        if group == "/":
            continue
        datasets = [t[group].ds for t in traces]
        if "chain" in datasets[0].dims:
            groups[group] = xr.concat([d.assign_coords(chain=[i]) for i, d in enumerate(datasets)],
                                      dim="chain", join="exact", compat="equals")
        else:
            if any(not d.equals(datasets[0]) for d in datasets[1:]):
                raise ValueError("HSSM constant/observed data differ")
            groups[group] = datasets[0]
    return xr.DataTree.from_dict(groups)


PREDICTIVE_TARGET = "heldout-existing-subjects:joint-choice-rt-density-seconds:untruncated:no-lapse:v1"


def save_exploratory_scores(table, root, model, train, test, diagnostics_pass=False,
                            validation_pass=False, limitations=""):
    """Publish one self-contained score file atomically; never publish partial scores."""
    root = Path(root)
    if not table.event_id.is_unique or set(table.event_id) != set(test.event_id):
        raise ValueError("Score rows must cover the exact held-out events once")
    if not table.model.eq(model).all():
        raise ValueError("Wrong model label in score table")
    if set(train.event_id) & set(test.event_id):
        raise ValueError("Train/test leakage")
    if any(np.isnan(table[k]).any() or np.isposinf(table[k]).any()
           for k in ("log_predictive_density", "lpd_first_half", "lpd_second_half")):
        raise ValueError("Invalid predictive densities")
    stamp = posterior_stamp(root/model)
    table = table.copy()
    for k, v in dict(**stamp, train_hash=dataset_hash(train), test_hash=dataset_hash(test),
                     predictive_target=PREDICTIVE_TARGET, diagnostics_pass=bool(diagnostics_pass),
                     validation_pass=bool(validation_pass), limitations=limitations).items():
        table[k] = v
    import tempfile
    with tempfile.NamedTemporaryFile(dir=root, suffix=".csv.gz", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        table.to_csv(temporary, index=False, compression="gzip")
        temporary.replace(root/"heldout_scores.csv.gz")
    finally:
        temporary.unlink(missing_ok=True)
    return table


def collect_exploratory_comparison(roots, train, test):
    """Show five-model progress; rank only complete, compatible, current scores.

    Split-draw differences measure integration sensitivity, not posterior mixing.
    This exploratory table does not certify mechanisms even after adding chains.
    """
    tables, progress = [], []
    for model, root in roots.items():
        path = Path(root)/"heldout_scores.csv.gz"
        if not path.exists():
            progress.append(dict(model=model, status="awaiting score"))
            continue
        try:
            stamp = posterior_stamp(Path(root)/model)
            table = pd.read_csv(path, float_precision="round_trip")
            expected = dict(**stamp, train_hash=dataset_hash(train), test_hash=dataset_hash(test),
                            predictive_target=PREDICTIVE_TARGET, model=model)
            if any(k not in table or not table[k].eq(v).all() for k, v in expected.items()):
                raise ValueError("stale posterior, data or scoring target; rescore this model")
            if not table.event_id.is_unique or set(table.event_id) != set(test.event_id):
                raise ValueError("different held-out events")
            actual = table.set_index("event_id").subj_idx.sort_index()
            expected_subjects = test.set_index("event_id").subj_idx.sort_index()
            if not actual.equals(expected_subjects):
                raise ValueError("different subject mapping")
            tables.append(table)
            progress.append(dict(model=model, status="scored", n_chains=stamp["n_chains"]))
        except (ValueError, KeyError, FileNotFoundError) as error:
            progress.append(dict(model=model, status=str(error)))
    progress = pd.DataFrame(progress)
    if len(tables) != len(roots):
        return progress, None
    scores = pd.concat(tables, ignore_index=True)
    result = compare_scores(scores, diagnostics_pass=False, allow_unvalidated=True)
    finite = np.isfinite(scores[["log_predictive_density", "lpd_first_half", "lpd_second_half"]]).all().all()
    if finite:
        first = scores.pivot(index="event_id", columns="model", values="lpd_first_half")
        second = scores.pivot(index="event_id", columns="model", values="lpd_second_half")
        best = result.iloc[0].model
        result["split_mc_change"] = [abs((first[m]-first[best]).sum() -
                                         (second[m]-second[best]).sum()) for m in result.model]
        result["integration_stable"] = result.split_mc_change.le(.1*result.se_subject_cluster)
        result["status"] = "EXPLORATORY_COMPARISON"
    else:
        result["status"] = "NONFINITE_PREDICTIVE_SUPPORT_NO_RANKING"
        result["integration_stable"] = False
    metadata = scores.groupby("model")[["n_chains", "posterior_draws", "diagnostics_pass", "validation_pass", "limitations"]].first()
    result = result.merge(metadata, on="model", validate="one_to_one")
    # A good rank never overrides numerical/observation/recovery limitations.
    result["valid_for_inference"] = False
    return progress, result


def score_run(model_dir, test_data, model=None, **kwargs):
    model_dir = Path(model_dir)
    model = model or model_dir.name
    if model != model_dir.name: raise ValueError("Model directory mismatch")
    traces, train = load_run_traces(model_dir)
    _, configs=load_chains(model_dir)
    stored=configs[0].get('dynamic')
    if 'dynamic' in kwargs and kwargs['dynamic']!=stored:
        raise ValueError('Scoring numerical configuration differs from fit')
    if stored is not None: kwargs['dynamic']=stored
    scored = score_heldout(traces, train, test_data, model, **kwargs)
    scored["lpd"] = scored["log_predictive_density"]
    return scored


def dynamic_loglik(draw, data, model, settings):
    from common.baseline_v import dynamic_trial_parameters, dynamic_cell_loglik
    validate_data(data,require_design=False)
    params=dynamic_trial_parameters(draw,data,model).reset_index(drop=True)
    signed=data.rt.to_numpy()*(2*data.response.to_numpy()-1)
    result=np.empty(len(data))
    for (mu,a,z,t0,effect,contrast),ix in params.groupby(['mu','a','z','t','effect','contrast']).groups.items():
        ix=np.asarray(ix,dtype=int)
        result[ix]=dynamic_cell_loglik(signed[ix],model,mu,a,z,t0,effect,contrast,settings)
    if np.isnan(result).any() or np.isposinf(result).any(): raise FloatingPointError('Invalid dynamic density')
    return result


def dynamic_numerical_checks(model, settings=None):
    """Paired off-state and mesh-refinement checks; return every numerical failure.

    L1 tolerance .02 and choice-probability tolerance .005 are prefit screening
    thresholds, not a guarantee of full-data likelihood accuracy. Posterior-range
    grid sensitivity is still required before any empirical interpretation.
    """
    import hddm
    from common.baseline_v import dynamic_solution,dynamic_settings
    s=dynamic_settings(settings);fine=dict(s,dx=s['dx']/2,dt=s['dt']/2)
    rows=[]
    for case,(mu,a,z) in enumerate([(-3.,.7,.35),(0.,1.3,.55),(3.,2.,.7)]):
        for effect in (0.,.4,1.):
            row=dict(model=model,case=case,mu=mu,a=a,z=z,effect=effect,contrast=.5)
            try:
                sol,scale,mass=dynamic_solution(model,mu,a,z,effect,.5,s)
                ref,_,_=dynamic_solution(model,mu,a,z,effect,.5,fine)
                grid=sol.t_domain
                row['mass']=mass
                row['refinement_l1']=sum(np.abs(sol.pdf(c)-np.interp(grid,ref.t_domain,ref.pdf(c))).sum()*s['dt'] for c in ['upper','lower'])
                row['refinement_p_upper']=abs(sol.prob('upper')-ref.prob('upper'))
                if effect==0:
                    native_upper=hddm.wfpt.pdf_array(grid*scale,mu,0.,a,z,0.,0.,0.,err=1e-8)*scale
                    native_lower=hddm.wfpt.pdf_array(-grid*scale,mu,0.,a,z,0.,0.,0.,err=1e-8)*scale
                    row['off_l1']=float((abs(sol.pdf('upper')-native_upper).sum()+abs(sol.pdf('lower')-native_lower).sum())*s['dt'])
                    row['off_p_upper']=abs(sol.prob('upper')-native_upper.sum()*s['dt'])
                row['passed']=bool(row['refinement_l1']<.02 and row['refinement_p_upper']<.005
                    and row.get('off_l1',0)<.02 and row.get('off_p_upper',0)<.005)
            except Exception as exc:
                row.update(passed=False,error=repr(exc))
            rows.append(row)
    return pd.DataFrame(rows)


def dynamic_target_contrasts(run_dir, data, output=None):
    """Diagnose target contrasts themselves, including natural M3 rates."""
    import arviz as az
    from common.baseline_v import cell_key
    frames,configs=load_chains(run_dir);derived=[]
    for frame in frames:
        values={}
        for match in ('match','mismatch'):
            per_subject=[frame[cell_key('v','Self',match,'Good',s)]-frame[cell_key('v','Self',match,'Bad',s)]
                         for s in sorted(data.subj_idx.unique())]
            values['mean_delta_v_'+match]=np.asarray(per_subject).mean(axis=0)
        if configs[0]['model']=='M3':
            for subject in data.subj_idx.unique():
                values['kappa_a.'+str(int(subject))]=np.exp(frame['log_kappa_a.'+str(int(subject))])
        derived.append(pd.DataFrame(values))
    posterior={col:np.stack([f[col].to_numpy() for f in derived]) for col in derived[0]}
    idata=az.from_dict(posterior=posterior)
    summary=az.summary(idata,round_to=None)
    for column in summary.index:
        summary.loc[column,'r_hat']=float(az.rhat(idata)[column])
        summary.loc[column,'ess_bulk']=float(az.ess(idata,method='bulk')[column])
        summary.loc[column,'ess_tail']=float(az.ess(idata,method='tail')[column])
    if output is not None: summary.to_csv(output)
    return summary
