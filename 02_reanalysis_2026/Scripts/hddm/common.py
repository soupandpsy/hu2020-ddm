"""Audited single-chain execution and bounded task-batch scheduling."""
import os
for _key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[_key]="1"
import importlib.util
import argparse, hashlib, importlib, json, random, subprocess, sys, time, traceback
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
OUTPUT=Path(os.environ.get("HU2026_OUTPUT",str(ROOT/"Output"))).expanduser().resolve()
TASK=os.environ.get("HU2026_TASK", "matching")
TASK_DIR=HERE/TASK
FIGURES_OUTPUT=None
def settings():
    cfg=json.loads((TASK_DIR/(TASK+"_config.json")).read_text())
    if os.environ.get('HU2026_LONG_REFIT')=='1':
        cfg.update(samples=60000,burn=30000,seed_base=202611210,refit_version='independent_60000_30000')
    return cfg
def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path, obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(obj,indent=2,default=str),encoding="utf-8")
def update_status(stage,workers):
    p=OUTPUT/"hddm/execution_status.json"
    status=json.loads(p.read_text()) if p.exists() else {}
    cfg=settings()
    done=[str(p.relative_to(OUTPUT)) for p in (OUTPUT/"hddm").glob("*/*/screening/r01/complete.json")]
    final=[str(p.relative_to(OUTPUT)) for p in (OUTPUT/"hddm").glob("*/*/final/r*/complete.json")]
    status.update(stage=stage,actual_parallel_limit=workers,screening_completed=len(done),
        final_extra_chains_completed=len(final),screening_completion_files=done,
        final_completion_files=final,updated_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()))
    write(p,status)
def module(task):
    spec=importlib.util.spec_from_file_location(task+"_models",HERE/task/(task+"_models.py"))
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
def build_model(data,name,cfg,conditions):
    import numpy as np
    import hddm
    if name not in cfg["models"]:raise ValueError(name)
    if set(data.stim.unique())!={0,1}:raise ValueError("Both stimulus directions required")
    if not data.groupby(conditions).stim.nunique().eq(2).all():
        raise ValueError("v conditions must contain both directions: v/dc aliasing")
    dep={"v":conditions}
    if name=="m4":dep["dc"]=conditions
    inc=["a","v","t"]+(["z"] if name in ("m1","m3","m4") else [])
    return hddm.HDDMStimCoding(data,stim_col="stim",split_param="v",
        drift_criterion=name in ("m2","m3","m4"),include=inc,depends_on=dep,
        p_outlier=cfg["p_outlier"],informative=cfg["informative"])
def load_data(task,sensitivity=False):
    import pandas as pd
    import numpy as np
    cfg=settings();data=pd.read_csv(ROOT/"Data"/(task+"_input.csv"))
    assert len(data)==cfg["fit_rows"][task] and data.subj_idx.nunique()==cfg["subjects"][task]
    assert not data[["subj_idx","rt","stim","response","val","id"]].isna().any().any()
    assert np.isfinite(data.rt).all() and data.rt.gt(0).all()
    assert set(data.response.unique())=={0,1} and set(data.stim.unique())=={0,1}
    assert ((data.response==data.stim).astype(int)==data.accuracy).all()
    if sensitivity:data=data[data.rt>.2]
    return data.reset_index(drop=True)
def fingerprint(task,model,cfg,sensitivity=False):
    files=[HERE/"common.py",TASK_DIR/(TASK+"_config.json"),HERE/task/(task+"_models.py"),ROOT/"Data"/(task+"_input.csv")]
    return {"task":task,"model":model,"sensitivity":sensitivity,
      "files":{str(p.relative_to(ROOT)):digest(p) for p in files},"settings":cfg,"code_layout_version":4}
def chain_dir(task,model,phase,chain):
    base=OUTPUT/"hddm"/task
    if phase=="sensitivity":return base/"sensitivity"/model/("r%02d"%chain)
    return base/model/phase/("r%02d"%chain)

def sample(task,name,chain,phase,output=None):
    import numpy as np
    import pandas as pd

    cfg=settings();smoke=phase=="smoke";probe=phase=="probe";sensitivity=phase=="sensitivity"
    out=Path(output) if output else chain_dir(task,name,phase,chain)
    seed=cfg["seed_base"]+cfg["tasks"].index(task)*1000+cfg["models"].index(name)*100+chain
    np.random.seed(seed);random.seed(seed)
    fp=fingerprint(task,name,cfg,sensitivity)
    if phase=="final":
        prior=OUTPUT/"hddm"/task/name/"screening/r01"
        assert json.loads((prior/"fingerprint.json").read_text())==fp,"r01 settings mismatch"
        assert (prior/"complete.json").exists()
    out.mkdir(parents=True,exist_ok=False)
    write(out/"fingerprint.json",fp)
    write(out/"settings.json",dict(cfg,task=task,model=name,chain=chain,seed=seed,phase=phase))
    data=load_data(task,sensitivity)
    if smoke:
        data=data[data.subj_idx.isin(sorted(data.subj_idx.unique())[:2])]
        data=data.groupby(["subj_idx","val","id","stim"],group_keys=False).head(12).reset_index(drop=True)
    data.to_csv(out/"used_data.csv.gz",index=False)
    start=time.time()
    try:
        m=module(task).build(data,name,cfg)
        # Save actual stochastic node names, distributions and priors for audit.
        priors={}
        for label,row in m.iter_stochastics():
            node=row["node"]
            priors[str(label)]={"distribution":type(node).__name__,
              "parents":{str(k):str(v) for k,v in node.parents.items()}}
        write(out/"priors.json",priors)
        assert any(str(x).startswith("z") for x in priors)==(name in ("m1","m3","m4"))
        assert any(str(x).startswith("dc") for x in priors)==(name in ("m2","m3","m4"))
        m.find_starting_values()
        m.sample(60 if (smoke or probe) else cfg["samples"],burn=30 if (smoke or probe) else cfg["burn"],
                 thin=cfg["thin"],db="pickle",dbname=str(out/"posterior.db"))
        traces=natural_traces(m.get_traces())
        assert np.isfinite(traces.to_numpy()).all()
        traces.to_csv(out/"traces.csv.gz",index=False)
        m.gen_stats().to_csv(out/"summary.csv")
        m.save(str(out/"model.hddm"))
        single_diagnostics(traces,out)
        ppc(m,data,out,2 if (smoke or probe) else cfg["ppc_draws_screening"],seed+500000)
        import resource
        write(out/"complete.json",{"draws":len(traces),"dic":float(m.dic),"seed":seed,
          "elapsed_seconds":time.time()-start,"max_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
          "phase":phase,"eligible_as_r01":phase=="screening"})
        print("COMPLETE",out,flush=True)
    except Exception:
        (out/"error.txt").write_text(traceback.format_exc())
        raise
def schedule(jobs,workers):
    active=[];pending=list(jobs);failed=[]
    phase=jobs[0][3];update_status(phase,workers)
    while pending or active:
        while pending and len(active)<workers:
            task,name,chain,phase=pending.pop(0)
            out=chain_dir(task,name,phase,chain)
            if out.exists():raise FileExistsError(out)
            out.parent.mkdir(parents=True,exist_ok=True)
            log=out.parent/("r%02d.stdout.log"%chain);err=out.parent/("r%02d.stderr.log"%chain)
            stdout=log.open("x");stderr=err.open("x")
            cmd=[sys.executable,"-u",str(HERE/task/(task+"_fit.py")),"--model",name,"--chain",str(chain),"--phase",phase]
            proc=subprocess.Popen(cmd,stdout=stdout,stderr=stderr)
            active.append((proc,stdout,stderr,(task,name,chain,phase)))
            print("START",proc.pid,task,name,chain,phase,flush=True)
        for item in active[:]:
            proc,stdout,stderr,job=item
            if proc.poll() is not None:
                stdout.close();stderr.close();active.remove(item)
                print("EXIT",job,proc.returncode,flush=True)
                if proc.returncode:failed.append(job);pending=[]
                update_status(phase if not failed else "failed_"+phase,workers)
        if active:time.sleep(2)
    if failed:raise RuntimeError("Failed jobs: "+str(failed))

"""Chain diagnostics and PPC; never compute multi-chain Rhat from one chain."""
import json, argparse
from pathlib import Path


def natural_traces(traces):
    """Keep stochastic traces and add inverse-logit z on its interpretable scale."""
    from scipy.special import expit
    traces=traces.copy()
    for c in list(traces):
        if c=="z_trans" or c.startswith("z_trans(") or c.startswith("z_subj_trans"):
            traces[c.replace("_trans","")]=expit(traces[c].to_numpy())
    return traces

def single_diagnostics(traces,out):
    import arviz as az
    vals=traces.to_numpy()
    rows=[]
    for col in traces:
        x=traces[col].to_numpy()[None,:]
        rows.append({"parameter":col,"ess_bulk":float(az.ess(x,method="bulk")),
          "ess_tail":float(az.ess(x,method="tail")),"mcse_mean":float(az.mcse(x,method="mean"))})
    pd.DataFrame(rows).to_csv(out/"single_chain_diagnostics.csv",index=False)

def ppc(model,data,out,draws,seed):
    import hddm
    np.random.seed(seed)
    pred=hddm.utils.post_pred_gen(model,samples=draws,append_data=False,parallel=False)
    pred.to_csv(out/"ppc.csv.gz")
    rows=[]
    # The library returns signed RT; choice must be recovered from its sign.
    for (node,draw),part in pred.groupby(level=[0,1]):
        indices=part.index.get_level_values(-1)
        obs=data.loc[indices]
        signed=part["rt"].to_numpy()
        for stim in sorted(obs.stim.unique()):
            mask=obs.stim.to_numpy()==stim
            pr=signed[mask]; ob=obs.loc[mask]
            row={"node":node,"draw":draw,"stim":int(stim),"n":len(ob),
                "observed_choice":float(ob.response.mean()),"predicted_choice":float((pr>0).mean())}
            for q in (.1,.5,.9):
                row["observed_q"+str(q)]=float(ob.rt.quantile(q))
                row["predicted_q"+str(q)]=float(np.quantile(np.abs(pr),q))
            rows.append(row)
    checks=pd.DataFrame(rows)
    checks.to_csv(out/"ppc_summary.csv",index=False)
    write(out/"ppc_metrics.json",{
      "choice_rmse":float(np.sqrt(np.average((checks.observed_choice-checks.predicted_choice)**2,weights=checks.n))),
      "rt_median_rmse_seconds":float(np.sqrt(np.average((checks["observed_q0.5"]-checks["predicted_q0.5"])**2,weights=checks.n))),
      "draws":draws,"seed":seed,"scope":"node/subject/condition checks; inspection still required"})

def screening(task=None):
    rows=[]
    for t in ([task] if task else settings()["tasks"]):
        for m in settings()["models"]:
            p=OUTPUT/"hddm"/t/m/"screening/r01"
            row={"task":t,"model":m,"chains":0,"status":"not_run","rhat":np.nan}
            if (p/"complete.json").exists():
                c=json.loads((p/"complete.json").read_text())
                d=pd.read_csv(p/"single_chain_diagnostics.csv")
                pp=json.loads((p/"ppc_metrics.json").read_text())
                row.update(chains=1,status="screening_only",dic=c["dic"],ess_bulk_min=d.ess_bulk.min(),
                  ess_tail_min=d.ess_tail.min(),choice_rmse=pp["choice_rmse"],rt_median_rmse_seconds=pp["rt_median_rmse_seconds"])
            rows.append(row)
    for t in sorted(set(r['task'] for r in rows)):
        dest=OUTPUT/'hddm'/t;dest.mkdir(parents=True,exist_ok=True)
        frame=pd.DataFrame([r for r in rows if r['task']==t])
        target=dest/'model_comparison.csv'
        if target.exists():raise FileExistsError(target)
        frame.to_csv(target,index=False)
        for model,part in frame.groupby('model'):
            (dest/model).mkdir(exist_ok=True)
            part.to_csv(dest/model/'screening_summary.csv',index=False)
    print("Comparison updated. Select one provisional model per task after reviewing traces, PPC, recovery and theory.")

def final(task_filter=None,model_name="m4"):
    import arviz as az
    import hddm
    selection={task_filter:{"model":model_name}} if task_filter else json.loads((OUTPUT/"hddm/selection.json").read_text())
    summaries=[]
    for task,sel in selection.items():
        if task_filter and task!=task_filter:continue
        model=sel["model"];base=OUTPUT/"hddm"/task/model
        paths=[base/"screening/r01"]+[base/"final"/("r%02d"%r) for r in (2,3,4)]
        fp=json.loads((paths[0]/"fingerprint.json").read_text())
        assert all(json.loads((p/"fingerprint.json").read_text())==fp for p in paths)
        assert all((p/"complete.json").exists() for p in paths)
        ts=[natural_traces(pd.read_csv(p/"traces.csv.gz")) for p in paths]
        assert all(list(x.columns)==list(ts[0].columns) and len(x)==len(ts[0]) for x in ts)
        arr=np.stack([x.to_numpy() for x in ts])
        dataset=az.from_dict(posterior={c:arr[:,:,i] for i,c in enumerate(ts[0].columns)})
        dest=base/"final/diagnostics";dest.mkdir(exist_ok=False)
        diag=az.summary(dataset,round_to="none")
        diag.to_csv(dest/"diagnostics.csv")
        dataset.to_netcdf(dest/"posterior.nc")
        passed=bool((diag.r_hat<settings()["rhat_max"]).all() and (diag.ess_bulk>=400).all() and (diag.ess_tail>=400).all())
        write(dest/"chain_sources.json",{"paths":[str(p.relative_to(OUTPUT)) for p in paths],"diagnostics_pass":passed,"selection_provisional":True})
        # Report only parameters actually estimated by this model.
        for c in ts[0]:
            if "_subj" in c or "_std" in c or "trans" in c:continue
            if not c.startswith(("v","z","dc")):continue
            values=np.concatenate([x[c].to_numpy() for x in ts])
            baseline=.5 if c=="z" else 0.
            summaries.append({"task":task,"model":model,"parameter":c,"mean":values.mean(),
              "q025":np.quantile(values,.025),"q975":np.quantile(values,.975),
              "p_above_baseline":float((values>baseline).mean()),"baseline":baseline,
              "diagnostics_pass":passed,"interpretation":"exploratory; model selection uncertainty not propagated"})
        comparisons=[]
        if task=="matching":
            for par in ("v","dc"):
                comparisons.extend([(f"{par}(Self.Good)",f"{par}(Self.Bad)"),
                                    (f"{par}(Self.Good)",f"{par}(Other.Good)")])
        elif task=="identity":
            comparisons=[(f"{par}(Good)",f"{par}(Bad)") for par in ("v","dc")]
        else:
            comparisons=[(f"{par}(Self)",f"{par}(Other)") for par in ("v","dc")]
        for left,right in comparisons:
            if left not in ts[0] or right not in ts[0]:continue
            delta=np.concatenate([(x[left]-x[right]).to_numpy() for x in ts])
            summaries.append({"task":task,"model":model,"parameter":left+" - "+right,
              "mean":delta.mean(),"q025":np.quantile(delta,.025),"q975":np.quantile(delta,.975),
              "p_above_baseline":float((delta>0).mean()),"baseline":0,
              "diagnostics_pass":passed,"interpretation":"constrained-model contrast; not identical to 2020 unconstrained v contrast"})
    for t in sorted(set(r['task'] for r in summaries)):
        pd.DataFrame([r for r in summaries if r['task']==t]).to_csv(OUTPUT/'hddm'/t/'parameter_summary.csv',index=False)

def matching_review():
    """Precise diagnostics, contrasts and explicitly labelled PPC conventions."""
    import arviz as az

    base=OUTPUT/'hddm/matching/m4'
    dest=base/'final/diagnostics'
    target=dest/'matching_precise_diagnostics.csv'
    if target.exists():raise FileExistsError(target)
    paths=[base/'screening/r01']+[base/'final'/('r%02d'%i) for i in (2,3,4)]
    fps=[json.loads((p/'fingerprint.json').read_text()) for p in paths]
    assert all(f==fps[0] for f in fps)
    ts=[natural_traces(pd.read_csv(p/'traces.csv.gz')) for p in paths]
    assert all(len(t)==5000 and list(t)==list(ts[0]) for t in ts)
    arrays={c:np.stack([t[c].to_numpy() for t in ts]) for c in ts[0]}
    diag=az.summary(az.from_dict(posterior=arrays),round_to='none')
    diag['mcse_over_sd']=diag.mcse_mean/diag.sd
    diag.to_csv(target)
    passed=bool((diag.r_hat<1.01).all() and (diag.ess_bulk>=400).all() and (diag.ess_tail>=400).all())
    params={c:x for c,x in arrays.items() if c.startswith(('v(','dc(','z')) and not any(s in c for s in ('subj','std','trans'))}
    for par in ('v','dc'):
        sg=arrays[par+'(Self.Good)'];sb=arrays[par+'(Self.Bad)'];og=arrays[par+'(Other.Good)'];ob=arrays[par+'(Other.Bad)']
        params.update({par+':Self.Good-Self.Bad':sg-sb,par+':Other.Good-Other.Bad':og-ob,
            par+':Self.Good-Other.Good':sg-og,par+':Self.Bad-Other.Bad':sb-ob,
            par+':interaction':(sg-sb)-(og-ob)})
    rows=[]
    for name,x in params.items():
        ref=.5 if name=='z' else 0.
        rows.append(dict(task='matching',model='m4',parameter=name,mean=x.mean(),sd=x.std(),
          q025=np.quantile(x,.025),q975=np.quantile(x,.975),p_above_baseline=(x>ref).mean(),baseline=ref,
          r_hat=float(az.rhat(x)),ess_bulk=float(az.ess(x,method='bulk')),ess_tail=float(az.ess(x,method='tail')),
          mcse_mean=float(az.mcse(x,method='mean')),diagnostics_pass=passed,
          interpretation='Conditional on provisional M4 and shared-direction v; Match-choice bias, not Good-self preference'))
    pd.DataFrame(rows).to_csv(dest/'matching_parameter_contrasts.csv',index=False)
    data=load_data('matching');rng=np.random.default_rng(2026102142);checks=[]
    for chain in range(1,5):
        pred=pd.read_csv(dest/('ppc_r%02d'%chain)/'ppc.csv.gz',index_col=[0,1,2])
        for draw,part in pred.groupby(level=1):
            idx=part.index.get_level_values(-1).astype(int)
            obs=data.loc[idx].reset_index(drop=True);signed=part.rt.to_numpy()
            for convention,max_rt in [('ddm_component',None),('contamination_abs_uniform_0_5',5.),('native_contamination_abs_uniform_0_10',10.)]:
                sim=signed.copy()
                if max_rt is not None:
                    mask=rng.random(len(sim))<.05
                    sim[mask]=rng.uniform(0,max_rt,mask.sum())*rng.choice([-1,1],mask.sum())
                for (ident,val,stim),indices in obs.groupby(['id','val','stim']).groups.items():
                    ix=np.asarray(indices);ob=obs.loc[ix];pr=sim[ix]
                    row=dict(chain=chain,draw=draw,convention=convention,id=ident,val=val,stim=stim,n=len(ix),
                        observed_choice=ob.response.mean(),predicted_choice=(pr>0).mean(),
                        observed_above_1_2=(ob.rt>1.2).mean(),predicted_above_1_2=(abs(pr)>1.2).mean())
                    for q in (.1,.5,.9):
                        row['observed_q'+str(q)]=ob.rt.quantile(q)
                        row['predicted_q'+str(q)]=np.quantile(abs(pr),q)
                    checks.append(row)
    pd.DataFrame(checks).to_csv(dest/'matching_condition_ppc.csv',index=False)
    predictive=pd.DataFrame(checks);summary_rows=[]
    for (conv,ident,val,stim),part in predictive.groupby(['convention','id','val','stim']):
        for metric in ['choice','q0.1','q0.5','q0.9','above_1_2']:
            values=part['predicted_'+metric]
            summary_rows.append(dict(convention=conv,id=ident,val=val,stim=stim,metric=metric,
              observed=part['observed_'+metric].iloc[0],predicted=values.mean(),q025=values.quantile(.025),q975=values.quantile(.975)))
    pd.DataFrame(summary_rows).to_csv(dest/'matching_ppc_summary.csv',index=False)

    write(dest/'matching_review.json',dict(diagnostics_pass=passed,max_rhat=diag.r_hat.max(),
        min_bulk_ess=diag.ess_bulk.min(),min_tail_ess=diag.ess_tail.min(),
        max_mcse_over_sd=diag.mcse_over_sd.max(),draws_per_chain=5000,chains=4,
        ppc_seed=2026102142,ppc_draws_per_chain=50,
        ppc_caveat='Two normalized contamination conventions are predictive sensitivity checks, not exact full-likelihood PPC. Installed likelihood has constant outlier floor without finite support; native simulator differs in signed density. No fitted likelihood changed.',
        ess_reference='https://mc-stan.org/learn-stan/diagnostics-warnings.html'))
    print('Matching precise diagnostics',passed,diag.r_hat.max(),diag.ess_bulk.min(),diag.ess_tail.min(),flush=True)

def matching_sensitivity():
    """Compare completed RT>.2 chains with the unchanged main fit."""
    import arviz as az
    base=OUTPUT/'hddm/matching/m4'
    paths=[OUTPUT/'hddm/matching/sensitivity/m4'/('r%02d'%i) for i in range(1,5)]
    if not all((p/'complete.json').exists() for p in paths):
        raise RuntimeError('All four sensitivity chains must finish before analysis')
    fps=[json.loads((p/'fingerprint.json').read_text()) for p in paths]
    if not all(f==fps[0] for f in fps):raise ValueError('Sensitivity fingerprints differ')
    ts=[natural_traces(pd.read_csv(p/'traces.csv.gz')) for p in paths]
    if not all(len(t)==5000 and list(t)==list(ts[0]) for t in ts):
        raise ValueError('Unexpected trace dimensions')
    arrays={c:np.stack([t[c].to_numpy() for t in ts]) for c in ts[0]}
    if not all(np.isfinite(x).all() for x in arrays.values()):raise ValueError('Nonfinite traces')
    dest=OUTPUT/'hddm/matching/sensitivity/m4/diagnostics';dest.mkdir(exist_ok=False)
    diag=az.summary(az.from_dict(posterior=arrays),round_to='none')
    diag['mcse_over_sd']=diag.mcse_mean/diag.sd
    diag.to_csv(dest/'diagnostics.csv')
    passed=bool((diag.r_hat<1.01).all() and (diag.ess_bulk>=400).all() and (diag.ess_tail>=400).all())
    for par in ('v','dc'):
        sg=arrays[par+'(Self.Good)'];sb=arrays[par+'(Self.Bad)'];og=arrays[par+'(Other.Good)'];ob=arrays[par+'(Other.Bad)']
        arrays.update({par+':Self.Good-Self.Bad':sg-sb,par+':Other.Good-Other.Bad':og-ob,
          par+':Self.Good-Other.Good':sg-og,par+':Self.Bad-Other.Bad':sb-ob,par+':interaction':sg-sb-og+ob})
    main=pd.read_csv(base/'final/diagnostics/matching_parameter_contrasts.csv')
    rows=[]
    for row in main.to_dict('records'):
        x=arrays[row['parameter']];ref=row['baseline']
        rows.append(dict(parameter=row['parameter'],main_mean=row['mean'],main_q025=row['q025'],main_q975=row['q975'],
          sensitivity_mean=x.mean(),sensitivity_q025=np.quantile(x,.025),sensitivity_q975=np.quantile(x,.975),
          sensitivity_p_above_baseline=(x>ref).mean(),baseline=ref,mean_shift=x.mean()-row['mean'],
          sensitivity_rhat=float(az.rhat(x)),sensitivity_ess_bulk=float(az.ess(x,method='bulk')),
          sensitivity_ess_tail=float(az.ess(x,method='tail')),sensitivity_mcse=float(az.mcse(x,method='mean')),
          sensitivity_diagnostics_pass=passed))
    pd.DataFrame(rows).to_csv(dest/'main_comparison.csv',index=False)
    write(dest/'review.json',dict(diagnostics_pass=passed,max_rhat=diag.r_hat.max(),min_bulk_ess=diag.ess_bulk.min(),
      min_tail_ess=diag.ess_tail.min(),chain_sources=[str(p.relative_to(OUTPUT)) for p in paths],
      interpretation='Additional RT>.2 robustness check; original author HDDM preprocessing retained fast responses. No cross-dataset chain merging.'))
    print('Matching sensitivity diagnostics',passed,flush=True)


"""Plot saved diagnostics without rerunning inference."""
from pathlib import Path
import argparse


def trace_plot_main():
    p=argparse.ArgumentParser();p.add_argument("--task",required=True);p.add_argument("--model",required=True)
    a=p.parse_args();base=OUTPUT/"hddm"/a.task/a.model
    paths=[base/"screening/r01"]+[base/"final"/("r%02d"%r) for r in (2,3,4)]
    paths=[p for p in paths if (p/"traces.csv.gz").exists()]
    if not paths:raise ValueError("No real posterior traces")
    ts=[natural_traces(pd.read_csv(p/"traces.csv.gz")) for p in paths]
    cols=[c for c in ts[0] if c.startswith(("v","z","dc")) and not any(s in c for s in ("_subj","_std","trans"))]
    out=FIGURES_OUTPUT or OUTPUT/"hddm"/a.task/"figures";out.mkdir(parents=True,exist_ok=True)
    target=out/(a.task+"_"+a.model+"_traces.png")
    if target.exists():raise FileExistsError(target)
    fig,axs=plt.subplots(len(cols),2,figsize=(12,2.3*len(cols)),squeeze=False)
    for i,c in enumerate(cols):
        for j,t in enumerate(ts):
            axs[i,0].plot(t[c],alpha=.5,label="r%02d"%(j+1))
            axs[i,1].hist(t[c],bins=40,density=True,histtype="step")
        axs[i,0].set_ylabel(c)
    axs[0,0].legend()
    fig.suptitle(a.task+" "+a.model+" — "+str(len(paths))+" chain(s); exploratory")
    fig.tight_layout();fig.savefig(target,dpi=160);plt.close(fig)

def matching_summary_plot():
    base=OUTPUT/'hddm/matching/m4/final/diagnostics'
    out=FIGURES_OUTPUT or OUTPUT/'hddm/matching/figures';out.mkdir(parents=True,exist_ok=True)
    table=pd.read_csv(base/'matching_parameter_contrasts.csv')
    target=out/'matching_m4_parameter_contrasts.png'
    if target.exists():raise FileExistsError(target)
    fig,axes=plt.subplots(1,3,figsize=(16,6),gridspec_kw={'width_ratios':[2,2,1]})
    for ax,par in zip(axes,['v','dc','z']):
        d=table[table.parameter.str.startswith(par+'(')|table.parameter.str.startswith(par+':')] if par!='z' else table[table.parameter=='z']
        y=np.arange(len(d));ax.errorbar(d['mean'],y,xerr=np.vstack([d['mean']-d.q025,d.q975-d['mean']]),fmt='o',capsize=3)
        ax.set_yticks(y);ax.set_yticklabels(d.parameter);ax.axvline(.5 if par=='z' else 0,color='gray',ls='--');ax.set_title(par+' posterior: 95% interval');ax.invert_yaxis()
    fig.suptitle('Matching M4: provisional model; conditional estimates, not causal preference')
    fig.tight_layout();fig.savefig(target,dpi=160);plt.close(fig)
    summary=pd.read_csv(base/'matching_ppc_summary.csv')
    fig,axes=plt.subplots(2,2,figsize=(14,9))
    conventions=summary.convention.unique();colors=['#2171b5','#d95f02','#238b45']
    for ax,metric in zip(axes.flat,['choice','q0.1','q0.5','q0.9']):
        for j,conv in enumerate(conventions):
            d=summary[(summary.metric==metric)&(summary.convention==conv)].sort_values(['id','val','stim'])
            y=np.arange(len(d));ax.errorbar(y+(j-1)*.16,d.predicted,yerr=np.vstack([d.predicted-d.q025,d.q975-d.predicted]),fmt='o',ms=3,capsize=2,color=colors[j],label=conv)
        ax.scatter(y,d.observed,color='black',marker='x',label='observed')
        ax.set_xticks(y);ax.set_xticklabels([f'{r.id}/{r.val}/stim{r.stim}' for r in d.itertuples()],rotation=45,ha='right',fontsize=8)
        ax.set_title(metric+(' (seconds)' if metric.startswith('q') else ' (proportion)'))
    axes[0,0].legend(fontsize=7);fig.suptitle('Matching predictive checks: component and two explicit contamination conventions')
    fig.tight_layout();fig.savefig(out/'matching_m4_ppc.png',dpi=160);plt.close(fig)


"""Minimum synthetic parameter-recovery gate; independent of empirical r01."""
from pathlib import Path
import argparse, json, random


def preflight_main():
    p=argparse.ArgumentParser();p.add_argument("--samples",type=int,default=2000)
    p.add_argument("--burn",type=int,default=1000)
    p.add_argument("--resume",action="store_true")
    a=p.parse_args()
    out=OUTPUT/"recovery/preflight";out.mkdir(parents=True,exist_ok=a.resume)
    if (out/"complete.json").exists():raise FileExistsError(out/"complete.json")
    cfg=settings();all_results=[]
    for ti,task in enumerate(cfg["tasks"]):
        seed=202610900+ti;np.random.seed(seed);random.seed(seed)
        conditions=module(task).CONDITIONS
        cells=[("Bad","Other"),("Bad","Self"),("Good","Other"),("Good","Self")] if task=="matching" else [("Bad","Other"),("Good","Other")] if task=="identity" else [("Bad","Other"),("Bad","Self")]
        frames=[];truth={}
        for ci,(val,ident) in enumerate(cells):
            key=".".join({"val":val,"id":ident}[c] for c in sorted(conditions))
            # HDDM uses sorted dependency names for parameter labels.
            v=1.3+.2*ci;dc=(-.25 if ci%2==0 else .25)
            truth["v("+key+")"]=v;truth["dc("+key+")"]=dc
            for subj in range(6):
                z=.58+(.01*(subj-2.5));aa=1.3+.03*(subj-2.5);tt=.25+.005*(subj-2.5)
                for stim in (0,1):
                    pars={"a":aa,"t":tt,"v":v+dc if stim else v-dc,"z":z if stim else 1-z,"sv":0,"sz":0,"st":0}
                    df,_=hddm.generate.gen_rand_data(pars,size=120,subjs=1,subj_noise=0,seed=seed+ci*100+subj*2+stim)
                    if not stim:df["response"]=1-df["response"]
                    df["subj_idx"]=subj;df["stim"]=stim;df["val"]=val;df["id"]=ident
                    if task=="identity":df["id"]="Self" if stim else "Other"
                    if task=="valence":df["val"]="Good" if stim else "Bad"
                    frames.append(df)
        truth.update(a=1.3,t=.25,z=.58)
        data=pd.concat(frames,ignore_index=True)
        dest=out/task;dest.mkdir(exist_ok=a.resume)

        if a.resume and (dest/"traces.csv.gz").exists():
            assert json.loads((dest/"truth.json").read_text())==truth
            saved=json.loads((dest/"settings.json").read_text())
            assert saved["samples"]==a.samples and saved["burn"]==a.burn and saved["seed"]==seed
            traces=natural_traces(pd.read_csv(dest/"traces.csv.gz"))
        else:
            data.to_csv(dest/"synthetic.csv.gz",index=False);write(dest/"truth.json",truth)
            write(dest/"settings.json",dict(seed=seed,samples=a.samples,burn=a.burn,model="m4",scope="minimum recovery, not full simulation-based calibration"))
            m=module(task).build(data,"m4",dict(cfg,p_outlier=0.))
            m.find_starting_values()
            m.sample(a.samples,burn=a.burn,thin=1,db="pickle",dbname=str(dest/"posterior.db"))
            traces=natural_traces(m.get_traces());traces.to_csv(dest/"traces.csv.gz",index=False)
        if not (dest/"single_chain_diagnostics.csv").exists():single_diagnostics(traces,dest)
        rows=[]
        for name,value in truth.items():
            if name not in traces:raise ValueError("Parameter label mismatch: "+name+" vs "+str(list(traces)))
            x=traces[name];lo,hi=x.quantile([.025,.975])
            # A minimum gate detects gross coding/sign errors; wider recovery remains mandatory.
            tolerance=.1 if name in ("z","t") else .4
            rows.append(dict(task=task,parameter=name,true=value,mean=x.mean(),q025=lo,q975=hi,
                absolute_error=abs(x.mean()-value),tolerance=tolerance,
                passed=bool(np.isfinite(x).all() and abs(x.mean()-value)<tolerance)))
        pd.DataFrame(rows).to_csv(dest/"recovery.csv",index=False);all_results.extend(rows)
        print(task,rows,flush=True)
    pd.DataFrame(all_results).to_csv(out/"recovery.csv",index=False)
    write(out/"complete.json",{"passed":all(r["passed"] for r in all_results),
        "scope":"minimum gross-error parameter recovery; single-chain synthetic fit is not a convergence or identifiability proof",
        "samples":a.samples,"burn":a.burn,"seed_base":202610900})

def matching_recovery():
    """Four bias scenarios, four independent chains each; no empirical posterior reuse."""
    import sys,subprocess,time,hashlib
    import arviz as az

    out=OUTPUT/'hddm/matching/recovery'
    scenarios={'neutral':(.5,[0.,0.,0.,0.]),'z_only':(.58,[0.,0.,0.,0.]),
      'dc_only':(.5,[-.25,.25,-.15,.15]),'opposing':(.58,[-.3,-.2,-.25,-.15])}
    if len(sys.argv)==2:
        out.mkdir(parents=True,exist_ok=False)
        recovery_settings()
        pending=[(s,c) for s in scenarios for c in range(1,5)];active=[];failed=[]
        while pending or active:
            while pending and len(active)<2 and not failed:
                s,c=pending.pop(0);dest=out/s;dest.mkdir(exist_ok=True)
                so=(dest/('r%02d.stdout.log'%c)).open('x');se=(dest/('r%02d.stderr.log'%c)).open('x')
                proc=subprocess.Popen([sys.executable,'-u',str(HERE/'matching/matching_recovery.py'),'--scenario',s,'--chain',str(c)],stdout=so,stderr=se)
                active.append((proc,so,se,s,c));print('START',s,c,flush=True)
            for item in active[:]:
                proc,so,se,s,c=item
                if proc.poll() is not None:
                    so.close();se.close();active.remove(item);print('EXIT',s,c,proc.returncode,flush=True)
                    if proc.returncode:failed.append((s,c));pending=[]
            if active:time.sleep(2)
        if failed:raise RuntimeError(failed)
        summarize_recovery()
        return
    s=sys.argv[2];chain=int(sys.argv[3]);si=list(scenarios).index(s)
    seed=202612000+si*100;z0,dcs=scenarios[s];rng=np.random.default_rng(seed)
    dest=out/s/('r%02d'%chain);dest.mkdir(exist_ok=False)
    frames=[];truth={'a':1.3,'t':.25,'z':z0}
    cells=[('Bad','Other'),('Bad','Self'),('Good','Other'),('Good','Self')]
    for ci,(val,ident) in enumerate(cells):
        v=1.3+.2*ci;dc=dcs[ci];key=ident+'.'+val
        truth['v('+key+')']=v;truth['dc('+key+')']=dc
        for subject in range(12):
            offset=(subject-5.5)/5.5
            aa=1.3+.08*offset;tt=.25+.015*offset;zz=z0+.02*offset
            for stim in (0,1):
                pars=dict(a=aa,t=tt,z=zz if stim else 1-zz,v=v+dc if stim else v-dc,sv=0,sz=0,st=0)
                d,_=hddm.generate.gen_rand_data(pars,size=60,subjs=1,subj_noise=0,seed=seed+ci*1000+subject*2+stim)
                if not stim:d['response']=1-d.response
                mask=rng.random(len(d))<.05
                d.loc[mask,'rt']=rng.uniform(.0001,5.,mask.sum());d.loc[mask,'response']=rng.integers(0,2,mask.sum())
                d['subj_idx']=subject;d['val']=val;d['id']=ident;d['stim']=stim;frames.append(d)
    data=pd.concat(frames,ignore_index=True)
    write(dest/'truth.json',truth);write(dest/'settings.json',dict(seed=seed+chain,samples=6000,burn=3000,scenario=s))
    data.to_csv(dest/'synthetic.csv.gz',index=False)
    np.random.seed(seed+chain);random.seed(seed+chain)
    m=module('matching').build(data,'m4',settings());m.find_starting_values()
    m.sample(6000,burn=3000,thin=1,db='pickle',dbname=str(dest/'posterior.db'))
    tr=natural_traces(m.get_traces());tr.to_csv(dest/'traces.csv.gz',index=False)
    write(dest/'complete.json',dict(draws=len(tr)))


def numerical_libraries():
    global np,pd,plt,hddm
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import hddm

def four_chain_ppc(task,model):
    base=OUTPUT/'hddm'/task/model
    paths=[base/'screening/r01']+[base/'final'/('r%02d'%i) for i in (2,3,4)]
    fp=json.loads((paths[0]/'fingerprint.json').read_text())
    assert all(json.loads((p/'fingerprint.json').read_text())==fp for p in paths)
    assert all((p/'complete.json').exists() for p in paths)
    dest=base/'final/diagnostics'
    if not (dest/'diagnostics.csv').exists():raise FileNotFoundError('Run four-chain diagnostics before PPC')
    targets=[dest/('ppc_r%02d'%(i+1)) for i in range(4)]
    for p in targets:
        if p.exists():raise FileExistsError(p)
    for i,(source,target) in enumerate(zip(paths,targets)):
        target.mkdir()
        fitted=hddm.load(str(source/'model.hddm'))
        ppc(fitted,load_data(task),target,settings()['ppc_draws_final']//4,
            settings()['seed_base']+800000+settings()['tasks'].index(task)*100+i)

def entry(task,operation):
    parser=argparse.ArgumentParser(description=task+' '+operation)
    if operation=='run':parser.add_argument('phase',choices=['smoke','screening','final','preflight','long-refit'])
    if operation=='fit':
        parser.add_argument('--phase',choices=['smoke','screening','final','sensitivity','long_refit'],default='screening')
        parser.add_argument('--chain',type=int,default=1)
    if operation=='diagnostics':
        parser.add_argument('--review',action='store_true',help='Matching precise contrasts and predictive conventions, after PPC')
        parser.add_argument('--long-refit',action='store_true')
    if operation=='plot':
        parser.add_argument('--kind',choices=['summary','traces'],default='summary' if task=='matching' else 'traces')
        parser.add_argument('--figures-output',type=Path,help='Save new figures here without changing the result-table source')
    if operation=='sensitivity':parser.add_argument('--analyze',action='store_true')
    if operation=='recovery':
        parser.add_argument('--summarize',action='store_true')
        parser.add_argument('--scenario',choices=['neutral','z_only','dc_only','opposing'])
        parser.add_argument('--chain',type=int,choices=[1,2,3,4])
    parser.add_argument('--model',choices=settings()['models'],default='m4')
    parser.add_argument('--workers',type=int,choices=range(1,5),default=1)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--plan',action='store_true')
    args=parser.parse_args()
    global OUTPUT,FIGURES_OUTPUT
    if getattr(args,"figures_output",None):FIGURES_OUTPUT=args.figures_output.expanduser().resolve()
    if args.output:OUTPUT=args.output.expanduser().resolve()
    os.environ['HU2026_OUTPUT']=str(OUTPUT)
    if (operation in ('sensitivity','recovery') or getattr(args,'review',False)) and (task!='matching' or args.model!='m4'):
        parser.error('These checks are implemented for Matching M4 only')
    if operation=='plot' and args.kind=='summary' and (task!='matching' or args.model!='m4'):
        parser.error('Summary plotting is implemented for Matching M4; use --kind traces')
    if args.plan:
        assert (ROOT/'Data'/(task+'_input.csv')).is_file()
        print(json.dumps(dict(task=task,operation=operation,input=str(ROOT/'Data'/(task+'_input.csv')),
          output_root=str(OUTPUT),task_output=str(OUTPUT/'hddm'/task)),indent=2));return
    numerical_libraries()
    if operation=='fit':
        if args.phase=='long_refit':os.environ['HU2026_LONG_REFIT']='1'
        sample(task,args.model,args.chain,args.phase)
    elif operation=='run':
        if args.phase=='preflight':sys.argv=[sys.argv[0]];preflight_main();return
        if args.phase=='screening':
            gate=OUTPUT/'recovery/preflight/complete.json'
            if not json.loads(gate.read_text())['passed']:raise ValueError('Preflight recovery gate failed')
            for m in settings()['models']:
                if not (chain_dir(task,m,'smoke',1)/'complete.json').exists():raise FileNotFoundError('Missing smoke fit: '+m)
        if args.phase=='long-refit':
            os.environ['HU2026_LONG_REFIT']='1';schedule([(task,args.model,i,'long_refit') for i in range(1,5)],args.workers);return
        models=[args.model] if args.phase=='final' else settings()['models']
        chains=[2,3,4] if args.phase=='final' else [1]
        jobs=[(task,m,c,args.phase) for m in models for c in chains]
        for t,m,c,phase in jobs:
            if chain_dir(t,m,phase,c).exists():raise FileExistsError(chain_dir(t,m,phase,c))
            if phase=='final':
                prior=chain_dir(t,m,'screening',1)
                if not (prior/'complete.json').exists():raise FileNotFoundError(prior)
                if json.loads((prior/'fingerprint.json').read_text())!=fingerprint(t,m,settings()):raise ValueError('Fingerprint mismatch; refit all chains in a fresh output root')
        schedule(jobs,args.workers)
    elif operation=='compare':screening(task)
    elif operation=='diagnostics':
        if args.long_refit:long_refit_review(task,args.model)
        elif args.review:matching_review()
        else:final(task,args.model)
    elif operation=='ppc':four_chain_ppc(task,args.model)
    elif operation=='plot':
        if args.kind=='summary':matching_summary_plot()
        else:sys.argv=[sys.argv[0],'--task',task,'--model',args.model];trace_plot_main()
    elif operation=='sensitivity':
        if args.analyze:matching_sensitivity()
        else:schedule([(task,'m4',i,'sensitivity') for i in range(1,5)],args.workers)
    elif operation=='recovery':
        if args.summarize:summarize_recovery();return
        if bool(args.scenario)!=bool(args.chain):parser.error('scenario and chain must be supplied together')
        sys.argv=[sys.argv[0],'matching']+([args.scenario,str(args.chain)] if args.scenario else [])
        matching_recovery()

def summarize_recovery():
    recovery_settings()
    import arviz as az
    out=OUTPUT/'hddm/matching/recovery'
    scenarios=['neutral','z_only','dc_only','opposing']
    rows=[]
    for s in scenarios:
        paths=[out/s/('r%02d'%c) for c in range(1,5)]
        ts=[pd.read_csv(p/'traces.csv.gz') for p in paths]
        arrays={p:np.stack([t[p].values for t in ts]) for p in ts[0]}
        diag=az.summary(az.from_dict(posterior=arrays),round_to='none');diag.to_csv(out/s/'diagnostics.csv')
        passed=bool((diag.r_hat<1.01).all() and (diag.ess_bulk>=400).all() and (diag.ess_tail>=400).all())
        truth=json.loads((paths[0]/'truth.json').read_text())
        for p,truth_value in truth.items():
            x=arrays[p].ravel();lo,hi=np.quantile(x,[.025,.975])
            rows.append(dict(scenario=s,parameter=p,truth=truth_value,mean=x.mean(),q025=lo,q975=hi,
              error=x.mean()-truth_value,covered=bool(lo<=truth_value<=hi),diagnostics_pass=passed,
              rhat=float(diag.loc[p,'r_hat']),ess_bulk=float(diag.loc[p,'ess_bulk'])))
    pd.DataFrame(rows).to_csv(out/'recovery.csv',index=False)
    write(out/'complete.json',dict(completed_scenarios=4,chains=16,all_diagnostics_pass=all(r['diagnostics_pass'] for r in rows)))

def regenerate_data_audit(source_repository=None):
    import csv,collections,hashlib
    R=Path(source_repository) if source_repository else ROOT.parent
    def read(p):
        with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
    reports={}
    for task,n,ns in [("matching",24528,42),("identity",14493,41),("valence",14431,41)]:
        old=read(R/f"01_reproduction_2020/Data/{task}_input.csv")
        rawpath=R/("01_reproduction_2020/Data/raw/confirm/MS_rep_"+("matchingTask" if task=="matching" else "categTask")+"_raw.csv")
        raw=read(rawpath)
        subjects={r["subj_idx"] for r in old}
        pairs=[(i+2,r) for i,r in enumerate(raw) if r["Subject"] in subjects and r["ACC"] in ("0","1") and (task=="matching" or r["Task"]=={"identity":"Id","valence":"Val"}[task])]
        print(task,len(old),len(subjects),len(pairs),"taskvalues",dict(collections.Counter(r.get("Task","") for r in raw)),flush=True)
        assert len(old)==n and len(subjects)==ns and len(pairs)==n
        invalid=[]; valid=[]; mapping=collections.defaultdict(set)
        for a,(row,b) in zip(old,pairs):
            stim=int(b["Match"]=="match") if task=="matching" else int(b["Identity"]=="Self") if task=="identity" else int(b["Morality"]=="Good")
            expected={"subj_idx":b["Subject"],"val":b["Morality"],"id":b["Identity"],"stim":str(stim),"response":str(int(stim==int(b["ACC"])))}
            assert all(a[k]==v for k,v in expected.items()),(task,row,a,b)
            assert abs(float(a["rt"])-float(b["RT"]))<1e-12
            keys={"matching":("m","n"),"identity":("h","j"),"valence":("u","y")}[task]
            resp=b["Resp"].lower(); corr=b.get("CorrResp",b.get("corrResp","")).lower()
            assert corr in keys
            upper=corr if stim else next(k for k in keys if k!=corr)
            mapping[b["Subject"]].add(upper)
            if resp not in keys:
                invalid.append({"source_row":row,"response_key":resp,"subject":b["Subject"]})
                continue
            assert int(a["response"])==int(resp==upper)
            valid.append(dict(a,source_row=row,accuracy=int(b["ACC"]),response_key=resp,upper_key=upper))
        assert all(len(x)==1 for x in mapping.values())
        prepared=read(ROOT/'Data'/f'{task}_input.csv')
        assert prepared==[{k:str(v) for k,v in row.items()} for row in valid], 'Prepared input differs from reconstructed source: '+task
        reports[task]={"original_rows":n,"subjects":ns,"fit_rows":len(valid),"excluded_invalid_keys":invalid,
          "rt_seconds_min":min(float(r["rt"]) for r in valid),"rt_seconds_max":max(float(r["rt"]) for r in valid),
          "rt_le_0p2":sum(float(r["rt"])<=.2 for r in valid),"key_mapping_verified":True,
          "original_sha256":hashlib.sha256((R/f"01_reproduction_2020/Data/{task}_input.csv").read_bytes()).hexdigest(),
          "raw_sha256":hashlib.sha256(rawpath.read_bytes()).hexdigest()}
    write(OUTPUT/'data_audit.json',reports)
    return reports

def report_policy():
    return json.loads((HERE/'report_config.json').read_text())

def long_refit_review(task,model='m4'):
    import arviz as az
    base=OUTPUT/'hddm'/task/model/'long_refit'
    paths=[base/('r%02d'%i) for i in range(1,5)]
    fingerprints=[json.loads((p/'fingerprint.json').read_text()) for p in paths]
    if not all(f==fingerprints[0] for f in fingerprints):raise ValueError('Long-refit fingerprints differ')
    ts=[pd.read_csv(p/'traces.csv.gz') for p in paths]
    if not all(len(t)==30000 and list(t)==list(ts[0]) for t in ts):raise ValueError('Expected four independent 30000-draw long chains')
    arrays={name:np.stack([t[name].values for t in ts]) for name in ts[0]}
    diag=az.summary(az.from_dict(posterior=arrays),round_to='none')
    dest=base/'diagnostics';dest.mkdir(parents=True,exist_ok=True)
    diag.to_csv(dest/'diagnostics.csv')
    failed=diag[(diag.r_hat>=1.01)|(diag.ess_bulk<400)|(diag.ess_tail<400)|~np.isfinite(diag[['r_hat','ess_bulk','ess_tail']]).all(axis=1)]
    failed.to_csv(dest/'failed_parameters.csv')
    data=load_data(task);subjects=report_policy()['long_review_subjects'][task];details=[]
    for subject in subjects:
        d=data[data.subj_idx==subject]
        details.append(dict(subject=subject,n=len(d),rt_min=d.rt.min(),rt_q01=d.rt.quantile(.01),rt_median=d.rt.median(),
          fast_rt_le_02=int((d.rt<=.2).sum()),accuracy=d.accuracy.mean(),choice=d.response.mean()))
    write(dest/'review.json',{task:dict(max_rhat=diag.r_hat.max(),min_bulk_ess=diag.ess_bulk.min(),min_tail_ess=diag.ess_tail.min(),
      failed_parameters=len(failed),group_failed=[x for x in failed.index if '_subj' not in x],subject_data=details)})
    cols=[c for c in ts[0] if (c in ['a','t'] or c.startswith('v(') or any(c==prefix+str(sub) for sub in subjects for prefix in ['a_subj.','t_subj.']))]
    fig,axes=plt.subplots(len(cols),1,figsize=(13,max(4,2*len(cols))),squeeze=False)
    for ax,col in zip(axes[:,0],cols):
        for i,t in enumerate(ts):ax.plot(t[col].values[::10],alpha=.6,label='r%02d'%(i+1))
        ax.set_ylabel(col)
    axes[0,0].legend();fig.suptitle(task+' M4: independent long-chain mixing diagnostics')
    fig.tight_layout();out=OUTPUT/'hddm'/task/'figures';out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/(task+'_m4_long_refit_traces.png'),dpi=160);plt.close(fig)

def regenerate_reports():
    import arviz as az
    policy=report_policy();selection={}
    for task,decision in policy['selection'].items():
        base=OUTPUT/'hddm'/task;table=pd.read_csv(base/'model_comparison.csv',keep_default_na=False)
        model=decision['model'];row=table[table.model==model].iloc[0].to_dict()
        traces=natural_traces(pd.read_csv(base/model/'screening/r01/traces.csv.gz'))
        group=[c for c in traces if '_subj' not in c and '_std' not in c and 'trans' not in c]
        ess=[float(az.ess(traces[c].values[None,:],method='bulk')) for c in group]
        shifts=[float((traces[c].iloc[-len(traces)//4:].mean()-traces[c].iloc[:len(traces)//4].mean())/traces[c].std()) for c in group if traces[c].std()>0]
        selection[task]=dict(decision,screening_metrics={k:v for k,v in row.items() if k in ['task','model','chains','status','rhat','dic','ess_bulk_min','ess_tail_min','choice_rmse','rt_median_rmse_seconds']},
          group_ess_bulk_min=min(ess),max_group_quarter_shift_sd=max(map(abs,shifts)),
          limitation='Single-chain screening cannot establish convergence; model selection uncertainty is not propagated.')
        table['provisional_selected']=table.model.eq(model)
        table['selection_status']=np.where(table.model.eq(model),'provisional_selected','not_selected_for_four_chain')
        table['screening_scope']='single_chain_descriptive_not_decisive'
        if task=='matching':
            passed=json.loads((base/'m4/final/diagnostics/matching_review.json').read_text())['diagnostics_pass']
            status='main_diagnostics_pass_with_validation_limits' if passed else 'main_nonconverged'
        else:
            review=json.loads((base/'m4/long_refit/diagnostics/review.json').read_text())[task]
            status='long_refit_nonconverged' if review['failed_parameters'] else 'long_refit_diagnostics_pass'
        table['followup_status']=np.where(table.model.eq(model),status,'screening_only')
        table.to_csv(base/'model_comparison.csv',index=False)
        for m,part in table.groupby('model'):part.to_csv(base/m/'screening_summary.csv',index=False)
    write(OUTPUT/'hddm/selection.json',selection)
    rows=pd.read_csv(OUTPUT/'hddm/matching/m4/final/diagnostics/matching_parameter_contrasts.csv').to_dict('records')
    sensitivity=json.loads((OUTPUT/'hddm/matching/sensitivity/m4/diagnostics/review.json').read_text())
    recovery=pd.read_csv(OUTPUT/'hddm/matching/recovery/recovery.csv')
    flags=list(policy['predictive_limitations'])+['provisional_model_selection']
    if not sensitivity['diagnostics_pass']:flags.append('sensitivity_diagnostics_failed')
    if not recovery.diagnostics_pass.all():flags.append('recovery_diagnostics_failed')
    for row in rows:
        row['result_status']='exploratory_conditional_estimate' if row['diagnostics_pass'] else 'nonconverged_not_for_inference'
        row['validation_limits']='; '.join(flags)
    for task in ['identity','valence']:
        review=json.loads((OUTPUT/'hddm'/task/'m4/long_refit/diagnostics/review.json').read_text())[task]
        rows.append(dict(task=task,model='m4',parameter='NOT_REPORTED',diagnostics_pass=review['failed_parameters']==0,
          interpretation='No posterior inference in this report; diagnostic evidence retained',
          result_status='stopped_nonconverged' if review['failed_parameters'] else 'not_interpreted',validation_limits=''))
    pd.DataFrame(rows).to_csv(OUTPUT/'results_summary.csv',index=False)

def verify_output_manifest():
    manifest=json.loads((HERE/'output_manifest.json').read_text())
    missing=[row['path'] for row in manifest if not (OUTPUT/row['path']).is_file()]
    if missing:raise FileNotFoundError('Missing public outputs: '+', '.join(missing))
    for row in manifest:
        p=OUTPUT/row['path']
        if p.stat().st_size==0:raise ValueError('Empty output: '+str(p))
        if p.suffix=='.json':json.loads(p.read_text())
        if p.suffix=='.csv':pd.read_csv(p)
    print('Verified all %d public Output paths'%len(manifest),flush=True)

def recovery_settings():
    import hashlib
    out=OUTPUT/'hddm/matching/recovery'
    scenarios={'neutral':(.5,[0.,0.,0.,0.]),'z_only':(.58,[0.,0.,0.,0.]),'dc_only':(.5,[-.25,.25,-.15,.15]),'opposing':(.58,[-.3,-.2,-.25,-.15])}
    write(out/'settings.json',dict(samples=6000,burn=3000,thin=1,subjects=12,trials_per_stim_cell=60,
            scenarios=scenarios,workers=2,seed_base=202612000,p_outlier=.05,
            contamination='Operational normalized signed-RT uniform(-5,5), .05 mixture; matches default constant floor within support, not globally exact installed likelihood.',
            scope='Four fixed scenarios, one synthetic dataset each and four chains; not SBC, model recovery or population-wide coverage estimate',
            script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
