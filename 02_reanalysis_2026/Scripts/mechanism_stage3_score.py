"""Checkpointed conditional event LOO; no historical outputs are changed."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import argparse,json,time,warnings,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mechanism_continue as m

def run(fit,out):
    out.mkdir(parents=True,exist_ok=True);start=time.time()
    s=np.load(fit);assert bool(s['diagnostics_pass']);model=str(s['model']);draws=s['draws']
    diag=m.recovery_diagnostics(draws,model);assert m.recovery_gate(diag)
    diag.to_csv(out/'verified_all_parameter_diagnostics.csv')
    frame,audit=m.prepare_real_matching();m.write_json(out/'data_audit.json',audit)
    # 800 equally spaced iterations per chain; retain chain identity for ESS.
    ids=np.linspace(0,draws.shape[1]-1,800,dtype=int);selected=draws[:,ids];labels=m.recovery_labels(model)
    ess=m.az.ess(m.az.from_dict(posterior={n:selected[:,:,j] for j,n in enumerate(labels)}),method='bulk')
    reff=min(1.,min(float(ess[n]) for n in labels)/3200)
    spec=dict(model=model,fit_sha256=hashlib.sha256(fit.read_bytes()).hexdigest(),indices=ids.tolist(),chains=4,reff=reff,
        holdout_unit='one event from an already observed participant',uncertainty_unit='42 participant clusters',
        target='conditional event density p(y_i|y_-i,theta) with shared deadline integrated by joint likelihood ratio',
        limitation='Deadline proxy derived from response records and held fixed: conditional comparison, not leakage-free prospective CV.',
        threshold=min(.7,1-1/np.log10(3200)))
    manifest=out/'specification.json'
    if manifest.exists():assert json.loads(manifest.read_text())==spec
    else:m.write_json(manifest,spec)
    path=out/'pointwise_conditional_loglik.npy';shape=(4,800,len(frame))
    ll=np.lib.format.open_memmap(path,mode='r+' if path.exists() else 'w+',dtype='float64',shape=shape)
    progress=out/'loglik_progress.json';completed=json.loads(progress.read_text())['completed'] if progress.exists() else 0
    context=m.recovery_loo_context(frame)
    for flat in range(completed,3200):
        ch,it=divmod(flat,800);theta,lam,eps,eta=m.real_fit_parameters(selected[ch,it],model)
        ll[ch,it]=m.recovery_conditional_loglik(context,theta,lam,eps,eta)
        if (flat+1)%25==0:
            ll.flush();m.write_json(progress,dict(status='running',completed=flat+1,total=3200,elapsed_seconds=time.time()-start))
        if (flat+1)%200==0:print(model,'loglik',flat+1,'/3200',flush=True)
    ll.flush();frame[['Subject','source_row','condition','stim','relation','choice','rt_new']].to_csv(out/'observation_index.csv',index_label='observation')
    rows=[];warn=[]
    for lo in range(0,len(frame),128):
        hi=min(lo+128,len(frame));block=np.asarray(ll[:,:,lo:hi]).reshape(3200,-1).T.copy()
        with warnings.catch_warnings(record=True) as ww:
            warnings.simplefilter('always');lw,k=m.az.psislw(-block,reff=reff)
        warn.extend(str(w.message) for w in ww)
        score=m.logsumexp(block+lw,axis=1);neff=1/np.exp(2*lw).sum(axis=1)
        for j in range(hi-lo):rows.append(dict(observation=lo+j,Subject=int(frame.Subject.iloc[lo+j]),source_row=int(frame.source_row.iloc[lo+j]),elpd_loo=score[j],pareto_k=k[j],psis_ess=neff[j],pareto_threshold=spec['threshold'],method='PSIS'))
    tab=pd.DataFrame(rows);tab.to_csv(out/'loo_pointwise.csv',index=False)
    bad=tab.loc[~np.isfinite(tab.pareto_k)|tab.pareto_k.gt(spec['threshold'])]
    bad.to_csv(out/'requires_exact_refit.csv',index=False)
    fig,ax=plt.subplots(figsize=(12,4));ax.scatter(tab.observation,tab.pareto_k,s=3);ax.axhline(spec['threshold'],color='red',ls='--');ax.set(xlabel='Observed event index',ylabel='Pareto k',title=model+' conditional event LOO diagnostics');fig.savefig(out/'pareto_k.png',dpi=180,bbox_inches='tight');plt.close(fig)
    m.write_json(out/'status.json',dict(status='requires_exact_refits' if len(bad) else 'conditional_psis_pass',model=model,max_pareto_k=float(tab.pareto_k.max()),bad_count=len(bad),bad_indices=bad.observation.tolist(),elpd_unrepaired=float(tab.elpd_loo.sum()),elapsed_seconds=time.time()-start,warnings=sorted(set(warn)),ranking_permitted=False,reason='Four-model convergence and PSIS gates plus conditional interpretation required'))
    print((out/'status.json').read_text(),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--fit',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();run(a.fit,a.out)
