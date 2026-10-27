"""Reproducible read-only audit and figures of the initial recovery batch.
Does not load posterior arrays or overwrite historical analysis.
"""
from pathlib import Path
import json, itertools, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'Output/hddm/matching/mechanism_recovery_20261025'
OUT=ROOT/'Output/hddm/matching/continuation_20261026/recovery_summary'
OUT.mkdir(parents=True,exist_ok=True)
MODELS=['H0','H1','H2','H3']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220})
def save(fig,name):
    fig.savefig(OUT/(name+'.png'),bbox_inches='tight');fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
rec=pd.read_csv(SRC/'parameter_recovery.csv')
assert np.allclose(rec.error,rec['mean']-rec.truth)
assert np.array_equal(rec.covered,(rec.lower<=rec.truth)&(rec.truth<=rec.upper))
rec.to_csv(OUT/'parameter_recovery.csv',index=False)
summ=rec.groupby(['generator','parameter'],sort=False).agg(n_independent_datasets=('dataset','size'),bias=('error','mean'),rmse=('error',lambda x:np.sqrt(np.mean(x*x))),coverage_count=('covered','sum'),coverage_fraction=('covered','mean'),mean_interval_width=('upper',lambda x: np.nan)).reset_index()
width=rec.assign(width=rec.upper-rec.lower).groupby(['generator','parameter'],sort=False).width.mean()
summ['mean_interval_width']=[width.loc[(r.generator,r.parameter)] for r in summ.itertuples()]
summ['limitation']='Two independent datasets per generator; descriptive coverage only; fixed hyperparameters, not SBC.'
summ.to_csv(OUT/'parameter_bias_coverage.csv',index=False)
# Group each parameter separately so identity lines and intervals retain original units.
params=rec.parameter.drop_duplicates().tolist()
fig,axs=plt.subplots(7,4,figsize=(14,20),constrained_layout=True)
colors=dict(zip(MODELS,['#0072B2','#E69F00','#009E73','#CC79A7']))
for ax,param in zip(axs.flat,params):
    q=rec[rec.parameter.eq(param)]
    ax.xaxis.set_major_locator(MaxNLocator(3));ax.yaxis.set_major_locator(MaxNLocator(4))
    for r in q.itertuples():
        ax.errorbar(r.truth,r.mean,yerr=[[r.mean-r.lower],[r.upper-r.mean]],fmt='o',ms=4,c=colors[r.generator],alpha=.7)
    lo=min(q.truth.min(),q.lower.min());hi=max(q.truth.max(),q.upper.max());pad=max((hi-lo)*.08,.001)
    ax.plot([lo-pad,hi+pad],[lo-pad,hi+pad],'k--',lw=.8);ax.set(xlim=(lo-pad,hi+pad),ylim=(lo-pad,hi+pad),title=param,xlabel='Generating value',ylabel='Posterior mean, 95% ETI')
for ax in list(axs.flat)[len(params):]:ax.set_visible(False)
fig.suptitle('Initial parameter recovery: 2 independent datasets per generator\nFixed hyperparameters; descriptive recovery, not SBC',fontsize=14)
axs.flat[-1].set_visible(True);axs.flat[-1].axis('off');axs.flat[-1].legend(handles=[Line2D([0],[0],color=colors[m],marker='o',label=m) for m in MODELS],loc='center',title='Generating model')
save(fig,'parameter_recovery_all')
fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
for ax,param in zip(axs,['mean_beta_v','mean_beta_dc']):
    q=rec[rec.parameter.eq(param)].reset_index(drop=True)
    ax.errorbar(q['mean'],np.arange(len(q)),xerr=[q['mean']-q.lower,q.upper-q['mean']],fmt='o',color='#0072B2',label='Posterior mean and 95% ETI')
    ax.scatter(q.truth,np.arange(len(q)),marker='x',c='black',label='Generating value')
    ax.set(yticks=np.arange(len(q)),yticklabels=q.dataset,xlabel=param);ax.axvline(0,color='grey',lw=.7);ax.legend(fontsize=8)
fig.suptitle('Recovery of free target effects (fixed-zero effects are not estimated)')
save(fig,'target_parameter_recovery')
diags=[];paretos=[];pairs=[];scores=[];comparisons=[]
mat=pd.DataFrame(0,index=MODELS,columns=MODELS+['PSIS unreliable'])
for dataset in sorted(rec.dataset.unique()):
    tabs={m:pd.read_csv(SRC/(dataset+'_'+m+'_loo.csv')) for m in MODELS}
    assert all(len(t)==25200 for t in tabs.values())
    assert all(np.array_equal(t.Subject,tabs['H0'].Subject) for t in tabs.values())
    good=all(np.isfinite(t[['elpd_loo','pareto_k']]).all().all() and (t.pareto_k<=t.pareto_threshold).all() for t in tabs.values())
    totals={m:t.elpd_loo.sum() for m,t in tabs.items()};order=sorted(totals,key=totals.get,reverse=True)
    for m,t in tabs.items():
        f=SRC/(dataset+'_'+m+'_diagnostics.csv');dg=pd.read_csv(f,index_col=0)
        diags.append(dict(dataset=dataset,model=m,n_parameters=len(dg),max_rhat=dg.r_hat.max(),min_bulk_ess=dg.ess_bulk.min(),min_tail_ess=dg.ess_tail.min(),n_rhat_failed=int((dg.r_hat>=1.01).sum()),n_bulk_failed=int((dg.ess_bulk<400).sum()),n_tail_failed=int((dg.ess_tail<400).sum()),finite=bool(np.isfinite(dg[['r_hat','ess_bulk','ess_tail']]).all().all())))
        paretos.append(dict(dataset=dataset,model=m,max_pareto_k=t.pareto_k.max(),threshold=t.pareto_threshold.min(),n_failed=int((t.pareto_k>t.pareto_threshold).sum())))
        scores.append(dict(dataset=dataset,model=m,elpd=totals[m],comparison_reliable=good))
        for ix,r in t[t.pareto_k>t.pareto_threshold].iterrows():
            pd.DataFrame([dict(dataset=dataset,model=m,zero_based_trial_index=int(ix),Subject=r.Subject,pareto_k=r.pareto_k,threshold=r.pareto_threshold)]).to_csv(OUT/'problematic_observation.csv',index=False)
    for a,b in itertools.combinations(MODELS,2):
        delta=tabs[a].elpd_loo-tabs[b].elpd_loo;cluster=delta.groupby(tabs[a].Subject).sum()
        pairs.append(dict(dataset=dataset,generator=dataset[:2],model_a=a,model_b=b,delta_elpd_a_minus_b=delta.sum(),subject_cluster_se=np.sqrt(len(cluster)*cluster.var(ddof=1)),n_subjects=len(cluster),comparison_reliable=good))
    top,second=order[:2];delta=tabs[top].elpd_loo-tabs[second].elpd_loo;cl=delta.groupby(tabs[top].Subject).sum();se=np.sqrt(len(cl)*cl.var(ddof=1))
    comparisons.append(dict(dataset=dataset,generator=dataset[:2],numerical_maximum=top,runner_up=second,delta_elpd=delta.sum(),subject_cluster_se=se,status=('unreliable' if not good else ('gap_exceeds_2SE' if delta.sum()>2*se else 'ambiguous'))))
    mat.loc[dataset[:2],top if good else 'PSIS unreliable']+=1
pd.DataFrame(diags).to_csv(OUT/'mcmc_diagnostics.csv',index=False)
pd.DataFrame(paretos).to_csv(OUT/'pareto_diagnostics.csv',index=False)
pd.DataFrame(pairs).to_csv(OUT/'all_pairwise_elpd.csv',index=False)
pd.DataFrame(scores).to_csv(OUT/'model_elpd.csv',index=False)
pd.DataFrame(comparisons).to_csv(OUT/'model_comparisons.csv',index=False)
mat.to_csv(OUT/'model_recovery_counts.csv',index_label='generator')
fig,ax=plt.subplots(figsize=(8,4),constrained_layout=True);im=ax.imshow(mat.to_numpy(),cmap='Blues',vmin=0,vmax=2)
for i in range(4):
    for j in range(5):ax.text(j,i,str(mat.iloc[i,j]),ha='center',va='center',color='black')
ax.set(xticks=range(5),xticklabels=mat.columns,yticks=range(4),yticklabels=MODELS,xlabel='Numerical ELPD maximum / unavailable comparison',ylabel='Generating model',title='Initial recovery counts: only 2 datasets per generator\nNumerical maxima include ambiguous comparisons, not decisive classification')
save(fig,'model_recovery_confusion')
fig,axs=plt.subplots(2,4,figsize=(15,7),constrained_layout=True)
for ax,dataset in zip(axs.flat,sorted(rec.dataset.unique())):
    tab=pd.DataFrame(pairs);tab=tab[tab.dataset.eq(dataset)];v=tab.delta_elpd_a_minus_b;s=tab.subject_cluster_se
    reliable=bool(tab.comparison_reliable.all());ax.errorbar(v,np.arange(6),xerr=s,fmt='o',color='#0072B2' if reliable else '#999999')
    ax.axvline(0,c='black',lw=.8);ax.set(yticks=range(6),yticklabels=tab.model_a+' - '+tab.model_b,title=dataset+('' if reliable else ' (PSIS unreliable)'),xlabel='Delta ELPD +/- 1 subject-cluster SE')
save(fig,'model_delta_elpd')
fig,axs=plt.subplots(1,3,figsize=(13,4),constrained_layout=True);dd=pd.DataFrame(diags)
for ax,col,limit in zip(axs,['max_rhat','min_bulk_ess','min_tail_ess'],[1.01,400,400]):
    ax.plot(np.arange(32),dd[col],'o',ms=4);ax.axhline(limit,c='red',ls='--');ax.set(xlabel='Fit index (dataset, model order)',ylabel=col)
save(fig,'mcmc_diagnostics')
fig,ax=plt.subplots(figsize=(10,4),constrained_layout=True);pp=pd.DataFrame(paretos)
ax.scatter(np.arange(32),pp.max_pareto_k);ax.plot(np.arange(32),pp.threshold,'r--',label='Sample-size adjusted threshold');ax.set(xlabel='Fit index (dataset, model order)',ylabel='Maximum Pareto k',title='Initial PSIS diagnostics before exact-LOO correction');ax.legend()
save(fig,'pareto_diagnostics')
# File timestamps approximate sequential batch elapsed time, including scoring overhead.
files=sorted(SRC.glob('*_fit.npz'),key=lambda p:p.stat().st_mtime)
gaps=np.diff([p.stat().st_mtime for p in files])/60
meta=dict(scope='simulation only; initial batch before exact-LOO correction',replicates_per_generator=2,fit_count=32,subjects=42,trials_per_subject=600,max_rhat=float(dd.max_rhat.max()),min_bulk_ess=float(dd.min_bulk_ess.min()),min_tail_ess=float(dd.min_tail_ess.min()),median_between_fit_completion_minutes=float(np.median(gaps)),remaining_fits_to_20_replicates_each=288,estimated_serial_hours_to_20_each=float(288*np.median(gaps)/60),timing_caveat='Approximate historical sequential completion spacing includes scoring and extensions; not a benchmark.',interpretation='dc depends on zero-sum log relation-difficulty reference; predictive maximum is not decisive mechanism recovery.',sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [SRC/'parameter_recovery.csv',SRC/'execution_status.json']})
(OUT/'summary.json').write_text(json.dumps(meta,indent=2))
(OUT/'README.md').write_text('''# Initial recovery summary (simulation only)\n\nAll figures use saved completed results. No simulation or MCMC was rerun. Historical source files are unchanged.\n\nParameter intervals are equal-tailed 95% intervals from the existing summary. Coverage is descriptive across only two independently generated datasets per generating model; these fixed hyperparameters do not constitute SBC. Fixed-zero effects are omitted from estimated-parameter figures. No individual-subject coverage is claimed.\n\nModel comparisons use conditional trial-level predictive scores under the shared-deadline likelihood, with uncertainty calculated from 42 subject-level sums of paired score differences. The numerical-maximum confusion counts include ambiguous comparisons and must not be called precise recovery success rates. The H3_r02 row is withheld from classification pending exact LOO correction of H1 trial 8158. These outputs retain the original PSIS scores for provenance; consult the continuation's corrected outputs when available.\n\nThe dc interpretation is conditional on the zero-sum log relation-difficulty reference. These simulated results are not evidence about which mechanism explains the real data.\n''',encoding='utf-8')
print(json.dumps(meta,indent=2))

