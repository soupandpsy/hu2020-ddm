"""Aggregate completed conditional CV outputs without sampling or overwriting history.
Run: python Scripts/mechanism_compare_outputs.py
Missing simulation correction or real four-model outputs are reported and skipped.
"""
from pathlib import Path
import hashlib,itertools,json,time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'Output/hddm/matching/continuation_20261026'
OLD=ROOT/'Output/hddm/matching/mechanism_recovery_20261025'
MODELS=['H0','H1','H2','H3']
TARGET='Same-subject held-out event conditional on response-derived deadline proxy; not unconditional prospective or leakage-free CV.'

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def diagnostic_ok(path):
    if not path.exists():return False
    d=pd.read_csv(path,index_col=0)[['r_hat','ess_bulk','ess_tail']]
    return bool(np.isfinite(d.to_numpy()).all() and (d.r_hat<1.01).all() and (d.ess_bulk>=400).all() and (d.ess_tail>=400).all())
def checked_table(path):
    t=pd.read_csv(path)
    assert len(t)==25200 and t.Subject.nunique()==42 and t.groupby('Subject').size().eq(600).all(),path
    assert np.isfinite(t[['elpd_loo','pareto_k','pareto_threshold']].to_numpy()).all(),path
    t['method']='PSIS'
    return t

def resolve_output(name,inputs):
    sources={str(p.relative_to(ROOT)):sha(p) for p in sorted(set(inputs))}
    payload=json.dumps(sources,sort_keys=True).encode();fingerprint=hashlib.sha256(payload).hexdigest()
    dest=BASE/name
    if dest.exists():
        old=dest/'provenance.json'
        if old.exists() and json.loads(old.read_text()).get('fingerprint')==fingerprint:
            print(name+': identical completed aggregation exists; skipped');return None
        dest=BASE/(name+'_'+time.strftime('%Y%m%d_%H%M%S'))
    dest.mkdir(parents=True,exist_ok=False)
    (dest/'provenance.json').write_text(json.dumps(dict(fingerprint=fingerprint,sources=sources),indent=2))
    return dest

def aggregate(datasets,dest,kind):
    pairs=[];scores=[];comparisons=[];diag=[];mat=pd.DataFrame(0,index=MODELS,columns=MODELS+['Unreliable'])
    for name,entry in datasets.items():
        tabs=entry['tables'];fit_ok=entry['fit_ok'];tables_ok={}
        for m,t in tabs.items():
            assert np.array_equal(t.Subject,tabs['H0'].Subject)
            unresolved=(t.pareto_k>t.pareto_threshold)&t.method.ne('exact_refit')
            tables_ok[m]=bool(not unresolved.any() and fit_ok[m])
            diag.append(dict(dataset=name,model=m,mcmc_pass=fit_ok[m],n_unresolved_psis=int(unresolved.sum()),n_exact_refits=int(t.method.eq('exact_refit').sum()),max_original_k=float(t.pareto_k.max()),max_exact_log_mcse=entry['mcse'].get(m,0.)))
        reliable=all(tables_ok.values());total={m:float(t.elpd_loo.sum()) for m,t in tabs.items()}
        for m in MODELS:scores.append(dict(dataset=name,model=m,elpd=total[m],model_predictive_diagnostics_pass=tables_ok[m],four_model_comparison_reliable=reliable))
        for a,b in itertools.combinations(MODELS,2):
            delta=tabs[a].elpd_loo-tabs[b].elpd_loo;cl=delta.groupby(tabs[a].Subject).sum()
            pairs.append(dict(dataset=name,model_a=a,model_b=b,delta_elpd_a_minus_b=float(delta.sum()),subject_cluster_se=float(np.sqrt(len(cl)*cl.var(ddof=1))),pair_diagnostics_pass=tables_ok[a] and tables_ok[b],four_model_comparison_reliable=reliable))
        order=sorted(total,key=total.get,reverse=True);win,second=order[:2]
        delta=tabs[win].elpd_loo-tabs[second].elpd_loo;cl=delta.groupby(tabs[win].Subject).sum();se=float(np.sqrt(len(cl)*cl.var(ddof=1)));gap=float(delta.sum())
        # Do not publish a winner/runner-up if even one candidate has failed diagnostics.
        row=dict(dataset=name,comparison_reliable=reliable,numerical_maximum=win if reliable else None,runner_up=second if reliable else None,delta_to_runner_up=gap if reliable else None,subject_cluster_se=se if reliable else None,status=('unreliable_no_ranking' if not reliable else ('gap_exceeds_2SE' if gap>2*se else 'ambiguous')))
        comparisons.append(row)
        if kind=='simulation':mat.loc[name[:2],win if reliable else 'Unreliable']+=1
        for m,t in tabs.items():t.to_csv(dest/(name+'_'+m+'_pointwise.csv'),index=False)
    pd.DataFrame(pairs).to_csv(dest/'all_pairwise_elpd.csv',index=False)
    pd.DataFrame(scores).to_csv(dest/'model_elpd.csv',index=False)
    pd.DataFrame(comparisons).to_csv(dest/'model_comparisons.csv',index=False)
    pd.DataFrame(diag).to_csv(dest/'predictive_diagnostics.csv',index=False)
    if kind=='simulation':mat.to_csv(dest/'model_recovery_counts.csv',index_label='generator')
    fig,axs=plt.subplots(2 if len(datasets)>1 else 1,4 if len(datasets)>1 else 1,figsize=(16,7) if len(datasets)>1 else (8,5),squeeze=False,constrained_layout=True)
    pp=pd.DataFrame(pairs)
    for ax,(name,entry) in zip(axs.flat,datasets.items()):
        t=pp[pp.dataset.eq(name)]
        for j,r in enumerate(t.itertuples()):ax.errorbar(r.delta_elpd_a_minus_b,j,xerr=r.subject_cluster_se,fmt='o',color='#0072B2' if r.pair_diagnostics_pass else '#999999')
        ax.axvline(0,c='black',lw=.7);ax.set(yticks=range(6),yticklabels=t.model_a+' - '+t.model_b,title=name+('' if t.four_model_comparison_reliable.all() else ' (no reliable ranking)'),xlabel='Delta ELPD +/- 1 subject-cluster SE')
    fig.suptitle('Conditional predictive comparison'+ ('; response-derived deadline proxy' if kind=='real' else '; exact-refit correction applied'))
    for ext in ['png','pdf']:fig.savefig(dest/('model_delta_elpd.'+ext),dpi=220,bbox_inches='tight')
    plt.close(fig)
    if kind=='simulation':
        fig,ax=plt.subplots(figsize=(8,4),constrained_layout=True);ax.imshow(mat.to_numpy(),vmin=0,vmax=2,cmap='Blues')
        for i in range(4):
            for j in range(5):ax.text(j,i,str(mat.iloc[i,j]),ha='center',va='center')
        ax.set(xticks=range(5),xticklabels=mat.columns,yticks=range(4),yticklabels=MODELS,xlabel='Numerical maximum (includes ambiguous comparisons)',ylabel='Generating model',title='Corrected recovery counts; only 2 datasets per generator')
        for ext in ['png','pdf']:fig.savefig(dest/('model_recovery_confusion.'+ext),dpi=220,bbox_inches='tight')
        plt.close(fig)
    summary=dict(kind=kind,comparisons=comparisons,predictive_target=TARGET if kind=='real' else 'Same-subject held-out event conditional on fixed exogenous simulation deadline proxy',difficulty_reference='Zero-sum log relation difficulty; dc mechanism interpretation depends on this reference.',uncertainty='SE = sqrt(42 * sample variance of subject-level paired score sums); figure bars are +/-1 SE, not 95% intervals.',ranking_rule='Numerical maximum only when all four candidates pass MCMC and corrected predictive diagnostics; gap >2 SE is descriptive, not a model selection theorem or unequivocal mechanism claim.',exact_mcse='Exact-refit MCSE is reported separately; it is not included in subject-cluster sampling SE.',recovery_limit='Two independent datasets per generator; counts include ambiguous maxima and cannot estimate precise recovery rates.' if kind=='simulation' else 'Real inference remains conditional on the observation model and proxy assumptions.')
    (dest/'summary.json').write_text(json.dumps(summary,indent=2))
    (dest/'README.md').write_text('# '+kind.title()+' conditional predictive comparisons\n\n'+summary['predictive_target']+'\n\n'+summary['difficulty_reference']+'\n\n'+summary['uncertainty']+'\n\n'+summary['ranking_rule']+'\n\n'+summary['exact_mcse']+'\n\n'+summary['recovery_limit']+'\n',encoding='utf8')
    print(str(dest))

def recovery():
    correction=BASE/'exact_loo_H3_r02_H1_8158/corrected_loo.csv';res=correction.parent/'result.json'
    if not correction.exists() or not res.exists():print('recovery_corrected: exact-refit output incomplete; skipped');return
    result=json.loads(res.read_text());assert result.get('diagnostics_pass') is True and result.get('status')=='complete'
    ix=int(result['index']);assert ix==8158 and np.isfinite(result['exact_elpd'])
    datasets={};inputs=[correction,res]
    for name in [g+'_r'+str(r).zfill(2) for g in MODELS for r in [1,2]]:
        entry=dict(tables={},fit_ok={},mcse={})
        for m in MODELS:
            f=OLD/(name+'_'+m+'_loo.csv');dg=OLD/(name+'_'+m+'_diagnostics.csv');inputs.extend([f,dg]);t=checked_table(f)
            if name=='H3_r02' and m=='H1':
                corrected=pd.read_csv(correction);assert np.array_equal(t.Subject,corrected.Subject)
                unchanged=np.arange(len(t))!=ix
                assert np.allclose(t.loc[unchanged,'elpd_loo'],corrected.loc[unchanged,'elpd_loo'],rtol=0,atol=1e-10)
                assert np.isclose(corrected.loc[ix,'elpd_loo'],result['exact_elpd'])
                t.loc[ix,'elpd_loo']=result['exact_elpd'];t.loc[ix,'method']='exact_refit';entry['mcse'][m]=float(result['log_predictive_mcse'])
            entry['tables'][m]=t;entry['fit_ok'][m]=diagnostic_ok(dg)
        datasets[name]=entry
    dest=resolve_output('recovery_corrected',inputs)
    if dest:aggregate(datasets,dest,'simulation')

def real():
    files={m:BASE/('real_'+m)/('real_'+m+'_loo.csv') for m in MODELS}
    if not all(p.exists() for p in files.values()):print('real_comparison: four model score files incomplete; skipped');return
    entry=dict(tables={},fit_ok={},mcse={});inputs=[]
    for m,f in files.items():
        t=checked_table(f);dg=f.with_name('real_'+m+'_diagnostics.csv');inputs.extend([f,dg]);entry['fit_ok'][m]=diagnostic_ok(dg);mcse=[]
        for res in sorted(f.parent.glob('reloo_*/result.json')):
            r=json.loads(res.read_text());inputs.append(res)
            if r.get('diagnostics_pass') is not True:continue
            ix=int(r['index']);assert 0<=ix<len(t) and res.parent.name=='reloo_'+str(ix)
            assert r.get('conditional_on_record_derived_proxy') is True and np.isfinite(r['exact_elpd'])
            t.loc[ix,'elpd_loo']=r['exact_elpd'];t.loc[ix,'method']='exact_refit';mcse.append(float(r['log_predictive_mcse']))
        entry['tables'][m]=t;entry['mcse'][m]=max(mcse,default=0.)
    # Missing diagnostics are an incomplete output, not a license to rank.
    if not all(p.exists() for p in inputs):print('real_comparison: diagnostics incomplete; skipped');return
    dest=resolve_output('real_comparison',inputs)
    if dest:aggregate({'Real Matching':entry},dest,'real')

if __name__=='__main__':
    recovery();real()
