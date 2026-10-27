"""Individual parameter recovery, retaining one posterior array at a time.
Subjects are nested within eight simulated datasets, NOT independent recovery replicates.
"""
from pathlib import Path
import gc,json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'Output/hddm/matching/mechanism_recovery_20261025'
OUT=ROOT/'Output/hddm/matching/continuation_20261026/recovery_summary'
OUT.mkdir(parents=True,exist_ok=True)
names=['log_a','log_t_excess','logit_z','v_intercept','v_other','v_other_contrast','beta_v','dc_intercept','dc_other','dc_other_contrast','beta_dc']
rows=[];sources={}
for model in ['H0','H1','H2','H3']:
    active=[j for j in range(11) if not (j==6 and model in ['H0','H2']) and not (j==10 and model in ['H0','H1'])]
    k=len(active)
    for rep in [1,2]:
        dataset=model+'_r'+str(rep).zfill(2);dp=SRC/(dataset+'_data.npz');fp=SRC/(dataset+'_'+model+'_fit.npz')
        with np.load(dp) as data:truth=data['truth_theta'];subjects=data['subject_ids']
        with np.load(fp) as fit:
            assert bool(fit['diagnostics_pass'])
            labels=fit['labels'].tolist();draws=fit['draws']
        assert len(subjects)==42 and truth.shape==(42,11)
        for i,sub in enumerate(subjects):
            for j in active:
                label='subject_'+str(sub)+'_'+names[j];idx=labels.index(label)
                assert idx==2*k+4+i*k+active.index(j)
                vals=draws[:,:,idx].reshape(-1);lo,hi=np.quantile(vals,[.025,.975]);mean=float(vals.mean());tr=float(truth[i,j])
                rows.append(dict(dataset=dataset,generator=model,Subject=int(sub),parameter=names[j],truth=tr,mean=mean,lower=float(lo),upper=float(hi),bias=mean-tr,covered=bool(lo<=tr<=hi),posterior_sd=float(vals.std(ddof=1))))
        del vals,draws;gc.collect()
        sources[dp.name]=hashlib.sha256(dp.read_bytes()).hexdigest()
        print('summarized '+dataset,flush=True)
t=pd.DataFrame(rows);t.to_csv(OUT/'individual_parameter_recovery.csv',index=False)
s=t.assign(width=t.upper-t.lower).groupby(['dataset','generator','parameter'],sort=False).agg(n_subjects=('Subject','size'),mean_bias=('bias','mean'),rmse=('bias',lambda x:np.sqrt(np.mean(x*x))),covered_subject_count=('covered','sum'),descriptive_subject_coverage=('covered','mean'),mean_interval_width=('width','mean')).reset_index()
s['independent_recovery_replicates_per_generator']=2
s['interpretation']='Subjects within one dataset share hyperparameters and posterior information; subject coverage is descriptive, not independent recovery replication or SBC.'
s.to_csv(OUT/'individual_dataset_summary.csv',index=False)
colors={'H1':'#E69F00','H2':'#009E73','H3':'#CC79A7'}
fig,axs=plt.subplots(2,4,figsize=(16,8),constrained_layout=True)
for row,param in enumerate(['beta_v','beta_dc']):
    datasets=t[t.parameter.eq(param)].dataset.unique()
    for ax,dataset in zip(axs[row],datasets):
        q=t[t.parameter.eq(param)&t.dataset.eq(dataset)];color=colors[dataset[:2]]
        ax.errorbar(q.truth,q['mean'],yerr=[q['mean']-q.lower,q.upper-q['mean']],fmt='o',ms=3,alpha=.6,lw=.6,color=color)
        lo=min(q.truth.min(),q.lower.min());hi=max(q.truth.max(),q.upper.max());pad=(hi-lo)*.04
        ax.plot([lo-pad,hi+pad],[lo-pad,hi+pad],'k--',lw=.8)
        ax.set(xlim=(lo-pad,hi+pad),ylim=(lo-pad,hi+pad),xlabel='Individual generating '+param,ylabel='Posterior mean and 95% ETI',title=dataset+' ('+str(len(q))+' subjects)')
fig.suptitle('Individual target-effect recovery under the generating model\n42 dependent subject estimates per dataset; only 2 independent datasets per generator')
for ext in ['png','pdf']:fig.savefig(OUT/('individual_target_recovery.'+ext),dpi=220,bbox_inches='tight')
plt.close(fig)
meta=dict(rows=len(t),datasets=8,independent_datasets_per_generator=2,subjects_per_dataset=42,source_hashes=sources,estimated_parameters='All active individual coefficients under the correct generating model; structurally fixed-zero coefficients excluded.',interval='95% equal-tailed posterior interval',limitations='Subject estimates share group hyperparameters. 336 subject instances are not 336 independent simulation replicates. No SBC claim; no pooled success-rate precision claim.',difficulty_reference='dc interpretation conditional on zero-sum log relation difficulty.')
(OUT/'individual_recovery_metadata.json').write_text(json.dumps(meta,indent=2))
(OUT/'individual_recovery_README.md').write_text('# Individual recovery\n\n'+meta['estimated_parameters']+'\n\n'+meta['interval']+'.\n\n'+meta['limitations']+'\n\n'+meta['difficulty_reference']+'\n',encoding='utf8')
print('completed '+str(len(t))+' individual parameter rows',flush=True)
