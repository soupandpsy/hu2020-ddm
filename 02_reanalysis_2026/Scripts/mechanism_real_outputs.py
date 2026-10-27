"""Verified-fit summaries and checkpointed PPC for the frozen Matching model.

Run inside hu2026-continuation:
 python Scripts/mechanism_real_outputs.py --fit PATH --output DIRECTORY
No MCMC, original-data edits, historical experiment execution, or model selection.
"""
import os
os.environ.setdefault('MPLBACKEND','Agg')
import argparse, hashlib, json, time, traceback
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mechanism_continue as m

def savefig(fig,path):
    fig.savefig(path.with_suffix('.png'),dpi=240,bbox_inches='tight')
    fig.savefig(path.with_suffix('.pdf'),bbox_inches='tight')
    plt.close(fig)

def stats(frame,cuts):
    table=m.real_behavior_statistics(frame,cuts)
    seen=frame.loc[frame.choice.ge(0)&frame.condition.isin(['Good_Self','Bad_Self'])].copy()
    seen['correct']=seen.choice==seen.stim
    rows=[]
    for (cond,stim,rel),g in frame.groupby(['condition','stim','relation'],sort=True):
        rows.append(dict(group=f'{cond}|stim={stim}|{rel}',metric='correct_all_trials_rate',value=float((g.choice.ge(0)&g.choice.eq(g.stim)).mean())))
    for stim in [0,1]:
        fast=seen.loc[seen.stim.eq(stim)&seen.rt_new.le(cuts[stim])]
        for metric,values in [
            ('fast_accuracy_good_minus_bad',fast.groupby(['Subject','condition']).correct.mean().unstack()),
            ('fast_correct_rt_good_minus_bad',fast.loc[fast.correct].groupby(['Subject','condition']).rt_new.mean().unstack())]:
            paired=values.reindex(columns=['Good_Self','Bad_Self']).dropna()
            rows.append(dict(group=f'Self|stim={stim}',metric=metric,value=float((paired.Good_Self-paired.Bad_Self).mean())))
            rows.append(dict(group=f'Self|stim={stim}',metric=metric+'_paired_subjects',value=len(paired)))
    return pd.concat([table,pd.DataFrame(rows)],ignore_index=True)

def posterior_deadlines(frame,theta,lam,eps,rng,width):
    nodes=7 if width==.001 else 21
    gx,gw=m.quadrature_nodes(nodes);result=[]
    for th,item in zip(theta,m.original_prepare(frame)):
        lower=max(item['proxy']-width,item['maxrt']);upper=item['proxy']+width
        ds=lower+(gx+1)*(upper-lower)/2;lw=np.log(gw)
        for (cond,stim,rel),rt,choice,nt in item['groups']:
            if nt:
                a,t,z,mu=m.unpack(th,cond,stim,rel,lam)
                lw+=nt*np.log(m.mix_survival(ds,a,t,z,mu,eps))
        result.append(rng.choice(ds,p=np.exp(lw-m.logsumexp(lw))))
    return np.asarray(result)

def ppc(fit,draws,model,frame,out,n,seed,width):
    """Each posterior replicate commits independently; restarting never redraws a seed."""
    folder=out/'ppc_checkpoints';folder.mkdir(exist_ok=True)
    setting=out/'ppc_settings.json';rng=np.random.default_rng(seed)
    cuts={stim:float(frame.loc[frame.choice.ge(0)&frame.stim.eq(stim)&frame.condition.isin(['Good_Self','Bad_Self']),'rt_new'].quantile(.2)) for stim in [0,1]}
    spec=dict(seed=seed,draws_per_chain=n,fit_sha256=hashlib.sha256(fit.read_bytes()).hexdigest(),
        posterior_indices=[rng.choice(draws.shape[1],n,replace=False).tolist() for _ in range(4)],
        fast_cutoffs=cuts,fast_definition='Observed pooled Self legal-response q20 separately for match/mismatch; fixed before PPC; no effect-based threshold selection',
        predictive_scope='Conditional, same-participant in-sample PPC; shared posterior deadline per participant',
        deadline_halfwidth=width,deadline_quadrature_nodes=7 if width==.001 else 21)
    if setting.exists():
        previous=json.loads(setting.read_text())
        assert previous['fit_sha256']==spec['fit_sha256'] and previous['seed']==seed and previous['draws_per_chain']==n and previous['deadline_halfwidth']==width
        spec=previous
    else:m.write_json(setting,spec)
    observed=stats(frame,cuts).assign(source='observed',chain=-1,draw=-1)
    observed.to_csv(out/'ppc_observed.csv',index=False)
    for ch,indices in enumerate(spec['posterior_indices']):
        for it,idx in enumerate(indices):
            target=folder/f'chain{ch}_draw{idx}.csv'
            if target.exists():continue
            draw_rng=np.random.default_rng(np.random.SeedSequence([seed,ch,int(idx)]))
            theta,lam,eps,eta=m.real_fit_parameters(draws[ch,idx],model)
            ds=posterior_deadlines(frame,theta,lam,eps,draw_rng,width)
            old=m.deadlines
            try:
                m.deadlines=pd.Series(ds,index=m.subjects)
                replica,_=m.simulate(frame,theta,lam,eps,eta,draw_rng,width=0.)
            finally:m.deadlines=old
            result=stats(replica,cuts).assign(source='replicated',chain=ch,draw=int(idx))
            tmp=target.with_suffix('.tmp');result.to_csv(tmp,index=False);tmp.replace(target)
            m.write_json(out/'ppc_progress.json',dict(model=model,completed=len(list(folder.glob('*.csv'))),target=4*n,last_chain=ch,last_draw=int(idx)))
            print('PPC',model,ch,it+1,'/',n,flush=True)
    tables=[observed]+[pd.read_csv(folder/f'chain{ch}_draw{idx}.csv') for ch,indices in enumerate(spec['posterior_indices']) for idx in indices]
    allrows=pd.concat(tables,ignore_index=True);allrows.to_csv(out/'ppc_all.csv',index=False)
    reps=allrows.loc[allrows.source.eq('replicated')]
    summ=reps.groupby(['group','metric']).value.agg(mean='mean',median='median',lower=lambda x:x.quantile(.025),upper=lambda x:x.quantile(.975),finite_replicates='count').reset_index()
    summ=summ.merge(observed[['group','metric','value']].rename(columns={'value':'observed'}),on=['group','metric'],validate='one_to_one')
    summ['observed_in_95_predictive_interval']=summ.observed.between(summ.lower,summ.upper)
    summ.to_csv(out/'ppc_summary.csv',index=False)
    plot_ppc(summ,out,model)

def plot_ppc(table,out,model):
    panels=[('accuracy_choice_timeout_ppc',['accuracy_given_response','correct_all_trials_rate','match_choice_given_response','timeout_rate']),
        ('rt_quantile_ppc',['correct_rt_q10','correct_rt_q50','correct_rt_q90','error_rt_q10','error_rt_q50','error_rt_q90']),
        ('fast_good_self_ppc',['fast_match_choice_good_minus_bad','fast_accuracy_good_minus_bad','fast_correct_rt_good_minus_bad'])]
    for name,metrics in panels:
        cols=2 if len(metrics)==6 else 1;rows=int(np.ceil(len(metrics)/cols))
        fig,axes=plt.subplots(rows,cols,figsize=(15 if cols==2 else 10,3.4*rows),squeeze=False,constrained_layout=True)
        for ax,metric in zip(axes.flat,metrics):
            t=table.loc[table.metric.eq(metric)].sort_values('group');x=np.arange(len(t))
            ax.errorbar(x,t['median'],yerr=[t['median']-t.lower,t.upper-t['median']],fmt='o',color='#0072B2',capsize=3,label='Posterior predictive median + 95% interval')
            ax.scatter(x,t.observed,color='#D55E00',marker='x',s=40,zorder=4,label='Observed')
            ax.set_xticks(x);ax.set_xticklabels(t.group,rotation=70,ha='right',fontsize=7)
            ax.set_title(metric.replace('_',' '));ax.grid(axis='y',alpha=.2)
            if 'minus_bad' in metric:ax.axhline(0,color='gray',ls='--',lw=.8)
            if 'rt' in metric:ax.set_ylabel('Seconds')
        axes.flat[0].legend(fontsize=8)
        fig.suptitle(f'{model}: same-participant posterior predictive checks')
        savefig(fig,out/name)

def posterior(draws,labels,model,out,fit):
    records=[]
    names=['mean_beta_v','mean_beta_dc','sd_beta_v','sd_beta_dc','contamination','illegal_rate','relation_coef1','relation_coef2']
    for name in names:
        v=draws[:,:,labels.index(name)].ravel() if name in labels else np.zeros(draws.shape[0]*draws.shape[1])
        lo,median,hi=np.quantile(v,[.025,.5,.975])
        records.append(dict(parameter=name,mean=float(v.mean()),median=float(median),lower=float(lo),upper=float(hi),fixed_zero=name not in labels))
    pd.DataFrame(records).to_csv(out/'key_parameter_summary.csv',index=False)
    fig,axes=plt.subplots(2,2,figsize=(10,6),constrained_layout=True)
    for ax,name in zip(axes.flat,names[:4]):
        if name in labels:
            j=labels.index(name)
            for ch in range(4):ax.hist(draws[ch,:,j],bins=45,density=True,histtype='step',alpha=.8,label=f'Chain {ch+1}')
        else:ax.text(.5,.5,'Fixed at zero by model',transform=ax.transAxes,ha='center')
        ax.axvline(0,color='gray',ls='--');ax.set_title(name);ax.set_ylabel('Posterior density')
    axes.flat[0].legend(fontsize=8);fig.suptitle(f'{model}: target contrasts and between-participant SDs')
    savefig(fig,out/'key_parameter_posteriors')
    diagpath=fit.with_name(fit.stem.replace('_fit','')+'_diagnostics.csv')
    diag=pd.read_csv(diagpath,index_col=0)
    assert m.recovery_gate(diag),'Failed MCMC diagnostic gate'
    diag.to_csv(out/'all_parameter_diagnostics.csv',index_label='parameter')
    stage=fit.with_name(fit.stem.replace('_fit','')+'_segment1.npz')
    if stage.exists():
        with np.load(stage) as segment:
            acceptance=[]
            for key in ['acceptance','global_acceptance','collapsed_acceptance']:
                if key in segment:
                    values=segment[key]
                    for ch in range(values.shape[0]):
                        for j,value in enumerate(values[ch]):
                            acceptance.append(dict(update=key,chain=ch,parameter_or_block_index=j,acceptance_rate=float(value),zero_acceptance=bool(value==0)))
        if acceptance:
            pd.DataFrame(acceptance).to_csv(out/'proposal_acceptance.csv',index=False)
            m.write_json(out/'proposal_acceptance_notes.json',dict(scope='First production segment only; calibration and any precision extension are not included.',
                interpretation='Acceptance rates describe proposal behavior, not NUTS divergences. Zero acceptance requires inspection; there is no universal acceptance cutoff for this mixed Metropolis/Gibbs/slice sampler.',
                zero_acceptance_blocks=sum(x['zero_acceptance'] for x in acceptance)))
    fig,axes=plt.subplots(1,3,figsize=(13,3.5),constrained_layout=True)
    for ax,col,threshold in zip(axes,['r_hat','ess_bulk','ess_tail'],[1.01,400,400]):
        ax.plot(np.arange(len(diag)),diag[col],'.',ms=2);ax.axhline(threshold,color='#D55E00',ls='--');ax.set_title(col);ax.set_xlabel('Basic / derived parameter index')
    savefig(fig,out/'all_mcmc_diagnostics')
    fig,axes=plt.subplots(2,1,figsize=(11,5),constrained_layout=True)
    for ax,name in zip(axes,['mean_beta_v','mean_beta_dc']):
        if name in labels:
            for ch in range(4):ax.plot(draws[ch,:,labels.index(name)],lw=.4,alpha=.6,label=f'Chain {ch+1}')
        ax.set_ylabel(name);ax.set_xlabel('Retained iteration')
    savefig(fig,out/'target_traceplots')

def reference_sensitivity(draws,labels,model,out):
    active=m.allowed(model);klen=len(active)
    rows=[];max_error=0.;invalid=0
    # Equal-spacing downsample only for this deterministic derived-posterior display.
    for ch in range(4):
        for it in np.linspace(0,draws.shape[1]-1,min(1000,draws.shape[1]),dtype=int):
            row=draws[ch,it];theta,lam,eps,eta=m.real_fit_parameters(row,model)
            r=np.exp(lam);vg=np.exp(theta[:,3]+.5*theta[:,6]);vb=np.exp(theta[:,3]-.5*theta[:,6]);dcd=theta[:,10]
            means=np.zeros(11);sds=np.zeros(11);means[active]=row[:klen];sds[active]=row[klen:2*klen]
            evdiff=np.exp(means[3]+.5*(sds[3]**2+.25*sds[6]**2))*(np.exp(.5*means[6])-np.exp(-.5*means[6]))
            for ref,shift in [('zero_sum_log',0.),('mean_multiplier_one',(1-r.mean())/2),('identity_only_one',(1-r[0])/2),('both_different_one',(1-r[2])/2)]:
                ok=bool(1-shift>0 and np.all(r+shift>0))
                if not ok:invalid+=1
                else:
                    rp=(r+shift)/(1-shift)
                    for v,dc in [(vg,theta[:,7]+.5*theta[:,10]),(vb,theta[:,7]-.5*theta[:,10])]:
                        vp=(1-shift)*v;dp=dc+shift*v
                        max_error=max(max_error,float(np.max(abs(vp+dp-v-dc))),float(np.max(abs(dp[:,None]-vp[:,None]*rp-(dc[:,None]-v[:,None]*r)))))
                rows.append(dict(chain=ch,draw=int(it),reference=ref,valid_positive_multipliers=ok,
                    population_mean_logv_contrast=means[6],population_mean_dc_contrast=means[10]+shift*evdiff,
                    sample_mean_dc_contrast=float(np.mean(dcd+shift*(vg-vb))),shift=float(shift)))
    assert max_error<1e-10
    table=pd.DataFrame(rows);table.to_csv(out/'reference_transformed_draws.csv',index=False)
    summaries=[]
    for ref,g in table.groupby('reference',sort=False):
        for metric in ['population_mean_logv_contrast','population_mean_dc_contrast','sample_mean_dc_contrast']:
            lo,med,hi=g[metric].quantile([.025,.5,.975]);summaries.append(dict(reference=ref,metric=metric,mean=g[metric].mean(),median=med,lower=lo,upper=hi,invalid_draws=int((~g.valid_positive_multipliers).sum())))
    summary=pd.DataFrame(summaries);summary.to_csv(out/'reference_sensitivity_summary.csv',index=False)
    fig,axes=plt.subplots(1,3,figsize=(14,4),constrained_layout=True)
    for ax,metric in zip(axes,summary.metric.unique()):
        t=summary.loc[summary.metric.eq(metric)];y=np.arange(len(t))
        ax.errorbar(t['median'],y,xerr=[t['median']-t.lower,t.upper-t['median']],fmt='o',capsize=3)
        ax.set_yticks(y);ax.set_yticklabels(t.reference,fontsize=8);ax.axvline(0,color='gray',ls='--');ax.set_title(metric.replace('_',' '),fontsize=9)
    fig.suptitle(f'{model}: equivalent drift representations; not prior-refitted models')
    savefig(fig,out/'difficulty_reference_sensitivity')
    m.write_json(out/'reference_sensitivity_notes.json',dict(max_effective_drift_error=max_error,invalid_positive_multiplier_draws=invalid,
        interpretation='Deterministic posterior transformation with pushed-forward priors; no likelihood or predictions change. Original dc interpretation depends on zero-sum log relationship difficulty.',
        population_mean='Analytic expectation over independent normal participant coefficients under each group posterior draw; transformed dc distribution is generally non-normal.',
        sample_mean='Average contrast across the 42 fitted participant parameter values; not the group coefficient.',
        invalid_policy='All draws retained and flagged. Invalid transformed positive multipliers prevent interpreting that reference as a valid positive-v parameterization.'))

def main():
    p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--output',required=True);p.add_argument('--draws-per-chain',type=int,default=50);p.add_argument('--seed',type=int,default=2026160000);p.add_argument('--width',type=float,default=.001);p.add_argument('--skip-ppc',action='store_true');args=p.parse_args()
    assert args.width>0
    fit=Path(args.fit).resolve();out=Path(args.output).resolve();out.mkdir(exist_ok=True,parents=True)
    start=time.time();state=out/'output_status.json'
    m.write_json(state,dict(status='running',fit=str(fit),pid=os.getpid()))
    try:
        saved=np.load(fit);assert bool(saved['diagnostics_pass'])
        model=str(saved['model']);draws=saved['draws'];labels=saved['labels'].tolist()
        # Explicitly check observed-data provenance before publishing real-data outputs.
        dataset=fit.parent/'real_data.npz'
        assert dataset.exists() and str(np.load(dataset)['data_kind'])=='observed','Real data provenance missing'
        assert str(saved['data_sha256'])==hashlib.sha256(dataset.read_bytes()).hexdigest()
        frame,audit=m.prepare_real_matching();m.write_json(out/'observed_data_audit.json',audit)
        posterior(draws,labels,model,out,fit);reference_sensitivity(draws,labels,model,out)
        if not args.skip_ppc:ppc(fit,draws,model,frame,out,args.draws_per_chain,args.seed,args.width)
        m.write_json(state,dict(status='complete',model=model,fit=str(fit),elapsed_seconds=time.time()-start,ppc_completed=not args.skip_ppc,exit_code=0))
    except Exception:
        m.write_json(state,dict(status='failed',error=traceback.format_exc(),elapsed_seconds=time.time()-start,exit_code=1));raise

if __name__=='__main__':main()
