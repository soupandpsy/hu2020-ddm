"""Checkpointed extended real-data PPC; historical sources remain read-only."""
import os
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[name]='1'
import json, time, hashlib, traceback
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mechanism_continue as m
import mechanism_real_outputs as old

BASE=m.ROOT/'Output/hddm/matching/continuation_stage3_20261026'
OUT=BASE/'H2_ppc'
FIT=m.ROOT/'Output/hddm/matching/diagnostic_H2H3_20261026_1425/H2_repair_v2_precision/formal_fit.npz'
original_stats=old.stats
EDGES={}
GRID=np.arange(.2,1.80001,.025)

def stats(frame,cuts):
    base=original_stats(frame,cuts); rows=[]
    def add(group,metric,value):rows.append(dict(group=group,metric=metric,value=float(value)))
    seen=frame.loc[frame.choice.ge(0)].copy();seen['correct']=seen.choice.eq(seen.stim)
    for (cond,stim),allg in frame.groupby(['condition','stim']):
        g=seen.loc[seen.condition.eq(cond)&seen.stim.eq(stim)];group=f'{cond}|stim={stim}|pooled'
        add(group,'timeout_rate',allg.choice.eq(-1).mean());add(group,'illegal_rate',allg.choice.eq(-2).mean())
        add(group,'accuracy_given_response',g.correct.mean());add(group,'match_choice_given_response',g.choice.mean())
        add(group,'correct_all_trials_rate',(allg.choice.ge(0)&allg.choice.eq(allg.stim)).mean())
        for label,correct in [('correct',True),('error',False)]:
            rt=g.loc[g.correct.eq(correct),'rt_new']
            add(group,label+'_mean_rt',rt.mean())
            for q in [.1,.5,.9]:add(group,f'{label}_rt_q{int(q*100)}',rt.quantile(q))
            for x in GRID:add(group,f'{label}_cdf_{x:.3f}',np.mean(rt<=x) if len(rt) else np.nan)
    # Signal=Match, response=Match; loglinear correction (count+.5)/(n+1).
    dp=[]
    for (sub,cond),g in seen.groupby(['Subject','condition']):
        hit=g.loc[g.stim.eq(1),'choice'];fa=g.loc[g.stim.eq(0),'choice']
        if len(hit) and len(fa):dp.append(dict(Subject=sub,condition=cond,value=norm.ppf((hit.sum()+.5)/(len(hit)+1))-norm.ppf((fa.sum()+.5)/(len(fa)+1))))
    dt=pd.DataFrame(dp)
    for cond,g in dt.groupby('condition'):add(cond+'|dprime','dprime_subject_mean',g.value.mean())
    dt=dt.pivot(index='Subject',columns='condition',values='value')
    add('Self|dprime','dprime_good_minus_bad',(dt.Good_Self-dt.Bad_Self).mean())
    # Fixed observed pooled-Self quintile boundaries. Paired subject differences.
    ss=seen.loc[seen.condition.isin(['Good_Self','Bad_Self'])]
    for stim in [0,1]:
        g=ss.loc[ss.stim.eq(stim)].copy();g['bin']=np.digitize(g.rt_new,EDGES[stim],right=True)
        for b in range(5):
            gb=g.loc[g.bin.eq(b)];group=f'Self|stim={stim}|bin={b+1}'
            for metric,col in [('binned_accuracy_good_minus_bad','correct'),('binned_match_choice_good_minus_bad','choice')]:
                vals=gb.groupby(['Subject','condition'])[col].mean().unstack().reindex(columns=['Good_Self','Bad_Self']).dropna()
                add(group,metric,(vals.Good_Self-vals.Bad_Self).mean());add(group,metric+'_paired_subjects',len(vals))
            for cond in ['Good_Self','Bad_Self']:
                add(group,cond+'_bin_fraction',len(gb.loc[gb.condition.eq(cond)])/max(1,len(g.loc[g.condition.eq(cond)])))
    return pd.concat([base,pd.DataFrame(rows)],ignore_index=True)

def save(fig,name):old.savefig(fig,OUT/name)

def interval(ax,t,ylabel):
    x=np.arange(len(t));ax.errorbar(x,t['median'],yerr=[t['median']-t.lower,t.upper-t['median']],fmt='o',color='#0072B2',capsize=3,label='Predictive median / 95% interval')
    ax.scatter(x,t.observed,marker='x',color='#D55E00',s=45,label='Observed',zorder=5)
    ax.set_xticks(x);ax.set_xticklabels(t.group.str.replace('|pooled','',regex=False),rotation=40,ha='right');ax.set_ylabel(ylabel);ax.grid(axis='y',alpha=.2)

def plots(table,out,model):
    pooled=table.loc[table.group.str.endswith('|pooled')]
    for name,metrics in [('accuracy_timeout_ppc',['accuracy_given_response','match_choice_given_response','timeout_rate','illegal_rate']),('rt_means_quantiles_ppc',['correct_mean_rt','error_mean_rt','correct_rt_q10','correct_rt_q50','correct_rt_q90','error_rt_q90'])]:
        fig,axes=plt.subplots(len(metrics)//2,2,figsize=(14,4*len(metrics)//2),constrained_layout=True)
        for ax,metric in zip(axes.flat,metrics):
            t=pooled.loc[pooled.metric.eq(metric)].sort_values('group');interval(ax,t,'Seconds' if 'rt' in metric else 'Proportion');ax.set_title(metric.replace('_',' '))
        axes.flat[0].legend(fontsize=9);fig.suptitle('H2: same-participant observed versus posterior predictive');save(fig,name)
    fig,axes=plt.subplots(4,4,figsize=(16,13),sharex=True,sharey=True,constrained_layout=True)
    for row,cond in enumerate(m.conditions):
        for col,(stim,label) in enumerate([(1,'correct'),(1,'error'),(0,'correct'),(0,'error')]):
            ax=axes[row,col];t=table.loc[table.group.eq(f'{cond}|stim={stim}|pooled')&table.metric.str.startswith(label+'_cdf_')].copy();t['x']=t.metric.str.split('_').str[-1].astype(float);t=t.sort_values('x')
            ax.fill_between(t.x.to_numpy(),t.lower.to_numpy(),t.upper.to_numpy(),color='#0072B2',alpha=.2);ax.plot(t.x,t['median'],color='#0072B2');ax.plot(t.x,t.observed,color='#D55E00',ls='--');ax.set_title(f'{cond}: '+('Match' if stim else 'Mismatch')+' '+label,fontsize=10);ax.set_xlabel('RT (s)');ax.set_ylabel('Cumulative proportion');ax.grid(alpha=.15)
    fig.suptitle('H2 RT distributions: observed dashed orange; predictive blue median / 95% pointwise band');save(fig,'rt_distribution_ppc')
    fig,axes=plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
    metrics=['correct_rt_good_minus_bad','accuracy_good_minus_bad','fast_accuracy_good_minus_bad','fast_match_choice_good_minus_bad','fast_correct_rt_good_minus_bad','dprime_good_minus_bad']
    for ax,metric in zip(axes.flat,metrics):
        t=table.loc[table.metric.eq(metric)].sort_values('group');interval(ax,t,'Seconds' if 'rt' in metric else 'Difference');ax.axhline(0,color='gray',ls=':');ax.set_title(metric.replace('_',' '),fontsize=10)
    axes.flat[0].legend(fontsize=8);fig.suptitle('H2: paired-subject Good-self minus Bad-self contrasts');save(fig,'fast_good_self_ppc')
    fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for row,stim in enumerate([0,1]):
        for col,metric in enumerate(['binned_accuracy_good_minus_bad','binned_match_choice_good_minus_bad']):
            ax=axes[row,col];t=table.loc[table.group.str.startswith(f'Self|stim={stim}|bin=')&table.metric.eq(metric)].sort_values('group');x=np.arange(1,6)
            ax.fill_between(x,t.lower.to_numpy(),t.upper.to_numpy(),color='#0072B2',alpha=.2);ax.plot(x,t['median'],'o-',color='#0072B2',label='Predictive median / 95% interval');ax.plot(x,t.observed,'x--',color='#D55E00',label='Observed');ax.axhline(0,color='gray',ls=':');ax.set_xticks(x);ax.set_xlabel('Fixed observed pooled-Self RT quintile bin');ax.set_ylabel('Good-self minus Bad-self');ax.set_title(('Match' if stim else 'Mismatch')+' '+metric.replace('binned_','').replace('_good_minus_bad',''))
    axes.flat[0].legend(fontsize=8);save(fig,'fast_shape_ppc')
    fig,ax=plt.subplots(figsize=(9,4),constrained_layout=True);interval(ax,table.loc[table.metric.eq('dprime_subject_mean')].sort_values('group'),'d-prime');ax.legend(fontsize=8);ax.set_title('H2 d-prime, conditional on legal response; loglinear count correction');save(fig,'dprime_ppc')

def main():
    OUT.mkdir(parents=True,exist_ok=True);start=time.time();status=OUT/'output_status.json'
    m.write_json(status,dict(status='running',pid=os.getpid(),fit=str(FIT),seed=2026262000))
    try:
        saved=np.load(FIT);assert bool(saved['diagnostics_pass']);draws=saved['draws'];labels=saved['labels'].tolist();model=str(saved['model']);assert model=='H2'
        frame,audit=m.prepare_real_matching();dataset=np.load(FIT.parent/'real_data.npz');assert str(dataset['data_kind'])=='observed'
        for name,value in [('choices',frame.choice),('rt',frame.rt_new),('trial_subject',frame.Subject),('source_row',frame.source_row)]:
            assert np.array_equal(dataset[name],value.to_numpy(),equal_nan=True),name
        m.write_json(OUT/'observed_data_audit.json',audit)
        for stim in [0,1]:EDGES[stim]=frame.loc[frame.choice.ge(0)&frame.stim.eq(stim)&frame.condition.isin(['Good_Self','Bad_Self']),'rt_new'].quantile([.2,.4,.6,.8]).to_numpy()
        m.write_json(OUT/'extended_statistic_definitions.json',dict(fixed_RT_quintile_boundaries={str(k):v.tolist() for k,v in EDGES.items()},dprime='Mean across participants of norm.ppf((hits+.5)/(legal Match n+1))-norm.ppf((false alarms+.5)/(legal Mismatch n+1)); conditional on responding, pooled mismatch relation design',denominators='Accuracy/choice conditioned on legal response; correct_all_trials includes timeout/illegal as not correct; timeout and illegal divided by all design events.',intervals='200 posterior predictive replicates, 50 independently selected draws per chain; 95% pointwise predictive intervals; not simultaneous bands',source_sha256=hashlib.sha256(FIT.read_bytes()).hexdigest(),model_fidelity='Reused frozen generator and observation model: mixture WFPT + uniform contaminant + timeout at shared per-subject posterior deadline + illegal Bernoulli; existing per-replicate survival integration assertion retained.'))
        old.posterior(draws,labels,model,OUT,FIT);old.reference_sensitivity(draws,labels,model,OUT)
        tab=pd.read_csv(OUT/'key_parameter_summary.csv');t=tab.loc[tab.parameter.isin(['mean_beta_v','mean_beta_dc'])]
        fig,ax=plt.subplots(figsize=(7,3),constrained_layout=True);ax.errorbar(t['median'],np.arange(len(t)),xerr=[t['median']-t.lower,t.upper-t['median']],fmt='o',capsize=5);ax.set_yticks(np.arange(len(t)));ax.set_yticklabels(['beta log(v): fixed zero','beta dc: zero-sum reference']);ax.axvline(0,color='gray',ls='--');ax.set_title('H2 target coefficients: posterior median / 95% credible interval');save(fig,'target_forest')
        old.stats=stats;old.plot_ppc=plots;old.ppc(FIT,draws,model,frame,OUT,50,2026262000,.001)
        summary=pd.read_csv(OUT/'ppc_summary.csv');summary['predictive_bias']=summary['median']-summary.observed;summary['outside_distance']=np.maximum(summary.lower-summary.observed,0)+np.maximum(summary.observed-summary.upper,0);summary.to_csv(OUT/'ppc_summary.csv',index=False)
        flags=summary.loc[~summary.observed_in_95_predictive_interval & ~summary.metric.str.contains('_cdf_|paired_subjects|bin_fraction')].copy();flags.to_csv(OUT/'ppc_misfit_flags.csv',index=False)
        m.write_json(status,dict(status='complete',model=model,elapsed_seconds=time.time()-start,replicates=200,metrics=len(summary),flagged_non_cdf_metrics=len(flags),exit_code=0))
        print(flags.to_string(index=False),flush=True)
    except Exception:
        m.write_json(status,dict(status='failed',error=traceback.format_exc(),elapsed_seconds=time.time()-start));raise

if __name__=='__main__':main()
