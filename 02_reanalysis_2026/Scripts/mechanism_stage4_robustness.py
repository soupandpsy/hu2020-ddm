"""Stage 4: diagnosed conditional sensitivity, preserving stage-3 artifacts."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import argparse,time,json,hashlib,traceback
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import logsumexp,betaln
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mechanism_continue as m
from mechanism_stage3_robustness import time_integral,quantile

ROOT=m.ROOT/'Output/hddm/matching'
def run(model):
    fit=ROOT/('diagnostic_H2H3_20261026_1425/H2_repair_v2_precision/formal_fit.npz' if model=='H2' else 'continuation_stage3_20261026/H3/formal_fit.npz')
    out=ROOT/'continuation_stage4_20261027'/('robustness_H2_support' if model=='H2' else 'robustness_H3');out.mkdir(parents=True,exist_ok=True)
    if (out/'status.json').exists() and json.loads((out/'status.json').read_text()).get('status')=='complete':print('Already complete',model,flush=True);return
    start=time.time();m.write_json(out/'status.json',dict(status='running',pid=os.getpid(),model=model))
    saved=np.load(fit);assert bool(saved['diagnostics_pass']);dr=saved['draws'];labels=saved['labels'].tolist();ids=np.linspace(0,dr.shape[1]-1,400,dtype=int);selected=dr[:,ids];del dr
    frame,audit=m.prepare_real_matching();items=m.original_prepare(frame)
    maxrt=float(frame.loc[frame.choice.ge(0),'rt_new'].max());assert maxrt<1.5
    cols=['mean_beta_dc','sd_beta_dc','contamination'] if model=='H2' else ['mean_beta_v','mean_beta_dc','sd_beta_v','sd_beta_dc','contamination']
    values={n:selected[:,:,labels.index(n)] for n in cols}
    reff=min(1.,min(float(m.az.ess(m.az.from_dict(posterior={n:v}),method='bulk')[n]) for n,v in values.items())/selected.shape[0]/selected.shape[1])
    weights={'baseline':np.zeros(selected.shape[:2]),'contamination_upper_1.5':np.zeros(selected.shape[:2])}
    if model=='H3':
        ep=values['contamination'];weights['contamination_Beta_1_9']=-10*np.log1p(-ep)+betaln(1,19)-betaln(1,9)
        weights['target_prior_scales_doubled']=sum(-np.log(2)+.5*values[name]**2*(1/sd**2-1/(2*sd)**2) for name,sd in [('mean_beta_v',.3),('mean_beta_dc',.5)])
        for width in [.01,.025]:weights['deadline_width_'+str(width)]=np.zeros(selected.shape[:2])
    spec=dict(model=model,fit_sha256=hashlib.sha256(fit.read_bytes()).hexdigest(),draw_indices=ids.tolist(),max_observed_legal_rt=maxrt,alternative_upper=1.5,original_upper=2.,relative_efficiency=reff,cases=list(weights),
        upper_change='Uniform contamination RT [.1,2] to [.1,1.5], with response density normalization AND deadline survival changed; all priors unchanged.',
        other_changes='H3 only: deadline halfwidth .001 to .01/.025; contamination prior Beta(1,19) to Beta(1,9); jointly double population beta_log(v) prior SD .3 to .6 and beta_dc SD .5 to 1, random-effect SD priors unchanged.',
        scope='Posterior importance reweighting, no refit. Diagnose each target separately. Not predictive model comparison. Wider proxy uncertainty is not systematic proxy-shift sensitivity.',
        gate='Pareto k<0.7, raw weight ESS>=400, weight Rhat<1.01 and all weighted-influence Rhat<1.01/bulk ESS>=400/tail ESS>=400. Failed targets need refits and their numbers are diagnostic-only.')
    m.write_json(out/'specification.json',spec);m.write_json(out/'data_audit.json',audit)
    checks=[]
    for ch in range(4):
        cache=out/f'likelihood_ratios_chain{ch}.npz'
        if cache.exists():
            old=np.load(cache)
            for name in weights:
                if name in old:weights[name][ch]=old[name]
            continue
        for it,row in enumerate(selected[ch]):
            theta,lam,eps,eta=m.real_fit_parameters(row,model)
            before=sum(m.subject_ll(th,lam,eps,item) for th,item in zip(theta,items))
            after=sum(m.subject_ll(th,lam,eps,item,upper=1.5) for th,item in zip(theta,items))
            weights['contamination_upper_1.5'][ch,it]=after-before
            if model=='H3':
                base=sum(time_integral(th,lam,eps,item,.001,7) for th,item in zip(theta,items))
                for width in [.01,.025]:
                    changed=sum(time_integral(th,lam,eps,item,width,21) for th,item in zip(theta,items))
                    weights['deadline_width_'+str(width)][ch,it]=changed-base
                    if it in [0,133,266,399]:
                        diff=[time_integral(th,lam,eps,item,width,41)-time_integral(th,lam,eps,item,width,21) for th,item in zip(theta,items)]
                        checks.append(dict(chain=ch,original_draw=int(ids[it]),width=width,max_subject_error=float(np.max(np.abs(diff))),joint_error=float(sum(diff))))
            if (it+1)%50==0:
                print(model,ch,it+1,'/400',flush=True);m.write_json(out/'progress.json',dict(completed=ch*400+it+1,total=1600,elapsed_seconds=time.time()-start))
        np.savez_compressed(cache,**{name:arr[ch] for name,arr in weights.items()},draw_indices=ids)
        pd.DataFrame([x for x in checks if x['chain']==ch]).to_csv(out/f'quadrature_chain{ch}.csv',index=False)
    if model=='H3':
        checktab=pd.concat([pd.read_csv(out/f'quadrature_chain{ch}.csv') for ch in range(4)])
        checktab.to_csv(out/'new_width_quadrature_check.csv',index=False);assert checktab.max_subject_error.max()<1e-4
    diagnostics=[];records=[];chainrecords=[];weightrows=[]
    for case,lw in weights.items():
        norm=np.exp(lw-logsumexp(lw));ess=1/(norm**2).sum()
        if case=='baseline':pareto=0.;wrhat=1.
        else:
            _,pareto=m.az.psislw(lw.ravel(),reff=reff);pareto=float(pareto);wrhat=float(m.az.rhat(m.az.from_dict(posterior={'weight':np.exp(lw-lw.max())})).weight)
        infl={name:norm*norm.size*(x-np.sum(norm*x)) for name,x in values.items()}
        diag=m.az.summary(m.az.from_dict(posterior=infl),round_to='none');diag.to_csv(out/(case+'_influence_diagnostics.csv'))
        passed=bool(pareto<.7 and ess>=400 and wrhat<1.01 and m.recovery_gate(diag))
        diagnostics.append(dict(case=case,pareto_k=pareto,raw_weight_ess=float(ess),max_normalized_weight=float(norm.max()),weight_rhat=wrhat,max_influence_rhat=float(diag.r_hat.max()),min_influence_bulk=float(diag.ess_bulk.min()),min_influence_tail=float(diag.ess_tail.min()),reliable=passed,status='conditional_sensitivity_supported' if passed else 'requires_refit'))
        for name,x in values.items():
            mean=float(np.sum(norm*x));lo,med,hi=quantile(x.ravel(),norm.ravel());mcse=float(m.az.mcse(m.az.from_dict(posterior={'influence':infl[name]}),method='mean').influence)
            records.append(dict(case=case,parameter=name,mean=mean,median=med,lower=lo,upper=hi,mean_mcse=mcse,prob_positive=float(np.clip(norm[x>0].sum(),0,1)),reliable=passed,interpretation='posterior sensitivity' if passed else 'DIAGNOSTIC ONLY: failed weights'))
            for ch in range(4):chainrecords.append(dict(case=case,parameter=name,chain=ch,weight_mass=float(norm[ch].sum()),weighted_mean=float(np.sum(norm[ch]*x[ch])/norm[ch].sum())))
        for ch in range(4):
            for it in range(400):weightrows.append(dict(case=case,chain=ch,original_draw=int(ids[it]),log_ratio=float(lw[ch,it]),normalized_weight=float(norm[ch,it])))
    dt=pd.DataFrame(diagnostics);tab=pd.DataFrame(records);dt.to_csv(out/'importance_diagnostics.csv',index=False);tab.to_csv(out/'sensitivity_parameter_summary.csv',index=False);pd.DataFrame(chainrecords).to_csv(out/'chain_sensitivity.csv',index=False);pd.DataFrame(weightrows).to_csv(out/'importance_weights.csv',index=False)
    fig,axes=plt.subplots(1,len(cols),figsize=(5*len(cols),5),constrained_layout=True)
    for ax,name in zip(np.atleast_1d(axes),cols):
        t=tab.loc[tab.parameter.eq(name)].reset_index(drop=True)
        for i,row in t.iterrows():
            if row.reliable:ax.errorbar(row['median'],i,xerr=[[row['median']-row.lower],[row.upper-row['median']]],fmt='o',color='#0072B2',capsize=3)
            else:ax.text(.5,i,'REQUIRES REFIT',transform=ax.get_yaxis_transform(),ha='center',color='#D55E00')
        ax.set_yticks(range(len(t)));ax.set_yticklabels(t.case,fontsize=8);ax.set_title(name);ax.grid(axis='x',alpha=.2)
        if name.startswith('mean_beta'):ax.axvline(0,color='gray',ls='--')
    fig.suptitle(model+': diagnosed targeted importance sensitivity (failed targets withheld)')
    for ext in ['png','pdf']:fig.savefig(out/('targeted_robustness.'+ext),dpi=180,bbox_inches='tight')
    plt.close(fig)
    report=['# '+model+' targeted robustness','', 'Exact normalized importance ratios are evaluated on 400 evenly spaced draws from each of four converged original chains. Results are conditional on the original difficulty reference and deadline proxy. No model is refitted here.','', 'All observed legal RTs are below 1.5 seconds: maximum '+str(maxrt)+'. The contamination-support perturbation changes both density normalization and censored survival. It is a different observation model, not a sampler change.','', 'Gate results:']
    for row in diagnostics:report.append('- '+row['case']+': '+row['status']+'; k='+format(row['pareto_k'],'.3f')+', weight ESS='+format(row['raw_weight_ess'],'.1f')+'.')
    report+=['','Numbers for failed targets in CSV files are diagnostic only and cannot support robustness claims. A direct refit with the altered observation model is required.','', 'For H3, target group prior SDs are doubled jointly (.3 to .6 for beta_log(v), .5 to 1 for beta_dc); participant SD priors remain unchanged. Deadline halfwidth perturbations preserve the proxy center and proper Uniform normalization and use 21-node shared-deadline integration, checked at 41 nodes. These checks do not test systematic proxy bias.','', 'See sensitivity_parameter_summary.csv, importance_diagnostics.csv, chain_sensitivity.csv, importance_weights.csv, and targeted_robustness.png. No mechanism ranking follows from within-model sensitivity.']
    (out/'REPORT.md').write_text('\n'.join(report),encoding='utf-8');m.write_json(out/'status.json',dict(status='complete',model=model,elapsed_seconds=time.time()-start,diagnostics=diagnostics));print((out/'status.json').read_text(),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('model',choices=['H2','H3']);args=p.parse_args();run(args.model)
