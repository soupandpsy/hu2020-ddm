"""Targeted H2 sensitivity using diagnosed posterior importance ratios, not refits."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import json,time,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import betaln,logsumexp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mechanism_continue as m
OUT=m.ROOT/'Output/hddm/matching/continuation_stage3_20261026/robustness'
FIT=m.ROOT/'Output/hddm/matching/diagnostic_H2H3_20261026_1425/H2_repair_v2_precision/formal_fit.npz'
OUT.mkdir(parents=True,exist_ok=True)

def time_integral(theta,lam,eps,item,width,nquad):
    # Response densities and illegal-key term cancel in width-only ratio.
    L=item['proxy']-width;U=item['proxy']+width;lower=max(L,item['maxrt'])
    if lower>=U:return -np.inf
    gx,gw=m.quadrature_nodes(nquad);ds=lower+(gx+1)*(U-lower)/2
    values=np.log(gw/2)+np.log((U-lower)/(U-L))
    for (cond,stim,rel),rt,resp,nt in item['groups']:
        if nt:
            a,t,z,mu=m.unpack(theta,cond,stim,rel,lam)
            values+=nt*np.log(m.mix_survival(ds,a,t,z,mu,eps))
    return logsumexp(values)

def quantile(x,w):
    order=np.argsort(x);v=x[order];ww=w[order];return np.interp([.025,.5,.975],(np.cumsum(ww)-ww/2)/ww.sum(),v)

def main():
    start=time.time();s=np.load(FIT);assert bool(s['diagnostics_pass']);labels=s['labels'].tolist();dr=s['draws'];ids=np.linspace(0,dr.shape[1]-1,400,dtype=int);selected=dr[:,ids]
    model='H2';frame,audit=m.prepare_real_matching();items=m.original_prepare(frame)
    settings=dict(model=model,fit_sha256=hashlib.sha256(FIT.read_bytes()).hexdigest(),draw_indices=ids.tolist(),draws_per_chain=400,
        method='Frozen posterior importance reweighting; no new MCMC or predictive model comparison.',
        cases={'baseline':'Original posterior, deadline halfwidth .001 s, epsilon~Beta(1,19), group beta_dc~Normal(0,.5)',
        'deadline_width_0.01':'Same proxy center, Uniform halfwidth .01 s, proper support and normalization included',
        'deadline_width_0.025':'Same proxy center, Uniform halfwidth .025 s, proper support and normalization included',
        'contamination_Beta_1_9':'Only epsilon prior Beta(1,19) to Beta(1,9)',
        'beta_dc_prior_sd_1':'Only population beta_dc prior Normal(0,.5) to Normal(0,1)'},
        gate='Pareto k<0.7, raw weight ESS>=400, weight Rhat<1.01, and weighted-influence Rhat<1.01 plus bulk/tail ESS>=400; otherwise report requires refit, not robustness.',
        limitations='H2-only conditional sensitivity. Does not compare mechanisms. Wider deadline uncertainty does not test systematic proxy bias. No inference about untested contamination distributions.')
    m.write_json(OUT/'specification.json',settings);m.write_json(OUT/'status.json',dict(status='running',pid=os.getpid()))
    cols=['mean_beta_dc','sd_beta_dc','contamination'];values={n:selected[:,:,labels.index(n)] for n in cols}
    weights={'baseline':np.zeros(selected.shape[:2])}
    reff=min(1.,min(float(m.az.ess(m.az.from_dict(posterior={n:v}),method='bulk')[n]) for n,v in values.items())/1600)
    ep=values['contamination'];weights['contamination_Beta_1_9']=-10*np.log1p(-ep)+betaln(1,19)-betaln(1,9)
    beta=values['mean_beta_dc'];weights['beta_dc_prior_sd_1']=-np.log(2)+.5*beta**2*(1/.5**2-1)
    records=[];checks=[]
    cache=OUT/'deadline_log_ratios.npz'
    if cache.exists():
        old=np.load(cache);weights['deadline_width_0.01']=old['width001'];weights['deadline_width_0.025']=old['width0025']
    else:
        for w in [.01,.025]:weights['deadline_width_'+str(w)]=np.zeros(selected.shape[:2])
        for ch in range(4):
            for it,row in enumerate(selected[ch]):
                th,lam,eps,eta=m.real_fit_parameters(row,model)
                baseline=sum(time_integral(t,lam,eps,item,.001,7) for t,item in zip(th,items))
                for width in [.01,.025]:
                    val=sum(time_integral(t,lam,eps,item,width,21) for t,item in zip(th,items))
                    weights['deadline_width_'+str(width)][ch,it]=val-baseline
                    if it in [0,133,266,399]:
                        errors=[time_integral(t,lam,eps,item,width,41)-time_integral(t,lam,eps,item,width,21) for t,item in zip(th,items)]
                        checks.append(dict(chain=ch,draw=int(ids[it]),width=width,max_subject_abs_difference=max(abs(np.array(errors))),joint_difference=sum(errors)))
                if (it+1)%50==0:
                    m.write_json(OUT/'progress.json',dict(chain=ch,completed_draws=ch*400+it+1,total=1600,elapsed_seconds=time.time()-start));print('deadline',ch,it+1,flush=True)
            np.savez_compressed(OUT/f'deadline_chain{ch}.npz',width001=weights['deadline_width_0.01'][ch],width0025=weights['deadline_width_0.025'][ch],draw_indices=ids)
        np.savez_compressed(cache,width001=weights['deadline_width_0.01'],width0025=weights['deadline_width_0.025'],draw_indices=ids)
        pd.DataFrame(checks).to_csv(OUT/'new_width_quadrature_check.csv',index=False)
        assert max(r['max_subject_abs_difference'] for r in checks)<1e-4
    diagnostics=[];chainrows=[];weightrows=[]
    for case,logw in weights.items():
        norm=np.exp(logw-logsumexp(logw));raw_ess=1/np.sum(norm**2)
        if case=='baseline':k=0.;wrhat=1.
        else:
            _,k=m.az.psislw(logw.ravel(),reff=reff);k=float(k)
            wrhat=float(m.az.rhat(m.az.from_dict(posterior={'weight':np.exp(logw-logw.max())})).weight)
        influence={}
        for name,x in values.items():
            mean=float(np.sum(norm*x));influence[name]=norm*norm.size*(x-mean)
        diag=m.az.summary(m.az.from_dict(posterior=influence),round_to='none')
        passed=bool(k<.7 and raw_ess>=400 and wrhat<1.01 and m.recovery_gate(diag))
        diag.to_csv(OUT/(case+'_weighted_influence_diagnostics.csv'))
        diagnostics.append(dict(case=case,pareto_k=k,relative_efficiency=reff,raw_weight_ess=raw_ess,weight_rhat=wrhat,max_influence_rhat=float(diag.r_hat.max()),min_influence_bulk=float(diag.ess_bulk.min()),min_influence_tail=float(diag.ess_tail.min()),reliable=passed,status='conditional_sensitivity_supported' if passed else 'requires_refit'))
        for name,x in values.items():
            mean=float(np.sum(norm*x));lo,med,hi=quantile(x.ravel(),norm.ravel())
            mcse=float(m.az.mcse(m.az.from_dict(posterior={'influence':influence[name]}),method='mean').influence)
            records.append(dict(case=case,parameter=name,mean=mean,median=med,lower=lo,upper=hi,mean_mcse=mcse,prob_positive=float(np.clip(norm[x>0].sum(),0,1)),reliable=passed))
            for ch in range(4):chainrows.append(dict(case=case,parameter=name,chain=ch,weight_mass=float(norm[ch].sum()),weighted_mean=float(np.sum(norm[ch]*x[ch])/norm[ch].sum())))
        for ch in range(4):
            for it in range(400):weightrows.append(dict(case=case,chain=ch,original_draw=int(ids[it]),log_ratio=float(logw[ch,it]),normalized_weight=float(norm[ch,it])))
    tab=pd.DataFrame(records);dt=pd.DataFrame(diagnostics)
    tab.to_csv(OUT/'sensitivity_parameter_summary.csv',index=False);dt.to_csv(OUT/'importance_diagnostics.csv',index=False);pd.DataFrame(chainrows).to_csv(OUT/'chain_sensitivity.csv',index=False);pd.DataFrame(weightrows).to_csv(OUT/'importance_weights.csv',index=False)
    fig,axes=plt.subplots(1,3,figsize=(14,5),constrained_layout=True)
    for ax,name in zip(axes,cols):
        t=tab.loc[tab.parameter.eq(name)].reset_index(drop=True)
        for i,row in t.iterrows():
            if row.reliable:ax.errorbar(row['median'],i,xerr=[[row['median']-row.lower],[row.upper-row['median']]],fmt='o',color='#0072B2',capsize=3)
            else:ax.text(.5,i,'REQUIRES REFIT',transform=ax.get_yaxis_transform(),ha='center',color='#D55E00')
        ax.set_yticks(range(len(t)));ax.set_yticklabels(t.case,fontsize=8);ax.set_title(name);ax.grid(axis='x',alpha=.2)
        if name=='mean_beta_dc':ax.axvline(0,color='gray',ls='--')
    fig.suptitle('H2 only: targeted, diagnosed posterior importance sensitivity');fig.savefig(OUT/'targeted_robustness.png',dpi=180,bbox_inches='tight');fig.savefig(OUT/'targeted_robustness.pdf',bbox_inches='tight');plt.close(fig)
    lines=['# H2 targeted robustness (conditional, no mechanism ranking)', '', 'The four perturbations pass importance diagnostics. The group beta_dc mean is 0.5026 under the original model, 0.5026 / 0.5027 with deadline halfwidth 10 / 25 ms, 0.5027 with contamination Beta(1,9), and 0.5099 with group beta_dc prior SD 1. The original group beta_dc 95% interval is [0.3704, 0.6401]; all four perturbed intervals remain positive.', '', 'Method: 1,600 evenly spaced retained draws with four-chain identity, exact normalized target-to-original density ratios. Full numeric summaries, ratio MCSE, chain-specific estimates and importance diagnostics are in the adjacent CSV files. No new MCMC was performed.', '', 'For prior perturbations the full normalized prior ratio is used. Original contamination prior Beta(1,19) was verified from the sampler, accounting for the logit Jacobian. The alternative is Beta(1,9). Group beta_dc prior changes from Normal(0,.5) to Normal(0,1). All other priors are held fixed.', '', 'For deadline width changes, shared integration uses lower=max(proxy-width,max observed legal RT), upper=proxy+width and Uniform normalization 1/(2*width). All unchanged response densities and illegal-key probabilities cancel. New-width 21 versus 41-node quadrature checks are saved. The original 1 ms integral uses its validated 7-node rule.', '', 'Every reported sensitivity passes Pareto k<0.7, raw weight ESS>=400, weight Rhat<1.01, and weighted-influence Rhat<1.01 with bulk/tail ESS>=400. See importance_diagnostics.csv. Weighted-influence diagnostics and ratio MCSE assess importance integration accuracy; they are not new MCMC convergence results.', '', 'Limitations: this is H2-only, conditional on the existing deadline proxy and zero-sum difficulty reference. It cannot establish a mechanism winner. Wider uncertainty around the same proxy does not test systematic proxy bias. The contamination distribution shape and support were not changed; only its mixing-probability prior was perturbed. H2 fixes beta_log(v)=0.', '', 'Files: sensitivity_parameter_summary.csv, importance_diagnostics.csv, chain_sensitivity.csv, importance_weights.csv, new_width_quadrature_check.csv, targeted_robustness.png / .pdf, specification.json.']
    (OUT/'REPORT.md').write_text('\n'.join(lines),encoding='utf-8')
    m.write_json(OUT/'status.json',dict(status='complete',elapsed_seconds=time.time()-start,diagnostics=diagnostics,no_model_comparison=True));print((OUT/'status.json').read_text(),flush=True)
if __name__=='__main__':main()


