"""Focused conditional posterior slices for previously poorly mixing t/a blocks."""
import sys,json,time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import mechanism_continue as m
out=m.ROOT/'Output/hddm/matching/diagnostic_H2H3_20261026_1425/common_sampler_review'
out.mkdir(parents=True,exist_ok=True)
model='H1';folder=m.BASE/'real_H1';saved=np.load(folder/'real_H1_fit.npz');draws=saved['draws'];active=m.allowed(model);k=len(active)
center=draws.mean(axis=(0,1));theta,lam,eps,eta=m.real_fit_parameters(center,model)
frame,audit=m.prepare_real_matching();items=m.original_prepare(frame)
group=np.zeros(11);sd=np.zeros(11);group[active]=center[:k];sd[active]=center[k:2*k]
records=[];summaries=[]
fig,axes=plt.subplots(3,2,figsize=(12,12),constrained_layout=True)
for row,sub in enumerate([7304,7323,7343]):
    si=m.subjects.index(sub);j0=2*k+4+si*k
    cha=draws[:,:,j0+list(active).index(0)];cht=draws[:,:,j0+list(active).index(1)]
    aa=np.linspace(np.quantile(cha,.005)-.06,np.quantile(cha,.995)+.06,61)
    tt=np.linspace(np.quantile(cht,.005)-.04,np.quantile(cht,.995)+.04,81)
    values=np.empty((len(tt),len(aa)))
    for ti,lt in enumerate(tt):
        for ai,la in enumerate(aa):
            th=theta[si].copy();th[0]=la;th[1]=lt
            values[ti,ai]=m.subject_ll(th,lam,eps,items[si])-.5*np.sum(((th[active]-group[active])/sd[active])**2)
    values-=values.max();profile=values.max(axis=1);peaks,properties=find_peaks(profile,prominence=1.)
    for ti,lt in enumerate(tt):
        for ai,la in enumerate(aa):records.append(dict(subject=sub,log_a=la,log_t_excess=lt,relative_conditional_log_posterior=values[ti,ai]))
    ax=axes[row,0];ax.contourf(np.exp(aa),.1+np.exp(tt),values,levels=np.linspace(-20,0,21),extend='min')
    for ch in range(4):ax.scatter(np.exp(cha[ch,::100]),.1+np.exp(cht[ch,::100]),s=5,alpha=.3,label=f'Chain {ch+1}')
    ax.set(xlabel='a',ylabel='t (seconds)',title=f'Subject {sub}: nuisance parameters held fixed')
    axes[row,1].plot(.1+np.exp(tt),profile)
    legal=frame.loc[frame.Subject.eq(sub)&frame.choice.ge(0),'rt_new']
    for rt in sorted(legal[legal.between(.1+np.exp(tt.min()),.1+np.exp(tt.max()))].unique()):axes[row,1].axvline(rt,color='gray',alpha=.15,lw=.5)
    axes[row,1].set(xlabel='t (seconds)',ylabel='Profile conditional log posterior (max=0)',title='Gray lines: observed RTs; diagnostic slice only')
    summaries.append(dict(subject=sub,interior_peaks_prominence_gt_1=len(peaks),peak_t=(.1+np.exp(tt[peaks])).tolist(),peak_relative_log_posterior=profile[peaks].tolist(),chain_t_means=(.1+np.exp(cht)).mean(axis=1).tolist()))
axes[0,0].legend(fontsize=8);fig.suptitle('Historical H1 failed fit: targeted t/a mixing investigation, NOT posterior inference')
for ext in ['png','pdf']:fig.savefig(out/('conditional_t_a_surface.'+ext),dpi=180,bbox_inches='tight')
pd.DataFrame(records).to_csv(out/'conditional_t_a_surface.csv',index=False)
m.write_json(out/'surface_summary.json',dict(model=model,source='historical failed H1',subjects=summaries,interpretation='Conditional slices with all other parameters fixed to overall chain means; these are not marginalized posterior densities. Used to diagnose local proposal geometry, not infer effects.'))
print(json.dumps(summaries,indent=2))
