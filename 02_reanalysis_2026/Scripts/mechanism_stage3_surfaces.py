import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import numpy as np,pandas as pd
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mechanism_stage3_sampler as s
m=s.m; out=s.BASE/'diagnostics';out.mkdir(exist_ok=True)
frame,_=m.prepare_real_matching();items=m.original_prepare(frame);records=[];summ=[]
fig,axes=plt.subplots(3,2,figsize=(12,12),constrained_layout=True)
for row,(model,sub) in enumerate([('H0',7304),('H1',7323),('H3',7338)]):
    warm=np.load(s.source(model))['draws'];labels=m.recovery_labels(model);active=m.allowed(model);k=len(active);center=warm.mean(axis=(0,1));th,lam,eps,eta=m.real_fit_parameters(center,model);si=m.subjects.index(sub)
    aa=warm[:,:,labels.index(f'subject_{sub}_log_a')];tt=warm[:,:,labels.index(f'subject_{sub}_log_t_excess')]
    ga=np.linspace(np.quantile(aa,.001)-.08,np.quantile(aa,.999)+.08,51);gt=np.linspace(np.quantile(tt,.001)-.03,np.quantile(tt,.999)+.03,101)
    z=np.empty((len(gt),len(ga)))
    for i,t in enumerate(gt):
        for j,a in enumerate(ga):
            cand=th[si].copy();cand[0]=a;cand[1]=t
            z[i,j]=m.subject_ll(cand,lam,eps,items[si])-.5*np.sum(((cand[active]-center[:k])/center[k:2*k])**2)
    z-=z.max();profile=z.max(axis=1);peaks,_=find_peaks(profile,prominence=1)
    axes[row,0].contourf(np.exp(ga),.1+np.exp(gt),z,levels=np.linspace(-20,0,21),extend='min')
    axes[row,0].set(xlabel='a',ylabel='t (s)',title=f'{model} / {sub}: fixed-nuisance slice')
    axes[row,1].plot(.1+np.exp(gt),profile);axes[row,1].set(xlabel='t (s)',ylabel='Relative conditional log posterior')
    f=frame[frame.Subject.eq(sub)&frame.choice.ge(0)]
    for rt in f.rt_new:
        if .1+np.exp(gt[0])<rt<.1+np.exp(gt[-1]):axes[row,1].axvline(rt,alpha=.15,color='gray')
    for i,t in enumerate(gt):
        for j,a in enumerate(ga):records.append(dict(model=model,subject=sub,log_a=a,log_t_excess=t,relative_logpost=z[i,j]))
    summ.append(dict(model=model,subject=sub,interior_peaks_prominence1=len(peaks),peak_t=(.1+np.exp(gt[peaks])).tolist(),peak_relative_logpost=profile[peaks].tolist()))
fig.suptitle('Targeted conditional slices; not marginalized posteriors or mechanism evidence')
fig.savefig(out/'targeted_t_a_surface.png',dpi=170);pd.DataFrame(records).to_csv(out/'targeted_t_a_surface.csv',index=False)
m.write_json(out/'targeted_t_a_summary.json',dict(summaries=summ,caveat='All nuisance coordinates fixed at pooled historical means, with a coarsely profiled. Modes here establish only conditional structure; they are not marginalized posterior probabilities. Gray lines denote observed response times.'))
print(summ)
