import sys
import numpy as np,pandas as pd
from scipy.stats import rankdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mechanism_stage4_sampler as s
m=s.m
model=sys.argv[1];stage=sys.argv[2] if len(sys.argv)>2 else 'diagnostic';folder=s.BASE/model
path=folder/('formal_fit.npz' if stage=='formal' else 'diagnostic_combined.npz');dr=np.load(path)['draws'];labels=m.recovery_labels(model)
selected=[f'subject_7338_{p}' for p in ['log_a','log_t_excess','v_intercept','dc_intercept']];rows=[]
fig,axes=plt.subplots(4,3,figsize=(15,10),constrained_layout=True)
for r,name in enumerate(selected):
    x=dr[:,:,labels.index(name)];rank=rankdata(x).reshape(x.shape)
    for c in range(4):
        axes[r,0].plot(x[c],lw=.4,alpha=.6,label=f'chain{c}');axes[r,1].hist(x[c],bins=50,histtype='step',density=True);axes[r,2].hist(rank[c],bins=20,histtype='step',density=True)
    axes[r,0].set_ylabel(name)
fig.suptitle(f'{model} {stage}: subject7338 trace / marginal density / pooled ranks');axes[0,0].legend();fig.savefig(folder/(stage+'_7338_rank.png'),dpi=160);plt.close(fig)
t=.1+np.exp(dr[:,:,labels.index('subject_7338_log_t_excess')]);a=np.exp(dr[:,:,labels.index('subject_7338_log_a')])
fig,ax=plt.subplots(figsize=(7,5),constrained_layout=True)
for c in range(4):
    mask=t[c]>.465;rows.append(dict(model=model,stage=stage,chain=c,n=len(mask),high_t_count=int(mask.sum()),high_t_fraction=float(mask.mean()),mode_changes=int(np.sum(np.diff(mask)!=0)),mean_t=float(t[c].mean()),mean_a=float(a[c].mean())))
    ax.scatter(a[c,::10],t[c,::10],s=3,alpha=.3,label=f'chain{c}')
ax.set(xlabel='a',ylabel='t (s)',title=f'{model} {stage}: subject7338');ax.legend();fig.savefig(folder/(stage+'_7338_joint.png'),dpi=160);pd.DataFrame(rows).to_csv(folder/(stage+'_7338_mode_occupancy.csv'),index=False)
