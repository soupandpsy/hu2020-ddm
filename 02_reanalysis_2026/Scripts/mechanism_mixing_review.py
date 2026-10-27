"""Diagnostics only for failed real fits; no posterior inference or model ranking."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE=Path(__file__).resolve().parents[1]/'Output/hddm/matching/continuation_20261026'
out=BASE/'real_mixing_review';out.mkdir(exist_ok=True)
summaries=[];chains=[]
for model in ['H0','H1']:
    folder=BASE/('real_'+model)
    diag=pd.read_csv(folder/('real_'+model+'_diagnostics.csv'),index_col=0)
    bad=(diag.r_hat>=1.01)|(diag.ess_bulk<400)|(diag.ess_tail<400)
    diag.loc[bad].to_csv(out/(model+'_failed_parameters.csv'))
    fit=np.load(folder/('real_'+model+'_fit.npz'));draws=fit['draws'];labels=fit['labels'].tolist()
    targets=diag.sort_values('r_hat',ascending=False).head(6).index.tolist()
    fig,axes=plt.subplots(len(targets),2,figsize=(13,2.2*len(targets)),constrained_layout=True)
    for row,name in enumerate(targets):
        if name not in labels:continue
        j=labels.index(name)
        for ch in range(4):
            v=draws[ch,:,j]
            axes[row,0].plot(v,lw=.35,alpha=.7,label=f'Chain {ch+1}')
            axes[row,1].hist(v,bins=45,histtype='step',density=True,label=f'Chain {ch+1}')
            chains.append(dict(model=model,parameter=name,chain=ch,mean=v.mean(),sd=v.std(),q025=np.quantile(v,.025),q975=np.quantile(v,.975)))
        axes[row,0].set_ylabel(name);axes[row,0].set_xlabel('Retained draw')
    axes[0,1].legend();fig.suptitle(model+': FAILED convergence — diagnostic plots only')
    for ext in ['png','pdf']:fig.savefig(out/(model+'_mixing.'+ext),dpi=180,bbox_inches='tight')
    plt.close(fig)
    summaries.append(dict(model=model,chains=draws.shape[0],retained_per_chain=draws.shape[1],max_rhat=float(diag.r_hat.max()),min_bulk=float(diag.ess_bulk.min()),min_tail=float(diag.ess_tail.min()),bad_rhat=int((diag.r_hat>=1.01).sum()),bad_bulk=int((diag.ess_bulk<400).sum()),bad_tail=int((diag.ess_tail<400).sum()),failed_any=int(bad.sum())))
    del draws,fit
pd.DataFrame(chains).to_csv(out/'chain_locations.csv',index=False)
pd.DataFrame(summaries).to_csv(out/'diagnostic_summary.csv',index=False)
(out/'summary.json').write_text(json.dumps(dict(status='real_fit_mixing_failure',models=summaries,note='Saved fits are diagnostic evidence only. Do not report posterior intervals, PPC or model ranking as reliable.'),indent=2))
print(json.dumps(summaries,indent=2))
