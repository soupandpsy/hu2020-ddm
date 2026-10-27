"""Render stage-4 sensitivity summaries with explicit failed-target labels."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]/'Output/hddm/matching/continuation_stage4_20261027'
for folder in [root/'robustness_H2_support',root/'robustness_H3']:
    source=folder/'sensitivity_parameter_summary.csv'
    if not source.exists():continue
    tab=pd.read_csv(source);cols=tab.parameter.unique();ncols=min(3,len(cols));nrows=int(np.ceil(len(cols)/ncols))
    fig,axes=plt.subplots(nrows,ncols,figsize=(6*ncols,3.7*nrows),squeeze=False,constrained_layout=True)
    for ax,name in zip(axes.flat,cols):
        t=tab.loc[tab.parameter.eq(name)].reset_index(drop=True)
        for i,row in t.iterrows():
            if row.reliable:ax.errorbar(row['median'],i,xerr=[[row['median']-row.lower],[row.upper-row['median']]],fmt='o',color='#0072B2',capsize=3)
            else:ax.text(.5,i,'REQUIRES REFIT',transform=ax.get_yaxis_transform(),ha='center',va='center',color='#D55E00')
        ax.set_ylim(-.5,len(t)-.5);ax.set_yticks(range(len(t)));ax.set_yticklabels(t.case,fontsize=8);ax.set_title(name);ax.grid(axis='x',alpha=.2)
        if name.startswith('mean_beta'):ax.axvline(0,color='gray',ls='--')
    for ax in list(axes.flat)[len(cols):]:ax.axis('off')
    fig.suptitle(folder.name.replace('robustness_','')+': diagnosed targeted importance sensitivity')
    for ext in ['png','pdf']:fig.savefig(folder/('targeted_robustness.'+ext),dpi=180,bbox_inches='tight')
    plt.close(fig)
