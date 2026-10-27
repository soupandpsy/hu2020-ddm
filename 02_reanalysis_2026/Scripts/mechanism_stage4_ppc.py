"""Same predeclared stage-three PPC statistics for an additional passed fit."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import argparse, inspect, json, time, traceback, hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mechanism_stage3_ppc as ext
m=ext.m;old=ext.old

def run():
    p=argparse.ArgumentParser();p.add_argument('--fit',required=True);p.add_argument('--output',required=True);p.add_argument('--seed',type=int,default=2026273000);a=p.parse_args()
    fit=Path(a.fit).resolve();out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=True);ext.OUT=out
    start=time.time();status=out/'output_status.json';m.write_json(status,dict(status='running',pid=os.getpid(),fit=str(fit),seed=a.seed))
    try:
        saved=np.load(fit);assert bool(saved['diagnostics_pass']);draws=saved['draws'];labels=saved['labels'].tolist();model=str(saved['model'])
        frame,audit=m.prepare_real_matching();sourceaudit=json.loads((fit.parent/'data_audit.json').read_text());assert sourceaudit==audit
        # Reconstruct observed events against immutable historical array reference.
        ref=m.ROOT/'Output/hddm/matching/diagnostic_H2H3_20261026_1425/H2_repair_v2_precision/real_data.npz';data=np.load(ref);assert str(data['data_kind'])=='observed'
        for name,value in [('choices',frame.choice),('rt',frame.rt_new),('trial_subject',frame.Subject),('source_row',frame.source_row)]:assert np.array_equal(data[name],value.to_numpy(),equal_nan=True),name
        m.write_json(out/'observed_data_audit.json',audit)
        m.write_json(out/'provenance.json',dict(fit_sha256=hashlib.sha256(fit.read_bytes()).hexdigest(),historical_observed_reference=str(ref),reference_sha256=hashlib.sha256(ref.read_bytes()).hexdigest(),verification='Current raw-event arrays match immutable observed-data reference; fit data_audit equals current source hashes and coding audit.'))
        for stim in [0,1]:ext.EDGES[stim]=frame.loc[frame.choice.ge(0)&frame.stim.eq(stim)&frame.condition.isin(['Good_Self','Bad_Self']),'rt_new'].quantile([.2,.4,.6,.8]).to_numpy()
        m.write_json(out/'extended_statistic_definitions.json',dict(fixed_RT_quintile_boundaries={str(k):v.tolist() for k,v in ext.EDGES.items()},statistics='Same mechanism_stage3_ppc.stats definitions as H2, no data-dependent revision; d-prime conditional on legal response with loglinear count correction.',replicates=200,intervals='95% pointwise predictive intervals, not simultaneous bands',sampler_scope='Same participants and full design; conditional shared posterior deadline.'))
        old.posterior(draws,labels,model,out,fit);old.reference_sensitivity(draws,labels,model,out)
        # Reuse plotted statistics, change model labels only in the isolated namespace.
        exec(compile(inspect.getsource(ext.plots).replace('H2',model),'<stage4_model_plot_labels>','exec'),ext.__dict__)
        old.stats=ext.stats;old.plot_ppc=ext.plots;old.ppc(fit,draws,model,frame,out,50,a.seed,.001)
        summary=pd.read_csv(out/'ppc_summary.csv');summary['predictive_bias']=summary['median']-summary.observed;summary['outside_distance']=np.maximum(summary.lower-summary.observed,0)+np.maximum(summary.observed-summary.upper,0);summary.to_csv(out/'ppc_summary.csv',index=False)
        flags=summary.loc[~summary.observed_in_95_predictive_interval & ~summary.metric.str.contains('_cdf_|paired_subjects|bin_fraction')];flags.to_csv(out/'ppc_misfit_flags.csv',index=False)
        h2=m.ROOT/'Output/hddm/matching/continuation_stage3_20261026/H2_ppc';both=[]
        for mod,folder in [('H2',h2),(model,out)]:
            t=pd.read_csv(folder/'key_parameter_summary.csv');both.append(t.loc[t.parameter.isin(['mean_beta_v','mean_beta_dc'])].assign(model=mod))
        params=pd.concat(both,ignore_index=True);params.to_csv(out/'H2_H3_target_forest.csv',index=False)
        fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
        for ax,name in zip(axes,['mean_beta_v','mean_beta_dc']):
            t=params.loc[params.parameter.eq(name)];ax.errorbar(t['median'],np.arange(len(t)),xerr=[t['median']-t.lower,t.upper-t['median']],fmt='o',capsize=4);ax.set_yticks(np.arange(len(t)));ax.set_yticklabels([str(r.model)+(' (fixed zero)' if r.fixed_zero else '') for r in t.itertuples()]);ax.axvline(0,color='gray',ls='--');ax.set_title(name.replace('mean_',''));ax.set_xlabel('Posterior median / 95% credible interval')
        fig.suptitle('Target contrasts under zero-sum log difficulty reference');old.savefig(fig,out/'H2_H3_target_forest')
        other=pd.read_csv(h2/'ppc_summary.csv');comp=other.merge(summary,on=['group','metric'],suffixes=('_H2','_H3'),validate='one_to_one');assert np.allclose(comp.observed_H2,comp.observed_H3,equal_nan=True);comp['both_outside']=(~comp.observed_in_95_predictive_interval_H2)&(~comp.observed_in_95_predictive_interval_H3);comp.to_csv(out/'H2_H3_ppc_comparison.csv',index=False)
        m.write_json(status,dict(status='complete',model=model,elapsed_seconds=time.time()-start,replicates=200,metrics=len(summary),flagged_non_cdf_metrics=len(flags),exit_code=0))
    except Exception:
        m.write_json(status,dict(status='failed',error=traceback.format_exc(),elapsed_seconds=time.time()-start));raise
if __name__=='__main__':run()
