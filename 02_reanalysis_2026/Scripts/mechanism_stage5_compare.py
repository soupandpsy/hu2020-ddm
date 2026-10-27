"""Aggregate only twenty passed, complete, identically split K-fold fits."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import logsumexp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'Output/hddm/matching/continuation_stage5_20261027';CV=B/'kfold5_repair';OUT=B/'predictive_comparison';OUT.mkdir(parents=True,exist_ok=True)
missing=[];tables={};mcse={};foldrows=[]
for model in ['H0','H1','H2','H3']:
    tt=[];variance=0.
    for fold in range(5):
        d=CV/f'{model}_fold{fold}';p=d/'status.json'
        if not p.exists():missing.append(dict(model=model,fold=fold,status='not_started'));continue
        s=json.loads(p.read_text());g=s.get('diagnostics',{})
        if s['status']!='complete' or not g.get('all_parameters_pass') or not (d/'heldout_scores.csv').exists():missing.append(dict(model=model,fold=fold,status=s['status']));continue
        assert g['max_rhat']<1.01 and g['min_bulk_ess']>=400 and g['min_tail_ess']>=400
        t=pd.read_csv(d/'heldout_scores.csv');assert len(t)==5040 and t.groupby('Subject').size().eq(120).all();t['fold']=fold;tt.append(t)
        ll=np.load(d/'heldout_pointwise_loglik.npy',mmap_mode='r');assert ll.shape==(4,800,5040)
        chain_scores=np.array([np.sum(logsumexp(ll[ch],axis=0)-np.log(800)) for ch in range(4)])
        approx_mcse=float(np.std(chain_scores,ddof=1)/2);variance+=approx_mcse**2
        foldrows.append(dict(model=model,fold=fold,elpd=float(t.elpd.sum()),approx_mcse_from_four_chains=approx_mcse,max_rhat=g['max_rhat'],min_bulk=g['min_bulk_ess'],min_tail=g['min_tail_ess']))
    if len(tt)==5:
        t=pd.concat(tt).sort_values('observation').reset_index(drop=True);assert np.array_equal(t.observation,np.arange(25200));tables[model]=t;mcse[model]=np.sqrt(variance)
if missing:
    pd.DataFrame(missing).to_csv(OUT/'incomplete_fits.csv',index=False)
    (OUT/'status.json').write_text(json.dumps(dict(status='comparison_blocked',ranking_permitted=False,missing_or_failed=missing),indent=2))
    print('Comparison blocked:',len(missing),'incomplete/failed fits');raise SystemExit(0)
ref=tables['H0'];assert all(np.array_equal(t.source_row,ref.source_row) and np.array_equal(t.fold,ref.fold) for t in tables.values())
for model,t in tables.items():t.assign(model=model).to_csv(OUT/(model+'_heldout_pointwise.csv'),index=False)
total={m:t.elpd.sum() for m,t in tables.items()};best=max(total,key=total.get);rows=[];pairs=[]
for model,t in tables.items():
    delta=t.elpd-tables[best].elpd;clusters=delta.groupby(t.Subject).sum();se=float(np.sqrt(42*clusters.var(ddof=1)))
    rows.append(dict(model=model,elpd=float(total[model]),delta_elpd_to_numerical_best=float(delta.sum()),subject_cluster_se=se,approx_elpd_mcse=mcse[model],approx_difference_mcse=float(np.hypot(mcse[model],mcse[best])) if model!=best else 0))
for a,b in itertools.combinations(tables,2):
    d=tables[a].elpd-tables[b].elpd;cl=d.groupby(ref.Subject).sum();pairs.append(dict(model_a=a,model_b=b,delta_a_minus_b=float(d.sum()),subject_cluster_se=float(np.sqrt(42*cl.var(ddof=1))),approx_mcse=float(np.hypot(mcse[a],mcse[b]))))
tab=pd.DataFrame(rows).sort_values('elpd',ascending=False);tab.to_csv(OUT/'four_model_comparison.csv',index=False);pd.DataFrame(pairs).to_csv(OUT/'pairwise_differences.csv',index=False);pd.DataFrame(foldrows).to_csv(OUT/'fold_diagnostics_scores.csv',index=False)
fig,ax=plt.subplots(figsize=(9,4),constrained_layout=True);y=np.arange(4);ax.errorbar(tab.delta_elpd_to_numerical_best,y,xerr=tab.subject_cluster_se,fmt='o',capsize=5);ax.set_yticks(y);ax.set_yticklabels(tab.model);ax.axvline(0,color='gray',ls='--');ax.set_xlabel('Delta K-fold ELPD relative to numerical best (bars: +/-1 participant-cluster SE)');ax.set_title('Same-participant conditional five-fold prediction');ax.invert_yaxis();fig.savefig(OUT/'delta_elpd_se.png',dpi=180);fig.savefig(OUT/'delta_elpd_se.pdf');plt.close(fig)
state=dict(status='conditional_comparison_complete',ranking_permitted=True,numerical_best=best,models=rows,folds=5,heldout_events_per_model=25200,training_diagnostics_pass=20,interpretation='Conditional on full-record-derived deadline proxy and same-participant design; numerical predictive advantage is not a unique psychological mechanism or evidence of absolute model adequacy.',se='Paired score differences summed within 42 participants. Ignores extra uncertainty from choosing a single fold partition.',mcse='Approximate delta-method diagnostic from variation in four independent chain-specific aggregate log predictive scores, summed in quadrature across independent fold fits; only four chain estimates per fold, not a high-precision MCSE guarantee.')
(OUT/'status.json').write_text(json.dumps(state,indent=2))
lines=['# 四模型条件五折预测比较','20/20训练拟合及全部留出评分通过。每模型25200事件恰好预测一次，四模型同fold同source_row；没有删除高Pareto-k数据。','|模型|ELPD|相对最高分差|被试聚类SE|近似MCSE|','|---|---:|---:|---:|---:|']
for r in tab.itertuples():lines.append(f'|{r.model}|{r.elpd:.2f}|{r.delta_elpd_to_numerical_best:.2f}|{r.subject_cluster_se:.2f}|{r.approx_elpd_mcse:.2f}|')
lines+=['上述排名仅为当前固定截止代理、同被试事件预测目标下的相对预测表现。不能以微小分差宣布机制胜者；所有模型仍须结合完整PPC的RT、超时、快速反应形状失配解释。共享先验、编码、关系难度参照与观察模型未改变。','误差棒为±1 SE，不是后验可信区间。SE以42被试为聚类单元；单一分折方案的不确定性未额外估计。MCSE是四链变化的近似诊断。','![ELPD difference](delta_elpd_se.png)']
(OUT/'REPORT.md').write_text('\n\n'.join(lines),encoding='utf-8');print(json.dumps(state,indent=2))
