"""Four-model PPC synthesis; no fitting, ranking, or historical writes."""
import os
os.environ['MPLBACKEND']='Agg'
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'Output/hddm/matching/continuation_stage4_20261027';OUT=BASE/'four_model_ppc';OUT.mkdir(exist_ok=True)
FOLDERS={m:BASE/(m+'_ppc') for m in ['H0','H1','H3']};FOLDERS['H2']=ROOT/'Output/hddm/matching/continuation_stage3_20261026/H2_ppc'
MODELS=['H0','H1','H2','H3'];COLORS={'H0':'#777777','H1':'#CC79A7','H2':'#0072B2','H3':'#009E73'};MARKERS={'H0':'v','H1':'^','H2':'o','H3':'s'}
tables=[];params=[];settings={}
for model in MODELS:
    folder=FOLDERS[model];state=json.loads((folder/'output_status.json').read_text());assert state['status']=='complete' and state['replicates']==200
    tables.append(pd.read_csv(folder/'ppc_summary.csv').assign(model=model));params.append(pd.read_csv(folder/'key_parameter_summary.csv').assign(model=model));settings[model]=json.loads((folder/'ppc_settings.json').read_text())
assert all(settings[m]['fast_cutoffs']==settings['H2']['fast_cutoffs'] for m in MODELS)
df=pd.concat(tables,ignore_index=True);par=pd.concat(params,ignore_index=True)
obs=df.pivot(index=['group','metric'],columns='model',values='observed');assert all(np.allclose(obs[m],obs['H2'],equal_nan=True) for m in MODELS)
df.to_csv(OUT/'four_model_ppc_summary.csv',index=False);par.to_csv(OUT/'four_model_key_parameter_summary.csv',index=False)
counts=[]
for model in MODELS:
    for metric in ['correct_mean_rt','error_mean_rt','timeout_rate','accuracy_given_response','match_choice_given_response','illegal_rate']:
        t=df.loc[df.model.eq(model)&df.group.str.endswith('|pooled')&df.metric.eq(metric)]
        counts.append(dict(model=model,metric=metric,groups=len(t),outside=int((~t.observed_in_95_predictive_interval).sum()),median_bias_min=t.predictive_bias.min(),median_bias_max=t.predictive_bias.max()))
count=pd.DataFrame(counts);count.to_csv(OUT/'four_model_misfit_counts.csv',index=False)
focus=df.loc[df.metric.isin(['dprime_good_minus_bad','correct_rt_good_minus_bad','accuracy_good_minus_bad','fast_accuracy_good_minus_bad','fast_match_choice_good_minus_bad','fast_correct_rt_good_minus_bad'])].copy();focus.to_csv(OUT/'four_model_target_behavior.csv',index=False)

def save(fig,name):
    for suffix in ['png','pdf']:fig.savefig(OUT/(name+'.'+suffix),dpi=220,bbox_inches='tight')
    plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(11,5),constrained_layout=True)
for ax,name in zip(axes,['mean_beta_v','mean_beta_dc']):
    for y,model in enumerate(MODELS):
        r=par.loc[par.model.eq(model)&par.parameter.eq(name)].iloc[0];ax.errorbar(r['median'],y,xerr=[[r['median']-r.lower],[r.upper-r['median']]],fmt=MARKERS[model],color=COLORS[model],capsize=4)
        if r.fixed_zero:ax.annotate('fixed zero',(r['median'],y),xytext=(7,7),textcoords='offset points',fontsize=9)
    ax.set_yticks(range(4));ax.set_yticklabels(MODELS);ax.axvline(0,color='gray',ls=':');ax.set_title('beta log(v)' if name=='mean_beta_v' else 'beta dc');ax.set_xlabel('Posterior median / 95% credible interval');ax.set_ylim(-.5,3.5)
fig.suptitle('Four-model target contrasts under zero-sum log difficulty reference');save(fig,'four_model_target_forest')
fig,axes=plt.subplots(2,2,figsize=(16,11),constrained_layout=True)
for ax,metric in zip(axes.flat,['correct_mean_rt','timeout_rate','accuracy_given_response','error_mean_rt']):
    groups=sorted(df.loc[df.group.str.endswith('|pooled'),'group'].unique());x=np.arange(len(groups))
    for j,model in enumerate(MODELS):
        t=df.loc[df.model.eq(model)&df.metric.eq(metric)].set_index('group').loc[groups]
        ax.errorbar(x+(j-1.5)*.12,t['median'],yerr=[t['median']-t.lower,t.upper-t['median']],fmt=MARKERS[model],color=COLORS[model],capsize=2,label=model,ms=4)
    ax.scatter(x,t.observed,color='#D55E00',marker='x',s=40,label='Observed',zorder=6);ax.set_xticks(x);ax.set_xticklabels([g.replace('|pooled','') for g in groups],rotation=40,ha='right');ax.set_ylabel('Seconds' if 'rt' in metric else 'Proportion');ax.set_title(metric.replace('_',' '));ax.grid(axis='y',alpha=.15)
axes.flat[0].legend(ncol=3,fontsize=9,loc='lower left');fig.suptitle('Four-model same-participant PPC: median / 95% pointwise predictive interval');save(fig,'four_model_absolute_fit_ppc')
fig,axes=plt.subplots(2,3,figsize=(15,9),constrained_layout=True)
for ax,(group,metric) in zip(axes.flat,[('Self|stim=1','correct_rt_good_minus_bad'),('Self|stim=0','correct_rt_good_minus_bad'),('Self|dprime','dprime_good_minus_bad'),('Self|stim=1','accuracy_good_minus_bad'),('Self|stim=1','fast_accuracy_good_minus_bad'),('Self|stim=0','fast_accuracy_good_minus_bad')]):
    for x,model in enumerate(MODELS):
        r=df.loc[df.model.eq(model)&df.group.eq(group)&df.metric.eq(metric)].iloc[0];ax.errorbar(x,r['median'],yerr=[[r['median']-r.lower],[r.upper-r['median']]],fmt=MARKERS[model],color=COLORS[model],capsize=4)
    ax.axhline(r.observed,color='#D55E00',ls='--',label='Observed');ax.axhline(0,color='gray',ls=':',alpha=.5);ax.set_xticks(range(4));ax.set_xticklabels(MODELS);ax.set_xlim(-.5,3.5);ax.set_ylabel('Seconds' if 'rt' in metric else 'Difference');ax.set_title(group.replace('|',' / ')+'\n'+metric.replace('_',' '),fontsize=10);ax.legend(fontsize=8)
fig.suptitle('Good-self minus Bad-self: paired-subject behavior / predictive 95% intervals');save(fig,'four_model_target_behavior_ppc')
fig,axes=plt.subplots(2,2,figsize=(13,9),constrained_layout=True)
for row,stim in enumerate([0,1]):
    for col,metric in enumerate(['binned_accuracy_good_minus_bad','binned_match_choice_good_minus_bad']):
        ax=axes[row,col]
        for j,model in enumerate(MODELS):
            t=df.loc[df.model.eq(model)&df.group.str.startswith(f'Self|stim={stim}|bin=')&df.metric.eq(metric)].sort_values('group');x=np.arange(1,6)+(j-1.5)*.04
            ax.errorbar(x,t['median'],yerr=[t['median']-t.lower,t.upper-t['median']],fmt=MARKERS[model]+'-',color=COLORS[model],capsize=2,label=model,ms=4)
        ax.plot(np.arange(1,6),t.observed,'x--',color='#D55E00',label='Observed');ax.axhline(0,color='gray',ls=':',alpha=.5);ax.set_xticks(range(1,6));ax.set_xlabel('Fixed pooled-Self observed RT quintile bin');ax.set_ylabel('Good-self minus Bad-self');ax.set_title(('Match' if stim else 'Mismatch')+' '+metric.replace('binned_','').replace('_good_minus_bad',''))
axes.flat[0].legend(ncol=3,fontsize=8);fig.suptitle('Fast-response shape: fixed bins and 95% pointwise predictive intervals');save(fig,'four_model_fast_shape_ppc')

lines=['|模型|指标|越界分组 / 总分组|','|---|---|---|']+[f'|{r.model}|{r.metric}|{r.outside}/{r.groups}|' for r in count.itertuples()]
targetlines=['|模型|分组|指标|观测|预测中位数 [95%区间]|','|---|---|---|---:|---:|']
for r in focus.itertuples():targetlines.append(f'|{r.model}|{r.group.replace("|"," / ")}|{r.metric}|{r.observed:.5f}|{r.median:.5f} [{r.lower:.5f}, {r.upper:.5f}]|')
intro='''# H0–H3 四模型完整后验预测检验

四模型均使用通过全参数MCMC诊断的真实数据正式后验。每个模型四链各随机抽50个后验样本，共200份完整42×600实验设计复制，模型间独立随机种子。真实观测表及固定快速反应阈值已逐模型核对完全一致。没有对失败链作推断，没有重跑恢复分析，没有修改历史模型、原始事件或先验。

PPC涵盖条件×Match/Mismatch×关系的选择/准确率、全部事件正确率、超时和非法率、正确及错误RT均值/分位数；另保存各条件跨关系汇总RT CDF、条件均值差、合法响应条件下的d-prime及快速反应形状。Good-self减Bad-self对比先按被试配对再取平均，每个复本保留实际配对人数。d-prime采用Match为信号、Match选择为肯定反应和(count+.5)/(n+1)校正，仅在合法响应内计算。快速q20及五分位边界沿用生成前固定的观测Self池化分布，分Match状态。

这些结果检验同被试、条件于原截止代理的样本内模型充分性，不是留出预测评分。95%区间是逐统计量预测区间，不是同时置信带。多个相依指标的越界数仅为描述，不是全局显著性检验。200复本的尾分位数存在Monte Carlo不确定性。PPC差异不作为正式模型排名，正式预测比较须使用另行完成且诊断可靠的留出评分。

群体β_log(v)、β_dc森林图展示各模型自己的后验，固定为零的参数明确标识。dc数值及心理解释依赖零和log关系难度参照；H3参照变换结果见H3_ppc。更灵活的H3也未必足以解释平均RT、超时和快速反应形状，不应把较好预测或正参数区间表述为唯一机制已经识别。

## 按预设汇总指标的失配

'''
images=['four_model_target_forest','four_model_absolute_fit_ppc','four_model_target_behavior_ppc','four_model_fast_shape_ppc']
report=intro+'\n'.join(lines)+'\n\n## Good-self目标行为\n\n'+'\n'.join(targetlines)+'\n\n## 图形\n\n'+'\n\n'.join(f'![{name}]({name}.png)' for name in images)
report+='\n\n所有模型的1437个分组统计量详见four_model_ppc_summary.csv；逐复本结果及RT分布图在各模型ppc目录。研究事实以CSV/JSON为准，既有图片和结果未覆盖。\n'
(OUT/'REPORT.md').write_text(report,encoding='utf-8')
for model in ['H0','H1']:
    folder=FOLDERS[model];sub=df.loc[df.model.eq(model)];ps=par.loc[par.model.eq(model)];local='# '+model+' 完整PPC\n\n200/200复本已完成。统计定义、预测范围与限制见../four_model_ppc/REPORT.md。源拟合通过全参数MCMC门槛；事件来源已核对。\n\n'
    local+='\n'.join(f'- {r.metric}: {r.outside}/{r.groups}个汇总分组的观测值位于95%逐统计量预测区间外。' for r in count.loc[count.model.eq(model)].itertuples())
    local+='\n\n'+ '\n\n'.join(f'![{name}]({name}.png)' for name in ['accuracy_timeout_ppc','rt_means_quantiles_ppc','rt_distribution_ppc','fast_good_self_ppc','fast_shape_ppc','dprime_ppc'])
    (folder/'REPORT.md').write_text(local,encoding='utf-8')
cells=[dict(cell_type='markdown',metadata={},source=intro.splitlines(True)),dict(cell_type='code',execution_count=None,metadata={},outputs=[],source=['from pathlib import Path\n','import pandas as pd\n','from IPython.display import display, Image\n','root=Path.cwd()\n',"if not (root/'four_model_ppc_summary.csv').exists():\n","    root=root/'Output/hddm/matching/continuation_stage4_20261027/four_model_ppc'\n","display(pd.read_csv(root/'four_model_misfit_counts.csv'))\n","display(pd.read_csv(root/'four_model_target_behavior.csv'))\n"])]
for name in images:cells.append(dict(cell_type='code',execution_count=None,metadata={},outputs=[],source=[f"display(Image(filename=str(root/'{name}.png')))\n"]))
(OUT/'four_model_ppc_review.ipynb').write_text(json.dumps(dict(cells=cells,metadata=dict(kernelspec=dict(display_name='Python 3',language='python',name='python3')),nbformat=4,nbformat_minor=5),ensure_ascii=False,indent=2),encoding='utf-8')
(OUT/'summary.json').write_text(json.dumps(dict(status='complete',models=MODELS,replicates_per_model=200,metric_rows=len(df),fixed_fast_cutoffs=settings['H2']['fast_cutoffs'],source_folders={k:str(v) for k,v in FOLDERS.items()},scope='PPC only, not predictive model ranking'),indent=2),encoding='utf-8')
print(count.to_string(index=False));print(focus.to_string(index=False))
