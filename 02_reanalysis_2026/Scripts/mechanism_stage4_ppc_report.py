"""H2/H3 PPC comparison report from completed immutable PPC tables."""
import os
os.environ['MPLBACKEND']='Agg'
import json
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'Output/hddm/matching/continuation_stage4_20261027/H3_ppc'
assert json.loads((OUT/'output_status.json').read_text())['status']=='complete'
df=pd.read_csv(OUT/'H2_H3_ppc_comparison.csv');h3=pd.read_csv(OUT/'ppc_summary.csv')
sel=(df.group.str.endswith('|pooled')&df.metric.isin(['correct_mean_rt','error_mean_rt','timeout_rate','accuracy_given_response']))|df.metric.isin(['dprime_good_minus_bad','correct_rt_good_minus_bad','accuracy_good_minus_bad','fast_accuracy_good_minus_bad','fast_match_choice_good_minus_bad','fast_correct_rt_good_minus_bad'])
key=df.loc[sel].copy();key.to_csv(OUT/'H2_H3_key_ppc_comparison.csv',index=False)
counts=[]
for metric in ['correct_mean_rt','error_mean_rt','timeout_rate','accuracy_given_response']:
    t=df.loc[df.group.str.endswith('|pooled')&df.metric.eq(metric)]
    counts.append(dict(metric=metric,groups=len(t),H2_outside=int((~t.observed_in_95_predictive_interval_H2).sum()),H3_outside=int((~t.observed_in_95_predictive_interval_H3).sum()),both_outside=int(t.both_outside.sum())))
pd.DataFrame(counts).to_csv(OUT/'shared_misfit_counts.csv',index=False)
def save(fig,name):
    for suffix in ['png','pdf']:fig.savefig(OUT/(name+'.'+suffix),dpi=220,bbox_inches='tight')
    plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(15,5),constrained_layout=True)
for ax,metric in zip(axes,['timeout_rate','correct_mean_rt']):
    t=df.loc[df.group.str.endswith('|pooled')&df.metric.eq(metric)].sort_values('group');x=np.arange(len(t))
    for model,offset,color,marker in [('H2',-.12,'#0072B2','o'),('H3',.12,'#009E73','s')]:
        ax.errorbar(x+offset,t['median_'+model],yerr=[t['median_'+model]-t['lower_'+model],t['upper_'+model]-t['median_'+model]],fmt=marker,color=color,capsize=3,label=model+' predictive 95% interval')
    ax.scatter(x,t.observed_H2,marker='x',color='#D55E00',label='Observed',zorder=5);ax.set_xticks(x);ax.set_xticklabels(t.group.str.replace('|pooled','',regex=False),rotation=40,ha='right');ax.set_ylabel('Seconds' if metric.endswith('rt') else 'Proportion');ax.set_title(metric.replace('_',' '));ax.grid(axis='y',alpha=.15)
axes[0].legend(fontsize=9);fig.suptitle('H2 and H3: common absolute-fit checks');save(fig,'H2_H3_common_misfit')
fig,axes=plt.subplots(2,2,figsize=(10,8),constrained_layout=True)
for ax,(group,metric) in zip(axes.flat,[('Self|stim=1','correct_rt_good_minus_bad'),('Self|stim=1','accuracy_good_minus_bad'),('Self|stim=1','fast_accuracy_good_minus_bad'),('Self|dprime','dprime_good_minus_bad')]):
    t=df.loc[df.group.eq(group)&df.metric.eq(metric)].iloc[0]
    for x,model,color in [(0,'H2','#0072B2'),(1,'H3','#009E73')]:
        med=t['median_'+model];ax.errorbar(x,med,yerr=[[med-t['lower_'+model]],[t['upper_'+model]-med]],fmt='o',color=color,capsize=5)
    ax.axhline(t.observed_H2,color='#D55E00',ls='--',label='Observed');ax.axhline(0,color='gray',ls=':',alpha=.5);ax.set_xticks([0,1]);ax.set_xticklabels(['H2','H3']);ax.set_xlim(-.5,1.5);ax.set_title(metric.replace('_',' '),fontsize=10);ax.set_ylabel('Seconds' if 'rt' in metric else 'Difference');ax.legend(fontsize=8)
fig.suptitle('Good-self minus Bad-self: paired-subject contrasts / predictive 95% intervals');save(fig,'H2_H3_target_behavior_ppc')
lines=['|分组|指标|观测|H2预测中位数 [95%区间]|H3预测中位数 [95%区间]|两者均越界|','|---|---|---:|---:|---:|---|']
for r in key.itertuples():lines.append(f'|{r.group.replace("|"," / ")}|{r.metric}|{r.observed_H2:.5f}|{r.median_H2:.5f} [{r.lower_H2:.5f}, {r.upper_H2:.5f}]|{r.median_H3:.5f} [{r.lower_H3:.5f}, {r.upper_H3:.5f}]|{r.both_outside}|')
intro='''# H3完整PPC及H2/H3共同失配

H3复用stage3新正式拟合，全部参数最大Rhat=1.001241、最低bulk ESS=3094.94、tail ESS=3445.16，每链保留8000样本。原始源文件哈希、fit data_audit以及25200真实事件的选择、RT、被试和来源行数组均核对通过；未修改数据。H2复用已完成的stage3 PPC。

两模型均为四链各50个后验抽样，共200个完整设计预测复制，采用完全相同的统计定义。保留原WFPT密度、污染、非法选择以及每个被试共享的截止时间后验。每个副本的数值密度与生存概率一致性检查保留；没有重跑旧恢复实验。每个副本、后验抽样索引、种子、观测值、分组分位数、预测区间、偏差及越界标记均保留。

这里检验的是条件化、同被试、样本内PPC充分性；不能把哪个PPC更接近观测直接解释为正式留出模型排名。95%区间是逐统计量区间，不是同时置信带。多个相依统计量的越界数仅作诊断，不是全局显著性检验。200复制的尾分位数有Monte Carlo不确定性。

H3群体β_log(v)中位数0.38523 [0.26240, 0.50586]，β_dc中位数0.63795 [0.47673, 0.79828]。这些后验在当前模型与零和log关系参照下为正；共同PPC失配仍限制机制解释，不可宣称已充分解释实验。H2的β_log(v)固定为零，在森林图明确标识。

快速阈值沿用预先固定的观测Self合并q20规则，分Match/Mismatch。五分位边界亦在生成前从同一池化观测分布固定，没有按效应选择。Good−Bad差值为被试内配对后平均，保存每次实际配对人数。d-prime以Match为信号、Match选择为肯定反应，仅在合法响应内计算，使用(count+.5)/(n+1)校正后跨被试平均，不含超时作为选择。

## 共同失配计数

'''
intro+='\n'.join(f'- {r["metric"]}: H2 {r["H2_outside"]}/{r["groups"]}、H3 {r["H3_outside"]}/{r["groups"]}，共同越界 {r["both_outside"]}/{r["groups"]}。' for r in counts)
intro+='''

## 参照解释稳健性

H3等预测重参数化时，dc群体对比中位数随参照从约0.540（both-different=1）到0.714（identity-only=1），原零和log参照约0.638。所展示参照下区间均为正，但这是同一预测模型的坐标变换及推前先验，不是重新拟合的先验敏感性检验，也不能穷尽所有可能参照或识别唯一心理解释。β_log(v)对比在这些变换下不变。详见reference_sensitivity_summary.csv及notes中的预测等价误差。

## 关键数值

'''
images=['H2_H3_common_misfit','H2_H3_target_behavior_ppc','H2_H3_target_forest','difficulty_reference_sensitivity','rt_distribution_ppc','fast_shape_ppc','dprime_ppc']
report=intro+'\n'.join(lines)+'\n\n## 图形\n\n'+'\n\n'.join(f'![{name}]({name}.png)' for name in images)
(OUT/'REPORT.md').write_text(report,encoding='utf-8')
cells=[dict(cell_type='markdown',metadata={},source=intro.splitlines(True)),dict(cell_type='code',execution_count=None,metadata={},outputs=[],source=['from pathlib import Path\n','import pandas as pd\n','from IPython.display import display, Image\n','root=Path.cwd()\n',"if not (root/'H2_H3_key_ppc_comparison.csv').exists():\n","    root=root/'Output/hddm/matching/continuation_stage4_20261027/H3_ppc'\n","display(pd.read_csv(root/'H2_H3_key_ppc_comparison.csv'))\n"])]
for name in images:cells.append(dict(cell_type='code',execution_count=None,metadata={},outputs=[],source=[f"display(Image(filename=str(root/'{name}.png')))\n"]))
(OUT/'H2_H3_ppc_review.ipynb').write_text(json.dumps(dict(cells=cells,metadata=dict(kernelspec=dict(display_name='Python 3',language='python',name='python3')),nbformat=4,nbformat_minor=5),ensure_ascii=False,indent=2),encoding='utf-8')
print(pd.DataFrame(counts).to_string(index=False));print('\n'.join(lines))
