"""Evidence-only report. Never ranks incomplete or failed CV fits."""
import json,hashlib
from pathlib import Path
from datetime import datetime,timezone,timedelta
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
O=ROOT/'Output/hddm/matching';S3=O/'continuation_stage3_20261026';S4=O/'continuation_stage4_20261027';OUT=O/'continuation_stage5_20261027';OUT.mkdir(exist_ok=True,parents=True)
def js(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def md(df):
    fmt=lambda x: (format(x,'.5f') if isinstance(x,float) else str(x)).replace('|',r'\|')
    return '\n'.join(['|'+'|'.join(map(str,df.columns))+'|','|'+'|'.join(['---']*len(df.columns))+'|']+['|'+'|'.join(fmt(x) for x in row)+'|' for row in df.itertuples(index=False,name=None)])
rows=[]
for model in ['H0','H1','H2','H3']:
    path=(S4/model/'formal_summary.json') if model in ['H0','H1'] else (O/'diagnostic_H2H3_20261026_1425/H2_repair_v2_precision/formal_summary.json') if model=='H2' else S3/'H3/formal_summary.json'
    x=js(path);loo=js((S3 if model=='H2' else S4)/(model+'_loo/status.json'))
    rows.append(dict(model=model,retained_per_chain=x['retained_draws_per_chain'],max_rhat=x['max_rhat'],min_bulk_ess=x['min_bulk_ess'],min_tail_ess=x['min_tail_ess'],mcmc_pass=x['all_parameters_pass'],high_pareto_k_events=loo['bad_count'],max_pareto_k=loo['max_pareto_k']))
tab=pd.DataFrame(rows);tab.to_csv(OUT/'four_model_diagnostics.csv',index=False)
folds=[]
for folder in sorted((S4/'kfold5').glob('H?_fold?')):
    x=js(folder/'status.json');d=x.get('diagnostics',{});folds.append(dict(fit=folder.name,status=x['status'],max_rhat=d.get('max_rhat'),bulk=d.get('min_bulk_ess'),tail=d.get('min_tail_ess')))
ft=pd.DataFrame(folds);ft.to_csv(OUT/'initial_kfold_diagnostics.csv',index=False)
checks=[]
for x in js(S3/'historical_integrity_before.json'):
    p=ROOT/x['path'];h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    checks.append(dict(path=x['path'],sha256=h.hexdigest(),unchanged=h.hexdigest().lower()==x['sha256'].lower()))
assert all(x['unchanged'] for x in checks)
(OUT/'historical_integrity_check.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
ppcpaths={'H2':S3/'H2_ppc',**{model:S4/(model+'_ppc') for model in ['H0','H1','H3']}}
ppc=[]
for model,p in ppcpaths.items():
    t=pd.read_csv(p/'ppc_summary.csv')
    pick=t.loc[t.group.isin(['Self|stim=1','Self|stim=0','Self|dprime']) & t.metric.isin(['correct_rt_good_minus_bad','fast_accuracy_good_minus_bad','dprime_good_minus_bad','accuracy_good_minus_bad']),['group','metric','observed','median','lower','upper','observed_in_95_predictive_interval']].copy();pick.insert(0,'model',model);ppc.append(pick)
ppt=pd.concat(ppc).sort_values(['metric','group','model']);ppt.to_csv(OUT/'key_ppc_results.csv',index=False)
when=datetime.now(timezone(timedelta(hours=8))).isoformat(timespec='seconds')
lines=['# Hu et al. (2020) Matching 再分析：已执行结果与当前阻碍',f'更新时间：{when}（北京时间）。',
'## 当前结论',
'真实数据H0–H3四模型正式MCMC已经全部通过。完整四模型PPC也已完成，每模型200个复制实验。正式可靠的四模型预测比较尚未完成：原PSIS各模型都有高Pareto-k事件，第一轮统一五折的20个训练拟合均未通过短链诊断；没有用这些失败链产生模型排名或ΔELPD图。后续修复目录与日志保留真实进度。',
'所有模型都不能充分再现当前Matching行为：四模型各有8/8条件正确RT均值和8/8条件超时率位于各自95%逐统计量预测区间外；快速Match准确率优势亦全部低估。因此，不能因为H2先收敛或某个模型数值分数较高，就宣称原论文的v解释被推翻。',
'## 数据、比较范围与模型可比性',
'保留42名被试的全部25200事件：24527合法响应、649超时、24非法按键。上边界编码Match；全部来源行号可追溯。没有删除被试或事件。零和log关系难度参照、原始先验与观察似然在H0–H3基线拟合中不变；仅改采样转移，所以复用已通过的H2。污染支持上界变化另列为敏感性，绝不与原模型混合比较。',
'五折预测目标是同一批被试中留出的事件：每折每被试480个训练事件、120个测试事件，四模型完全共用划分。共享截止值只用训练响应条件化，逐个预测测试事件，不以其他测试响应为条件。截止代理本身仍来自完整原始记录；这是给定代理的条件预测，不是无信息泄漏的前瞻预测，也不能外推新被试。原PSIS保存的是共享截止积分后的条件逐事件密度，不能把它当作普通独立因子似然。',
'## 全数据正式诊断',md(tab),
'四模型均四链、每链8000保留样本，全参数（包括条件派生量）统一要求Rhat<1.01、bulk/tail ESS>=400。H0/H1对7338增加冻结的条件t–a混合MH后，正式每链跨峰约938–1068次；原失败短链为0–5次。全后验、chain日志、警告、rank/trace、t–a联合图及源代码哈希均保存。高总体接受率本身不能排除跨峰困难。',
'## 原论文与模拟证据',
'[Hu et al. (2020)](https://doi.org/10.1525/collabra.301)报告Good-self加工优势，并以原HDDM分析将其关联到信息积累率。这里是改变参数化及观察规则后的探索性再分析，不是原模型直接复现。行为现象与其机制解释须分开评价。',
'既有32/32模拟恢复拟合通过；H3两份模拟数据经精确留一修正后，H3对H1的ΔELPD为29.63±8.02和35.82±7.91（被试聚类SE）。每类只有两份独立数据，证据仅适用于这些生成设定，不能给出可靠恢复率、SBC或真实机制结论。本轮未重复已通过恢复、旧数值审计或网格搜索。',
'## 真实数据完整PPC',
'每模型4链各50个样本，沿用完整42×600实验设计与相同观察过程。正确率/选择以合法反应为分母，超时以所有事件为分母；d-prime在合法响应内按被试用loglinear校正。快速反应阈值和五分位边界固定为观测Self两种效价合并后、分别按Match/Mismatch计算，所有模型及复制保持一致。区间为95%逐统计量预测区间，非同时区间；越界数仅为失配诊断，不是多重检验后的显著性计数。',
'快速Match Good-self准确率优势观测为0.204；H0/H1/H2/H3预测中位数分别约−0.002/0.027/0.022/0.065，全部区间未覆盖。H3改善了部分整体准确率和RT差异，但仍不能再现快速形状；Match RT优势观测−82.9ms，H3约−43.2ms，H2约−15.8ms。',
md(ppt),
'[完整四模型PPC报告、森林图与Notebook](../continuation_stage4_20261027/four_model_ppc/REPORT.md)',
'## 稳健性与尚不能排除的解释',
'H2/H3对截止半宽10/25ms、污染比例先验及目标群体先验尺度的有限扰动，后验重要性诊断通过；目标效应在各自模型内保持正向。这些不是额外MCMC拟合，也不是模型选择证据。等预测关系参照变换保持预测，但改变H3中dc对比的数值；推前先验不等于在新参照下设独立先验再拟合。',
'改变污染RT均匀分布上界2→1.5秒时重要性重加权失效：H2 k=3.126、权重ESS=2.59；H3 k=2.767、ESS=1.02。已拒绝以这些伪后验区间声称稳健或不稳健；需要直接重拟合。放宽同一截止代理的局部范围不能排除代理系统偏差。所有模型共同PPC失配提示当前观察/认知模型族仍有局限，不能把某一系数唯一归结为特定心理阶段。',
'## 预测比较的实际阻碍',md(ft),
'第一轮五折使用全数据后验训练固定提议。减去训练事件后，多个被试出现不同的t/a/z/v联动模式；正在以各折自身短链重新训练精确校正的冻结条件提议。不能只延长旧链或放宽验收门槛。全部训练拟合与预测评分通过后，才允许汇总逐事件分差、42被试聚类SE及ΔELPD图。',
'## 产物与复现',
'- 基线H0/H1：../continuation_stage4_20261027/H0 与 H1；H3：../continuation_stage3_20261026/H3；H2：../diagnostic_H2H3_20261026_1425/H2_repair_v2_precision。\n- 四模型逐观测conditional loglik、Pareto-k及异常清单：stage3/H2_loo及stage4/H0_loo、H1_loo、H3_loo。\n- 五折分割、日志、20份失败短链及完整诊断：../continuation_stage4_20261027/kfold5。\n- 新修复与污染支持直接重拟合：本stage5目录，以各status.json为准。\n- four_model_diagnostics.csv、key_ppc_results.csv、initial_kfold_diagnostics.csv和historical_integrity_check.json与本报告同目录。\n- 原始数据、既有Notebook、历史拟合及核心脚本哈希核对全部未变；不推送GitHub。',
'方法参考：[非因子分解条件LOO](https://mc-stan.org/loo/articles/loo2-non-factorized.html)、[Pareto-k诊断及重拟合/K-fold处理](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)。']
comparison=js(OUT/'predictive_comparison/status.json') if (OUT/'predictive_comparison/status.json').exists() else {}
comparison_done=comparison.get('status')=='conditional_comparison_complete'
if comparison_done:
    lines[3]='真实数据四模型正式MCMC及完整PPC均已完成；修复后的20个五折训练拟合与留出评分均通过，条件预测比较已完成。请结合预测分差、被试聚类SE及全部模型的PPC失配解释，不能直接宣称唯一心理机制。'
    lines.insert(4,'[可靠条件五折比较及ΔELPD/SE](predictive_comparison/REPORT.md)')
lines.append('## 最新执行状态')
for sub in ['kfold5_repair','support_refits']:
    path=OUT/sub/'queue_status.json'
    if path.exists():lines.append(sub+': '+js(path).get('status','unknown')+'。详情见该目录状态、日志与报告。')
report='\n\n'.join(lines)+'\n';(OUT/'REPORT.md').write_text(report,encoding='utf-8')
nb=dict(nbformat=4,nbformat_minor=5,metadata=dict(kernelspec=dict(name='python3',display_name='Python 3',language='python')),cells=[dict(cell_type='markdown',metadata={},source=report.splitlines(True)),dict(cell_type='code',metadata={},execution_count=1,source=['import pandas as pd\n','pd.read_csv("four_model_diagnostics.csv")\n'],outputs=[dict(output_type='display_data',data={'text/plain':[tab.to_string(index=False)],'text/html':[tab.to_html(index=False)]},metadata={})])])
(OUT/'matching_reanalysis_current.ipynb').write_text(json.dumps(nb,ensure_ascii=False,indent=1),encoding='utf-8')
(OUT/'analysis_status.json').write_text(json.dumps(dict(updated=when,full_data_mcmc_pass=4,ppc_completed=4,reliable_four_model_comparison_complete=comparison_done,ranking_permitted=comparison_done,initial_kfold_pass=0,initial_kfold_total=20,historical_files_unchanged=True),indent=2),encoding='utf-8')
print(OUT/'REPORT.md')



