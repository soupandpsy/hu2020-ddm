"""Render a factual progress/results receipt and append it to the original notebook."""
from pathlib import Path
import json,datetime,hashlib,shutil
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'Output/hddm/matching/continuation_20261026'
def read(p):return json.loads(p.read_text()) if p.exists() else None
def main():
    now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(timespec='seconds')
    queue=read(OUT/'queue_status.json') or {}
    original=read(OUT/'recovery_summary/summary.json') or {}
    exact=read(OUT/'exact_loo_H3_r02_H1_8158/result.json')
    real={m:read(OUT/('fit_'+m+'_status.json')) for m in ['H0','H1','H2','H3']}
    lines=['# Matching H0–H3 后续分析：运行结果与进度',f'更新时间：{now}（北京时间）。','',
      '研究问题：现有 Matching 数据能否区分 Good-self 优先效应中的 log(v) 条件对比与 dc 对比。固定 H0 无目标效应、H1 仅 log(v)、H2 仅 dc、H3 两者都有；个体层面同时固定不允许的目标效应及其随机方差。',
      '', '**解释约束：dc 依赖零和 log 关系难度参照。等预测的重参数表示不能由这些数据本身选定心理学解释。**',
      '', '## 实际执行与保存',
      '- 核对32份历史拟合、原始数据哈希、Notebook函数、Windows/WSL/Docker进程；没有重跑300组数值核验、320份前向模拟或1296点网格。',
      '- 复用原 Python 3.9.13 / HDDM 1.0.1RC / NumPy 1.23.3 环境；新独立容器 hu2026-continuation，不恢复历史暂停容器。',
      '- 新任务每模型四条独立链、BLAS单线程，最多两个模型任务并发；容器上限8 CPU/4 GiB。初始Windows可用内存约3 GB，双拟合实测容器占用约0.7 GiB（共享页计数与各进程RSS不同）。',
      '- 每任务独立目录、固定seed、进度日志、退出码及阶段状态。历史结果与原始数据保持原样，不推送GitHub。',
      '- 每种机制20份需额外288拟合；历史中位间隔约38.13分钟，单队列约183小时。此次优先修复与真实数据，不启动这一多日扩充；每模型仍仅2个独立重复。',
      '', '## 模拟证据',
      f"32/32通过原设定全参数诊断。最差 Rhat={original.get('max_rhat','待读取')}，bulk ESS={original.get('min_bulk_ess','待读取')}，tail ESS={original.get('min_tail_ess','待读取')}。",
      'H1生成的两份数据，H1−H2 ΔELPD分别64.27（SE13.67）、51.96（SE12.53）；H2生成时分别−43.55（SE10.59）、−41.36（SE9.01）。SE由42名被试内逐试次差值求和后计算。与H3接近时不能据最大分数宣布机制被识别。',
      '参数恢复区间和覆盖只描述这些固定超参数下的有限重复，并非SBC，也不是精确覆盖率估计。混淆矩阵区分数值最高者和证据明确程度。',
      '', '## 异常点与精确留一',
      '异常为H3_r02数据、H1候选模型、零基索引8158（Subject7317），Pareto-k=0.7274236，样本量相关阈值0.676439。模拟choice=1、rt_new=0.411121秒；模板内原始RT/ACC/Resp不是该模拟观测，不能混用。',
      '80个跨链后验样本的条件log likelihood与完整联合似然减留一联合似然最大误差2.61e−14。精确重拟合只把该事件从训练集暂留出，保留42名被试及25199个训练事件，再预测留出事件；没有删除异常数据。',
      ('精确重拟合已完成：\n```json\n'+json.dumps(exact,ensure_ascii=False,indent=2)+'\n```') if exact else '精确重拟合尚在运行；原始H3_r02排名继续标记为不可靠，不能提前宣称已修复。',
      '', '## 真实数据状态']
    for m,s in real.items():lines.append(f"- {m}: {(s or {}).get('status','尚未启动')}" )
    lines += ['', '真实数据保留42×600=25200事件：合法响应24527、超时649、非法按键24。真实数据与模拟数据分别保存。尚未完成的拟合不输出后验结论。',
      '预测目标为同一批被试、给定记录派生截止代理的留出事件预测。截止代理取自无反应记录，因此该条件化评分不是完全无信息泄漏的前瞻交叉验证，不支持无条件的正式机制胜者声明。超阈值点安排精确重拟合；未修复前阻止正式排名。',
      '', '## PPC与敏感性预定义',
      'PPC每模型四链各50个后验样本，保留完整设计；检查条件×Match×关系下的选择比例、准确率、超时及非法率、正确/错误RT分位数。快速反应截止为观察数据中Self两种效价合并、分Match状态的RT第20百分位数，固定用于所有模型与复制数据，不按Good-self结果挑选。',
      '零和参照敏感性使用逐后验样本的等预测重参数变换；它是参数解释检查，不是不同先验下重新拟合的模型比较。另对H3将共享截止半宽由1ms放宽到10ms和25ms并重新拟合，属于条件截止不确定性检查，不能校正未知系统性代理偏差。',
      '', '## 输出索引',
      '- recovery_summary/：原始32拟合的参数恢复、偏差/描述性覆盖、混淆矩阵、ΔELPD及SE、Rhat/ESS/Pareto-k，CSV/JSON及PNG/PDF。',
      '- problem_point.json、problem_point_loglik_audit.csv：异常模拟观测与联合似然差核对。',
      '- exact_loo_H3_r02_H1_8158/：精确重拟合样本、诊断、heldout_loglik、修正评分（完成后）。',
      '- real_H0/ ... real_H3/：真实数据拟合、全参数诊断、LOO、精确重拟合（如需）、posterior_predictive/图表（完成后）。',
      '- sensitivity_H3_width_0.01/、sensitivity_H3_width_0.025/：截止敏感性（完成后）。',
      '- queue_status.json、各stage_status.json、*_console.log：实时进度、seed对应脚本及退出码。',
      '', '## 方法依据',
      '- [非因子分解模型的条件LOO](https://mc-stan.org/loo/articles/loo2-non-factorized.html)：使用联合似然比获得条件逐点预测密度。',
      '- [ArviZ精确重拟合](https://python.arviz.org/projects/stats/en/latest/api/generated/arviz_stats.reloo.html)：对超过样本量相关Pareto阈值的点重新拟合。',
      '', '## 当前队列', '```json',json.dumps(queue,ensure_ascii=False,indent=2),'```']
    comparisons=read(OUT/'real_comparison/summary.json')
    corrected=read(OUT/'recovery_corrected/summary.json')
    if corrected:lines += ['## 已完成的修正恢复比较','```json',json.dumps(corrected['comparisons'],ensure_ascii=False,indent=2),'```']
    if comparisons:lines += ['## 已完成的真实条件预测比较','```json',json.dumps(comparisons['comparisons'],ensure_ascii=False,indent=2),'```']
    estimates=[];ppc_checks=[]
    for folder in ['real_H0','real_H1','real_H2','real_H3','sensitivity_H3_width_0.01','sensitivity_H3_width_0.025']:
        path=OUT/folder/'posterior_predictive'
        if (path/'key_parameter_summary.csv').exists():
            tab=pd.read_csv(path/'key_parameter_summary.csv');tab['analysis']=folder;estimates.append(tab)
        if (path/'ppc_summary.csv').exists():
            tab=pd.read_csv(path/'ppc_summary.csv')
            flagged=tab.loc[~tab.observed_in_95_predictive_interval.astype(bool)].copy();flagged['analysis']=folder;ppc_checks.append(flagged)
    if estimates:
        table=pd.concat(estimates,ignore_index=True);table.to_csv(OUT/'all_key_posterior_estimates.csv',index=False)
        shown=table.loc[table.parameter.isin(['mean_beta_v','mean_beta_dc']),['analysis','parameter','mean','lower','upper','fixed_zero']]
        markdown='| '+' | '.join(shown.columns)+' |\n| '+' | '.join(['---']*len(shown.columns))+' |\n'
        markdown+='\n'.join('| '+' | '.join(str(v) if not isinstance(v,float) else f'{v:.5f}' for v in row)+' |' for row in shown.itertuples(index=False,name=None))
        lines += ['## 已完成的关键后验估计（95%等尾区间）',markdown]
    if ppc_checks:
        flags=pd.concat(ppc_checks,ignore_index=True);flags.to_csv(OUT/'ppc_outside_interval.csv',index=False)
        lines += ['## PPC检查',f'已完成的PPC中，有{len(flags)}项观察统计量落在逐项95%预测区间之外，详见ppc_outside_interval.csv。这是逐项描述性模型失配检查，不是校正多重比较后的拒绝检验。']
    if estimates and any(t.analysis.str.startswith('sensitivity').any() for t in estimates):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        combined=pd.concat(estimates);selected=combined.loc[combined.analysis.isin(['real_H3','sensitivity_H3_width_0.01','sensitivity_H3_width_0.025'])]
        fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
        for ax,p in zip(axes,['mean_beta_v','mean_beta_dc']):
            t=selected.loc[selected.parameter.eq(p)].reset_index(drop=True)
            ax.errorbar(t['median'],range(len(t)),xerr=[t['median']-t.lower,t.upper-t['median']],fmt='o',capsize=4)
            ax.set_yticks(range(len(t)));ax.set_yticklabels(t.analysis);ax.axvline(0,color='gray',ls='--');ax.set_title(p+'; 95% posterior interval')
        fig.suptitle('Shared-deadline uncertainty: H3 sensitivity under zero-sum reference')
        for ext in ['png','pdf']:fig.savefig(OUT/('deadline_sensitivity.'+ext),dpi=240,bbox_inches='tight')
        plt.close(fig)
    text='\n\n'.join(lines)+'\n';(OUT/'REPORT.md').write_text(text,encoding='utf-8')
    receipt=dict(updated_beijing=now,queue=queue,real_fit_status=real,exact=exact,scope='simulation evidence and actual real fitting progress; incomplete outputs are not conclusions')
    (OUT/'progress_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
    nbpath=ROOT/'matching_observation_recovery.ipynb';backup=OUT/'matching_observation_recovery_before_continuation.ipynb'
    if not backup.exists():shutil.copy2(nbpath,backup)
    nb=json.loads(nbpath.read_text(encoding='utf-8'));marker='continuation_20261026_receipt'
    cell=dict(cell_type='markdown',metadata={'tags':[marker]},source=text.splitlines(keepends=True))
    existing=[i for i,c in enumerate(nb['cells']) if marker in c.get('metadata',{}).get('tags',[])]
    if existing:nb['cells'][existing[0]]=cell
    else:nb['cells'].append(cell)
    nbpath.write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n',encoding='utf-8')
    print(str(OUT/'REPORT.md'))
if __name__=='__main__':main()
