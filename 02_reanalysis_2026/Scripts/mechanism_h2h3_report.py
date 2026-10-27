from pathlib import Path
import json,datetime
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'Output/hddm/matching/diagnostic_H2H3_20261026_1425'
def get(path):return json.loads(path.read_text()) if path.exists() else None
now=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(timespec='seconds')
rows=[];status={};focus=[];integrity=[]
for name in ['H2','H3','H2_repair_v1','H3_repair_v1','H2_repair_v2','H3_repair_v2','H2_repair_v2_precision']:
    folder=BASE/name;status[name]=get(folder/'status.json')
    for stage in ['diagnostic','formal']:
        s=get(folder/(stage+'_summary.json'))
        if s:rows.append(dict(run=name,stage=stage,**s))
        dg=folder/(stage+'_a_t_diagnostics.csv')
        if dg.exists():
            tab=pd.read_csv(dg);group=tab.loc[tab.parameter.str.startswith(('mean_','sd_'))]
            indiv=tab.loc[tab.parameter.str.startswith('subject_')]
            focus.append(dict(run=name,stage=stage,group_a_t_max_rhat=group.r_hat.max(),group_a_t_min_bulk=group.ess_bulk.min(),group_a_t_min_tail=group.ess_tail.min(),
                individual_a_t_max_rhat=indiv.r_hat.max(),individual_a_t_min_bulk=indiv.ess_bulk.min(),individual_a_t_min_tail=indiv.ess_tail.min(),worst_parameter=indiv.sort_values('r_hat',ascending=False).iloc[0].parameter))
    acc=[];shifts=[]
    for file in folder.glob('diagnostic_chain[0-3].npz'):
        saved=np.load(file)
        integrity.append(dict(run=name,chain=file.stem,seed=int(saved['seed']),draws=saved['draws'].shape[0],all_draws_finite=bool(np.isfinite(saved['draws']).all())))
        if 'empirical_jump_counts' in saved:
            counts=saved['empirical_jump_counts']
            for i,(accepted,tried) in enumerate(counts):
                acc.append(dict(chain=file.stem,subject_index=i,accepted=int(accepted),proposed=int(tried),rate=accepted/tried if tried else np.nan))
        if 'mode_shift_counts' in saved:
            for i,(accepted,tried) in enumerate(saved['mode_shift_counts']):
                shifts.append(dict(chain=file.stem,subject_index=i,accepted=int(accepted),proposed=int(tried),rate=accepted/tried if tried else np.nan))
    if acc:pd.DataFrame(acc).to_csv(folder/'joint_proposal_acceptance.csv',index=False)
    if shifts:pd.DataFrame(shifts).to_csv(folder/'mode_shift_acceptance.csv',index=False)
pd.DataFrame(rows).to_csv(BASE/'diagnostic_comparison.csv',index=False)
pd.DataFrame(focus).to_csv(BASE/'a_t_focus.csv',index=False)
pd.DataFrame(integrity).to_csv(BASE/'chain_integrity.csv',index=False)
resources=get(BASE/'resource_summary.json')
text=f'''# 真实Matching H2/H3：并行短诊断与采样器修正

更新时间：{now}（北京时间）

## 执行设计

18逻辑核心，启动前CPU负载约16%，Windows可用内存约2.8 GiB，无活动研究采样。最多同时2个模型×4链，BLAS单线程；容器限制8核心、4 GiB。H0/H1历史拟合、原始数据未修改，无GitHub推送。

首轮每链1000次校准（舍弃500），再2000次诊断（舍弃500），保留1500次。链分别记录PID、seed、日志、退出码和样本。首轮严重t/a混合问题阻止正式采样。

## 共性问题与修正

历史H1中7304、7323、7343的固定其他参数t–a条件后验切片均有两个明显峰，提示单个局部二次提议不足以有效跨峰。此切片仅诊断，不是边际后验或心理机制证据。

修正版保留全部原始更新，额外加入四成分Student-t独立联合提议（3个t区间成分，权重各0.3；宽尾成分权重0.1）。初轮短链仅用于拟合并冻结提议，不改变先验。每个提议仍按完整原始似然、原先验和正反提议密度精确MH校正；保留原更新确保探索。

数值验证包含独立SciPy密度核对和双向详细平衡代数检查。修正版再次仅运行2000次/链，舍弃500，使用全新种子。短诊断尚未达标时不报告正式后验或模型排名。

v1没有解决全部问题：H2最大Rhat仍1.308，H3为1.905，且H3某链某被试独立提议接受率仅2.2%。v2使用前两轮样本共同冻结提议，另加入与当前状态无关、正反等概率的峰间位移及小尺度对称噪声；因此提议比为1，按原始后验差精确校正。v2也是新的2000次短验证，不与旧样本拼接作推断。

v2短诊断中H3全部通过并进入正式采样。H2的t/a最大Rhat降为1.00639，全部参数仍有relation_coef1 Rhat=1.01475、contamination Rhat=1.01208，最低tail ESS=364，因此安排一次2000次/链、保持同一冻结转移核的有界精度续接。该续接仅与v2同核的1500次保留样本合并检查，不合并早期不同核的训练样本；门槛未降低。

## 实际诊断

'''
if rows:
    text+='|运行|阶段|最大Rhat|最低bulk ESS|最低tail ESS|t/a最大Rhat|全部通过|\n|---|---|---:|---:|---:|---:|---|\n'
    for r in rows:text+=f"|{r['run']}|{r['stage']}|{r['max_rhat']:.5f}|{r['min_bulk_ess']:.1f}|{r['min_tail_ess']:.1f}|{r['a_t_max_rhat']:.5f}|{r['all_parameters_pass']}|\n"
text+='\n## 当前状态\n\n```json\n'+json.dumps(status,ensure_ascii=False,indent=2)+'\n```\n'
text+='\n## 实测资源\n\n```json\n'+json.dumps(resources,ensure_ascii=False,indent=2)+'\n```\n'
text+='\n## 产出\n\n每个模型目录包括specification.json、逐链NPZ与日志、全部参数diagnostics.csv、a/t诊断子表、链轨迹/密度PNG/PDF。修正版另有proposal_components.npz、proposal_validation.json、joint_proposal_acceptance.csv。common_sampler_review保存共性多峰证据。诊断失败的样本只能用于采样器分析。dc解释仍依赖零和关系难度参照。\n'
(BASE/'REPORT.md').write_text(text,encoding='utf-8')
(BASE/'status_summary.json').write_text(json.dumps(dict(updated_beijing=now,runs=status,diagnostics=rows,resources=resources),indent=2),encoding='utf-8')
nbpath=ROOT/'matching_observation_recovery.ipynb';nb=json.loads(nbpath.read_text(encoding='utf-8'));tag='h2h3_short_diagnostics_20261026'
cell=dict(cell_type='markdown',metadata={'tags':[tag]},source=text.splitlines(keepends=True))
indices=[i for i,c in enumerate(nb['cells']) if tag in c.get('metadata',{}).get('tags',[])]
if indices:nb['cells'][indices[0]]=cell
else:nb['cells'].append(cell)
nbpath.write_text(json.dumps(nb,ensure_ascii=False,indent=1)+'\n',encoding='utf-8')
print(str(BASE/'REPORT.md'))
