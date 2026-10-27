# 第二部分：H0–H3 的 v / dc 竞争机制再分析

[English](README.md) · [项目总览](../README_zh.md) · [综合结果 Notebook](matching_h0_h3_results.ipynb)

状态日期：**2026-10-28**。全数据四模型正式拟合及完整 PPC 已完成；**可靠五折预测比较尚未完成，第三部分尚未正式分析**。本文是第二部分当前入口；带日期的旧报告及 Notebook 是过程记录，不与本文并列作为最新结论。

## 研究问题 → 方法选择

原论文行为优势可以复现，但行为优势是否必须由刺激相关漂移 v 解释？在共同观察规则和难度参照下，Good-self−Bad-self 差异是否需要 v、方向性漂移偏置 dc，或二者共同变化？本分析仅针对 Matching；不把 Identity/Valence 历史未收敛估计当作机制证据，也不是对原论文 HDDM 参数的直接数值复现。

42 名被试 × 600 个事件，共 25,200：24,527 合法响应、649 超时、24 非法按键。上边界为 Match；条件按形状关联编码。合法选择、超时、非法事件分别进入观察模型，不能把无反应反编码为选择。原论文行为筛选口径与本分析全事件口径不同。

| 模型 | Good-self−Bad-self 的 log(v) 对比 | dc 对比 |
|---|---|---|
| H0 | 每名被试固定为 0 | 每名被试固定为 0 |
| H1 | 自由 | 固定为 0 |
| H2 | 固定为 0 | 自由 |
| H3 | 自由 | 自由 |

所有模型共享其余结构：层级 a/t/z、其他条件效应、零和 log 关系难度、污染与非法事件过程、被试共享截止积分。β_log(v) 不是历史 M4 原始尺度的 v 差 0.666。详细定义见 [历史设计方案](matching_model_design.md)；方案中的大规模恢复预算不代表已完成工作。

## 诊断发现 → 方法调整 → 结果

| 逻辑步骤 | 发现及调整 | 阅读入口 |
|---|---|---|
| 行为与历史探索 | 原论文选定统计复现；M0–M4 探索的 PPC、恢复及跨任务收敛暴露限制 | [第一部分](../01_reproduction_2020/README_zh.md)、[历史 M0–M4](README_historical_m0_m4_zh.md)、[历史结果](results.ipynb) |
| 回到实验设计 | 区分合法响应/超时/非法按键，审查 Mismatch 标签组成和快速选择曲线 | [行为诊断](matching_behavior_review.ipynb) |
| 检查可区分性 | 固定参数模拟不能替代恢复；难度参照影响 dc 解释，旧参数代入不等于新模型 PPC | [前向模拟](matching_mechanism_simulation.ipynb)、[结构核查](matching_mechanism_audit.ipynb) |
| 修正观察模型与采样 | 使用一致生成—似然、共同截止积分；改善转移与精度，不靠放宽门槛 | [观察模型与恢复日志](matching_observation_recovery.ipynb) |
| 检查真实拟合 | 四链每链 8,000 保留样本，四模型均通过全参数 R̂<1.01、bulk/tail ESS≥400 | [诊断表](Output/hddm/matching/continuation_stage5_20261027/four_model_diagnostics.csv) |
| 检查模型充分性 | 每模型 200 PPC 复本；四模型各有 8/8 条件正确 RT 均值和超时率失配，快速优势均低估 | [PPC 统计表](Output/hddm/matching/continuation_stage4_20261027/four_model_ppc/four_model_ppc_summary.csv)、[综合 Notebook](matching_h0_h3_results.ipynb) |
| 检查预测比较 | 各模型 PSIS 均有高 k；首轮 20 套五折短链均失败，转为逐折修复 | [首轮诊断](Output/hddm/matching/continuation_stage5_20261027/initial_kfold_diagnostics.csv)、[修复状态](Output/hddm/matching/continuation_stage5_20261027/kfold5_repair/queue_status.json) |

全数据最大 R̂：H0=1.001248、H1=1.001703、H2=1.002674、H3=1.001241。收敛是计算诊断，不能保证心理模型充分。快速 Match 准确率优势观测 0.204；H0/H1/H2/H3 预测中位数约 −0.002/0.027/0.022/0.065，95% 逐统计量预测区间均未覆盖观测。

修复后仅 H2 的 `fold4`（从 0 编号的第五折）完成正式拟合与留出预测；H0 同折 R̂=1.013752、tail ESS=63.26，阻止推广。**五折 × 四模型 = 20 套训练拟合，不是“五十个比较”。目前无可靠 ΔELPD 排名。**

## 科学结论边界与下一步

- 有限模拟恢复已有 32/32 拟合通过，但每个生成模型只有两份独立数据；不是可靠恢复率、覆盖率或 SBC 验证。修正后的模拟比较不等于真实数据模型排名。
- H3 改污染 RT 上界 2→1.5 秒的直接重拟合通过，两目标效应方向保持；H2 tail ESS=210.63 尚不足。新上界下的 H3 PPC 未生成，不能直接复用基线 PPC。
- 截止代理来自完整原始记录。即使五折完成，也仅支持给定代理、同被试留出事件的条件预测，不能声称无信息泄漏的前瞻预测或外推新被试。
- dc 的解释依赖难度参照；全模型共同 PPC 失配。不能宣称 H2 推翻原论文的 v 解释，不能宣布 H3 是唯一机制。
- 下一步需解决逐折后验混合、完成全部评分与不确定性汇总、补齐敏感性和观察机制核查。**第三部分目前只作待开展事项，没有正式分析或结论；进入第三部分不意味着第二部分已完成可靠机制选择。**

## Notebook 与图表索引

主线入口是 [matching_h0_h3_results.ipynb](matching_h0_h3_results.ipynb)，保存表格与内嵌图片，可直接在 GitHub 阅读。

| 图 | 内容 | 解释边界 |
|---|---|---|
| II-1 | 有限模拟目标参数恢复 | 两份/生成模型；不是稳定覆盖率 |
| II-2 | 四模型群体目标参数森林图 | 模型内估计；固定为零项不是估计结果 |
| II-3 | 绝对 PPC：准确率、RT、超时 | 样本内充分性，不是预测排名 |
| II-4 | Good-self−Bad-self 行为对比 PPC | 指标和分母同观测对齐 |
| II-5 | 快速反应形状 PPC | 观测固定分箱；95% 逐统计量预测区间 |

第一部分图表及历史 M4 诊断见 [复现 Notebook](../01_reproduction_2020/reproduction_results.ipynb) 与 [历史结果 Notebook](results.ipynb)。[目录、归档与命名说明](docs/ARCHIVE_AND_NAMING.md)明确 B/M/H 三组标识；原始数据、历史结果和采样文件原位保留。

## 复现方法与发布范围

轻量阅读复现：从仓库根目录执行，输出写入新的 `.preview/`，不覆盖已存研究结果。

```bash
python -m pip install -r requirements-notebooks.txt
python -m jupyter nbconvert --execute --to notebook --output-dir .preview 01_reproduction_2020/reproduction_results.ipynb
python -m jupyter nbconvert --execute --to notebook --output-dir .preview 02_reanalysis_2026/matching_h0_h3_results.ipynb
python -m jupyter nbconvert --execute --to notebook --output-dir .preview 02_reanalysis_2026/results.ipynb
```

这三个展示 Notebook 读取已保存 CSV/JSON/PNG，不运行 MCMC。历史行为/模拟/恢复 Notebook 保留代码和旧输出；其中模拟与恢复依赖 HDDM、SciPy、ArviZ、Linux fork 及未发布的大型中间样本，不能用展示环境宣称全流程已复现。

H0–H3 的 `Scripts/mechanism_*.py` 保留原始实现；部分按 `matching_observation_recovery.ipynb` 的固定单元编号读取模型。不要重排单元。采样/报告脚本可能创建或覆盖阶段产物，不能作为轻量重读命令；完整研究重跑前必须在隔离副本恢复相应 NPZ、模型数据库与原运行环境，再核对各脚本依赖。本次不重跑采样。

GitHub 仅发布白名单中的源代码、已有公开输入、紧凑汇总/状态和图；不发布完整后验、逐事件 loglik、数据库或本地日志。图表可以重读，后验计算不能仅凭此克隆完整再生。[证据 SHA-256 清单](docs/evidence_manifest.json)记录当前发布的结果来源；[整理验证记录](docs/VALIDATION.md)说明执行与保留范围。
