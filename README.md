# Good Me Bad Me · 正式实验复现与决策机制分析

Hu 等人（2020）发现，与 Bad-self 相比，人们对 Good-self 的匹配判断更快。第一部分复现原文的三项任务。第二部分聚焦 Matching，检验起点是否随形状条件改变。第三部分沿用第二部分 M2 的条件结构，分别检验证据跳跃和边界收缩能否改善对选择和反应时的预测。[原论文](https://doi.org/10.1525/collabra.301) · [作者数据与分析](https://github.com/hcp4715/moralSelf_ddm)

| 分析 | 模型与主要结果 | 查看结果 |
| --- | --- | --- |
| [第一部分 · 原研究复现](01_reproduction_2020/README.md) | 分别拟合 Matching、Identity、Valence。Good-self 减去 Bad-self 的漂移率后验均值依次为 0.88、0.51、0.53，三个 95% 可信区间均在零以上。 | [行为图与统计表](01_reproduction_2020/outputs/behavior/) · [Matching 漂移率图](01_reproduction_2020/outputs/matching/Figure_4A.pdf) |
| [第二部分 · 条件起点](02_reanalysis_2026/README.md) | M1 共享起点，M2 按形状估计四个起点。M2 的 DIC 较低，Good-self 的起点更靠近 Match。M1 的两项留出误差略低。 | [起点后验](02_reanalysis_2026/outputs/z_posterior.csv) · [拟合比较](02_reanalysis_2026/outputs/fit_model_comparison.csv) · [留出比较](02_reanalysis_2026/outputs/heldout_prediction.csv) |
| [第三部分 · 证据跳跃](03_alternative_mechanisms/scripts/03_levy_flight_bayesflow.ipynb) | 两模型均固定边界。与 Gaussian 相比，Lévy 的平均留出联合对数预测密度提高 0.217，反应时第 90 百分位（P90）的平均绝对误差减少 10.55 ms。个体 α 的恢复误差与先验均值基准接近。 | [联合密度图](03_alternative_mechanisms/outputs/levy/individual_window_recovery_optimization/heldout_rt_choice_density.pdf) · [参数恢复](03_alternative_mechanisms/outputs/levy/individual_window_recovery_optimization/parameter_recovery.csv) |
| [第三部分 · 边界收缩](03_alternative_mechanisms/scripts/03_collapsing_boundary_hssm.ipynb) | 两模型均使用高斯增量。与固定边界相比，收缩边界使留出 P90 误差从 83.28 降至 56.70 ms，条件 DIC 从 −8828.37 降至 −10506.60。 | [联合密度图](03_alternative_mechanisms/outputs/collapsing/m2_centered_ndt/heldout_rt_density_by_condition.pdf) · [留出比较](03_alternative_mechanisms/outputs/collapsing/m2_centered_ndt/heldout_model_comparison.csv) · [链轨迹](03_alternative_mechanisms/outputs/collapsing/m2_centered_ndt/trace_diagnostics.pdf) |

第三部分的两条路线各自与同一路线的对照模型比较。路线二的桥采样缺少分布重叠，且部分留出试次的预测密度为零，因此没有可报告的后验模型概率或有限的平均对数预测密度差。[两条路线的设定、诊断与复算说明](03_alternative_mechanisms/README.md)

## 数据、代码与汇报

- [正式实验原始数据](01_reproduction_2020/data/raw/)和[三个任务的模型输入](01_reproduction_2020/data/derived/)保存在第一部分。第二、三部分直接读取其中的 Matching 数据。
- 各阶段的 `scripts/` 保存分析入口，`outputs/` 保存相应的结果表和原始图。第三部分两张联合密度图均采用 10 ms 分箱。同一刺激条件下，两种按键的曲线共用分母，单条曲线的面积对应该按键的选择率。
- [三阶段建模与结果汇报（PPTX，8 页）](presentation/good-self-models-and-results.pptx)
- [已拟合模型下载](https://github.com/soupandpsy/hu2020-ddm/releases/tag/analysis-2026-11-28)包含第一部分三项任务的 HDDM 采样链和第三部分固定／收缩边界的后验抽样。下载所需 ZIP 并解压到仓库根目录，即可还原 Notebook 使用的输出路径。
