# 第三部分 · 计算模型拓展

两条路线都问，同一批 Matching 数据中，第二部分 M2 对选择和慢反应的预测还能否改善？两者沿用 M2 的条件结构。漂移率 `v` 和非决策时间 `t` 按八种试次条件估计，相对起点 `z` 按四种形状估计，每人有一个跨条件的初始边界间距。两条路线各自拟合一个对照模型，在本路线内比较。[第二部分 M2](../02_reanalysis_2026/scripts/02_static_v_dc.ipynb)

两个 Notebook 均从 `03_alternative_mechanisms/scripts/` 运行。输入来自第一部分的 [Matching 建模表](../01_reproduction_2020/data/derived/matching_hddm.csv)和[原始 Matching 记录](../01_reproduction_2020/data/raw/MS_rep_matchingTask_raw.csv)。42 名被试共有 24,528 个有效作答。沿用第二部分的分层划分（种子 2026），19,599 个试次用于拟合，4,929 个试次用于留出评价。个人超时记录用于估计作答截止时间。两条路线的密度图均以 10 ms 分箱，分别展示观测与预测的选择和反应时联合分布。

## 路线一 · 固定边界下的 Gaussian 与 Lévy

两模型采用相同的固定边界和条件起点，改变证据每一步的随机增量分布。Gaussian 固定稳定性参数 α = 2，Lévy 为每名被试估计一个跨条件共享的 α，允许证据向任一方向出现较大跳跃。BayesFlow 在模拟数据上训练网络，再对真实被试逐人估计参数和模型概率。[已执行 Notebook](scripts/03_levy_flight_bayesflow.ipynb) · [Liss 等人（2026）](https://doi.org/10.1016/j.cogpsych.2026.101836)

| 留出指标 | Gaussian | Lévy |
| --- | ---: | ---: |
| 每试次平均联合对数预测密度 | −0.136 | 0.081 |
| 反应时第 90 百分位数平均绝对误差（ms） | 96.28 | 85.73 |
| Match 选择率平均绝对误差（百分点） | 6.76 | 4.24 |

Lévy 的总体留出评分更高，逐人比较时有 19/42 人的评分更高。模拟数据中，个体 α 的恢复误差为 0.1990，先验均值基准为 0.1987。128 个由 Gaussian 生成的模拟个体中有 63 个被判为 Lévy。总体预测改善与个体参数识别可分别从以下结果核查。

- **预测**　[总体联合评分](outputs/levy/individual_window_recovery_optimization/heldout_joint_lpd.csv) · [逐人评分](outputs/levy/individual_window_recovery_optimization/heldout_joint_lpd_by_subject.csv) · [预测误差](outputs/levy/individual_window_recovery_optimization/heldout_errors.csv) · [选择与反应时联合密度图](outputs/levy/individual_window_recovery_optimization/heldout_rt_choice_density.pdf) · [慢反应与选择率图](outputs/levy/individual_window_recovery_optimization/heldout_slow_rt_and_choice.pdf)
- **恢复与训练**　[模型恢复](outputs/levy/individual_window_recovery_optimization/model_recovery.csv) · [参数恢复](outputs/levy/individual_window_recovery_optimization/parameter_recovery.csv) · [参数恢复图](outputs/levy/individual_window_recovery_optimization/parameter_recovery.pdf) · [两档训练预算的配对比较](outputs/levy/individual_window_recovery_optimization/parameter_recovery_budget_comparison.csv) · [训练损失图](outputs/levy/individual_window_recovery_optimization/training_loss.pdf)

**复算。** 正式环境为 Python 3.12、BayesFlow 2.1.0、Keras 3.15.0 和 JAX 0.12.1，CPU 模拟使用 Numba 0.68.1。每批 32 名模拟被试，每轮 500 批，30 轮共更新 15,000 次。依次运行 Notebook 会读取固定验证集、训练记录和[第 30 轮检查点](outputs/levy/individual_window_recovery_optimization/checkpoints/epoch_30.weights.h5)，然后复算参数恢复和留出预测。小档预算的两份基准 CSV 位于 [individual_window_flow_check](outputs/levy/individual_window_flow_check/)。

### 数值步长

模拟步长 `dt = 0.00003125` 秒。两个候选模型和两档训练预算共用这一设置。步长核对使用 Match 选择率与两侧反应时第 90 百分位数，对照半步长模拟和 Gaussian 解析参照。[Liss 等人的公开模拟代码](https://github.com/jvliss/time_to_jump/tree/01ee42cf72736de2800b64c846d7d9fdb5e73f6b/3_simulators)

## 路线二 · 高斯积累下的固定边界与收缩边界

两模型共用高斯积累和个人层级参数。固定模型的上下边界为 ±a，收缩模型为 ±[a − s·tan(θ)]，其中 `a` 为初始半间距，`s` 为扣除非决策时间后的积累时间。固定模型取 θ = 0，收缩模型为全组估计一个 θ，先验为 0 至 1.2 弧度的均匀分布。HSSM 接收 JAX 实现的首达时间似然，按个人有效作答窗口归一化。[已执行 Notebook](scripts/03_collapsing_boundary_hssm.ipynb)

| 指标 | 固定边界 | 收缩边界 |
| --- | ---: | ---: |
| 拟合试次条件 DIC | −8,828.37 | −10,506.60 |
| 留出反应时第 90 百分位数平均绝对误差（ms） | 83.28 | 56.70 |
| 留出 Match 选择率平均绝对误差（百分点） | 2.20 | 1.73 |

收缩边界模型在这三项指标上均更低。两种模型均无发散，最大 R̂ 分别为 1.0056 和 1.0096。桥采样的对数密度比范围不重叠，未得到可报告的后验模型概率。固定和收缩模型分别有 7 个和 3 个留出试次的预测密度为零，两者的平均留出对数预测密度均为 −∞，无法用该指标判断优劣。

- **比较与预测**　[条件 DIC](outputs/collapsing/m2_centered_ndt/conditional_dic.csv) · [留出预测误差](outputs/collapsing/m2_centered_ndt/heldout_model_comparison.csv) · [桥采样状态](outputs/collapsing/m2_centered_ndt/posterior_model_probabilities.csv) · [逐模型留出对数密度](outputs/collapsing/m2_centered_ndt/heldout_log_density.csv) · [选择与反应时联合密度图](outputs/collapsing/m2_centered_ndt/heldout_rt_density_by_condition.pdf) · [慢反应与选择率图](outputs/collapsing/m2_centered_ndt/heldout_slow_rt_and_rates.pdf)
- **参数与诊断**　[采样汇总](outputs/collapsing/m2_centered_ndt/sampling_summary.csv) · [链轨迹图](outputs/collapsing/m2_centered_ndt/trace_diagnostics.pdf) · [参数区间](outputs/collapsing/m2_centered_ndt/parameter_intervals.pdf) · [分链后验表](outputs/collapsing/m2_centered_ndt/parameter_posterior_by_chain.csv)

**复算。** 在 Python 3.12 的 HSSM、JAX 和 NumPyro 环境中运行。每个模型用 4 条链，每链预热 1,500 次，保留 1,000 次抽样，目标接受率为 0.90。完整运行会将两个模型的后验分别写入 `outputs/collapsing/m2_centered_ndt/fixed_fit.nc` 和 `collapsing_fit.nc`。如已有这两份文件，可先执行采样前的定义单元，跳过采样单元，接着从读取两份后验的单元继续复算诊断、模型比较及绘图。
