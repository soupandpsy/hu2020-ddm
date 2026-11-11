# M0 / M1：实现核验与有限恢复验证

本页整理 2026-11-07 至 2026-11-11 已保存的 M0/M1 证据。历史模拟和当前实证分析的数据版本、预算及评分方法分别记录，不混用结果。

## 模型实现

M0 调用 HDDMStimCoding 原生 Wiener 首达时间密度。M1 仅增加 Self 内 Good−Bad 的 logit 起点效应，采用非中心化层级。三组固定参数的中心化／非中心化观测 log-likelihood 差与先验 Jacobian 对照差均为零；这是指定参数点上的实现核验，不是全参数域证明。

完整输入上的 778 个共同随机节点先验一致，M1 另有 44 个随机节点（群体均值、群体标准差及 42 个标准正态个体变量）。保存轨迹的变换复核使用当时两条各 4,000 个保留样本的模拟链，最大重构误差小于 1.2e−16；该记录不是当前实证四链。

- [非中心化变换核验](comparison_M1_noncentered_equivalence.json)
- [共同先验核验](comparison_M1_nc_full_common_prior_check.json)
- [历史轨迹映射核验](comparison_M1_nc_saved_trace_mapping_check.json)
- [实际先验节点](comparison_actual_prior_nodes.json)

历史先验核验文件中的“完整四链诊断待完成”是生成该文件时的状态；当前实证诊断另见下方。

## 拟合前模拟

每种情形 42 名模拟参与者，训练集与独立测试集各 24,527 条，每情形仅一份数据实现。零效应情形的个体起点效应全为零；弱效应与较强效应的群体均值分别为 0.15、0.6，个体差异 SD 均为 0.2，均在 logit 起点尺度。这些数值是项目的模拟设计，不是 HDDM 官方规定的强度。

三个场景的已保存链通过采样诊断。较强效应场景的目标恢复与预测评分稳定性通过当时检查；零／弱情形的评分积分精度不足，弱效应情形的 Mismatch 漂移对比区间未覆盖真值。较强效应场景的 ΔELPD 为 216.364，按被试计算的 SE 为 27.017。单次模拟不能提供可靠的重复覆盖率或模型识别频率。

这些 ELPD 和两半样本积分检查属于历史恢复实验；当前实证比较采用 DIC＋PPC，不以此重新添加实证运行门槛。

- [生成设计与来源](comparison_recovery_provenance.json)
- [三个模拟场景](comparison_recovery_scenarios.csv)
- [目标参数恢复](comparison_m1_all_recovery_targets.csv)
- [模型区分记录](comparison_model_recovery_decisions.csv)
- [零／弱效应后续评分](pre_empirical_m1_rescore_20261108.json)

## 当前数据与实证后验

完整输入预检覆盖 19,580 条训练试次、42 人、336 个被试×条件单元；M0 和 M1 构建及有限初始化均通过。[输入预检记录](M0_M1_input_checks_20261108.json)。预检中的 `mcmc_draws: 0` 指该次构建检查，不描述后续实证采样。

当前两个模型各有四条完整实证链，每链保留 5,000 个样本。M0 的 828 项后验变量与派生对比、M1 的 957 项均通过 R-hat<1.01 和 ESS-bulk/ESS-tail≥400。MCSE、后验指纹、DIC 及 PPC 使用的输入版本均已保存。

| 模型 | 最大 R-hat | 最小 ESS-bulk | 最小 ESS-tail | DIC | PPC 复制数 |
|---|---:|---:|---:|---:|---:|
| M0 | 1.00138 | 3939 | 6724 | −3552.29 | 2000 |
| M1 | 1.00349 | 2050 | 2451 | −3602.35 | 2000 |

[M0 诊断](M0_20261111_diagnostics.csv) · [M1 诊断](../01_starting_point_z/M1_20261111_diagnostics.csv) · [DIC](M0_M1_20261111_dic.csv)

PPC 显示两模型均高估慢反应尾部，八条件 RT 第 90 百分位平均高估约 0.150 秒。当前未截断似然尚未处理 RT 下限筛选和 Deadline。诊断通过支持后验样本可靠性，不代表观察模型完整、恢复验证充分或心理机制已获因果证明。
