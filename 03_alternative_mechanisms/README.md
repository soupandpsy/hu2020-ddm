# 第三部分 · 计算模型拓展

## 路线一 · 固定边界下的 Gaussian 与 Lévy

[打开已执行的单被试 Notebook](scripts/03_levy_flight_bayesflow.ipynb)。使用 Hu 等人（2020）正式实验的 Matching 试次。两个模型共用八种刺激条件、个人作答窗口和分层留出划分。Gaussian 固定 α=2，Lévy 为每位被试单独估计 α。正负编码表示实际选择。反应时模型以窗口内的有效作答为条件。[Liss 等人（2026）的参考论文](https://doi.org/10.1016/j.cogpsych.2026.101836)

### 已保存的结果

- [独立模拟的参数恢复](outputs/levy/individual_window_recovery_optimization/parameter_recovery.csv)和[小档与中档的配对对照](outputs/levy/individual_window_recovery_optimization/parameter_recovery_budget_comparison.csv)。两档使用同一批模拟真值。增加训练量后，漂移差的恢复误差下降，个人 α 的恢复误差仍接近先验均值基准。
- [真人留出联合评分](outputs/levy/individual_window_recovery_optimization/heldout_joint_lpd.csv)和[慢反应及选择率误差](outputs/levy/individual_window_recovery_optimization/heldout_errors.csv)。4,929 次留出作答的平均联合评分为 Gaussian −0.136、Lévy 0.081。Lévy 总体得分较高，但逐人得分更高的被试为 19/42。对个体跳跃机制的解释仍受 α 恢复能力限制。
- [训练损失](outputs/levy/individual_window_recovery_optimization/training_loss.pdf)、[参数恢复图](outputs/levy/individual_window_recovery_optimization/parameter_recovery.pdf)、[慢反应与选择率](outputs/levy/individual_window_recovery_optimization/heldout_slow_rt_and_choice.pdf)和[选择与反应时联合密度](outputs/levy/individual_window_recovery_optimization/heldout_rt_choice_density.pdf)。

### 复算

从 `03_alternative_mechanisms/scripts/` 运行 Notebook。原始与派生 Matching 表已位于第一部分的 `data/`。训练环境使用 Python 3.12、BayesFlow 2.1.0、Keras 3.15.0 和 JAX 0.12.1，CPU 模拟使用 Numba 0.68.1。选择对应的 Jupyter Python 内核。

每批 32 名模拟被试，每轮 500 批，共 30 轮和 15,000 次更新。仓库提供固定验证集、训练记录和[第 30 轮检查点](outputs/levy/individual_window_recovery_optimization/checkpoints/epoch_30.weights.h5)。从这些文件载入即可复算参数恢复和留出预测。小档的两份基准 CSV 用于配对恢复比较。

### 数值步长

模拟器使用 `dt = 0.00003125` 秒。步长在固定边界模拟中按 Match 选择率与两侧 P90 核对，并对照半步长及 Gaussian 解析参照。两个候选和两档训练预算共用这个步长。[Liss 等人的公开模拟代码](https://github.com/jvliss/time_to_jump/tree/01ee42cf72736de2800b64c846d7d9fdb5e73f6b/3_simulators)
