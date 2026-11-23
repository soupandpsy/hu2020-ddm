# 第二部分 · 刺激编码再分析

使用 Hu 等人（2020）正式实验的 Matching 数据，检验起点偏差是否随刺激条件变化。

## Notebook 与结果

[Matching 条件起点 Notebook](scripts/02_static_v_dc.ipynb) 比较共享起点的 M1 与起点随形状条件变化的 M2。结果包括参数恢复、采样诊断、留出预测和拟合试次的 DIC 与 PPC-MSE。

[图表目录](outputs/) · [留出预测比较](outputs/heldout_prediction.csv) · [拟合试次模型比较](outputs/fit_model_comparison.csv)

## 运行

在 CPU DockerHDDM `hcp4715/hddm:1.0.1` 中使用 `python3` 内核，从 `scripts/` 运行 Notebook。输入是第一部分生成的 [Matching 试次表](../01_reproduction_2020/data/derived/matching_hddm.csv)。
