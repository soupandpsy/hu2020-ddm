# 第二部分 · 刺激编码再分析

Hu 等人（2020）正式实验的 Matching 任务中，好／坏 × 自我／他人这四种形状条件是否具有不同的决策起点？条件起点能否提高对同一批被试新试次的预测精度？

## 模型与数据

[Matching 条件起点 Notebook](scripts/02_static_v_dc.ipynb) 沿用第一部分 Matching 的刺激编码。实际选择 Match 编为 1，Mismatch 编为 0。两套模型均按实际匹配状态 × 效价 × 身份估计漂移率 `v` 和非决策时间 `t`。每位被试跨条件共享边界间距 `a`。模型之间只改变相对起点 `z` 的条件结构。

| 模型 | 同一被试的起点设定 |
| --- | --- |
| M1 | 四种形状条件共用一个 `z` |
| M2 | Good-self、Bad-self、Good-other、Bad-other 各有一个 `z` |

42 名被试按被试及试次条件分层抽取约 20% 的试次留出，两个模型共用 19,599 个拟合试次和 4,929 个留出试次。每个模型采样四条链，各 10,000 次，舍弃前 1,000 次。

## 结果

M2 中，Good-self 和 Bad-self 的组层起点均值分别为 0.5845 和 0.5390，二者的后验差值为 0.0455，95% 可信区间为 [0.0201, 0.0710]。Good-self 的起点更靠近 Match 边界。[四条件起点图](outputs/z_by_condition.pdf) · [自我条件差值图](outputs/z_self_contrast.pdf) · [后验表](outputs/z_posterior.csv)

| 指标（越低越好） | M1 | M2 |
| --- | ---: | ---: |
| 拟合试次 DIC | −2821.88 | −3040.34 |
| 拟合试次 PPC-MSE | 0.004987 | 0.005074 |
| 留出选择比例平均绝对误差（百分点） | 3.24 | 3.34 |
| 留出反应时分位数平均绝对误差（ms） | 87.68 | 89.04 |

DIC 支持 M2。就这次留出预测而言，M1 的两项误差均略低。DIC 和 PPC-MSE 使用各模型的第一条拟合链。PPC-MSE 在同一被试、同一试次条件内比较模拟与实测的 Match 选择率，再对模拟抽样及观测节点取平均。[拟合比较表](outputs/fit_model_comparison.csv) · [留出比较表](outputs/heldout_prediction.csv)

| 检查 | 图表 |
| --- | --- |
| 参数恢复与采样 | [模拟起点恢复](outputs/z_recovery.csv) · [链轨迹](outputs/trace_diagnostics.pdf) · [R̂](outputs/rhat.csv) |
| 漂移率与留出预测 | [漂移率后验](outputs/matching_drift_posterior.pdf) · [观测与预测反应时密度](outputs/heldout_rt_density.pdf) · [各条件预测区间](outputs/heldout_ppc.csv) · [选择比例](outputs/heldout_choice_proportions.csv) · [反应时分位数](outputs/heldout_rt_quantiles.csv) |

## 运行

从**仓库根目录**在 PowerShell 运行以下命令。容器会同时读到第一部分的 [Matching 试次表](../01_reproduction_2020/data/derived/matching_hddm.csv)，并将结果写入本部分的 [outputs/](outputs/)。

```powershell
docker run --rm --mount "type=bind,source=$($PWD.Path),target=/work" `
  --workdir /work/02_reanalysis_2026/scripts --env CUDA_VISIBLE_DEVICES=-1 `
  hcp4715/hddm:1.0.1 jupyter nbconvert --to notebook --execute --inplace `
  --ExecutePreprocessor.kernel_name=python3 --ExecutePreprocessor.timeout=-1 02_static_v_dc.ipynb
```
