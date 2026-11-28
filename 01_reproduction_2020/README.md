# 第一部分 · 正式实验复现

Hu 等人（2020）的正式实验中，Good-self 的行为优势和漂移率差异能否复现？R 脚本复现 Figure 3A–C，三份 DockerHDDM Notebook 分别拟合 Matching、Identity 和 Valence，复现 Figure 4A–C。

| 分析入口 | 图及结果 |
| --- | --- |
| [01_behavior.R](scripts/01_behavior.R) | [Figure 3A](outputs/behavior/Figure_3A.pdf) · [3B](outputs/behavior/Figure_3B.pdf) · [3C](outputs/behavior/Figure_3C.pdf) · [ANOVA](outputs/behavior/anova.csv) · [Bayes factor](outputs/behavior/bayes_models.csv) · [计划对比](outputs/behavior/planned_contrasts.csv) |
| [02_matching_hddm.ipynb](scripts/02_matching_hddm.ipynb) | [Figure 4A](outputs/matching/Figure_4A.pdf) · [后验对比](outputs/matching/posterior_contrasts.csv) |
| [03_identity_hddm.ipynb](scripts/03_identity_hddm.ipynb) | [Figure 4B](outputs/identity/Figure_4B.pdf) · [后验对比](outputs/identity/posterior_contrasts.csv) |
| [04_valence_hddm.ipynb](scripts/04_valence_hddm.ipynb) | [Figure 4C](outputs/valence/Figure_4C.pdf) · [后验对比](outputs/valence/posterior_contrasts.csv) |

正式实验的两份[原始数据](data/raw/)经 R 脚本处理后生成三份[模型输入](data/derived/)。模型输入中的反应时以秒计，Figure 3 中的反应时以毫秒计。原始 `ACC` 中，1 为正确，0 为错误，−1 为无响应，2 为无效键。

## 模型与复现结果

三个任务分别拟合刺激编码模型。上边界在 Matching、Identity 和 Valence 中分别代表 Match、Self 和 Good。`split_param="v"` 根据刺激翻转漂移率方向，`drift_criterion=False` 关闭附加漂移偏置。漂移率 `v` 和非决策时间 `t` 在 Matching 中按匹配状态 × 效价 × 身份估计，在两个分类任务中按效价 × 身份估计。各任务中，同一被试跨条件共享边界间距 `a` 和相对起点 `z`。三项最终模型的 `p_outlier=0` 已与作者归档核对。[作者的 Matching Notebook](https://github.com/hcp4715/moralSelf_ddm/blob/master/2_confirm_study/Results/3_hddm/HDDM_MS_Rep_match_stim.ipynb) · [作者的分类 Notebook](https://github.com/hcp4715/moralSelf_ddm/blob/master/2_confirm_study/Results/3_hddm/HDDM_MS_Rep_categ_stim.ipynb)

Matching 中，Good-self 比 Bad-self 平均快 82.56 ms，辨别力 d′ 高 0.526。下表为三个独立模型中 Good-self 减去 Bad-self 的漂移率差值。[行为计划对比](outputs/behavior/planned_contrasts.csv)

| 任务 | 后验均值差 | 95% 可信区间 | 链诊断与后验预测检查 |
| --- | ---: | --- | --- |
| Matching | 0.88 | [0.57, 1.19] | [轨迹](outputs/matching/trace_diagnostics.pdf) · [R̂](outputs/matching/rhat.csv) · [PPC](outputs/matching/ppc.csv) |
| Identity | 0.51 | [0.14, 0.88] | [轨迹](outputs/identity/trace_diagnostics.pdf) · [R̂](outputs/identity/rhat.csv) · [PPC](outputs/identity/ppc.csv) |
| Valence | 0.53 | [0.12, 0.93] | [轨迹](outputs/valence/trace_diagnostics.pdf) · [R̂](outputs/valence/rhat.csv) · [PPC](outputs/valence/ppc.csv) |

## 运行

在本部分目录运行 `Rscript scripts/01_behavior.R`。脚本使用 `dplyr`、`tidyr`、`ggplot2`、`tibble`、`afex` 和 `BayesFactor`。本次行为分析使用 R 4.6.0。

HDDM 使用 CPU 镜像 `hcp4715/hddm:1.0.1` 的 `python3` 内核。实际镜像包含 Python 3.9.13，HDDM 报告版本 `1.0.1RC`。从本部分目录在 PowerShell 运行三份 Notebook。

```powershell
docker run --rm --mount "type=bind,source=$($PWD.Path),target=/work" `
  --workdir /work/scripts --env CUDA_VISIBLE_DEVICES=-1 `
  hcp4715/hddm:1.0.1 jupyter nbconvert --to notebook --execute --inplace `
  --ExecutePreprocessor.timeout=-1 02_matching_hddm.ipynb 03_identity_hddm.ipynb 04_valence_hddm.ipynb
```

每项模型使用四条独立链，每链采样 10,000 次并丢弃前 1,000 次。运行会在各任务输出目录生成完整采样链和模型文件。

[论文](https://doi.org/10.1525/collabra.301) · [原作者正式实验分析](https://github.com/hcp4715/moralSelf_ddm/tree/master/2_confirm_study/Results) · [DockerHDDM 镜像](https://hub.docker.com/r/hcp4715/hddm)
