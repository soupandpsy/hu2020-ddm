# 第一部分 · 正式实验复现

使用 R 复现 Figure 3A–C，使用 DockerHDDM 复现 Figure 4A–C。

| 分析入口 | 图及结果 |
| --- | --- |
| [01_behavior.R](scripts/01_behavior.R) | [行为图及统计表](outputs/behavior) |
| [02_matching_hddm.ipynb](scripts/02_matching_hddm.ipynb) | [Figure 4A](outputs/matching/Figure_4A.pdf) · [后验对比](outputs/matching/posterior_contrasts.csv) |
| [03_identity_hddm.ipynb](scripts/03_identity_hddm.ipynb) | [Figure 4B](outputs/identity/Figure_4B.pdf) · [后验对比](outputs/identity/posterior_contrasts.csv) |
| [04_valence_hddm.ipynb](scripts/04_valence_hddm.ipynb) | [Figure 4C](outputs/valence/Figure_4C.pdf) · [后验对比](outputs/valence/posterior_contrasts.csv) |

正式实验的原始数据在 `data/raw/`，三个 HDDM 输入在 `data/derived/`，图表和诊断在 `outputs/`。反应时在模型输入中以秒计，Figure 3 以毫秒计。原始 `ACC` 中，1 为正确，0 为错误，−1 为无响应，2 为无效键。

## 运行

在本部分目录运行 `Rscript scripts/01_behavior.R`。脚本使用 `dplyr`、`tidyr`、`ggplot2`、`tibble`、`afex` 和 `BayesFactor`。本次行为分析使用 R 4.6.0。

HDDM 使用 CPU 镜像 `hcp4715/hddm:1.0.1` 的 `python3` 内核。实际镜像包含 Python 3.9.13，HDDM 报告版本 `1.0.1RC`。从本部分目录在 PowerShell 运行三份 Notebook。

```powershell
docker run --rm --mount "type=bind,source=$($PWD.Path),target=/work" `
  --workdir /work/scripts --env CUDA_VISIBLE_DEVICES=-1 `
  hcp4715/hddm:1.0.1 jupyter nbconvert --to notebook --execute --inplace `
  --ExecutePreprocessor.timeout=-1 02_matching_hddm.ipynb 03_identity_hddm.ipynb 04_valence_hddm.ipynb
```

每项模型使用四条独立链，每链采样 10,000 次并丢弃前 1,000 次。三项刺激编码模型均估计 `v`、`a`、`t` 和 `z`。`p_outlier=0` 取自作者的三个归档最终模型。运行会在各任务输出目录生成完整采样链和模型文件。

[论文](https://doi.org/10.1525/collabra.301) · [原作者正式实验分析](https://github.com/hcp4715/moralSelf_ddm/tree/master/2_confirm_study/Results) · [DockerHDDM 镜像](https://hub.docker.com/r/hcp4715/hddm)
