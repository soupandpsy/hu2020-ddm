# Hu et al.（2020）：计算复现与 HDDM 再分析

[English](README.md) · 中文

本仓库包含《Good me, bad me: Prioritization of the good-self during perceptual decision-making》的计算复现和探索性再分析。两部分的研究目的及模型设定不同，分别报告。

| 部分 | 内容 | 主要结果与限制 |
|---|---|---|
| [第一部分：2020 研究复现](01_reproduction_2020/README_zh.md) | Pilot 与 Confirm 行为统计、Pilot ex-Gaussian 拟合、Confirm HDDM | 复现了部分选定统计量和模型分析；文档列明分析设置差异及尚未完成的补充分析。 |
| [第二部分：2026 方法再分析](02_reanalysis_2026/README_zh.md) | 采用刺激编码 HDDM，区分起点偏置与漂移偏置 | Matching M4 达到主分析采样诊断标准，但预测拟合与参数恢复仍有限制；Identity 和 Valence 因未收敛，不报告推断性参数结果。 |

Matching 再分析中，Self-Good 与 Self-Bad 漂移幅度之差的后验均值为 0.666，95% 可信区间为 [0.387, 0.942]。该结果是给定 M4 的探索性估计。由于模型约束不同，它不构成对 2020 模型参数的直接数值复现，也不作因果解释。

两部分均按 `Data`、`Scripts`、`Output` 组织，并分别提供中英文 README，说明方法、结果、限制和软件要求。第二部分提供汇总表、诊断证据及结果图，不分发完整模型数据库和后验样本；重新计算依赖后验的结果需要相应拟合输出。具体复现范围见各部分说明。

数据和原始分析材料来自[作者仓库](https://github.com/hcp4715/moralSelf_ddm)。再分析采用[方法仓库](https://github.com/Chuan-Peng-Lab/ddm_stim_coding_posts)中的刺激编码框架。本仓库为独立复现与再分析，不是原研究的官方更新。
