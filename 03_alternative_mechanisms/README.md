# Part 3：VTZ 基线与条件起点效应

Matching 任务的层级扩散决策分析，比较条件性漂移率解释与起点效应。数据来自 Hu et al. (2020)，沿用仓库已有的被试筛选和反应编码。

| 分析入口 | 模型 | 四链 DIC | 后验诊断 |
|---|---|---:|---|
| [M0 · VTZ 基线](Scripts/notebooks/00_baseline_vtz.ipynb) | 八条件 v、t₀；被试内共享 a、z | −3552.29 | 通过 |
| [M1 · 条件起点效应](Scripts/notebooks/01_starting_point_z.ipynb) | M0 加 Self 内 Good−Bad 的 logit 起点效应 | −3602.35 | 通过 |

42 名被试、24,499 个 RT≥0.2 秒的有效响应中，固定划分的 19,580 个用于拟合、DIC 和 PPC。每个模型完成四条独立链，每链 10,000 次迭代、丢弃前 5,000 次，共保留 20,000 个样本。全部后验变量及目标对比达到 R-hat<1.01、ESS-bulk/ESS-tail≥400；MCSE 见完整诊断表。

DIC 差值 50.06 偏向 M1。M1 的平均起点差 Δz 为 0.042，95% HDI [0.025, 0.059]。Match 条件的 Self 内漂移差从 M0 的 0.948 变为 M1 的 0.754，两者的 95% HDI 均高于零。DIC 差值不是显著性检验或 Bayes factor。

每模型完成 2,000 个 PPC 复制数据集。两模型均高估慢反应尾部：八条件 RT 第 90 百分位平均高估约 0.150 秒。采样收敛和相对 DIC 优势不等于充分解释行为。当前使用有效响应子集的未截断 RT–choice 似然，尚未建模 RT 下限选择及 Deadline；恢复验证的范围也有限。

## 阅读与复现

两份 Notebook 已保存结果表和图，可直接在 GitHub 阅读。章节顺序为 Model Specification → Model Fitting → Model Diagnosis → Model Comparison & PPC → Statistical Inference。

重新运行所需环境见 [ENVIRONMENT.md](ENVIRONMENT.md)，输入说明见 [Data/README.md](Data/README.md)。所有计算开关默认关闭。Git 仅保存代码、输入数据、小型汇总和 Notebook 图表；原始链数据库、NetCDF 和逐次预测记录留在本地 `local_bulk/`，未随本分支发布。新克隆中重新执行 Notebook 会提示缺少本地后验；恢复原链或显式重新拟合后，才能重算轨迹、DIC 和 PPC。直接阅读已保存输出不需要原始链。

## 结果文件

- [DIC 比较](Output/tables/M0_M1_20261111_dic.csv)
- M0：[全参数诊断](Output/tables/M0_20261111_diagnostics.csv)、[目标参数](Output/tables/M0_20261111_posterior_summary.csv)、[分条件 PPC](Output/tables/M0_20261111_ppc_summary.csv)、[来源与指纹](Output/tables/M0_20261111_provenance.json)
- M1：[全参数诊断](Output/01_starting_point_z/M1_20261111_diagnostics.csv)、[目标参数](Output/01_starting_point_z/M1_20261111_posterior_summary.csv)、[分条件 PPC](Output/01_starting_point_z/M1_20261111_ppc_summary.csv)、[来源与指纹](Output/01_starting_point_z/M1_20261111_provenance.json)
- [实现与有限恢复验证](Output/tables/M0_M1_validation_20261111.md)

结果日期为 2026-11-11。来源文件中的 `mcmc_run: false` 表示该次汇总没有重新采样；四链数量及后验指纹记录的是已完成的采样。

本分支发布 M0、M1。M2、M4 尚无完整正式后验，M3 保持暂停，均不纳入此处的模型排名。

## 参考

- [Hu et al. (2020)：Good Me Bad Me](https://doi.org/10.1525/collabra.301)
- [原作者分析仓库](https://github.com/hcp4715/moralSelf_ddm)
- [dockerHDDM 官方工作流](https://github.com/hcp4715/dockerHDDM/blob/master/dockerHDDM_Workflow.ipynb)
