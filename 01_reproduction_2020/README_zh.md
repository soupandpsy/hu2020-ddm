# 第一部分：Hu et al.（2020）复现

> **阅读导航 / Reading path (2026-10-28):** [复现结果 Notebook](reproduction_results.ipynb) → [历史 M0–M4 探索](reproduction_results.ipynb) → [第二部分 H0–H3](../02_reanalysis_2026/README_zh.md)。历史结果已合入第一部分总 Notebook，不与 H0–H3 结果混用；五折预测比较尚未可靠完成，第三部分尚未正式分析。

## 逻辑链与图表导航 / Research sequence

研究问题是先确认 Good-self 行为现象能否重现；方法沿用 Pilot/Confirm 筛选、被试内统计和历史 HDDM。核查发现选定 t 统计吻合，但 ANOVA 设置、ex-Gaussian 与部分模型诊断有差异。因此先分别报告可复现的现象和模型限制，再由 M0–M4 的预测/恢复问题引出第二部分的观察模型与 H0–H3 比较，不用行为复现直接证明机制。

[复现 Notebook](reproduction_results.ipynb)依次呈现 I-1 Pilot Matching、I-2 Pilot Categorization、I-3 Confirm Matching、I-4/5 Confirm Identity/Valence、I-6 历史漂移参数图及 I-7–I-10 的 M0–M4 路径示意、参数对比、轨迹与PPC。I 系列编号仅作本仓库导航；原论文 Figure 2/3/4 编号保持不变。各图保留区间和样本口径说明。


[English](README.md) · 中文 · [项目首页](../README_zh.md)

## 研究范围

本部分复现《Good me, bad me: Prioritization of the good-self during perceptual decision-making》的主要计算分析，涵盖 Pilot 与 Confirm 行为统计、Pilot ex-Gaussian 和 Confirm HDDM，并单独保存敏感性分析。结果总体上复现了部分选定统计量和模型分析，但并非所有统计量、补充分析或模型设置均与论文完全一致。

原始材料见[作者仓库](https://github.com/hcp4715/moralSelf_ddm)。

## 复现结果

- Pilot 保留 29 名被试；Confirm Matching 保留 42 名，Categorization 保留 41 名。R 导出的 Confirm 输入与仓库内输入一致。15 项选定的 t 统计量在指定容差内一致（Pilot < .01；Confirm < .02），不代表严格符合论文打印值的舍入结果。
- Matching/Identity/Valence vtz p0 和 Identity vz p0 通过了配置中的四次运行诊断标准。Matching/Valence vz 只有一次运行，仅用于描述。**Identity vz p0.05 未通过诊断，不应合并其各次样本。**
- Figure 3C 将作者误用的 Identity 汇总改为 Valence 数据；误差线仍为均值 ± 1.96 SE。Pilot 个体线在 Subject × Identity 内连接 Good/Bad，Confirm 则连接每名被试的四个条件。
## 差异与局限

- Pilot Categorization RT 与作者 JASP 的条件列分配不同，对应七项 F/omega 不能视为等价比较。其余 33 项可比 ANOVA F 值在 1e-9 内一致。Omega 定义、包含贝叶斯因子及 Pilot Categorization Accuracy 的干扰因素设置也存在差异，因此不能声称所有统计量与论文相同。
- ex-Gaussian 使用 1,000 次 bootstrap，而论文报告 2,000 次。一个单元（7014，Val/Other/Bad，仅 3 个试次）未拟合；231 个完成拟合的单元筛选后得到 27 名符合条件的被试。与可取得的作者估计相比，mu/sigma/tau 最大差异约为 .0283/.0567/.0208 秒；不声称 ANOVA 等价。详见[各单元结果](Output/behavior/pilot/exgaussian/parameters_all.csv)及保留的警告记录。
- 补充 ex-Gaussian 贝叶斯与简单效应分析、Figure S3、完整 Table S2/S3 汇总和贝叶斯因子稳健性图尚未全部完成。论文未报告 Pilot DDM，本项目也未实现该部分。
- 本部分记录 2020 复现工作。上述模型通过的是当时配置中的诊断标准，不能自动等同于 2026 再分析采用的秩归一化 R-hat 和 ESS 标准。2026 再分析及其独立结论见[第二部分](../02_reanalysis_2026/README_zh.md)。

## 文件与复现方法

### 文件结构与软件环境

|目录|内容|
|---|---|
|`Scripts/pilot_behavior/`、`Scripts/confirm_behavior/`|R 分析和独立的行为绘图入口|
|`Scripts/hddm/`|Matching、Identity、Valence 三个任务；每个模型有 Fit、Run 和 JSON 配置，任务级 Plot/PPC 独立|
|`Data/`|原始 CSV、三份 Confirm 模型输入和精简历史后验样本|
|`Output/`|行为统计、模型汇总与诊断、结果图|
|`Output/sensitivity/`|Identity vz 的 p_outlier=0/.05 扩展分析|

`vtz` 允许 v/t 随条件变化并估计 z；`vz` 的 t 在条件间共享。`p0`/`p0p05` 表示 p_outlier=0/.05；`r01` 为运行编号。拟合和绘图通过 JSON 文件读取模型配置。

使用版本：**R 4.6.0；Python 3.7.12；HDDM 0.8.0；NumPy 1.19.5；SciPy 1.5.3；Matplotlib 3.3.2**。R 包包括 Rcpp、dplyr、tidyr、ggplot2、afex、BayesFactor、effectsize 和 retimes（仅 Pilot 使用，已被 CRAN 归档）。具体包版本见各研究结果中的 `dependencies.csv`、`sessionInfo.txt`。运行前需安装依赖；仓库不包含本地包库和完整模型数据库。

以下命令均在本目录执行。大写目录占位符需替换为实际输出路径；在 WSL 中使用 Linux 路径。

### 查看与重画已有结果

本节 Figure 2、3、4 及补充图沿用原论文图号。

```bash
Rscript --vanilla Scripts/pilot_behavior/pilot_plot.R --analysis-dir Output/behavior/pilot
Rscript --vanilla Scripts/confirm_behavior/confirm_plot.R --analysis-dir Output/behavior/confirm
python -B Scripts/hddm/figure04_plot.py
```

上述命令读取保存的表格或样本，不重新拟合。Figure 2 对应 Pilot 行为结果；Figure 3 对应 Confirm 行为结果；[Figure 4](Output/figures/figure04_vtz_p0_r01_drift.pdf) 对应 Confirm vtz p0 r01 的漂移率。Matching 的 mismatch 漂移率对应 Figure S6。Figure 4 使用随仓库提供的 35 个群体参数样本（每个 9,000 次抽样），不是多链合并结果；这些精简样本不能替代完整数据库用于诊断或 PPC。

### 从数据重新分析

**1. 处理原始数据。** 每个入口会显示新建的分析目录；Pilot 包含 ex-Gaussian 拟合。

```bash
Rscript --vanilla Scripts/pilot_behavior/pilot_analysis.R
Rscript --vanilla Scripts/confirm_behavior/confirm_analysis.R
```

绘制本次行为结果时，使用上面的 Plot 命令，将分析目录换成本次输出。R 入口支持可选的 `--data-dir`、`--output-dir`，不覆盖非空输出目录。筛选规则、反应编码和单位见 [Data/README.md](Data/README.md)。

**2. 拟合三个 Confirm 模型。** `CONFIRM_DIR` 为已完成的 Confirm R 分析目录；代码自动选择相应的导出表。

```bash
python -B Scripts/hddm/matching/vtz/matching_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
python -B Scripts/hddm/identity/vtz/identity_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
python -B Scripts/hddm/valence/vtz/valence_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
```

省略 `--behavior-dir` 时读取仓库内的 Confirm 输入。仅运行一次时使用 `--runs r01`。拟合样本和分析结果分别存入配置指定模型目录下的 `runs/`、`analysis_new_*`；已有运行目录不会被覆盖或自动续跑。单独分析已完成的新拟合使用 `--mode analyze --source new`。历史分析（`--source existing`）需要外部完整数据库；在其他电脑上使用时需更新配置路径。

**3. 检查诊断并绘制新结果。** 填入上一步显示的各模型分析目录。

```bash
python -B Scripts/hddm/matching/matching_plot.py --analysis-dir MATCHING_ANALYSIS_DIR
python -B Scripts/hddm/identity/identity_plot.py --analysis-dir IDENTITY_ANALYSIS_DIR
python -B Scripts/hddm/valence/valence_plot.py --analysis-dir VALENCE_ANALYSIS_DIR
python -B Scripts/hddm/figure04_plot.py --source new --matching-analysis MATCHING_ANALYSIS_DIR --identity-analysis IDENTITY_ANALYSIS_DIR --valence-analysis VALENCE_ANALYSIS_DIR
```

新 Figure 4 读取各分析对应的 r01 完整样本。默认 `--source bundled` 读取仓库内的历史精简样本；`--source existing` 读取配置指定的历史数据库。重绘会新建输出目录。解释结果前请检查 `diagnostics.csv` 和 `status.json`。

**可选 PPC（后验预测检查）**，需要完整的已拟合模型：

```bash
python -B Scripts/hddm/matching/matching_ppc.py --config Scripts/hddm/matching/vtz/matching_vtz_p0_config.json --source new --run r01 --samples 500
```

Identity 和 Valence 有对应任务入口。省略 `--samples` 时只检查文件引用，不进行模拟。

**Identity vz 敏感性分析：** 必须明确选择同一 `vz/` 目录中的 `identity_vz_p0_config.json` 或 `identity_vz_p0p05_config.json`，没有默认设置。

```bash
python -B Scripts/hddm/identity/vz/identity_vz_run.py --config Scripts/hddm/identity/vz/identity_vz_p0_config.json --mode all --runs r01 r02 r03 r04
```

### 目录组织

本项目按 [TIER Protocol 4.0](https://www.projecttier.org/tier-protocol/protocol-4-0/) 做最小化整理：`Data` 保存输入，`Scripts` 保存分析代码，`Output` 保存结果，其中敏感性分析位于 `Output/sensitivity`。运行顺序见上文命令。目前未提供完整的 TIER 数据附录及全流程主脚本，因此不宣称完全符合 TIER 全部要求。

## 参考资料

- [Hu et al. (2020), *Good me, bad me: Prioritization of the good-self during perceptual decision-making*](https://doi.org/10.1525/collabra.301)
- [原始数据与分析代码](https://github.com/hcp4715/moralSelf_ddm)
- [TIER Protocol 4.0](https://www.projecttier.org/tier-protocol/protocol-4-0/)

## 历史 M0–M4：方法与复现 / Historical methods

`historical` 表示历史探索，M0–M4 是五个候选模型；并非两个新的分析部分。结果已并入本部分唯一 Notebook。以下历史命令均从 `02_reanalysis_2026/` 运行，勿与本部分原论文 `vtz/vz` 脚本混用。全拟合昂贵，阅读结果只需运行展示 Notebook。

### Historical reproduction commands

| Directory or file        | Contents                                                                                 |
| ------------------------ | ---------------------------------------------------------------------------------------- |
| `Data/`                  | Model inputs for Matching, Identity, and Valence                                         |
| `Scripts/hddm/common.py` | Shared HDDM model construction, fitting, and analysis functions                          |
| `Scripts/hddm/matching/` | Matching model fitting, analysis, plotting, parameter recovery, and sensitivity analysis |
| `Scripts/hddm/identity/` | Identity model fitting, analysis, and plotting                                           |
| `Scripts/hddm/valence/`  | Valence model fitting, analysis, and plotting                                            |
| `Output/`                | Model comparisons, parameter estimates, diagnostics, and figures                         |
| `results.ipynb`          | Summary and presentation of saved analysis results                                       |

Public outputs and their sources are documented in [output_manifest.json](../02_reanalysis_2026/Scripts/hddm/output_manifest.json).

Fitting uses dockerHDDM (HDDM 1.0.1RC). The container image and sampling settings are specified in the task configuration.

From the `02_reanalysis_2026/` directory, inspect the reproduction plan:

```bash
python Scripts/hddm/reproduce.py plan --output Output_new
```

To run the full fitting pipeline in the configured HDDM environment:

```bash
python Scripts/hddm/reproduce.py run --output Output_new --workers 2
```

A complete run requires the adjacent `01_reproduction_2020/Data` directory and a new output directory that does not already exist. Tasks run sequentially, and `--workers` supports 1–4 concurrent jobs, defaulting to 1. Lower concurrency may be necessary on memory-constrained systems.

To verify the existing public output files:

```bash
python Scripts/hddm/reproduce.py verify --output Output
```

This command checks file integrity and basic structure without recomputing parameter estimates or assessing convergence.

The current complete fitting pipeline has not yet been rerun end-to-end after code reorganization. The repository provides model inputs, analysis code, and summarized results, but does not distribute full MCMC databases or draw-level posterior samples.

To rebuild Figure 1 from public summaries (requires NumPy and Matplotlib, without HDDM or local posterior databases):

```bash
python Scripts/hddm/matching/matching_plot.py --kind schematic --output Output --figures-output Figures_check
```

The figure and simulation metadata are written to `Figures_check/`, preserving the published outputs.

The effective-drift table is produced by `python Scripts/hddm/matching/matching_analyze.py diagnostics --effective-drift --output Output`, which requires posterior samples from the four formal chains and NumPy, pandas, and ArviZ; plotting reads public summaries only. The seed, time step, source hashes, and choices and RTs of simulated paths are recorded in `Output/hddm/matching/figures/matching_ddm_simulation.json`.

#### Figure I-7 Simulation Settings

Fixed parameters are a = 1.224, z = 0.568, t = 0.464 s, v = 2.214, and dc = 0.321. The first two paths per stimulus are shown without selection by outcome or appearance. Paths start at the product of the group means and stop at the first detected boundary contact. Time is measured from accumulation onset; RT adds t.

Paths approximate a unit-diffusion Wiener process using Euler–Maruyama steps ($`\Delta t=0.0001`$ s, seed `2026102301`), with interpolation within the final step to locate boundary contact. Panel A does not propagate posterior parameter uncertainty or between-participant variability. It illustrates only the DDM component, excluding the outlier mixture and experimental response windows, and is not a full posterior predictive check. The plug-in starting point $`\bar a\bar z`$ is not $`E[az]`$. Existing PPC and parameter-recovery limitations still apply.
