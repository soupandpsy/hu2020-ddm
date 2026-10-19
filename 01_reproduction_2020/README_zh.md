# 2020 年复现

[English](README.md) | 中文

包括 Pilot（预实验）和 Confirm（正式实验）的行为分析、Pilot ex-Gaussian、Confirm HDDM，以及单独列出的敏感性分析。[原作者仓库](https://github.com/hcp4715/moralSelf_ddm)。

## 文件与软件

|目录|内容|
|---|---|
|`code/pilot_behavior/`、`code/confirm_behavior/`|R 分析和独立的行为绘图入口|
|`code/hddm/`|Matching、Identity、Valence 三个任务；每个模型有 Fit、Run 和 JSON 配置，任务级 Plot/PPC 独立|
|`data/`|原始 CSV、三份 Confirm 模型输入和精简历史后验样本|
|`results/`|行为统计、模型汇总与诊断、结果图|
|`sensitivity/results/`|Identity vz 的 p_outlier=0/.05 扩展分析|

`vtz` 允许 v/t 随条件变化并估计 z；`vz` 的 t 在条件间共享。`p0`/`p0p05` 表示 p_outlier=0/.05；`r01` 为运行编号。JSON 文件参与拟合和绘图，请保留。

使用版本：**R 4.6.0；Python 3.7.12；HDDM 0.8.0；NumPy 1.19.5；SciPy 1.5.3；Matplotlib 3.3.2**。R 包包括 Rcpp、dplyr、tidyr、ggplot2、afex、BayesFactor、effectsize 和 retimes（仅 Pilot 使用，已被 CRAN 归档）。具体包版本见各研究结果中的 `dependencies.csv`、`sessionInfo.txt`。运行前需安装依赖；仓库不包含本地包库和完整模型数据库。

以下命令均在本目录执行。大写目录占位符需替换为实际输出路径；在 WSL 中使用 Linux 路径。

## 重画已有结果

```text
Rscript --vanilla code/pilot_behavior/pilot_plot.R --analysis-dir results/behavior/pilot
Rscript --vanilla code/confirm_behavior/confirm_plot.R --analysis-dir results/behavior/confirm
python -B code/hddm/figure04_plot.py
```

上述命令读取保存的表格或样本，不重新拟合。Figure 2 对应 Pilot 行为结果；Figure 3 对应 Confirm 行为结果；[Figure 4](results/figures/figure04_vtz_p0_r01_drift.pdf) 对应 Confirm vtz p0 r01 的漂移率。Matching 的 mismatch 漂移率对应 Figure S6。Figure 4 使用随仓库提供的 35 个群体参数样本（每个 9,000 次抽样），不是多链合并结果；这些精简样本不能替代完整数据库用于诊断或 PPC。

## 重新分析与拟合

**1. 处理原始数据。** 每个入口会显示新建的分析目录；Pilot 包含 ex-Gaussian 拟合。

```text
Rscript --vanilla code/pilot_behavior/pilot_analysis.R
Rscript --vanilla code/confirm_behavior/confirm_analysis.R
```

绘制本次行为结果时，使用上面的 Plot 命令，将分析目录换成本次输出。R 入口支持可选的 `--data-dir`、`--output-dir`，不覆盖非空输出目录。筛选规则、反应编码和单位见 [data/README.md](data/README.md)。

**2. 拟合三个 Confirm 模型。** `CONFIRM_DIR` 为已完成的 Confirm R 分析目录；代码自动选择相应的导出表。

```text
python -B code/hddm/matching/vtz/matching_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
python -B code/hddm/identity/vtz/identity_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
python -B code/hddm/valence/vtz/valence_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
```

省略 `--behavior-dir` 时读取仓库内的 Confirm 输入。仅运行一次时使用 `--runs r01`。拟合样本和分析结果分别存入配置指定模型目录下的 `runs/`、`analysis_new_*`；已有运行目录不会被覆盖或自动续跑。单独分析已完成的新拟合使用 `--mode analyze --source new`。历史分析（`--source existing`）需要外部完整数据库；在其他电脑上使用时需更新配置路径。

**3. 检查诊断并绘制新结果。** 填入上一步显示的各模型分析目录。

```text
python -B code/hddm/matching/matching_plot.py --analysis-dir MATCHING_ANALYSIS_DIR
python -B code/hddm/identity/identity_plot.py --analysis-dir IDENTITY_ANALYSIS_DIR
python -B code/hddm/valence/valence_plot.py --analysis-dir VALENCE_ANALYSIS_DIR
python -B code/hddm/figure04_plot.py --source new --matching-analysis MATCHING_ANALYSIS_DIR --identity-analysis IDENTITY_ANALYSIS_DIR --valence-analysis VALENCE_ANALYSIS_DIR
```

新 Figure 4 读取各分析对应的 r01 完整样本。默认 `--source bundled` 读取仓库内的历史精简样本；`--source existing` 读取配置指定的历史数据库。重绘会新建输出目录。解释结果前请检查 `diagnostics.csv` 和 `status.json`。

**可选 PPC（后验预测检查）**，需要完整的已拟合模型：

```text
python -B code/hddm/matching/matching_ppc.py --config code/hddm/matching/vtz/matching_vtz_p0_config.json --source new --run r01 --samples 500
```

Identity 和 Valence 有对应任务入口。省略 `--samples` 时只检查文件引用，不进行模拟。

**Identity vz 敏感性分析：** 必须明确选择同一 `vz/` 目录中的 `identity_vz_p0_config.json` 或 `identity_vz_p0p05_config.json`，没有默认设置。

```text
python -B code/hddm/identity/vz/identity_vz_run.py --config code/hddm/identity/vz/identity_vz_p0_config.json --mode all --runs r01 r02 r03 r04
```

## 结果状态与差异

- Pilot 保留 29 名被试；Confirm Matching 保留 42 名，Categorization 保留 41 名。R 导出的 Confirm 输入与仓库内输入一致。15 项选定的 t 统计量在指定容差内一致（Pilot < .01；Confirm < .02），不代表严格符合论文打印值的舍入结果。
- Matching/Identity/Valence vtz p0 和 Identity vz p0 通过了配置中的四次运行诊断标准。Matching/Valence vz 只有一次运行，仅用于描述。**Identity vz p0.05 未通过诊断，不应合并其各次样本。**
- Figure 3C 将作者误用的 Identity 汇总改为 Valence 数据；误差线仍为均值 ± 1.96 SE。Pilot 个体线在 Subject × Identity 内连接 Good/Bad，Confirm 则连接每名被试的四个条件。
- Pilot Categorization RT 与作者 JASP 的条件列分配不同，对应七项 F/omega 不能视为等价比较。其余 33 项可比 ANOVA F 值在 1e-9 内一致。Omega 定义、包含贝叶斯因子及 Pilot Categorization Accuracy 的干扰因素设置也存在差异，因此不能声称所有统计量与论文相同。
- ex-Gaussian 使用 1,000 次 bootstrap，而论文报告 2,000 次。一个单元（7014，Val/Other/Bad，仅 3 个试次）未拟合；231 个完成拟合的单元筛选后得到 27 名符合条件的被试。与可取得的作者估计相比，mu/sigma/tau 最大差异约为 .0283/.0567/.0208 秒；不声称 ANOVA 等价。详见[各单元结果](results/behavior/pilot/exgaussian/parameters_all.csv)及保留的警告记录。
- 补充 ex-Gaussian 贝叶斯与简单效应分析、Figure S3、完整 Table S2/S3 汇总和贝叶斯因子稳健性图尚未全部完成。论文未报告 Pilot DDM，本项目也未实现该部分。
- 整理期间核验了已有模型分析和重绘流程，没有重新执行完整 MCMC 或 PPC。`.gitignore` 排除包库、完整样本和重复运行输出。2026 年方法再分析尚未开展。
