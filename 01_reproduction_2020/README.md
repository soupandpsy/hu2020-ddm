# Part I: Hu et al. (2020) reproduction

> **阅读导航 / Reading path (2026-10-28):** [复现结果 Notebook](reproduction_results.ipynb) → [历史 M0–M4 探索](reproduction_results.ipynb) → [第二部分 H0–H3](../02_reanalysis_2026/README_zh.md)。历史结果已合入第一部分总 Notebook，不与 H0–H3 结果混用；五折预测比较尚未可靠完成，第三部分尚未正式分析。

## 逻辑链与图表导航 / Research sequence

研究问题是先确认 Good-self 行为现象能否重现；方法沿用 Pilot/Confirm 筛选、被试内统计和历史 HDDM。核查发现选定 t 统计吻合，但 ANOVA 设置、ex-Gaussian 与部分模型诊断有差异。因此先分别报告可复现的现象和模型限制，再由 M0–M4 的预测/恢复问题引出第二部分的观察模型与 H0–H3 比较，不用行为复现直接证明机制。

[复现 Notebook](reproduction_results.ipynb)依次呈现 I-1 Pilot Matching、I-2 Pilot Categorization、I-3 Confirm Matching、I-4/5 Confirm Identity/Valence、I-6 历史漂移参数图及 I-7–I-10 的 M0–M4 路径示意、参数对比、轨迹与PPC。I 系列编号仅作本仓库导航；原论文 Figure 2/3/4 编号保持不变。各图保留区间和样本口径说明。


English · [中文](README_zh.md) · [Project Home](../README.md)

## Research Scope

This section reproduces the principal computational analyses of *Good me, bad me: Prioritization of the good-self during perceptual decision-making*: Pilot and Confirm behavioral statistics, Pilot ex-Gaussian fits, and Confirm HDDM, with separate sensitivity analyses. Selected statistics and model analyses were reproduced, but not every statistic, supplementary analysis or model specification is identical to the paper.

Original materials: [author repository](https://github.com/hcp4715/moralSelf_ddm).

## Reproduction Results

- Pilot retained 29 participants; Confirm retained 42 for Matching and 41 for Categorization. Exported Confirm inputs matched the included tables. Fifteen selected t statistics agreed within tolerances (Pilot < .01; Confirm < .02), not necessarily exact printed rounding.
- Matching/Identity/Valence vtz p0 and Identity vz p0 passed the configured four-run diagnostic gate. Matching/Valence vz have one run and are descriptive only. **Identity vz p0.05 failed the gate; do not pool its runs.**
- Figure 3C corrects the author's Identity-summary reference to Valence data. Error bars remain mean +/- 1.96 SE; Pilot lines connect Good/Bad within Subject x Identity, while Confirm lines connect each subject's four conditions.
## Differences and Limitations

- Pilot Categorization RT differs from the author's JASP condition-column assignment; its seven F/omega comparisons are not equivalent. The other 33 comparable ANOVA F values agreed within 1e-9. Omega definitions, inclusion Bayes factors, and Pilot Categorization Accuracy nuisance-factor settings also differ; statistics are not all identical to the paper.
- Ex-Gaussian uses 1,000 bootstrap iterations rather than the reported 2,000. One cell (7014, Val/Other/Bad; 3 trials) was not fitted; 231 completed cells yield 27 eligible subjects. Maximum mu/sigma/tau differences from available author estimates were about .0283/.0567/.0208 seconds; ANOVA equivalence is not claimed. See [cell results](Output/behavior/pilot/exgaussian/parameters_all.csv) and retained warnings.
- Supplementary ex-Gaussian Bayesian/simple-effect analyses and Figure S3, full Table S2/S3 summaries, and Bayes-factor robustness figures remain incomplete. Pilot DDM was not reported in the paper and is not implemented here.
- These diagnostic passes refer to the historical configured criteria; they must not be equated automatically with the rank-normalized Rhat/ESS criteria of the 2026 analysis. See the independent [2026 reanalysis](../02_reanalysis_2026/README.md).

## Files and Reproducibility

### Files and Software

|Folder|Contents|
|---|---|
|`Scripts/pilot_behavior/`, `Scripts/confirm_behavior/`|R analysis and independent behavior plots|
|`Scripts/hddm/`|Matching, Identity and Valence; each model has Fit, Run and JSON settings; task Plot/PPC entries are separate|
|`Data/`|Original CSVs, three Confirm model inputs and compact historical posterior samples|
|`Output/`|Behavior statistics, model summaries/diagnostics and figures|
|`Output/sensitivity/`|Identity vz p_outlier=0/.05 extensions|

`vtz` allows condition-dependent v/t and estimates z; `vz` shares t across conditions. `p0`/`p0p05` mean p_outlier=0/.05; `r01` is a run ID. Fitting and plotting read the model settings from JSON files.

Runtime used: **R 4.6.0; Python 3.7.12; HDDM 0.8.0; NumPy 1.19.5; SciPy 1.5.3; Matplotlib 3.3.2**. R packages: Rcpp, dplyr, tidyr, ggplot2, afex, BayesFactor, effectsize and retimes (Pilot only; archived on CRAN). Exact versions are recorded in each study's `dependencies.csv` and `sessionInfo.txt`. Install dependencies before running; local libraries and full model databases are not included.

Run commands from this directory. Replace uppercase directory placeholders with actual output paths; use Linux paths when running in WSL.

### Redraw saved results

Figure 2, 3, 4 and supplementary figure numbers in this section follow the original paper.

```bash
Rscript --vanilla Scripts/pilot_behavior/pilot_plot.R --analysis-dir Output/behavior/pilot
Rscript --vanilla Scripts/confirm_behavior/confirm_plot.R --analysis-dir Output/behavior/confirm
python -B Scripts/hddm/figure04_plot.py
```

These use saved tables/samples, without refitting. Figure 2: Pilot behavior; Figure 3: Confirm behavior; [Figure 4](Output/figures/figure04_vtz_p0_r01_drift.pdf): Confirm vtz p0 r01 drift. Matching mismatch drift corresponds to Figure S6. Figure 4 uses 35 bundled group parameters (9,000 draws each), not pooled chains; these samples cannot replace full databases for diagnostics or PPC.

### Reanalyze and refit

**1. Process raw data.** Each entry prints a new analysis directory. Pilot includes ex-Gaussian fitting.

```bash
Rscript --vanilla Scripts/pilot_behavior/pilot_analysis.R
Rscript --vanilla Scripts/confirm_behavior/confirm_analysis.R
```

To draw the new behavior results, rerun the Plot commands above with the new analysis directories. R accepts optional `--data-dir` and `--output-dir`; nonempty outputs are not overwritten. Filtering, coding and units are explained in [Data/README.md](Data/README.md).

**2. Fit the three Confirm models.** `CONFIRM_DIR` is the completed Confirm R output; the correct exported table is selected automatically.

```bash
python -B Scripts/hddm/matching/vtz/matching_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
python -B Scripts/hddm/identity/vtz/identity_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
python -B Scripts/hddm/valence/vtz/valence_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
```

Omit `--behavior-dir` to use included Confirm inputs. Use `--runs r01` for one run. Fits create `runs/` and analyses create `analysis_new_*` under the configured model output; existing runs are not overwritten or resumed. To analyze completed new runs separately, use `--mode analyze --source new`. Historical analysis (`--source existing`) requires external databases and updated configuration paths on another computer.

**3. Check diagnostics and draw new results.** Use the analysis directories printed above.

```bash
python -B Scripts/hddm/matching/matching_plot.py --analysis-dir MATCHING_ANALYSIS_DIR
python -B Scripts/hddm/identity/identity_plot.py --analysis-dir IDENTITY_ANALYSIS_DIR
python -B Scripts/hddm/valence/valence_plot.py --analysis-dir VALENCE_ANALYSIS_DIR
python -B Scripts/hddm/figure04_plot.py --source new --matching-analysis MATCHING_ANALYSIS_DIR --identity-analysis IDENTITY_ANALYSIS_DIR --valence-analysis VALENCE_ANALYSIS_DIR
```

New Figure 4 uses each analysis's r01 full samples. The default `--source bundled` uses the included historical samples; `--source existing` uses configured historical databases. Redraws create fresh directories. Check `diagnostics.csv` and `status.json` before interpreting results.

**Optional PPC**, requiring a complete fitted model:

```bash
python -B Scripts/hddm/matching/matching_ppc.py --config Scripts/hddm/matching/vtz/matching_vtz_p0_config.json --source new --run r01 --samples 500
```

Identity and Valence have equivalent task entries. Omit `--samples` for file-reference checks only.

**Identity vz sensitivity:** explicitly choose either `identity_vz_p0_config.json` or `identity_vz_p0p05_config.json` in the same `vz/` folder; there is no default setting.

```bash
python -B Scripts/hddm/identity/vz/identity_vz_run.py --config Scripts/hddm/identity/vz/identity_vz_p0_config.json --mode all --runs r01 r02 r03 r04
```

### Organization

This is a minimal adaptation of [TIER Protocol 4.0](https://www.projecttier.org/tier-protocol/protocol-4-0/): `Data` holds inputs, `Scripts` holds executable analyses, and `Output` holds results, including `Output/sensitivity`. The commands above specify the execution order. A complete TIER data appendix and full-pipeline master script are not yet provided; this is not a claim of full TIER compliance.

## References

- [Hu et al. (2020), *Good me, bad me: Prioritization of the good-self during perceptual decision-making*](https://doi.org/10.1525/collabra.301)
- [Original data and analysis code](https://github.com/hcp4715/moralSelf_ddm)
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
