# 2020 reproduction

English | [中文](README_zh.md)

Pilot and Confirm behavioral analyses, Pilot ex-Gaussian, Confirm HDDM, and separate sensitivity analyses. [Original author repository](https://github.com/hcp4715/moralSelf_ddm).

## Files and software

|Folder|Contents|
|---|---|
|`Scripts/pilot_behavior/`, `Scripts/confirm_behavior/`|R analysis and independent behavior plots|
|`Scripts/hddm/`|Matching, Identity and Valence; each model has Fit, Run and JSON settings; task Plot/PPC entries are separate|
|`Data/`|Original CSVs, three Confirm model inputs and compact historical posterior samples|
|`Output/`|Behavior statistics, model summaries/diagnostics and figures|
|`Output/sensitivity/`|Identity vz p_outlier=0/.05 extensions|

`vtz` allows condition-dependent v/t and estimates z; `vz` shares t across conditions. `p0`/`p0p05` mean p_outlier=0/.05; `r01` is a run ID. Keep the JSON files: fitting and plotting use them.

Runtime used: **R 4.6.0; Python 3.7.12; HDDM 0.8.0; NumPy 1.19.5; SciPy 1.5.3; Matplotlib 3.3.2**. R packages: Rcpp, dplyr, tidyr, ggplot2, afex, BayesFactor, effectsize and retimes (Pilot only; archived on CRAN). Exact versions are recorded in each study's `dependencies.csv` and `sessionInfo.txt`. Install dependencies before running; local libraries and full model databases are not included.

Run commands from this directory. Replace uppercase directory placeholders with actual output paths; use Linux paths when running in WSL.

## Redraw saved results

```text
Rscript --vanilla Scripts/pilot_behavior/pilot_plot.R --analysis-dir Output/behavior/pilot
Rscript --vanilla Scripts/confirm_behavior/confirm_plot.R --analysis-dir Output/behavior/confirm
python -B Scripts/hddm/figure04_plot.py
```

These use saved tables/samples, without refitting. Figure 2: Pilot behavior; Figure 3: Confirm behavior; [Figure 4](Output/figures/figure04_vtz_p0_r01_drift.pdf): Confirm vtz p0 r01 drift. Matching mismatch drift corresponds to Figure S6. Figure 4 uses 35 bundled group parameters (9,000 draws each), not pooled chains; these samples cannot replace full databases for diagnostics or PPC.

## Reanalyze and refit

**1. Process raw data.** Each entry prints a new analysis directory. Pilot includes ex-Gaussian fitting.

```text
Rscript --vanilla Scripts/pilot_behavior/pilot_analysis.R
Rscript --vanilla Scripts/confirm_behavior/confirm_analysis.R
```

To draw the new behavior results, rerun the Plot commands above with the new analysis directories. R accepts optional `--data-dir` and `--output-dir`; nonempty outputs are not overwritten. Filtering, coding and units are explained in [Data/README.md](Data/README.md).

**2. Fit the three Confirm models.** `CONFIRM_DIR` is the completed Confirm R output; the correct exported table is selected automatically.

```text
python -B Scripts/hddm/matching/vtz/matching_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
python -B Scripts/hddm/identity/vtz/identity_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
python -B Scripts/hddm/valence/vtz/valence_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
```

Omit `--behavior-dir` to use included Confirm inputs. Use `--runs r01` for one run. Fits create `runs/` and analyses create `analysis_new_*` under the configured model output; existing runs are not overwritten or resumed. To analyze completed new runs separately, use `--mode analyze --source new`. Historical analysis (`--source existing`) requires external databases and updated configuration paths on another computer.

**3. Check diagnostics and draw new results.** Use the analysis directories printed above.

```text
python -B Scripts/hddm/matching/matching_plot.py --analysis-dir MATCHING_ANALYSIS_DIR
python -B Scripts/hddm/identity/identity_plot.py --analysis-dir IDENTITY_ANALYSIS_DIR
python -B Scripts/hddm/valence/valence_plot.py --analysis-dir VALENCE_ANALYSIS_DIR
python -B Scripts/hddm/figure04_plot.py --source new --matching-analysis MATCHING_ANALYSIS_DIR --identity-analysis IDENTITY_ANALYSIS_DIR --valence-analysis VALENCE_ANALYSIS_DIR
```

New Figure 4 uses each analysis's r01 full samples. The default `--source bundled` uses the included historical samples; `--source existing` uses configured historical databases. Redraws create fresh directories. Check `diagnostics.csv` and `status.json` before interpreting results.

**Optional PPC**, requiring a complete fitted model:

```text
python -B Scripts/hddm/matching/matching_ppc.py --config Scripts/hddm/matching/vtz/matching_vtz_p0_config.json --source new --run r01 --samples 500
```

Identity and Valence have equivalent task entries. Omit `--samples` for file-reference checks only.

**Identity vz sensitivity:** explicitly choose either `identity_vz_p0_config.json` or `identity_vz_p0p05_config.json` in the same `vz/` folder; there is no default setting.

```text
python -B Scripts/hddm/identity/vz/identity_vz_run.py --config Scripts/hddm/identity/vz/identity_vz_p0_config.json --mode all --runs r01 r02 r03 r04
```

## Result status and differences

- Pilot retained 29 participants; Confirm retained 42 for Matching and 41 for Categorization. Exported Confirm inputs matched the included tables. Fifteen selected t statistics agreed within tolerances (Pilot < .01; Confirm < .02), not necessarily exact printed rounding.
- Matching/Identity/Valence vtz p0 and Identity vz p0 passed the configured four-run diagnostic gate. Matching/Valence vz have one run and are descriptive only. **Identity vz p0.05 failed the gate; do not pool its runs.**
- Figure 3C corrects the author's Identity-summary reference to Valence data. Error bars remain mean +/- 1.96 SE; Pilot lines connect Good/Bad within Subject x Identity, while Confirm lines connect each subject's four conditions.
- Pilot Categorization RT differs from the author's JASP condition-column assignment; its seven F/omega comparisons are not equivalent. The other 33 comparable ANOVA F values agreed within 1e-9. Omega definitions, inclusion Bayes factors, and Pilot Categorization Accuracy nuisance-factor settings also differ; statistics are not all identical to the paper.
- Ex-Gaussian uses 1,000 bootstrap iterations rather than the reported 2,000. One cell (7014, Val/Other/Bad; 3 trials) was not fitted; 231 completed cells yield 27 eligible subjects. Maximum mu/sigma/tau differences from available author estimates were about .0283/.0567/.0208 seconds; ANOVA equivalence is not claimed. See [cell results](Output/behavior/pilot/exgaussian/parameters_all.csv) and retained warnings.
- Supplementary ex-Gaussian Bayesian/simple-effect analyses and Figure S3, full Table S2/S3 summaries, and Bayes-factor robustness figures remain incomplete. Pilot DDM was not reported in the paper and is not implemented here.
- Existing model analyses and redraw routes were validated during organization; full MCMC and PPC were not rerun. `.gitignore` excludes libraries, full samples and repeated outputs. The 2026 reanalysis has not started.

## Organization

This is a minimal adaptation of [TIER Protocol 4.0](https://www.projecttier.org/tier-protocol/protocol-4-0/): `Data` holds inputs, `Scripts` holds executable analyses, and `Output` holds results, including `Output/sensitivity`. The commands above specify the execution order. A complete TIER data appendix and full-pipeline master script are not yet provided; this is not a claim of full TIER compliance.
