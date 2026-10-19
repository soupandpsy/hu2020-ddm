# Hu et al. (2020): reproduction

This folder contains the code, necessary data and selected results for the Pilot and Confirm studies. It is the `01_reproduction_2020` workline; the 2026 reanalysis has not started.

## Files

```text
code/                       Analysis, per-model Fit/Run, and separate plotting/PPC
data/raw/pilot/              Original Pilot CSVs
data/raw/confirm/            Original Confirm CSVs
data/*_input.csv             Three Confirm HDDM input tables
data/posterior/              Saved group samples for Figure 4
results/behavior/pilot/      Pilot statistics, ex-Gaussian and figures
results/behavior/confirm/    Confirm statistics and figures
results/hddm/               Selected main-model summaries, diagnostics and plots
results/figures/             Figure 4 and labeled extensions
sensitivity/results/        Identity vz p_outlier=0/.05 extension
```

The current code is the organized reproduction version. This upload package contains only code, data, selected results, sensitivity results, README.md and .gitignore. Author materials are available at [hcp4715/moralSelf_ddm](https://github.com/hcp4715/moralSelf_ddm). In this local desktop workspace, author originals are under `../../../reference/moralSelf_ddm-master/`, historical scripts and full models under `../../../history/01_original_hddm/` and `../../../history/02_reproduction_workspace/`, and earlier organization materials under `../../../history/local_history/`. These external folders are not part of the GitHub package.

Names follow **study/task + necessary distinction + purpose**. Model labels are `vtz` (condition-dependent v/t, estimated z) and `vz` (condition-dependent v, shared t, estimated z); other parameters such as a still exist. `p0` means p_outlier=0, `p0p05` means .05, and `r01` identifies a run. Directory context supplies fields omitted from short filenames such as `summary.csv`.

## Redraw existing figures

Run from this folder, using your installed R and Python:

```text
Rscript --vanilla code/pilot_behavior/pilot_plot.R --analysis-dir results/behavior/pilot
Rscript --vanilla code/confirm_behavior/confirm_plot.R --analysis-dir results/behavior/confirm
python -B code/hddm/figure04_plot.py
```

The R plots read saved tables; they do not repeat statistics or ex-Gaussian fitting. Figure 4 reads the bundled NPZ group samples; it does not require HDDM or a desktop-specific path. Every redraw creates a new timestamped directory, leaving selected results in place.

|Paper figure|Saved output|Entry|
|---|---|---|
|Figure 2A|[Pilot Matching](results/behavior/pilot/figures/pilot_matching_behavior.pdf)|pilot_plot.R|
|Figure 2B|[Pilot Categorization](results/behavior/pilot/figures/pilot_categorization_behavior.pdf)|pilot_plot.R|
|Figure 3A|[Confirm Matching](results/behavior/confirm/figures/confirm_matching_behavior.pdf)|confirm_plot.R|
|Figure 3B|[Confirm Identity](results/behavior/confirm/figures/confirm_identity_behavior.pdf)|confirm_plot.R|
|Figure 3C|[Confirm Valence](results/behavior/confirm/figures/confirm_valence_behavior.pdf)|confirm_plot.R|
|Figure 4A/B/C|[Confirm drift posteriors](results/figures/figure04_vtz_p0_r01_drift.pdf)|figure04_plot.py|

PNG versions accompany the PDFs. Figure 4 uses vtz, p_outlier=0, r01, with 9,000 retained draws per parameter. Matching mismatch drift corresponds to Supplementary Figure S6; nondecision time and starting bias plots provide additional visualizations. The 35 group parameters in `data/posterior` cannot replace full-model databases for subject-level diagnostics, multi-chain analysis or PPC. The historical 50-bin histogram / 300-point spline curve method is retained.

## Analyze raw data

```text
Rscript --vanilla code/pilot_behavior/pilot_analysis.R
Rscript --vanilla code/confirm_behavior/confirm_analysis.R
```

Both entries read project-local raw CSVs and print a new output directory under `results/behavior/`. Pass that directory to the respective plotting entry. `--data-dir PATH` and `--output-dir PATH` are optional; an existing nonempty output directory is not overwritten.

Pilot analysis includes ex-Gaussian fitting with `retimes::timefit(iter=1000)`. Behavioral RT uses milliseconds; HDDM/ex-Gaussian uses seconds. Different trial filters are intentionally retained. See [data/README.md](data/README.md). Candidate HDDM inputs are exported separately and do not replace the three established Confirm inputs.

## Reproduce from raw data through new model figures

Run all commands from this project directory. Use R for the first two stages and your HDDM Python environment for fitting. In WSL, use Linux paths (for example `/mnt/c/...`), not `C:\...` paths.

1. Run `Rscript --vanilla code/pilot_behavior/pilot_analysis.R` and `Rscript --vanilla code/confirm_behavior/confirm_analysis.R`. Record the output directory printed by each command.
2. Run `Rscript --vanilla code/pilot_behavior/pilot_plot.R --analysis-dir PILOT_DIR` and `Rscript --vanilla code/confirm_behavior/confirm_plot.R --analysis-dir CONFIRM_DIR`, replacing the directory placeholders with the actual outputs.
3. Pass the completed Confirm directory directly to each model. No copying or renaming of CSVs is required:

```text
python -B code/hddm/matching/vtz/matching_vtz_run.py --mode check --source new --behavior-dir CONFIRM_DIR
python -B code/hddm/matching/vtz/matching_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
python -B code/hddm/identity/vtz/identity_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
python -B code/hddm/valence/vtz/valence_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
```

`--behavior-dir` selects `hddm_input/MS_match_hddm_stim.csv`, `MS_categ_id_hddm_stim.csv`, or `MS_categ_val_hddm_stim.csv` according to the model task. The same option works with all six model Run entries, including vz. It rejects Pilot directories. The selected input path is saved with the fit configuration, and each fit saves its own `input.csv`. Without this option, the included `data/*_input.csv` files are used. The original inputs are never overwritten. For a separate analysis of already fitted samples, omit `--behavior-dir` and use `--mode analyze --source new`.

4. Record each model's new `analysis_new_*` directory. Inspect its diagnostics before interpreting results. Pass it to the task Plot entry:

```text
python -B code/hddm/matching/matching_plot.py --analysis-dir MATCHING_ANALYSIS_DIR
python -B code/hddm/identity/identity_plot.py --analysis-dir IDENTITY_ANALYSIS_DIR
python -B code/hddm/valence/valence_plot.py --analysis-dir VALENCE_ANALYSIS_DIR
```

5. Draw Figure 4 from those analyses, using their r01 samples:

```text
python -B code/hddm/figure04_plot.py --source new --matching-analysis MATCHING_ANALYSIS_DIR --identity-analysis IDENTITY_ANALYSIS_DIR --valence-analysis VALENCE_ANALYSIS_DIR
```

For one task, use `--task matching --source new --matching-analysis MATCHING_ANALYSIS_DIR` (or the corresponding Identity/Valence options). Add `--check` to validate sample loading and curve calculation without saving figures. New-source plotting requires a vtz p0 analysis containing r01 and its full sample database. It reads the database recorded in that analysis; it never falls back to bundled samples. It plots r01, not pooled chains, and does not certify convergence. Output names retain task/model/setting/run, inside a fresh `results/figures/figure04_*` directory; `sources.json` records the selected source.

`--source bundled` (default) redraws the included historical group samples; `--source existing` reads the configured historical full databases. Neither chooses newly fitted samples automatically.

These commands complete the input-to-figure routing. They do not add the outstanding supplementary analyses or claim full paper replication. Full fitting and PPC were not rerun for this routing change.

## Fit or analyze a model

Each model has its own Fit and Run file. Fit defines one fit; Run selects configured runs and analyzes the samples. Plotting and PPC remain separate. JSON stores the model settings, seeds, input/output paths and historical sample locations.

Example: fit Matching vtz four times and analyze the new samples:

```text
python -B code/hddm/matching/vtz/matching_vtz_run.py --mode all --runs r01 r02 r03 r04
```

Use `--runs r01` for one run, or `--runs r02 r03 r04` to run the remaining configured IDs. Existing run folders are not overwritten or automatically resumed. New fits go under the configured model folder in `runs/`; analysis goes into a new `analysis_*` folder. Choose a different `output_root` in the config if you intentionally need another set of fits.

Analyze full historical samples already available locally:

```text
python -B code/hddm/matching/vtz/matching_vtz_run.py --mode analyze --source existing
```

The config's historical paths point to the local `history/02_reproduction_workspace/` outside this package. On another computer, supply the corresponding full models/databases and update those paths, or fit from the included inputs. Full databases are not part of this GitHub package.

Draw diagnostic plots from a selected analysis directory, for example:

```text
python -B code/hddm/matching/matching_plot.py --analysis-dir results/hddm/matching_vtz_p0
```

These plots require full samples. PPC is an independent optional step; for example `matching_ppc.py --config code/hddm/matching/vtz/matching_vtz_p0_config.json --run r01 --samples 500`. Add `--source new` to use a newly fitted model. Without `--samples`, PPC checks that the selected run has a model, input CSV and full sample database; it does not load the model or simulate. Only the selected run is required. New-run PPC uses that run's saved `input.csv`, not a potentially changed original input path. PPC was not rerun during this organization.

|Model setting|Selected results|Interpretation|
|---|---|---|
|Matching vtz p0|[Summary](results/hddm/matching_vtz_p0/summary.csv)|Four runs; diagnostic gate passed|
|Matching vz p0|[Summary](results/hddm/matching_vz_p0/summary.csv)|One run; descriptive only|
|Identity vtz p0|[Summary](results/hddm/identity_vtz_p0/summary.csv)|Four runs; diagnostic gate passed|
|Valence vtz p0|[Summary](results/hddm/valence_vtz_p0/summary.csv)|Four runs; diagnostic gate passed|
|Valence vz p0|[Summary](results/hddm/valence_vz_p0/summary.csv)|One run; descriptive only|
|Identity vz p0|[Summary](sensitivity/results/identity_vz_p0/summary.csv)|Extension; four runs; diagnostic gate passed|
|Identity vz p0p05|[Summary](sensitivity/results/identity_vz_p0p05/summary.csv)|Extension; diagnostic gate failed; do not pool|

Identity vz requires an explicit configuration; there is no default outlier setting:

```text
python -B code/hddm/identity/vz/identity_vz_run.py --config code/hddm/identity/vz/identity_vz_p0_config.json --mode all --runs r01 r02 r03 r04
python -B code/hddm/identity/vz/identity_vz_run.py --config code/hddm/identity/vz/identity_vz_p0p05_config.json --mode all --runs r01 r02 r03 r04
```

Add `--behavior-dir CONFIRM_DIR` to fit directly from a new Confirm R analysis. Both settings write to their separate directories under `sensitivity/results/`.

Each model folder contains `summary.csv`, `diagnostics.csv`, `contrasts.csv`, `status.json`, `config.json`, `sources.csv`, and two diagnostic PDFs. The model source table records run/seed/database paths, without integrity-check machinery.

## Runtime used

R 4.6.0; Python 3.7.12; HDDM 0.8.0; NumPy 1.19.5; SciPy 1.5.3; Matplotlib 3.3.2. R packages: Rcpp 1.1.2, dplyr 1.2.1, tidyr 1.3.2, ggplot2 4.0.3, afex 1.5.1, BayesFactor 0.9.12.4.8, effectsize 1.0.3 and retimes 0.1-2. The latter is archived on CRAN. The scripts support an installed project-local `.rlib`, which is excluded from Git. Actual R session/package records accompany the selected behavior results. No environment-reconstruction files are provided.

## Scientific limitations and differences

- Pilot retained 29 participants; Confirm retained 42 for Matching and 41 for Categorization. The three Confirm model-input tables remain unchanged. The prior input comparisons agreed; 15 paper t targets agreed within explicit tolerances (Pilot < .01, Confirm < .02), not necessarily exact printed rounding.
- Figure 3C uses Valence means/error bars, correcting the author's Identity-summary reference. Error bars remain mean +/- 1.96 SE. Pilot participant lines connect Good/Bad within Subject x Identity; Confirm connects each Subject's four ordered conditions. The different supplied plotting definitions are retained.
- The author Pilot Categorization RT JASP analysis assigns the condition columns in a different order, starting with `RT_Val_Bad_Other` in an Id/Good/Self cell. Its seven F/omega comparisons are not treated as equivalent to the label-based R analysis. The other 33 comparable behavioral ANOVA F values agreed within 1e-9.
- Current R omega definitions differ from older JASP. ANOVA inclusion Bayes factors have numerical differences; Pilot Categorization Accuracy JASP additionally treats Val and Id as nuisance terms, unlike the supplied R script. Paired-test directions must also be matched. These are not claimed as identical paper results.
- Pilot ex-Gaussian has 231 completed cells and one unfit cell: Subject 7014, Val/Other/Bad, only 3 valid trials. The downstream n>36/eight-valid-cell rule retains 27 subjects. Inspect [all parameters](results/behavior/pilot/exgaussian/parameters_all.csv). Compared with available author values, maximum mu/sigma/tau differences are about .0283/.0567/.0208 seconds. Refit ANOVA is not claimed identical; no comparable frequentist tau ANOVA was found in the saved author JASP file. Internal bootstrap optimizer convergence is not exposed by the legacy package.
- Noncentral-t precision warnings for effect-size intervals and data-frame coercion warnings are retained with the results. Completion alone is not a claim that all statistics match the paper.
- Existing HDDM numerical results and figure redraws were checked before and after organization. No new full MCMC or PPC was run. Pilot DDM was not reported in the published paper; it is not implemented here. **2026 reanalysis is not started.** Identity vz p0/.05 is an extension, separate from the main reproduction.

The supplementary analyses remain partial: the ex-Gaussian bootstrap uses 1,000 rather than the paper's 2,000 iterations; supplementary ex-Gaussian Bayesian/simple-effect analyses and Figure S3, full Table S2/S3 summaries, and Bayes-factor robustness figures are not all implemented.

## Upload scope

Upload this folder's code, data, selected results and documentation. `.gitignore` excludes the local package library, logs, caches, full model databases and newly generated duplicate runs. Full historical fitted models, original author files and earlier organization materials remain in the local reference/history folders outside this package; they have not been deleted. No upload has been performed.

## Complete current file directory and functions

This inventory lists every file in this upload package. The local package library, reference archive and caches are not included. New runtime directories are described below.

```text
01_reproduction_2020/
|-- code/  # Executable analysis, fitting and plotting code
|   |-- confirm_behavior/  # Confirm behavioral analysis
|   |   |-- confirm_analysis.R  # Raw-data processing, statistics and model-input export
|   |   `-- confirm_plot.R  # Behavior figures from saved analysis tables
|   |-- hddm/  # Confirm HDDM analyses
|   |   |-- identity/  # Identity task
|   |   |   |-- vtz/  # Condition-dependent v/t, estimated z
|   |   |   |   |-- identity_vtz_fit.py  # Define and fit one identity_vtz model run
|   |   |   |   |-- identity_vtz_p0_config.json  # Model settings, seeds and input/output locations
|   |   |   |   `-- identity_vtz_run.py  # Select inputs/runs, call Fit and analyze identity_vtz samples
|   |   |   |-- vz/  # Condition-dependent v, shared t, estimated z
|   |   |   |   |-- identity_vz_fit.py  # Define and fit one identity_vz model run
|   |   |   |   |-- identity_vz_p0_config.json  # Model settings, seeds and input/output locations
|   |   |   |   |-- identity_vz_p0p05_config.json  # Model settings, seeds and input/output locations
|   |   |   |   `-- identity_vz_run.py  # Select inputs/runs, call Fit and analyze identity_vz samples
|   |   |   |-- identity_plot.py  # Posterior and trace diagnostic plots for identity
|   |   |   `-- identity_ppc.py  # Posterior predictive checks for identity
|   |   |-- matching/  # Matching task
|   |   |   |-- vtz/  # Condition-dependent v/t, estimated z
|   |   |   |   |-- matching_vtz_fit.py  # Define and fit one matching_vtz model run
|   |   |   |   |-- matching_vtz_p0_config.json  # Model settings, seeds and input/output locations
|   |   |   |   `-- matching_vtz_run.py  # Select inputs/runs, call Fit and analyze matching_vtz samples
|   |   |   |-- vz/  # Condition-dependent v, shared t, estimated z
|   |   |   |   |-- matching_vz_fit.py  # Define and fit one matching_vz model run
|   |   |   |   |-- matching_vz_p0_config.json  # Model settings, seeds and input/output locations
|   |   |   |   `-- matching_vz_run.py  # Select inputs/runs, call Fit and analyze matching_vz samples
|   |   |   |-- matching_plot.py  # Posterior and trace diagnostic plots for matching
|   |   |   `-- matching_ppc.py  # Posterior predictive checks for matching
|   |   |-- valence/  # Valence task
|   |   |   |-- vtz/  # Condition-dependent v/t, estimated z
|   |   |   |   |-- valence_vtz_fit.py  # Define and fit one valence_vtz model run
|   |   |   |   |-- valence_vtz_p0_config.json  # Model settings, seeds and input/output locations
|   |   |   |   `-- valence_vtz_run.py  # Select inputs/runs, call Fit and analyze valence_vtz samples
|   |   |   |-- vz/  # Condition-dependent v, shared t, estimated z
|   |   |   |   |-- valence_vz_fit.py  # Define and fit one valence_vz model run
|   |   |   |   |-- valence_vz_p0_config.json  # Model settings, seeds and input/output locations
|   |   |   |   `-- valence_vz_run.py  # Select inputs/runs, call Fit and analyze valence_vz samples
|   |   |   |-- valence_plot.py  # Posterior and trace diagnostic plots for valence
|   |   |   `-- valence_ppc.py  # Posterior predictive checks for valence
|   |   `-- figure04_plot.py  # Figure 4 and related plots from bundled, historical or new samples
|   `-- pilot_behavior/  # Pilot behavioral analysis
|       |-- pilot_analysis.R  # Raw-data processing, statistics and model-input export; ex-Gaussian fitting
|       `-- pilot_plot.R  # Behavior figures from saved analysis tables
|-- data/  # Inputs and portable historical samples
|   |-- posterior/  # Group posterior samples for portable redraw
|   |   |-- README.md  # Scope, instructions and file explanations
|   |   |-- identity_vtz_p0_r01_samples.npz  # Saved historical r01 group posterior samples
|   |   |-- manifest.json  # Task/run/parameter mapping for NPZ files
|   |   |-- matching_vtz_p0_r01_samples.npz  # Saved historical r01 group posterior samples
|   |   `-- valence_vtz_p0_r01_samples.npz  # Saved historical r01 group posterior samples
|   |-- raw/  # Unmodified raw CSV inputs
|   |   |-- confirm/  # Confirmatory study
|   |   |   |-- MS_rep_categTask_raw.csv  # Original categorization trial data
|   |   |   `-- MS_rep_matchingTask_raw.csv  # Original matching trial data
|   |   `-- pilot/  # Pilot study
|   |       |-- MS_categTask_raw.csv  # Original categorization trial data
|   |       `-- MS_matchTask_raw.csv  # Original matching trial data
|   |-- README.md  # Scope, instructions and file explanations
|   |-- identity_input.csv  # Prepared Confirm identity model input
|   |-- matching_input.csv  # Prepared Confirm matching model input
|   `-- valence_input.csv  # Prepared Confirm valence model input
|-- results/  # Selected reference results
|   |-- behavior/  # Behavioral analyses
|   |   |-- confirm/  # Confirmatory study
|   |   |   |-- figures/  # Saved result figures
|   |   |   |   |-- confirm_identity_behavior.pdf  # Behavior figure (PDF)
|   |   |   |   |-- confirm_identity_behavior.png  # Behavior figure (PNG)
|   |   |   |   |-- confirm_matching_behavior.pdf  # Behavior figure (PDF)
|   |   |   |   |-- confirm_matching_behavior.png  # Behavior figure (PNG)
|   |   |   |   |-- confirm_valence_behavior.pdf  # Behavior figure (PDF)
|   |   |   |   `-- confirm_valence_behavior.png  # Behavior figure (PNG)
|   |   |   |-- MS_categ_behav_noTask_wide.csv  # Categorization summaries collapsed across tasks
|   |   |   |-- MS_categ_behav_wide.csv  # Participant-level categorization summaries by task
|   |   |   |-- MS_cross_taskeffect_wide.csv  # Participant-level effects joined across tasks
|   |   |   |-- MS_match_behav_wide.csv  # Participant-level matching summaries
|   |   |   |-- completed.txt  # Successful R entry completion marker
|   |   |   |-- confirm_anova_bayes_inclusion.csv  # ANOVA inclusion Bayes factors
|   |   |   |-- confirm_anova_bayes_models.csv  # Bayesian ANOVA model comparisons
|   |   |   |-- confirm_anova_frequentist.csv  # Frequentist behavioral ANOVA
|   |   |   |-- confirm_anova_omega_current_R.csv  # Current R omega effect sizes
|   |   |   |-- confirm_categ_ACC_long.csv  # Accuracy long-format analysis/plot table
|   |   |   |-- confirm_categ_ACC_noTask_long.csv  # Accuracy long-format analysis/plot table; tasks collapsed
|   |   |   |-- confirm_categ_RT_long.csv  # RT long-format analysis/plot table
|   |   |   |-- confirm_categ_RT_noTask_long.csv  # RT long-format analysis/plot table; tasks collapsed
|   |   |   |-- confirm_cross_task_correlations.csv  # Cross-task effect correlations
|   |   |   |-- confirm_exclusions.csv  # Excluded participants
|   |   |   |-- confirm_key_reproduction_checks.csv  # Selected paper-statistic comparisons
|   |   |   |-- confirm_match_RT_long.csv  # RT long-format analysis/plot table
|   |   |   |-- confirm_match_d_long.csv  # Sensitivity d-prime long-format analysis/plot table
|   |   |   |-- confirm_planned_tests.csv  # Planned paired comparisons
|   |   |   |-- confirm_sample_flow.csv  # Participant and trial counts through screening
|   |   |   |-- dependencies.csv  # R package versions
|   |   |   |-- reproduction_notes.txt  # Analysis-specific reproduction notes
|   |   |   |-- sessionInfo.txt  # R runtime and session information
|   |   |   `-- warnings.txt  # Warnings recorded during analysis
|   |   `-- pilot/  # Pilot study
|   |       |-- exgaussian/  # Pilot ex-Gaussian estimates and statistics
|   |       |   |-- exgaussian_anova.csv  # Frequentist ex-Gaussian ANOVA
|   |       |   |-- parameters_all.csv  # All ex-Gaussian cell estimates and fit status
|   |       |   `-- parameters_valid.csv  # Ex-Gaussian estimates retained for statistics
|   |       |-- figures/  # Saved result figures
|   |       |   |-- pilot_categorization_behavior.pdf  # Behavior figure (PDF)
|   |       |   |-- pilot_categorization_behavior.png  # Behavior figure (PNG)
|   |       |   |-- pilot_matching_behavior.pdf  # Behavior figure (PDF)
|   |       |   `-- pilot_matching_behavior.png  # Behavior figure (PNG)
|   |       |-- MS_categ_behav_noTask_wide.csv  # Categorization summaries collapsed across tasks
|   |       |-- MS_categ_behav_wide.csv  # Participant-level categorization summaries by task
|   |       |-- MS_cross_taskeffect_wide.csv  # Participant-level effects joined across tasks
|   |       |-- MS_match_behav_wide.csv  # Participant-level matching summaries
|   |       |-- completed.txt  # Successful R entry completion marker
|   |       |-- dependencies.csv  # R package versions
|   |       |-- pilot_anova_bayes_inclusion.csv  # ANOVA inclusion Bayes factors
|   |       |-- pilot_anova_bayes_models.csv  # Bayesian ANOVA model comparisons
|   |       |-- pilot_anova_frequentist.csv  # Frequentist behavioral ANOVA
|   |       |-- pilot_anova_omega_current_R.csv  # Current R omega effect sizes
|   |       |-- pilot_categ_ACC_long.csv  # Accuracy long-format analysis/plot table
|   |       |-- pilot_categ_ACC_noTask_long.csv  # Accuracy long-format analysis/plot table; tasks collapsed
|   |       |-- pilot_categ_RT_long.csv  # RT long-format analysis/plot table
|   |       |-- pilot_categ_RT_noTask_long.csv  # RT long-format analysis/plot table; tasks collapsed
|   |       |-- pilot_cross_task_correlations.csv  # Cross-task effect correlations
|   |       |-- pilot_exclusions.csv  # Excluded participants
|   |       |-- pilot_key_reproduction_checks.csv  # Selected paper-statistic comparisons
|   |       |-- pilot_match_RT_long.csv  # RT long-format analysis/plot table
|   |       |-- pilot_match_d_long.csv  # Sensitivity d-prime long-format analysis/plot table
|   |       |-- pilot_planned_tests.csv  # Planned paired comparisons
|   |       |-- pilot_sample_flow.csv  # Participant and trial counts through screening
|   |       |-- sessionInfo.txt  # R runtime and session information
|   |       `-- warnings.txt  # Warnings recorded during analysis
|   |-- figures/  # Saved result figures
|   |   |-- figure04_vtz_p0_r01_drift.pdf  # Combined Figure 4 (PDF)
|   |   |-- figure04_vtz_p0_r01_drift.png  # Combined Figure 4 (PNG)
|   |   |-- identity_vtz_p0_r01_drift.pdf  # Drift-rate posterior (PDF)
|   |   |-- identity_vtz_p0_r01_drift.png  # Drift-rate posterior (PNG)
|   |   |-- identity_vtz_p0_r01_nondecision.pdf  # Nondecision-time posterior (PDF)
|   |   |-- identity_vtz_p0_r01_nondecision.png  # Nondecision-time posterior (PNG)
|   |   |-- identity_vtz_p0_r01_starting_bias.pdf  # Starting bias posterior (PDF)
|   |   |-- identity_vtz_p0_r01_starting_bias.png  # Starting bias posterior (PNG)
|   |   |-- matching_vtz_p0_r01_match_drift.pdf  # Drift-rate posterior (PDF)
|   |   |-- matching_vtz_p0_r01_match_drift.png  # Drift-rate posterior (PNG)
|   |   |-- matching_vtz_p0_r01_match_nondecision.pdf  # Nondecision-time posterior (PDF)
|   |   |-- matching_vtz_p0_r01_match_nondecision.png  # Nondecision-time posterior (PNG)
|   |   |-- matching_vtz_p0_r01_mismatch_drift.pdf  # Drift-rate posterior (PDF)
|   |   |-- matching_vtz_p0_r01_mismatch_drift.png  # Drift-rate posterior (PNG)
|   |   |-- matching_vtz_p0_r01_mismatch_nondecision.pdf  # Nondecision-time posterior (PDF)
|   |   |-- matching_vtz_p0_r01_mismatch_nondecision.png  # Nondecision-time posterior (PNG)
|   |   |-- matching_vtz_p0_r01_starting_bias.pdf  # Starting bias posterior (PDF)
|   |   |-- matching_vtz_p0_r01_starting_bias.png  # Starting bias posterior (PNG)
|   |   |-- valence_vtz_p0_r01_drift.pdf  # Drift-rate posterior (PDF)
|   |   |-- valence_vtz_p0_r01_drift.png  # Drift-rate posterior (PNG)
|   |   |-- valence_vtz_p0_r01_nondecision.pdf  # Nondecision-time posterior (PDF)
|   |   |-- valence_vtz_p0_r01_nondecision.png  # Nondecision-time posterior (PNG)
|   |   |-- valence_vtz_p0_r01_starting_bias.pdf  # Starting bias posterior (PDF)
|   |   `-- valence_vtz_p0_r01_starting_bias.png  # Starting bias posterior (PNG)
|   `-- hddm/  # HDDM model results
|       |-- identity_vtz_p0/  # Model results for identity_vtz_p0
|       |   |-- config.json  # Configuration used for these results
|       |   |-- contrasts.csv  # Implemented posterior contrasts
|       |   |-- diagnostics.csv  # R-hat and effective sample size diagnostics
|       |   |-- posterior.pdf  # Group posterior diagnostic plots
|       |   |-- sources.csv  # Run, seed and full database locations
|       |   |-- status.json  # Analysis source, diagnostic status and pooling decision
|       |   |-- summary.csv  # Posterior parameter summaries by run and eligible pooled scope
|       |   `-- traceplots.pdf  # Sampling trace plots
|       |-- matching_vtz_p0/  # Model results for matching_vtz_p0
|       |   |-- config.json  # Configuration used for these results
|       |   |-- contrasts.csv  # Implemented posterior contrasts
|       |   |-- diagnostics.csv  # R-hat and effective sample size diagnostics
|       |   |-- posterior.pdf  # Group posterior diagnostic plots
|       |   |-- sources.csv  # Run, seed and full database locations
|       |   |-- status.json  # Analysis source, diagnostic status and pooling decision
|       |   |-- summary.csv  # Posterior parameter summaries by run and eligible pooled scope
|       |   `-- traceplots.pdf  # Sampling trace plots
|       |-- matching_vz_p0/  # Model results for matching_vz_p0
|       |   |-- config.json  # Configuration used for these results
|       |   |-- contrasts.csv  # Implemented posterior contrasts
|       |   |-- diagnostics.csv  # R-hat and effective sample size diagnostics
|       |   |-- posterior.pdf  # Group posterior diagnostic plots
|       |   |-- sources.csv  # Run, seed and full database locations
|       |   |-- status.json  # Analysis source, diagnostic status and pooling decision
|       |   |-- summary.csv  # Posterior parameter summaries by run and eligible pooled scope
|       |   `-- traceplots.pdf  # Sampling trace plots
|       |-- valence_vtz_p0/  # Model results for valence_vtz_p0
|       |   |-- config.json  # Configuration used for these results
|       |   |-- contrasts.csv  # Implemented posterior contrasts
|       |   |-- diagnostics.csv  # R-hat and effective sample size diagnostics
|       |   |-- posterior.pdf  # Group posterior diagnostic plots
|       |   |-- sources.csv  # Run, seed and full database locations
|       |   |-- status.json  # Analysis source, diagnostic status and pooling decision
|       |   |-- summary.csv  # Posterior parameter summaries by run and eligible pooled scope
|       |   `-- traceplots.pdf  # Sampling trace plots
|       `-- valence_vz_p0/  # Model results for valence_vz_p0
|           |-- config.json  # Configuration used for these results
|           |-- contrasts.csv  # Implemented posterior contrasts
|           |-- diagnostics.csv  # R-hat and effective sample size diagnostics
|           |-- posterior.pdf  # Group posterior diagnostic plots
|           |-- sources.csv  # Run, seed and full database locations
|           |-- status.json  # Analysis source, diagnostic status and pooling decision
|           |-- summary.csv  # Posterior parameter summaries by run and eligible pooled scope
|           `-- traceplots.pdf  # Sampling trace plots
|-- sensitivity/  # Identity vz sensitivity extension
|   `-- results/  # Selected reference results
|       |-- identity_vz_p0/  # Model results for identity_vz_p0
|       |   |-- config.json  # Configuration used for these results
|       |   |-- contrasts.csv  # Implemented posterior contrasts
|       |   |-- diagnostics.csv  # R-hat and effective sample size diagnostics
|       |   |-- posterior.pdf  # Group posterior diagnostic plots
|       |   |-- sources.csv  # Run, seed and full database locations
|       |   |-- status.json  # Analysis source, diagnostic status and pooling decision
|       |   |-- summary.csv  # Posterior parameter summaries by run and eligible pooled scope
|       |   `-- traceplots.pdf  # Sampling trace plots
|       `-- identity_vz_p0p05/  # Model results for identity_vz_p0p05
|           |-- config.json  # Configuration used for these results
|           |-- contrasts.csv  # Implemented posterior contrasts
|           |-- diagnostics.csv  # R-hat and effective sample size diagnostics
|           |-- posterior.pdf  # Group posterior diagnostic plots
|           |-- sources.csv  # Run, seed and full database locations
|           |-- status.json  # Analysis source, diagnostic status and pooling decision
|           |-- summary.csv  # Posterior parameter summaries by run and eligible pooled scope
|           `-- traceplots.pdf  # Sampling trace plots
|-- .gitignore  # Exclude local libraries, full samples, caches and repeated runs
`-- README.md  # Scope, instructions and file explanations
```

New runs create `results/behavior/<study>/analysis_*/` (R tables and `hddm_input/`), `results/hddm/<model>/runs/<run>/` (`input.csv`, `config.json`, `fit.log`, `posterior.db`, saved `model` files and `complete.json`), and `analysis_new_*/` (summaries, contrasts, diagnostics, config, sources and status). Model diagnostic plots go in the selected analysis's `plots_*/`; Figure 4 outputs go in `results/figures/figure04_*/` with `figure_index.csv` and `sources.json`. These generated directories are excluded by `.gitignore`. Sensitivity fits use the corresponding path under `sensitivity/results/`.
