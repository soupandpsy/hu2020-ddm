# Moral self-associations and perceptual decisions: HDDM reanalysis

[中文](README_zh.md) · [Full results and calculations](results.ipynb)

## Question and methods

Using the confirmatory data from Hu et al. (2020), this analysis examines whether stimuli associated with the good-self show faster evidence accumulation than those associated with the bad-self or good-other, while distinguishing evidence accumulation from response bias.

Identity (Self, Other) and moral valence (Good, Bad) define four conditions. Matching requires a shape–label match judgment; Identity and Valence require classification by identity and moral valence, respectively. The tasks include 42, 41 and 41 participants and 24,527, 14,493 and 14,431 valid trials. One invalid keypress was excluded from Matching.

Five models (M0–M4) were screened per task, followed by four-chain evaluation of candidate M4. M4 estimates condition-dependent drift rate v and drift bias dc, with shared starting point z, boundary separation a and nondecision time t. Matching and mismatching stimuli share the magnitude of v. Each screening fit used one chain, leaving uncertainty in model selection.

## Matching results

Drift rate was higher for good-self than for bad-self and good-other. The Good−Bad contrast was larger for Self than for Other, indicating an identity-by-valence difference in evidence accumulation.

|Drift-rate contrast|Mean difference|95% credible interval|
|---|---:|---|
|Good-self − Bad-self|0.666|[0.387, 0.942]|
|Good-self − Good-other|0.438|[0.160, 0.719]|
|Good-other − Bad-other|−0.009|[−0.278, 0.267]|
|Identity × valence interaction¹|0.675|[0.284, 1.059]|

¹ (Good-self−Bad-self)−(Good-other−Bad-other). Intervals are equal-tailed posterior intervals.

The starting point z was 0.568 [0.558, 0.578], indicating a bias toward the Match response. The good-self minus bad-self dc contrast was 0.600 [0.427, 0.774]. These biases describe response tendencies in matching judgments, rather than directly measuring a preference for the good-self.

![Matching parameter estimates and contrasts](Output/hddm/matching/figures/matching_m4_parameter_contrasts.png)

**Figure 1. Matching parameter estimates and condition contrasts.** Points are posterior means and lines are 95% credible intervals. v denotes evidence accumulation, dc drift bias, and z the relative starting point.

## Estimation stability

The four Matching chains produced similar parameter distributions and passed the convergence checks. Identity and Valence retained substantial between-chain differences after longer sampling, so their parameters are not interpreted.

|Task|Maximum Rhat|Minimum bulk ESS|Minimum tail ESS|Assessment|
|---|---:|---:|---:|---|
|Matching|1.0022|917|472|Passed|
|Identity (long refit)|1.3076|10|17|Not passed|
|Valence (long refit)|1.5282|7|27|Not passed|

Rhat assesses agreement between chains; ESS measures effective sampling information. The criteria were Rhat<1.01 and bulk/tail ESS≥400, supplemented by trace inspection.

![Matching four-chain traces and distributions](Output/hddm/matching/figures/matching_m4_traces.png)

**Figure 2. Matching traces and parameter distributions.** Colors distinguish chains; the left column shows sampling traces and the right column shows parameter distributions.

## Model checks and limitations

The model overpredicts the slow-RT tail (Figure 3). The experiment imposed a response window of approximately 0.9–1.2 seconds, whereas predictions were not subject to the same window. Trial-specific windows have not been recovered, so their contribution to the discrepancy remains unresolved.

![Matching posterior predictive checks](Output/hddm/matching/figures/matching_m4_ppc.png)

**Figure 3. Observed and predicted Matching outcomes.** Black crosses show observations; colored points and intervals show predictive means and 95% intervals of predictive statistics. Panels show choice proportions and RT quantiles. Two additional outlier-generation conventions assess sensitivity of predictions to these assumptions; neither exactly matches the full fitted outlier likelihood.

Excluding RT≤0.2 seconds changed the main group estimates little, but one individual parameter had tail ESS=342, below the criterion. All four parameter-recovery scenarios also contained convergence failures, leaving recovery performance unresolved.

The Matching results are consistent in direction with a good-self processing advantage, but model fit and recovery limit their interpretation. Parameter constraints differ from the 2020 model, preventing direct numerical comparison. Identity and Valence do not support an assessment of the original cross-task conclusion.

## Files and execution

The [notebook](results.ipynb) contains model comparisons, condition estimates and checks. Inputs are in `Data/`, calculations in `Scripts/hddm/`, and results in `Output/`. The notebook reads saved outputs without fitting models.

Within the dockerHDDM environment specified in the task configurations, run from this directory:

```sh
python Scripts/hddm/reproduce.py plan --output Output_new
python Scripts/hddm/reproduce.py run --output Output_new --workers 2
python Scripts/hddm/reproduce.py verify --output Output_new
```

A complete run requires the adjacent `01_reproduction_2020/Data` and a new output directory. The [output manifest](Scripts/hddm/output_manifest.json) identifies file producers. Postprocessing was checked using historical draws; the reorganized full fitting pipeline has not been rerun from scratch. Databases and draw-level samples are not distributed.

### Task entries and provenance

Each task retains its own configuration, models, run, analyze and plot entry points. Matching additionally has recovery and sensitivity entries. `run fit` fits one chain; the other Run phases schedule smoke, screening, final or independent long-refit jobs. Analyze operations execute independently:

```sh
python Scripts/hddm/matching/matching_run.py fit --model m4 --chain 1 --phase screening --output Output_new --plan
python Scripts/hddm/matching/matching_analyze.py compare --output Output_new --plan
python Scripts/hddm/matching/matching_analyze.py diagnostics --output Output_new --plan
python Scripts/hddm/matching/matching_analyze.py ppc --output Output_new --plan
```

Replace `matching` with `identity` or `valence` for their entries. Remove `--plan` to execute; `--help` lists each operation's options. Formal fitting appends r02–r04 to a compatible screening r01. Long refits use four independent chains. Compare reads existing screening outputs and does not rerun PPC.

In comparison tables, `status` and `chains` describe screening only. `selection_status` records provisional selection; `followup_status` distinguishes `screening_only`, `diagnostics_pass` and `nonconverged`. Passing convergence does not remove predictive or recovery limitations.

The pre-refactor public snapshot is commit `639aa66af934db39b9e2962a6c5154350266758e`, tagged `reanalysis-before-entry-refactor-20261023`. This identifies the archived public report, not the code that generated every historical fit. Original fit fingerprints remain authoritative and unchanged. The manifest's `producer` describes the current entry point; `historical_snapshot` identifies the earlier published artifact. Chain-source paths are relative to `Output/`; databases are not distributed. New fits record the changed `common.py` hash and cannot append to historical chains with different fingerprints. Use a fresh output directory for new fits, or the original fingerprint-matched code for historical continuation.

## References

[Hu et al. (2020)](https://doi.org/10.1525/collabra.301); [Pan et al. (2025), dockerHDDM](https://doi.org/10.1177/25152459241298700); [stimulus-coding code](https://github.com/Chuan-Peng-Lab/ddm_stim_coding_posts), commit `468965c8d4c991bb7633f3f38257a6bac978bce5`.
