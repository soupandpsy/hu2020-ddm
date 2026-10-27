# Part II: H0–H3 competition between v and dc

English · [完整中文说明](README_zh.md) · [Project home](../README.md) · [Current notebook](matching_h0_h3_results.ipynb)

Status: **28 October 2026**. Full-data fits and PPC are complete. **Reliable five-fold predictive comparison is incomplete; Part III has not formally started.** This is the current entry point. [M0–M4 guide](README_historical_m0_m4.md) and [results.ipynb](results.ipynb) are historical exploration, now treated as Part I methodological background.

## Question and method

Does the Good-self versus Bad-self advantage require stimulus-related drift v, directional drift criterion dc, or both? H0 fixes both target contrasts to zero for each subject; H1 frees log(v), H2 frees dc, and H3 frees both. All retain common nuisance structure, a/t/z, other condition effects, zero-sum log relation difficulty and the same observation model. This is an exploratory change in parameterization, not direct replication of the paper's model parameters.

The Matching model retains 42 × 600 = 25,200 events: 24,527 valid responses, 649 timeouts and 24 invalid keys. Match is the upper boundary. Response density, censoring and contamination use a consistent generative model with integration over a shared subject deadline. A log(v) contrast is not the historical raw-scale M4 contrast of 0.666.

## Diagnostic findings → adjustments → results

1. [Behavior review](matching_behavior_review.ipynb) distinguishes response coding, timeout proxies, mismatch-label composition and fast-choice patterns. Its historical behavior filtering differs from the full-event analysis.
2. [Forward simulation](matching_mechanism_simulation.ipynb) and [structural audit](matching_mechanism_audit.ipynb) show why parameter plug-in and finite grids cannot establish recovery or unique mechanisms.
3. [Observation/recovery log](matching_observation_recovery.ipynb) records consistent likelihood development and sampler changes. The old [design](matching_model_design.md) retains stage-specific plans, not current completion claims. Thirty-two simulation fits passed; only two datasets per generator were used, insufficient for precise recovery rates or coverage.
4. [Full-data diagnostics](Output/hddm/matching/continuation_stage5_20261027/four_model_diagnostics.csv) pass Rhat<1.01 and bulk/tail ESS≥400 for all four models (four chains, 8,000 retained draws each). Max Rhat: H0 1.001248, H1 1.001703, H2 1.002674, H3 1.001241.
5. Each model has 200 PPC replicates. All miss correct mean RT and timeout rate in all eight cells, and underestimate the fast Match advantage (observed .204; predicted medians −.002/.027/.022/.065). Sampling convergence does not establish model adequacy.
6. PSIS has high-k events for every model. All 20 initial five-fold short fits failed. Only repaired H2 fold4 (the fifth, zero-based fold) passed fitting and prediction; H0 fold4 still has Rhat=1.013752 and tail ESS=63.26. There is **no reliable ΔELPD ranking**. Five folds × four models means 20 training fits, not fifty comparisons.

## Limits and next steps

H3's direct refit with contamination upper support 1.5 instead of 2 seconds passed and retained positive target effects; its new-support PPC is absent. H2's corresponding short-fit tail ESS is 210.63. dc interpretation depends on the difficulty reference. The deadline proxy comes from the complete original records: even completed CV would be conditional same-subject event prediction, not leakage-free prospective or new-subject prediction. H2 has not displaced the original v account; H3 is not an established unique mechanism. Finish fold diagnostics and scoring, sensitivity and observation checks before stronger claims. Part III remains unanalysed.

## Reading, figures and reproduction

[Current notebook](matching_h0_h3_results.ipynb) contains saved figures II-1 (limited recovery), II-2 (target posteriors), II-3 (absolute PPC), II-4 (target behavioral PPC) and II-5 (fast shape), each with interpretation. [Part I notebook](../01_reproduction_2020/reproduction_results.ipynb) and [historical M4 notebook](results.ipynb) explain the earlier evidence. See [archive/naming map](docs/ARCHIVE_AND_NAMING.md), [evidence hashes](docs/evidence_manifest.json) and [validation](docs/VALIDATION.md).

From the repository root:

```bash
python -m pip install -r requirements-notebooks.txt
python -m jupyter nbconvert --execute --to notebook --output-dir .preview 01_reproduction_2020/reproduction_results.ipynb
python -m jupyter nbconvert --execute --to notebook --output-dir .preview 02_reanalysis_2026/matching_h0_h3_results.ipynb
python -m jupyter nbconvert --execute --to notebook --output-dir .preview 02_reanalysis_2026/results.ipynb
```

These display notebooks read saved CSV/JSON/PNG, without MCMC. Historical simulation/recovery logs require HDDM/SciPy/ArviZ, Linux fork and local intermediate samples; do not Run All. Mechanism scripts are preserved for inspection, some extract fixed notebook cell indices and some write results on import/execution. Full reruns require an isolated copy, missing posterior/intermediate files and the original environment. The Git whitelist publishes compact evidence and source, not databases, posterior arrays or per-event log likelihoods. No claim of a fully reproducible sampling pipeline is made.
