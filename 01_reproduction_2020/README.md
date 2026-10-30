# Part I: reproduction and historical HDDM

[中文](README_zh.md) · [Project home](../README.md) · [Part I notebook](01_reproduction_2020.ipynb)

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


## Historical exploration and reading order

The notebook proceeds from behavior (I-1–I-5), original HDDM (I-6), and historical M0–M4 paths, parameters, traces and PPC (I-7–I-10), to limitations. Historical code and outputs remain under `../02_reanalysis_2026/` for path compatibility. Matching M4 passed diagnostics but overpredicted slow RTs; Identity/Valence long chains and four recovery scenarios had convergence limitations. Behavioral reproduction does not identify a unique mechanism.

As of 28 October 2026, Part II H0–H3 full-data fits pass current sampling diagnostics, but all show RT, timeout and fast Good-self advantage PPC misfit. Reliable five-fold comparison is incomplete. H2 has not overturned the v account; H3 is not established as best or unique. Part III has not formally started.

## Essential reproduction

Install `requirements-notebooks.txt` and run from the repository root:

```bash
python -m jupyter nbconvert --execute --to notebook --output-dir .preview 01_reproduction_2020/01_reproduction_2020.ipynb
```

See [data definitions](Data/README.md). Run original-paper analysis and plotting from `01_reproduction_2020/`; analyses create new output directories and these plotting commands read saved results:

```bash
Rscript --vanilla Scripts/pilot_behavior/pilot_analysis.R
Rscript --vanilla Scripts/confirm_behavior/confirm_analysis.R
Rscript --vanilla Scripts/pilot_behavior/pilot_plot.R --analysis-dir Output/behavior/pilot
Rscript --vanilla Scripts/confirm_behavior/confirm_plot.R --analysis-dir Output/behavior/confirm
python -B Scripts/hddm/figure04_plot.py
```

Historical environment: R 4.6.0, Python 3.7.12, HDDM 0.8.0; R package versions are recorded in output `dependencies.csv` and `sessionInfo.txt`. Full HDDM fitting uses Run entries and JSON configurations under `Scripts/hddm/<task>/<model>/`. Matching example:

```bash
python -B Scripts/hddm/matching/vtz/matching_vtz_run.py --mode all --behavior-dir CONFIRM_DIR --runs r01 r02 r03 r04
```

Replace `CONFIRM_DIR` with the Confirm output; omit it to use bundled inputs. `vtz` allows condition-dependent v/t and estimates z; `vz` shares t; `p0/p0p05` specifies contamination and `r01` a run. Full databases are unpublished. Figure 4 uses 35 bundled group parameters with 9,000 samples each, not pooled multichain diagnostics.

For historical M0–M4, use the dockerHDDM environment in task configurations and run from `02_reanalysis_2026/`:

```bash
python Scripts/hddm/reproduce.py plan --output Output_new
python Scripts/hddm/reproduce.py run --output Output_new --workers 2
python Scripts/hddm/reproduce.py verify --output Output_new
```

Full sampling is expensive and was not rerun in this cleanup. The notebook only reads saved results. Continue to [Part II](../02_reanalysis_2026/README.md).
