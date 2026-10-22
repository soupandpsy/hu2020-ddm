# Hu et al. (2020): reproduction and HDDM reanalysis

English | [中文](README_zh.md)

This repository contains a computational reproduction and an exploratory reanalysis of *Good me, bad me: Prioritization of the good-self during perceptual decision-making*. The two analyses address different questions and use different model specifications.

| Section | Scope | Main findings and limitations |
|---|---|---|
| [Part I: 2020 reproduction](01_reproduction_2020/README.md) | Pilot and Confirm behavioral statistics, Pilot ex-Gaussian fits, and Confirm HDDM | Selected statistics and model analyses were reproduced. Differences in analysis settings and incomplete supplementary analyses are documented. |
| [Part II: 2026 reanalysis](02_reanalysis_2026/README.md) | Stimulus-coded HDDM distinguishing starting-point bias from drift criterion | Matching M4 meets the main sampling diagnostics, but predictive fit and parameter recovery remain unresolved. Identity and Valence estimates are withheld because their chains did not converge. |

In the Matching reanalysis, the posterior mean difference between Self-Good and Self-Bad drift magnitude is 0.666 (95% credible interval [0.387, 0.942]). This is an exploratory estimate conditional on M4, whose model constraints differ from the 2020 specification. It is not a direct numerical replication or a causal estimate.

Each section contains `Data`, `Scripts`, and `Output`, with a bilingual README describing methods, results, limitations, and software requirements. The second section includes summary tables, diagnostic evidence, and figures; complete fitted databases and posterior draws are not distributed. Recomputing posterior-based results requires the corresponding fit outputs. Detailed reproducibility limits are described in the section READMEs.

Data and original analysis materials derive from the [author repository](https://github.com/hcp4715/moralSelf_ddm). The reanalysis follows the stimulus-coding framework in the [method repository](https://github.com/Chuan-Peng-Lab/ddm_stim_coding_posts). This repository is an independent reproduction and reanalysis, not an official update of the original study.
