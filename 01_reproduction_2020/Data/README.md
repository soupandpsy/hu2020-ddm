# Data

`raw/pilot` and `raw/confirm` retain the author's CSV names and contents. Source: hcp4715/moralSelf_ddm, respectively 1_pilot_study/Results/1_preproc and 2_confirm_study/Results/1_preproc.

The three `*_input.csv` files are Confirm HDDM inputs: Matching 24,528 trials / 42 participants; Identity 14,493 and Valence 14,431 trials / 41 participants. R exports its candidate inputs to its own output directory; these do not overwrite the formal inputs. Pilot inputs must not replace Confirm inputs.

Behavioral RT is in milliseconds; HDDM/ex-Gaussian RT is in seconds. Pilot traditional analysis removes all RT <= .2 seconds; Confirm removes correct trials with RT <= .2 seconds. HDDM inputs retain valid responses without that RT exclusion. `response` represents choice, not correctness. `stim=1` means match, Self or Good for the respective tasks.

`posterior` contains only the saved r01 group samples needed for Figure 4 and its extensions, not full model databases. Its manifest records model, setting, run, seed, parameter names and draw counts.
