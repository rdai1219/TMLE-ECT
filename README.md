# TMLE-ECT

Simulation code for the balancing-weight TMLE framework described in
*Targeted Learning for Externally Controlled Trials and Beyond*
(Fang, He, and Dai). Implements a general targeted maximum likelihood
estimator for the class of tilted average treatment effect estimands
$\tau_h$ (Section 3), and reproduces the Section 5.2 / Appendix A.4
simulation study comparing MLE, IPW, AIPW, and TMLE under four
nuisance-model scenarios, for five choices of the tilting function
$h(x)$: a prespecified tilting function, and the ATE, ATT, ATC, and
ATO.

## Repository structure

```
TMLE-ECT/
├── R/
│   ├── tmle_h_methods.R          Estimator/nuisance-fitting library (DGP-agnostic)
│   └── run_moderateoverlap_trunc.R   Simulation driver (DGP-specific; sources the file above)
├── results/                      Output CSVs from the driver (included so results can be
│                                  inspected without re-running the simulation)
├── LICENSE
└── README.md
```

## Files

**`R/tmle_h_methods.R`** is a general-purpose module with no
simulation-specific code: nuisance fitting (`fit_Q`, `fit_g`, either by
GLM or Super Learner), the pure estimators (`mle_h`, `ipw_h`, `aipw_h`,
`tmle_h`), the target-step bootstrap (`tmle_bootstrap_h`), and a single
user-facing entry point, `estimate_tau_h()` (or
`estimate_tau_h_given_fits()` if you already have fitted nuisances),
that ties them together. It can be reused for any $\tau_h$ estimation
problem, not just this simulation.

**`R/run_moderateoverlap_trunc.R`** sources `tmle_h_methods.R` and adds
everything specific to this simulation study: the data-generating
mechanism (Section 5.2.1, moderate-overlap design), the correctly
specified and misspecified working models for $Q$ and $g$ (Section
5.2.3), the five $h(x)$ choices, the Monte Carlo truth computation
(including the per-replication truth for ATT/ATC/ATO described in
Remark 2, since their $h(x) = \eta\{\widehat g^0(x)\}$ depends on the
estimated propensity score), and the replication driver that produces
Tables 3 and 4 (and their Appendix A.4 analogs) for each $h$-choice.

Propensity scores are truncated to $[0.025, 0.975]$ throughout,
matching the default in the `tmle` R package (Gruber and van der Laan,
2012).

## Requirements

- R (tested under R \>= 4.x)
- No packages beyond base R are required to run
  `run_moderateoverlap_trunc.R` as shipped — it fits all nuisance
  models via `glm()`.
- `tmle_h_methods.R`'s `fit_Q()`/`fit_g()` also support
  `method = "SL"` (Super Learner) as an alternative to `method = "glm"`.
  This branch is not exercised by `run_moderateoverlap_trunc.R`, but if
  you use it, install the `SuperLearner` package (and any learner
  packages referenced in your `SL.library`, e.g. `gam` for `SL.gam`)
  and attach it with `library(SuperLearner)` before calling `fit_Q()`/
  `fit_g()` with `method = "SL"`. `SuperLearner::SuperLearner()` relies
  on some internal lookups that require the package to be attached,
  not just referenced via `::`.

## Running the simulation

```r
setwd("R")            # or: Session > Set Working Directory > To Source File Location, in RStudio
source("run_moderateoverlap_trunc.R")
```

This runs $B = 1{,}000$ replications for each of the five $h$-choices
(tilted, ATE, ATT, ATC, ATO), with $n_1 = 200$ treated and $n_0 = 300$
external-control units per replication, and writes 10 CSVs to
`results/`:

- `table3_<h>_trunc_moderateoverlap.csv` -- point-estimation
  performance (truth, mean, absolute bias, empirical SD, RMSE) for
  each of MLE/IPW/AIPW/TMLE, under each of the four nuisance-model
  scenarios.
- `table4_<h>_trunc_moderateoverlap.csv` -- TMLE inference performance
  (mean SE and coverage) for the EIF-Wald interval and the two
  target-step bootstrap intervals (normal and percentile), under each
  scenario.

Runtime is on the order of tens of minutes per $h$-choice on a typical
laptop, dominated by the 300-replicate target-step bootstrap within
each of the 1,000 simulation replications.

## Results included in this repository

The `results/` folder contains the CSVs as reported in the paper
(Tilted and ATT in the main text; ATE, ATC, and ATO in Appendix A.4).
Re-running `run_moderateoverlap_trunc.R` will regenerate these files;
because the driver seeds each replication explicitly
(`set.seed(seed_offset + b)`), results should be exactly reproducible
on the same R version.

## Correspondence

Questions about the methodology should be directed to the
corresponding author of the paper. Questions about this code may be
opened as an issue on this repository.
