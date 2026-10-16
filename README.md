# Individual Differences in The Intergenerational Transmission of Attention-Deficit/Hyperactivity Disorder: Genetic and Environmental Factors

This repository contains the R analysis, visualization, and sensitivity refit scripts for the study:

> **Pak M, Lee J, Lee E, Kim B-G, Kim G, Posner J, Joo YY, Cha J.**  
> *Individual Differences in The Intergenerational Transmission of Attention-Deficit/Hyperactivity Disorder: Genetic and Environmental Factors.*  
> *Translational Psychiatry* (under revision; preprint available at Research Square: https://doi.org/10.21203/rs.3.rs-7405635/v1).

---

## Overview

The scripts in this repository implement generalized random forest (GRF; Athey, Tibshirani, & Wager, 2019; R package `grf`) models in the Adolescent Brain Cognitive Development (ABCD) Study (Release 4.0) to:
1. Estimate the covariate-adjusted average association (average treatment effect, ATE) between parental ADHD symptoms (`asr_scr_adhd_t >= 65`) and child ADHD symptoms (`cbcl_scr_dsm5_adhd_t`).
2. Test for omnibus heterogeneity in the conditional association (`test_calibration` slope on out-of-bag causal forest predictions) across four nested covariate configurations (environmental, genetic [30 polygenic scores], and structural brain volume covariates) at baseline and 1-, 2-, and 3-year follow-ups.
3. Identify candidate covariates along which the fitted conditional association varies using a 500-iteration permutation-based variable importance test (after Gong et al., 2021) and characterize their partial dependence curves.
4. Execute comprehensive revision sensitivity analyses (`revision_r1/`), including:
   - Measured baseline child BMI sensitivity refits (`Table S14`: substituting measured child BMI for BMI PGS, and including both simultaneously).
   - Longitudinal causal forest refits adjusting for baseline child ADHD symptoms at 1-, 2-, and 3-year follow-ups (`Table S15`).
   - Non-oversampled European-ancestry sensitivity refits under both sex-matched ($N = 726$) and full unexposed ($N = 6,106$) comparison groups (`Table S16`).
   - Shared-informant exclusion sensitivity refits (excluding parental psychopathology indicators, child sleep problems, or all parent-reported covariates).

---

## Repository Structure

| Path | Description |
|---|---|
| `Analysis2.R` | Complete primary pipeline: ABCD Release 4.0 table loading (`baseline_year_1_arm_1` and longitudinal waves), variable construction, near-zero-variance filtering (`caret` v6.0-93), $k$-nearest-neighbor imputation (`impute` v1.72.3, $k = 10$), sex-matched control downsampling (`set.seed(1)`), European SMOTE oversampling (`DMwR` v0.4.1), 25,000-tree causal forest fitting (`5 seeds × 5,000 trees`, clustered on `site_id_l`), permutation importance testing, partial dependence computation, and regional brain volume correlation tests (`Tables S4–S11`). |
| `Visualization.R` | Script generating Figure 1 (ATE, calibration test, and longitudinal trajectory plots) and Figure 2 (top 20 variable importance bar plots and partial dependence curves). |
| `Code.R` | Standalone template script for fitting the 25,000-tree causal forest, 500-permutation variable importance threshold, and partial dependence curves on a single prepared analytic data frame (`data_multi_env_gen`). |
| `revision_r1/extract_summary_diagnostics.R` | Reproducible extraction of raw/analytic CBCL $t$-score standard deviations, out-of-bag individual CATE 95% confidence intervals, partial dependence quartile contrasts (`Q1–Q3`, `Q0–Q4`), PGS discovery sample size rank correlations, and sample verification. |
| `revision_r1/run_refit_grf_models.R` | Reproducible causal forest sensitivity refits for non-oversampled European samples (`Table S16`), baseline-adjusted longitudinal outcomes (`Table S15`), and shared-informant covariate exclusions. |
| `revision_r1/run_bmi_sensitivity_refit.R` | Reproducible causal forest sensitivity refits evaluating measured baseline child BMI alongside or in place of the BMI polygenic score (`Table S14`). |
| `revision_r1/results/` | Aggregate machine-readable JSON summaries and Markdown tables (`summary_diagnostics.json`, `sensitivity_refits.json`, `table_s14_bmi.md`, `table_s16_european.md`). |

---

## Data Availability

The data analyzed in this study are from the **Adolescent Brain Cognitive Development (ABCD) Study Data Release 4.0** (`doi:10.15154/1523041`), available to qualified researchers with an approved Data Use Certification (DUC) through the NIMH Data Archive (https://nda.nih.gov).

In strict compliance with the ABCD Data Use Certification, **no individual-level participant data, subject identifiers (`subjectkey`), or derived participant-level `.RData` workspaces are distributed in this repository.** To reproduce the analyses, obtain ABCD Release 4.0 tabular files from the NIMH Data Archive and place them in the working directory expected by `Analysis2.R`.

---

## Software Environment

As documented in the manuscript and Supplementary Material:
- **R**: v4.2.1 (RStudio 2024.4.2.764; revision sensitivity refits verified under R 4.4)
- **`grf`**: v2.3.2 (generalized random forests; Athey et al., 2019)
- **`caret`**: v6.0-93 (near-zero-variance filtering of structural brain features)
- **`impute`**: v1.72.3 ($k$-nearest-neighbor missing data imputation, `impute.knn`)
- **`DMwR`**: v0.4.1 (synthetic minority oversampling in the secondary European analysis; archived on CRAN on 2021-03-16 and installable from https://cran.r-project.org/src/contrib/Archive/DMwR/DMwR_0.4.1.tar.gz)

Note on `grf::test_calibration`: `grf::test_calibration` returns one-sided $p$-values (`pnorm(t.value, lower.tail = FALSE)`). Throughout the revised manuscript and Supplementary Material, all calibration test $p$-values are reported as two-sided (`2 * pnorm(-abs(t.value))`).

---

## License

This code is released under the [MIT License](LICENSE).

## Contact

Corresponding author: **Jiook Cha, PhD**, Department of Psychology, Seoul National University (`connectome@snu.ac.kr`).
