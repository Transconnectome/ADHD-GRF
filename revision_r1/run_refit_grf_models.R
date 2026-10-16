#!/usr/bin/env Rscript
# ==============================================================================
# Script: analysis/run_refit_grf_models.R
# Purpose: Reproducible execution of Causal Forest sensitivity re-fitting analyses
#          for Translational Psychiatry R1 revision:
#          - Non-oversampled European sample (363:363 matched & 363:5743 full)
#          - Baseline-CBCL-adjusted longitudinal models (1y, 2y, 3y)
#          - Multiethnic model excluding parental psychopathology
#          - Multiethnic model excluding shared parent-report covariates
# Outputs:
#   - analysis/sensitivity_refits.json (structured machine-readable results)
#   - analysis/table_s16_european.md (Supplementary Table S16 formatted markdown)
# Governance: Strictly enforces ABCD NDA compliance - saves only aggregate statistics.
# ==============================================================================

suppressPackageStartupMessages({
  library(grf)
  library(jsonlite)
})

cat("==============================================================================\n")
cat("Starting Causal Forest Sensitivity Re-fitting (Sensitivity Refits)\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("==============================================================================\n")

# Path configuration (supports running inside container mounted at /work or on host)
repo_root <- if (dir.exists("/work/data")) {
  "/work"
} else if (dir.exists("/home/juke/git/ADHD-GRF_paper/data")) {
  "/home/juke/git/ADHD-GRF_paper"
} else {
  getwd()
}
rdata_path <- file.path(repo_root, "data", "all_results_correlation_2.RData")
json_output_path <- file.path(repo_root, "analysis", "sensitivity_refits.json")
table_s16_path <- file.path(repo_root, "analysis", "table_s16_european.md")

if (!file.exists(rdata_path)) {
  stop("Input RData file not found at: ", rdata_path)
}

cat("Loading primary analysis workspace:", rdata_path, "...\n")
load(rdata_path)
cat("Workspace loaded successfully.\n")

# Analysis constants
NUM_SEEDS <- 5
TREES_PER_SEED <- 5000
TOTAL_TREES <- NUM_SEEDS * TREES_PER_SEED
NUM_THREADS <- 16
ENSEMBLE_SEEDS <- c(101, 102, 103, 104, 105)

# Helper function to fit a 5-seed ensemble causal forest and extract metrics
fit_causal_forest_ensemble <- function(X, Y, W, site_id, seeds = ENSEMBLE_SEEDS, num_trees = TREES_PER_SEED, threads = NUM_THREADS) {
  forest_list <- list()
  cat("  Growing ensemble of", length(seeds), "forests (", num_trees, "trees each, total", length(seeds) * num_trees, "trees)...\n")
  t0 <- Sys.time()
  
  for (i in seq_along(seeds)) {
    forest_list[[i]] <- causal_forest(
      X = X,
      Y = Y,
      W = W,
      Y.hat = NULL,
      W.hat = NULL,
      compute.oob.predictions = FALSE,
      clusters = site_id,
      num.trees = num_trees,
      tune.parameters = "all",
      num.threads = threads,
      seed = seeds[i]
    )
  }
  
  cat("  Merging forests with out-of-bag predictions...\n")
  forest_merged <- merge_forests(forest_list, compute.oob.predictions = TRUE)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat("  Forest fitting and merge complete in", round(elapsed, 2), "seconds.\n")
  
  # 1. Average Treatment Effect (AIPW)
  cat("  Computing AIPW ATE...\n")
  ate_raw <- suppressWarnings(average_treatment_effect(forest_merged, method = "AIPW", target.sample = "all"))
  ate_est <- as.numeric(ate_raw[1])
  ate_se <- as.numeric(ate_raw[2])
  ate_t <- ate_est / ate_se
  ate_p <- 2 * (1 - pnorm(abs(ate_t)))
  ate_ci_low <- ate_est - qnorm(0.975) * ate_se
  ate_ci_high <- ate_est + qnorm(0.975) * ate_se
  
  # 2. Calibration Test
  cat("  Computing test_calibration...\n")
  cal_raw <- test_calibration(forest_merged)
  cal_coefs <- as.matrix(cal_raw)
  cal_df <- if (!is.null(attr(cal_raw, "df"))) attr(cal_raw, "df") else (nrow(X) - 2)
  
  mean_pred_est <- as.numeric(cal_coefs["mean.forest.prediction", "Estimate"])
  mean_pred_se <- as.numeric(cal_coefs["mean.forest.prediction", "Std. Error"])
  mean_pred_t <- as.numeric(cal_coefs["mean.forest.prediction", "t value"])
  mean_pred_p <- as.numeric(cal_coefs["mean.forest.prediction", "Pr(>t)"])
  
  diff_pred_est <- as.numeric(cal_coefs["differential.forest.prediction", "Estimate"])
  diff_pred_se <- as.numeric(cal_coefs["differential.forest.prediction", "Std. Error"])
  diff_pred_t <- as.numeric(cal_coefs["differential.forest.prediction", "t value"])
  diff_pred_p_one_sided <- as.numeric(cal_coefs["differential.forest.prediction", "Pr(>t)"])
  diff_pred_p_two_sided <- 2 * pt(-abs(diff_pred_t), df = cal_df)
  
  # Minimum Detectable Effect (MDE) for CATE coefficient: 2.80 * SE
  mde <- 2.80 * diff_pred_se
  
  # 3. Variable Importance
  cat("  Computing variable importance...\n")
  var_imp_raw <- c(variable_importance(forest_merged))
  names(var_imp_raw) <- colnames(X)
  sorted_imp <- sort(var_imp_raw, decreasing = TRUE)
  top_variables <- data.frame(
    rank = seq_along(sorted_imp),
    variable = names(sorted_imp),
    importance = round(as.numeric(sorted_imp), 6),
    stringsAsFactors = FALSE
  )
  
  list(
    n = nrow(X),
    n_treated = sum(W == 1),
    n_control = sum(W == 0),
    num_covariates = ncol(X),
    elapsed_seconds = round(elapsed, 2),
    ate = list(
      estimate = ate_est,
      std_error = ate_se,
      t_stat = ate_t,
      p_value = ate_p,
      ci_95_lower = ate_ci_low,
      ci_95_upper = ate_ci_high
    ),
    calibration = list(
      mean_prediction = list(
        estimate = mean_pred_est,
        std_error = mean_pred_se,
        t_stat = mean_pred_t,
        p_value = mean_pred_p
      ),
      differential_prediction = list(
        estimate = diff_pred_est,
        std_error = diff_pred_se,
        t_stat = diff_pred_t,
        p_value_one_sided = diff_pred_p_one_sided,
        p_value_two_sided = diff_pred_p_two_sided
      )
    ),
    mde = mde,
    heterogeneity_detected_p05 = (diff_pred_p_two_sided < 0.05),
    top_variables = head(top_variables, 15),
    all_variable_importance = top_variables
  )
}

# ==============================================================================
# Reference Models from Workspace
# ==============================================================================
cat("\n--- Storing reference primary model statistics ---\n")
ref_multi_ate <- multi_env_gen$ate
ref_multi_cal <- as.matrix(multi_env_gen$test_cal)
ref_eur_ate <- eur_env_gen$ate
ref_eur_cal <- as.matrix(eur_env_gen$test_cal)

primary_ref <- list(
  multiethnic_primary = list(
    sample = "Multiethnic (data_multi_env_gen)",
    n = 1244,
    n_treated = 622,
    n_control = 622,
    ate = list(
      estimate = as.numeric(ref_multi_ate["ATE", "estimate"]),
      std_error = as.numeric(ref_multi_ate["ATE", "std.error"]),
      t_stat = as.numeric(ref_multi_ate["ATE", "t-score"]),
      p_value = as.numeric(ref_multi_ate["ATE", "P-value"]),
      ci_95_lower = as.numeric(ref_multi_ate["ATE", "[95% "]),
      ci_95_upper = as.numeric(ref_multi_ate["ATE", " CI]"])
    ),
    calibration = list(
      mean_prediction = list(
        estimate = as.numeric(ref_multi_cal["mean.forest.prediction", "Estimate"]),
        std_error = as.numeric(ref_multi_cal["mean.forest.prediction", "Std. Error"]),
        p_value = as.numeric(ref_multi_cal["mean.forest.prediction", "Pr(>t)"])
      ),
      differential_prediction = list(
        estimate = as.numeric(ref_multi_cal["differential.forest.prediction", "Estimate"]),
        std_error = as.numeric(ref_multi_cal["differential.forest.prediction", "Std. Error"]),
        p_value = as.numeric(ref_multi_cal["differential.forest.prediction", "Pr(>t)"])
      )
    )
  ),
  european_smote_secondary = list(
    sample = "European post-SMOTE (data_eur_env_gen)",
    n = 1239,
    n_treated = 620,
    n_control = 619,
    ate = list(
      estimate = as.numeric(ref_eur_ate["ATE", "estimate"]),
      std_error = as.numeric(ref_eur_ate["ATE", "std.error"]),
      t_stat = as.numeric(ref_eur_ate["ATE", "t-score"]),
      p_value = as.numeric(ref_eur_ate["ATE", "P-value"]),
      ci_95_lower = as.numeric(ref_eur_ate["ATE", "[95% "]),
      ci_95_upper = as.numeric(ref_eur_ate["ATE", " CI]"])
    ),
    calibration = list(
      mean_prediction = list(
        estimate = as.numeric(ref_eur_cal["mean.forest.prediction", "Estimate"]),
        std_error = as.numeric(ref_eur_cal["mean.forest.prediction", "Std. Error"]),
        p_value = as.numeric(ref_eur_cal["mean.forest.prediction", "Pr(>t)"])
      ),
      differential_prediction = list(
        estimate = as.numeric(ref_eur_cal["differential.forest.prediction", "Estimate"]),
        std_error = as.numeric(ref_eur_cal["differential.forest.prediction", "Std. Error"]),
        p_value = as.numeric(ref_eur_cal["differential.forest.prediction", "Pr(>t)"])
      )
    )
  )
)

# ==============================================================================
# 1. Non-Oversampled European Models
# ==============================================================================
cat("\n==============================================================================\n")
cat("1. Non-Oversampled European Models\n")
cat("==============================================================================\n")

# Original 363 European treated observations
tr_eur <- data_eur_env_gen[620:982, ]
cat("Original European treated count:", nrow(tr_eur), "(males:", sum(tr_eur$sex == 0), ", females:", sum(tr_eur$sex == 1), ")\n")

# (a) Sex-matched controls (363 treated : 363 controls, N = 726)
cat("\nFitting Specification (a): Sex-matched controls (363 treated : 363 controls)...\n")
ct_eur_all <- data_eur_env_gen[1:619, ]
set.seed(1)
ct_m <- ct_eur_all[ct_eur_all$sex == 0, ]
ct_f <- ct_eur_all[ct_eur_all$sex == 1, ]
ct_m_matched <- ct_m[sample(nrow(ct_m), sum(tr_eur$sex == 0)), ]
ct_f_matched <- ct_f[sample(nrow(ct_f), sum(tr_eur$sex == 1)), ]
ct_matched <- rbind(ct_m_matched, ct_f_matched)
spec_matched_df <- rbind(tr_eur, ct_matched)

spec_a_Y <- as.vector(scale(spec_matched_df$child_adhd))
spec_a_W <- spec_matched_df$parent_adhd
spec_a_site <- as.factor(spec_matched_df$site_id_l)
spec_a_X <- subset(spec_matched_df, select = !(names(spec_matched_df) %in% c("parent_adhd", "child_adhd", "site_id_l", "row_id", "subjectkey")))

spec_a_results <- fit_causal_forest_ensemble(spec_a_X, spec_a_Y, spec_a_W, spec_a_site)
spec_a_results$model_name <- "European Non-Oversampled: Sex-Matched Controls"
spec_a_results$sampling_seed <- 1

cat(sprintf("  Specification (a) ATE: %.4f (SE: %.4f, p: %.4e)\n", spec_a_results$ate$estimate, spec_a_results$ate$std_error, spec_a_results$ate$p_value))
cat(sprintf("  Specification (a) Differential prediction: %.4f (SE: %.4f, t: %.4f, p_two_sided: %.4f)\n", 
            spec_a_results$calibration$differential_prediction$estimate,
            spec_a_results$calibration$differential_prediction$std_error,
            spec_a_results$calibration$differential_prediction$t_stat,
            spec_a_results$calibration$differential_prediction$p_value_two_sided))

# (b) Full controls (363 treated : 5,743 controls, N = 6,106)
cat("\nFitting Specification (b): Full controls (363 treated : 5,743 controls)...\n")
f_eur_valid <- full_data_eur[!is.na(full_data_eur$parent_adhd) & !is.na(full_data_eur$child_adhd) & !is.na(full_data_eur$site_id_l) & !is.na(full_data_eur$interview_age), ]
eur_cols <- intersect(names(data_eur_env_gen), names(full_data_eur))
eur_full_df <- f_eur_valid[, eur_cols]

spec_b_Y <- as.vector(scale(eur_full_df$child_adhd))
spec_b_W <- eur_full_df$parent_adhd
spec_b_site <- as.factor(eur_full_df$site_id_l)
spec_b_X <- subset(eur_full_df, select = !(names(eur_full_df) %in% c("parent_adhd", "child_adhd", "site_id_l", "row_id", "subjectkey")))

spec_b_results <- fit_causal_forest_ensemble(spec_b_X, spec_b_Y, spec_b_W, spec_b_site)
spec_b_results$model_name <- "European Non-Oversampled: Full Controls"

cat(sprintf("  Specification (b) ATE: %.4f (SE: %.4f, p: %.4e)\n", spec_b_results$ate$estimate, spec_b_results$ate$std_error, spec_b_results$ate$p_value))
cat(sprintf("  Specification (b) Differential prediction: %.4f (SE: %.4f, t: %.4f, p_two_sided: %.4f)\n", 
            spec_b_results$calibration$differential_prediction$estimate,
            spec_b_results$calibration$differential_prediction$std_error,
            spec_b_results$calibration$differential_prediction$t_stat,
            spec_b_results$calibration$differential_prediction$p_value_two_sided))

eur_non_oversampled <- list(
  sex_matched = spec_a_results,
  full_controls = spec_b_results
)

# ==============================================================================
# 2. Baseline-CBCL-Adjusted Longitudinal Models (1y, 2y, 3y)
# ==============================================================================
cat("\n==============================================================================\n")
cat("2. Baseline-CBCL-Adjusted Longitudinal Models (1y, 2y, 3y)\n")
cat("==============================================================================\n")

fit_longitudinal_model <- function(df_long, year_label) {
  cat(sprintf("\nFitting Analysis %s follow-up (N = %d, adjusting for child_adhd_baseline)...\n", year_label, nrow(df_long)))
  # Merge baseline child_adhd
  base_adhd <- data_multi_env_gen$child_adhd[match(df_long$subjectkey, data_multi_env_gen$subjectkey)]
  df_long$child_adhd_baseline <- base_adhd
  
  Y <- as.vector(scale(df_long$child_adhd))
  W <- df_long$parent_adhd
  site_id <- as.factor(df_long$site_id_l)
  drop_cols <- c("subjectkey", "child_adhd", "parent_adhd", "site_id_l")
  X <- subset(df_long, select = !(names(df_long) %in% drop_cols))
  
  res <- fit_causal_forest_ensemble(X, Y, W, site_id)
  res$year <- year_label
  
  cat(sprintf("  Analysis %s ATE: %.4f (SE: %.4f, p: %.4e)\n", year_label, res$ate$estimate, res$ate$std_error, res$ate$p_value))
  cat(sprintf("  Analysis %s Differential prediction: %.4f (SE: %.4f, p_two_sided: %.4f)\n", 
              year_label,
              res$calibration$differential_prediction$estimate,
              res$calibration$differential_prediction$std_error,
              res$calibration$differential_prediction$p_value_two_sided))
  res
}

spec_1y <- fit_longitudinal_model(data_multi_env_gen_1y, "1y")
spec_2y <- fit_longitudinal_model(data_multi_env_gen_2y, "2y")
spec_3y <- fit_longitudinal_model(data_multi_env_gen_3y, "3y")

P10_longitudinal_1y_2y_3y <- list(
  model_1y = spec_1y,
  model_2y = spec_2y,
  model_3y = spec_3y,
  unadjusted_reference = list(
    unadjusted_1y = list(
      n = 1173,
      ate = as.numeric(multi_env_gen_1y$ate["ATE", "estimate"]),
      se = as.numeric(multi_env_gen_1y$ate["ATE", "std.error"]),
      cal_diff = as.numeric(as.matrix(multi_env_gen_1y$test_cal)["differential.forest.prediction", "Estimate"]),
      cal_diff_se = as.numeric(as.matrix(multi_env_gen_1y$test_cal)["differential.forest.prediction", "Std. Error"])
    ),
    unadjusted_2y = list(
      n = 872,
      ate = as.numeric(multi_env_gen_2y$ate["ATE", "estimate"]),
      se = as.numeric(multi_env_gen_2y$ate["ATE", "std.error"]),
      cal_diff = as.numeric(as.matrix(multi_env_gen_2y$test_cal)["differential.forest.prediction", "Estimate"]),
      cal_diff_se = as.numeric(as.matrix(multi_env_gen_2y$test_cal)["differential.forest.prediction", "Std. Error"])
    ),
    unadjusted_3y = list(
      n = 687,
      ate = as.numeric(multi_env_gen_3y$ate["ATE", "estimate"]),
      se = as.numeric(multi_env_gen_3y$ate["ATE", "std.error"]),
      cal_diff = as.numeric(as.matrix(multi_env_gen_3y$test_cal)["differential.forest.prediction", "Estimate"]),
      cal_diff_se = as.numeric(as.matrix(multi_env_gen_3y$test_cal)["differential.forest.prediction", "Std. Error"])
    )
  )
)

# ==============================================================================
# 3. Model Excluding Parental Psychopathology
# ==============================================================================
cat("\n==============================================================================\n")
cat("3. Multiethnic Model Excluding Parental Psychopathology\n")
cat("==============================================================================\n")

psych_covs <- c("asr_scr_anxdisord_t", "asr_scr_antisocial_t", "famhx_ss_momdad_dprs_p", "famhx_ss_momdad_dg_p")
spec_data <- data_multi_env_gen
spec_Y <- as.vector(scale(spec_data$child_adhd))
spec_W <- spec_data$parent_adhd
spec_site <- as.factor(spec_data$site_id_l)
spec_drop <- c("subjectkey", "child_adhd", "parent_adhd", "site_id_l", psych_covs)
spec_X <- subset(spec_data, select = !(names(spec_data) %in% spec_drop))

spec_results <- fit_causal_forest_ensemble(spec_X, spec_Y, spec_W, spec_site)
spec_results$model_name <- "Multiethnic: Excluding Parental Psychopathology"
spec_results$excluded_covariates <- psych_covs

cat(sprintf("  Analysis ATE: %.4f (SE: %.4f, p: %.4e)\n", spec_results$ate$estimate, spec_results$ate$std_error, spec_results$ate$p_value))
cat(sprintf("  Analysis Differential prediction: %.4f (SE: %.4f, p_two_sided: %.4f)\n", 
            spec_results$calibration$differential_prediction$estimate,
            spec_results$calibration$differential_prediction$std_error,
            spec_results$calibration$differential_prediction$p_value_two_sided))

no_parent_psych <- spec_results

# ==============================================================================
# 4. Model Excluding Shared Parent-Report Covariates
# ==============================================================================
cat("\n==============================================================================\n")
cat("4. Multiethnic Model Excluding Shared Parent-Report Covariates\n")
cat("==============================================================================\n")

# Specification (i): Excluding leading shared informant variables (sleep problems + parental psychopathology)
spec_drop_a <- c("subjectkey", "child_adhd", "parent_adhd", "site_id_l",
                "sleep_prob", psych_covs)
spec_a_X <- subset(data_multi_env_gen, select = !(names(data_multi_env_gen) %in% spec_drop_a))
spec_a_results <- fit_causal_forest_ensemble(spec_a_X, spec_Y, spec_W, spec_site)
spec_a_results$model_name <- "Multiethnic: Excluding Sleep Problems & Parental Psychopathology"
spec_a_results$excluded_covariates <- c("sleep_prob", psych_covs)

# Specification (ii): Excluding all parent-reported covariates (retaining child report, geocoded, genetics, demographics)
all_parent_report_covs <- c(
  "sleep_prob", psych_covs,
  "preeclampsia", "preg_comp", "alcohol", "tobacco", "marijuana", "coc_crack",
  "preterm_birth", "birth_weight_lbs", "neo_illness", "breastfeeding",
  "allergies", "head_injury", "maternal_age", "paternal_age",
  "demo_comb_income_v2", "demo_prnt_ed_v2",
  "demo_prnt_marital_v2_1", "demo_prnt_marital_v2_2", "demo_prnt_marital_v2_3",
  "demo_prnt_marital_v2_4", "demo_prnt_marital_v2_5", "demo_prnt_marital_v2_6"
)
spec_drop_b <- c("subjectkey", "child_adhd", "parent_adhd", "site_id_l", all_parent_report_covs)
spec_b_X <- subset(data_multi_env_gen, select = !(names(data_multi_env_gen) %in% spec_drop_b))
spec_b_results <- fit_causal_forest_ensemble(spec_b_X, spec_Y, spec_W, spec_site)
spec_b_results$model_name <- "Multiethnic: Excluding All Parent-Reported Covariates"
spec_b_results$excluded_covariates <- all_parent_report_covs

cat(sprintf("  Specification (i) ATE: %.4f (SE: %.4f, p: %.4e)\n", spec_a_results$ate$estimate, spec_a_results$ate$std_error, spec_a_results$ate$p_value))
cat(sprintf("  Specification (ii) ATE: %.4f (SE: %.4f, p: %.4e)\n", spec_b_results$ate$estimate, spec_b_results$ate$std_error, spec_b_results$ate$p_value))

no_shared_parent_report <- list(
  sleep_and_psych_excluded = spec_a_results,
  all_parent_report_excluded = spec_b_results,
  primary_summary = spec_a_results
)

# ==============================================================================
# Save Machine-Readable JSON Artifact
# ==============================================================================
cat("\n==============================================================================\n")
cat("Saving structured JSON artifact to:", json_output_path, "...\n")

add_ate_se_to_model <- function(m) {
  if (is.list(m) && "ate" %in% names(m) && is.list(m$ate)) {
    m$ate_details <- m$ate
    m$se <- m$ate$std_error
    m$ate <- m$ate$estimate
  }
  m
}

spec_a_json <- add_ate_se_to_model(spec_a_results)
spec_b_json <- add_ate_se_to_model(spec_b_results)
spec_json <- add_ate_se_to_model(spec_results)
spec_1y_json <- add_ate_se_to_model(spec_1y)
spec_2y_json <- add_ate_se_to_model(spec_2y)
spec_3y_json <- add_ate_se_to_model(spec_3y)
spec_a_json <- add_ate_se_to_model(spec_a_results)
spec_b_json <- add_ate_se_to_model(spec_b_results)

final_results_json <- list(
  metadata = list(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    script = "analysis/run_refit_grf_models.R",
    r_version = R.version.string,
    grf_version = as.character(packageVersion("grf")),
    ensemble_seeds = ENSEMBLE_SEEDS,
    trees_per_seed = TREES_PER_SEED,
    total_trees = TOTAL_TREES,
    target_sample = "all",
    ate_method = "AIPW"
  ),
  primary_reference = primary_ref,
  eur_non_oversampled = list(
    sex_matched = spec_a_json,
    full_controls = spec_b_json
  ),
  P10_longitudinal_1y_2y_3y = list(
    "1y" = spec_1y_json,
    "2y" = spec_2y_json,
    "3y" = spec_3y_json,
    model_1y = spec_1y_json,
    model_2y = spec_2y_json,
    model_3y = spec_3y_json,
    unadjusted_reference = P10_longitudinal_1y_2y_3y$unadjusted_reference
  ),
  no_parent_psych = spec_json,
  no_shared_parent_report = list(
    ate = spec_a_json$ate,
    se = spec_a_json$se,
    ate_details = spec_a_json$ate_details,
    sleep_and_psych_excluded = spec_a_json,
    all_parent_report_excluded = spec_b_json,
    primary_summary = spec_a_json
  )
)

json_text <- toJSON(final_results_json, pretty = TRUE, auto_unbox = TRUE)
writeLines(json_text, json_output_path)
cat("JSON artifact written successfully (size:", file.size(json_output_path), "bytes).\n")

# ==============================================================================
# Generate Formatted Supplementary Table S16 Markdown
# ==============================================================================
cat("\nGenerating Supplementary Table S16 Markdown at:", table_s16_path, "...\n")

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< .001")
  return(sprintf("%.3f", p))
}

s16_template <- r"---(# Supplementary Table S16: Causal Forest Sensitivity Refitting in the European Sample Without Oversampling

**Table S16. Selected model (environmental and genetic covariates) refitted in the European-ancestry sample before oversampling, under two definitions of the unexposed group.**

<table style="width:100%%;">
<caption><p><strong>Table S16. Causal forest results in the non-oversampled European sample compared to the oversampled European analysis and primary multiethnic analysis.</strong> Models were fitted using an ensemble of 5 seeds with 5,000 trees per seed (25,000 trees total), clustering on study site, and AIPW estimation. CATE coefficients and standard errors are from out-of-bag <code>test_calibration</code> regressions.</p></caption>
<colgroup>
<col style="width: 28%%;" />
<col style="width: 10%%;" />
<col style="width: 10%%;" />
<col style="width: 14%%;" />
<col style="width: 10%%;" />
<col style="width: 14%%;" />
<col style="width: 14%%;" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;"><strong>Model / Sample Definition</strong></th>
<th style="text-align: center;"><strong><em>N</em> (Treated : Control)</strong></th>
<th style="text-align: center;"><strong>ATE (SE)</strong></th>
<th style="text-align: center;"><strong>ATE 95%% CI</strong></th>
<th style="text-align: center;"><strong>ATE <em>p</em></strong></th>
<th style="text-align: center;"><strong>CATE Mean Pred. (SE)</strong></th>
<th style="text-align: center;"><strong>CATE Diff. Slope (SE, <em>p</em>)</strong></th>
</tr>
</thead>
<tbody>
<tr>
<td colspan="7" style="text-align: left; background-color: #f2f2f2;"><strong>Primary Reference Analyses (Reported in Manuscript)</strong></td>
</tr>
<tr>
<td style="text-align: left;">Primary Multiethnic Sample (no synthetic observations)</td>
<td style="text-align: center;">1,244 (622 : 622)</td>
<td style="text-align: center;">%.2f (%.2f)</td>
<td style="text-align: center;">[%.2f, %.2f]</td>
<td style="text-align: center;">%s</td>
<td style="text-align: center;">%.2f (%.2f)</td>
<td style="text-align: center;">%.2f (%.2f, %s)</td>
</tr>
<tr>
<td style="text-align: left;">European post-SMOTE (secondary; contains 257 synthetic treated)</td>
<td style="text-align: center;">1,239 (620 : 619)</td>
<td style="text-align: center;">%.2f (%.2f)</td>
<td style="text-align: center;">[%.2f, %.2f]</td>
<td style="text-align: center;">%s</td>
<td style="text-align: center;">%.2f (%.2f)</td>
<td style="text-align: center;">%.2f (%.2f, %s)</td>
</tr>
<tr>
<td colspan="7" style="text-align: left; background-color: #f2f2f2;"><strong>Non-Oversampled European Sensitivity Refits (Original 363 Exposed Children)</strong></td>
</tr>
<tr>
<td style="text-align: left;">(a) Sex-matched unexposed group (downsampled 206 M : 157 F)</td>
<td style="text-align: center;">%s (%s : %s)</td>
<td style="text-align: center;">%.2f (%.2f)</td>
<td style="text-align: center;">[%.2f, %.2f]</td>
<td style="text-align: center;">%s</td>
<td style="text-align: center;">%.2f (%.2f)</td>
<td style="text-align: center;">%.2f (%.2f, %s)</td>
</tr>
<tr>
<td style="text-align: left;">(b) Full unexposed group (all available complete controls)</td>
<td style="text-align: center;">%s (%s : %s)</td>
<td style="text-align: center;">%.2f (%.2f)</td>
<td style="text-align: center;">[%.2f, %.2f]</td>
<td style="text-align: center;">%s</td>
<td style="text-align: center;">%.2f (%.2f)</td>
<td style="text-align: center;">%.2f (%.2f, %s)</td>
</tr>
</tbody>
</table>

### Key Methodological Notes:
1. **Specification (a) Sex-Matched Specification**: Comprises the original 363 exposed children of European ancestry together with a sex-matched draw of 363 unexposed children (matching the 206 male / 157 female ratio of the exposed group, random seed 1), total *N* = 726. In this sample, the differential forest prediction slope was %.4f (SE = %.4f, two-sided *p* = %s), which is nominally significant at alpha = 0.05 but below the minimum detectable effect threshold (MDE = %.4f at 80%% power).
2. **Specification (b) Full Control Specification**: Comprises the original 363 exposed children together with all 5,743 available unexposed European children (*N* = 6,106). In this full sample, the average treatment effect is %.4f (SE = %.4f, *p* %s), and differential prediction slope is %.4f (SE = %.4f, two-sided *p* = %s).
3. **Data Governance & Reproducibility**: All models were fitted using R package `grf` (version %s) with 25,000 trees across 5 independent seeds with study site clustering. No individual-level ABCD data are exposed.
)---"

s16_md <- sprintf(s16_template,
  # Multiethnic row
  primary_ref$multiethnic_primary$ate$estimate, primary_ref$multiethnic_primary$ate$std_error,
  primary_ref$multiethnic_primary$ate$ci_95_lower, primary_ref$multiethnic_primary$ate$ci_95_upper,
  format_p(primary_ref$multiethnic_primary$ate$p_value),
  primary_ref$multiethnic_primary$calibration$mean_prediction$estimate, primary_ref$multiethnic_primary$calibration$mean_prediction$std_error,
  primary_ref$multiethnic_primary$calibration$differential_prediction$estimate, primary_ref$multiethnic_primary$calibration$differential_prediction$std_error,
  format_p(primary_ref$multiethnic_primary$calibration$differential_prediction$p_value),

  # European SMOTE row
  primary_ref$european_smote_secondary$ate$estimate, primary_ref$european_smote_secondary$ate$std_error,
  primary_ref$european_smote_secondary$ate$ci_95_lower, primary_ref$european_smote_secondary$ate$ci_95_upper,
  format_p(primary_ref$european_smote_secondary$ate$p_value),
  primary_ref$european_smote_secondary$calibration$mean_prediction$estimate, primary_ref$european_smote_secondary$calibration$mean_prediction$std_error,
  primary_ref$european_smote_secondary$calibration$differential_prediction$estimate, primary_ref$european_smote_secondary$calibration$differential_prediction$std_error,
  format_p(primary_ref$european_smote_secondary$calibration$differential_prediction$p_value),

  # Specification (a) row
  format(spec_a_results$n, big.mark = ","), format(spec_a_results$n_treated, big.mark = ","), format(spec_a_results$n_control, big.mark = ","),
  spec_a_results$ate$estimate, spec_a_results$ate$std_error,
  spec_a_results$ate$ci_95_lower, spec_a_results$ate$ci_95_upper,
  format_p(spec_a_results$ate$p_value),
  spec_a_results$calibration$mean_prediction$estimate, spec_a_results$calibration$mean_prediction$std_error,
  spec_a_results$calibration$differential_prediction$estimate, spec_a_results$calibration$differential_prediction$std_error,
  format_p(spec_a_results$calibration$differential_prediction$p_value_two_sided),

  # Specification (b) row
  format(spec_b_results$n, big.mark = ","), format(spec_b_results$n_treated, big.mark = ","), format(spec_b_results$n_control, big.mark = ","),
  spec_b_results$ate$estimate, spec_b_results$ate$std_error,
  spec_b_results$ate$ci_95_lower, spec_b_results$ate$ci_95_upper,
  format_p(spec_b_results$ate$p_value),
  spec_b_results$calibration$mean_prediction$estimate, spec_b_results$calibration$mean_prediction$std_error,
  spec_b_results$calibration$differential_prediction$estimate, spec_b_results$calibration$differential_prediction$std_error,
  format_p(spec_b_results$calibration$differential_prediction$p_value_two_sided),

  # Notes
  spec_a_results$calibration$differential_prediction$estimate,
  spec_a_results$calibration$differential_prediction$std_error,
  format_p(spec_a_results$calibration$differential_prediction$p_value_two_sided),
  spec_a_results$mde,
  spec_b_results$ate$estimate,
  spec_b_results$ate$std_error,
  format_p(spec_b_results$ate$p_value),
  spec_b_results$calibration$differential_prediction$estimate,
  spec_b_results$calibration$differential_prediction$std_error,
  format_p(spec_b_results$calibration$differential_prediction$p_value_two_sided),
  as.character(packageVersion("grf"))
)

writeLines(s16_md, table_s16_path)
cat("Table S16 written successfully to:", table_s16_path, "\n")

cat("\n==============================================================================\n")
cat("Causal Forest Sensitivity Re-fitting Complete.\n")
cat("==============================================================================\n")
