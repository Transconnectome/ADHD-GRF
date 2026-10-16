#!/usr/bin/env Rscript
# ==============================================================================
# Script: analysis/run_bmi_sensitivity_refit.R
# Purpose: Reproducible execution of Table S14 Causal Forest
#          sensitivity re-fitting analyses for Translational Psychiatry R1 revision:
#          - Specification (a) (bmi_substituting_pgs): Replace BMI PGS (BMImulti) with measured
#                  baseline child BMI (66 covariates).
#          - Specification (b) (bmi_plus_pgs): Include both measured baseline child BMI and
#                  BMI PGS (BMImulti) simultaneously (67 covariates).
#          - Relocate parental psychopathology refit to shared-informant sensitivity analysis.
# Outputs:
#   - Updates analysis/sensitivity_refits.json with measured_child_bmi
#   - Generates analysis/table_s14_bmi.md (formatted Supplementary Table S14)
# Governance: Strictly enforces ABCD NDA compliance - saves only aggregate statistics.
# ==============================================================================

suppressPackageStartupMessages({
  library(grf)
  library(jsonlite)
})

cat("==============================================================================\n")
cat("Starting Measured Child BMI Sensitivity Re-fitting (Table S14)\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("==============================================================================\n")

# Path configuration (supports running inside Docker container mounted at /work or on host)
repo_root <- if (dir.exists("/work/data")) {
  "/work"
} else if (dir.exists("/home/juke/git/ADHD-GRF_paper/data")) {
  "/home/juke/git/ADHD-GRF_paper"
} else {
  getwd()
}

rdata_path <- file.path(repo_root, "data", "all_results_correlation_2.RData")
json_output_path <- file.path(repo_root, "analysis", "sensitivity_refits.json")
table_s14_path <- file.path(repo_root, "analysis", "table_s14_bmi.md")

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
    ate = ate_est,
    se = ate_se,
    ate_details = list(
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
# Extract and Verify Baseline Measured Child BMI
# ==============================================================================
cat("\n--- Extracting measured child BMI for Multiethnic analytical sample ---\n")
m_keys <- data_multi_env_gen$subjectkey

# Measured child BMI: from full_data (computed as anthroweight1lb / anthro_1_height_in^2 * 703; Analysis2.R:315)
# with KNN-imputed values for treated observations from data_tr_imp (Analysis2.R:486, 563)
bmi_vec <- full_data$bmi[match(m_keys, full_data$subjectkey)]
m_tr_imp <- match(m_keys, data_tr_imp$subjectkey)
has_tr_imp <- !is.na(m_tr_imp) & !is.na(data_tr_imp$bmi[m_tr_imp])
bmi_vec[has_tr_imp] <- data_tr_imp$bmi[m_tr_imp[has_tr_imp]]

cat(sprintf("Child BMI extracted for %d subjects: Min = %.2f, 1st Qu = %.2f, Median = %.2f, Mean = %.2f, 3rd Qu = %.2f, Max = %.2f, NAs = %d\n",
            length(bmi_vec), min(bmi_vec, na.rm=TRUE), quantile(bmi_vec, 0.25, na.rm=TRUE),
            median(bmi_vec, na.rm=TRUE), mean(bmi_vec, na.rm=TRUE), quantile(bmi_vec, 0.75, na.rm=TRUE),
            max(bmi_vec, na.rm=TRUE), sum(is.na(bmi_vec))))

# Outcome, treatment, and clustering
spec_Y <- as.vector(scale(data_multi_env_gen$child_adhd))
spec_W <- data_multi_env_gen$parent_adhd
spec_site <- as.factor(data_multi_env_gen$site_id_l)

drop_base <- c("subjectkey", "child_adhd", "parent_adhd", "site_id_l")
base_X <- subset(data_multi_env_gen, select = !(names(data_multi_env_gen) %in% drop_base))

# ==============================================================================
# 1. Specification (a): Substituting BMI PGS with Measured Child BMI
# ==============================================================================
cat("\n==============================================================================\n")
cat("1. Fitting Specification (a): Substituting BMImulti with Measured Child BMI\n")
cat("==============================================================================\n")

bmi_a_X <- subset(base_X, select = !(names(base_X) %in% "BMImulti"))
bmi_a_X$bmi <- bmi_vec
cat("Specification (a) covariate matrix:", nrow(bmi_a_X), "rows x", ncol(bmi_a_X), "columns (BMImulti excluded, bmi included)\n")

bmi_a_results <- fit_causal_forest_ensemble(bmi_a_X, spec_Y, spec_W, spec_site)
bmi_a_results$model_name <- "Multiethnic: Substituting BMI PGS with Measured Child BMI"
bmi_a_results$specification <- "bmi_substituting_pgs"
bmi_a_results$substituted_covariate <- "BMImulti"
bmi_a_results$added_covariate <- "bmi"

cat(sprintf("  Specification (a) ATE: %.4f (SE: %.4f, 95%% CI: [%.4f, %.4f], p: %.4e)\n",
            bmi_a_results$ate_details$estimate, bmi_a_results$ate_details$std_error,
            bmi_a_results$ate_details$ci_95_lower, bmi_a_results$ate_details$ci_95_upper,
            bmi_a_results$ate_details$p_value))
cat(sprintf("  Specification (a) CATE Diff slope: %.4f (SE: %.4f, t: %.4f, p_two_sided: %.4f)\n",
            bmi_a_results$calibration$differential_prediction$estimate,
            bmi_a_results$calibration$differential_prediction$std_error,
            bmi_a_results$calibration$differential_prediction$t_stat,
            bmi_a_results$calibration$differential_prediction$p_value_two_sided))

# ==============================================================================
# 2. Specification (b): Including Both BMI PGS and Measured Child BMI
# ==============================================================================
cat("\n==============================================================================\n")
cat("2. Fitting Specification (b): Including Both BMImulti and Measured Child BMI\n")
cat("==============================================================================\n")

bmi_b_X <- base_X
bmi_b_X$bmi <- bmi_vec
cat("Specification (b) covariate matrix:", nrow(bmi_b_X), "rows x", ncol(bmi_b_X), "columns (both BMImulti and bmi included)\n")

bmi_b_results <- fit_causal_forest_ensemble(bmi_b_X, spec_Y, spec_W, spec_site)
bmi_b_results$model_name <- "Multiethnic: Including Both BMI PGS and Measured Child BMI"
bmi_b_results$specification <- "bmi_plus_pgs"
bmi_b_results$added_covariate <- "bmi"

cat(sprintf("  Specification (b) ATE: %.4f (SE: %.4f, 95%% CI: [%.4f, %.4f], p: %.4e)\n",
            bmi_b_results$ate_details$estimate, bmi_b_results$ate_details$std_error,
            bmi_b_results$ate_details$ci_95_lower, bmi_b_results$ate_details$ci_95_upper,
            bmi_b_results$ate_details$p_value))
cat(sprintf("  Specification (b) CATE Diff slope: %.4f (SE: %.4f, t: %.4f, p_two_sided: %.4f)\n",
            bmi_b_results$calibration$differential_prediction$estimate,
            bmi_b_results$calibration$differential_prediction$std_error,
            bmi_b_results$calibration$differential_prediction$t_stat,
            bmi_b_results$calibration$differential_prediction$p_value_two_sided))

# ==============================================================================
# Update Machine-Readable JSON Artifact
# ==============================================================================
cat("\n==============================================================================\n")
cat("Updating structured JSON artifact at:", json_output_path, "...\n")

if (file.exists(json_output_path)) {
  existing_json <- fromJSON(json_output_path, simplifyVector = FALSE)
} else {
  stop("Existing refit JSON file not found at: ", json_output_path)
}

# 1. Record true measured child BMI refits under measured_child_bmi
existing_json$measured_child_bmi <- list(
  bmi_substituting_pgs = bmi_a_results,
  bmi_plus_pgs = bmi_b_results,
  primary_summary = bmi_a_results
)

# 2. Relocate parental psychopathology exclusion model to no_shared_parent_report
if (!is.null(existing_json$no_parent_psych)) {
  old_p8 <- existing_json$no_parent_psych
  old_p8$model_name <- "Multiethnic: Excluding Parental Psychopathology (Relocated to shared-informant sensitivity / Analysis)"
  existing_json$no_shared_parent_report$parent_psych_excluded <- old_p8
  cat("  Parental psychopathology refit successfully relocated to no_shared_parent_report$parent_psych_excluded.\n")
}

# 3. Update metadata
existing_json$metadata$bmi_refit_updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
existing_json$metadata$bmi_refit_script <- "analysis/run_bmi_sensitivity_refit.R"

# Write updated JSON
writeLines(toJSON(existing_json, pretty = TRUE, auto_unbox = TRUE), json_output_path)
cat("JSON artifact updated successfully (size:", file.size(json_output_path), "bytes).\n")

# ==============================================================================
# Generate Formatted Supplementary Table S14 Markdown
# ==============================================================================
cat("\nGenerating Supplementary Table S14 Markdown at:", table_s14_path, "...\n")

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("< .001")
  return(sprintf("%.3f", p))
}

ref_multi <- existing_json$primary_reference$multiethnic_primary

table_s14_template <- r"---(# Supplementary Table S14: Sensitivity Refits for Measured Child Body Mass Index

**Table S14. Generalized random forest results for the selected environmental and genetic model refitted with measured child body mass index.**

<table style="width:100%%;">
<caption><p><strong>Table S14. Generalized random forest results for the selected model refitted with measured child body mass index.</strong> Models were fitted in the primary multiethnic sample (<em>N</em> = 1,244; 622 exposed, 622 unexposed) using an ensemble of 5 seeds with 5,000 trees per seed (25,000 trees total), clustering on study site, and AIPW estimation. CATE coefficients and standard errors are from out-of-bag <code>test_calibration</code> regressions.</p></caption>
<colgroup>
<col style="width: 32%%;" />
<col style="width: 12%%;" />
<col style="width: 10%%;" />
<col style="width: 14%%;" />
<col style="width: 10%%;" />
<col style="width: 11%%;" />
<col style="width: 11%%;" />
</colgroup>
<thead>
<tr>
<th style="text-align: left;"><strong>Model / Specification</strong></th>
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
<td colspan="7" style="text-align: left; background-color: #f2f2f2;"><strong>Primary Reference Model (Reported in Manuscript)</strong></td>
</tr>
<tr>
<td style="text-align: left;">Primary Multiethnic Model (all 66 covariates; BMI PGS included, measured child BMI excluded)</td>
<td style="text-align: center;">1,244 (622 : 622)</td>
<td style="text-align: center;">%.2f (%.2f)</td>
<td style="text-align: center;">[%.2f, %.2f]</td>
<td style="text-align: center;">%s</td>
<td style="text-align: center;">%.2f (%.2f)</td>
<td style="text-align: center;">%.2f (%.2f, %s)</td>
</tr>
<tr>
<td colspan="7" style="text-align: left; background-color: #f2f2f2;"><strong>Measured Child BMI Sensitivity Refits</strong></td>
</tr>
<tr>
<td style="text-align: left;">(a) Substituting BMI PGS with measured baseline child BMI (66 covariates)</td>
<td style="text-align: center;">1,244 (622 : 622)</td>
<td style="text-align: center;">%.2f (%.2f)</td>
<td style="text-align: center;">[%.2f, %.2f]</td>
<td style="text-align: center;">%s</td>
<td style="text-align: center;">%.2f (%.2f)</td>
<td style="text-align: center;">%.2f (%.2f, %s)</td>
</tr>
<tr>
<td style="text-align: left;">(b) Including both BMI PGS and measured baseline child BMI (67 covariates)</td>
<td style="text-align: center;">1,244 (622 : 622)</td>
<td style="text-align: center;">%.2f (%.2f)</td>
<td style="text-align: center;">[%.2f, %.2f]</td>
<td style="text-align: center;">%s</td>
<td style="text-align: center;">%.2f (%.2f)</td>
<td style="text-align: center;">%.2f (%.2f, %s)</td>
</tr>
</tbody>
</table>

*Note.* Child body mass index (`bmi`) was calculated from measured baseline height and weight (`anthro_1_height_in` and `anthroweight1lb`) according to the standard formula (weight in lbs / height in inches^2 * 703; Table S3). In specification (a) (`bmi_substituting_pgs`), the body mass index polygenic score (`BMImulti`) was replaced by measured child `bmi` (66 covariates total). In specification (b) (`bmi_plus_pgs`), both measured child `bmi` and the BMI polygenic score (`BMImulti`) were included simultaneously (67 covariates total). In specification (a), the top five covariates by variable importance were %s. In specification (b), the top five covariates by variable importance were %s.
)---"

format_top5 <- function(top_vars) {
  items <- sapply(1:min(5, nrow(top_vars)), function(i) {
    sprintf("%s (`%s`, %.3f)", top_vars$variable[i], top_vars$variable[i], top_vars$importance[i])
  })
  paste(items, collapse = ", ")
}

s14_md <- sprintf(
  table_s14_template,
  # Primary reference
  ref_multi$ate$estimate, ref_multi$ate$std_error,
  ref_multi$ate$ci_95_lower, ref_multi$ate$ci_95_upper,
  format_p(ref_multi$ate$p_value),
  ref_multi$calibration$mean_prediction$estimate, ref_multi$calibration$mean_prediction$std_error,
  ref_multi$calibration$differential_prediction$estimate, ref_multi$calibration$differential_prediction$std_error,
  format_p(ref_multi$calibration$differential_prediction$p_value),
  # Specification (a)
  bmi_a_results$ate_details$estimate, bmi_a_results$ate_details$std_error,
  bmi_a_results$ate_details$ci_95_lower, bmi_a_results$ate_details$ci_95_upper,
  format_p(bmi_a_results$ate_details$p_value),
  bmi_a_results$calibration$mean_prediction$estimate, bmi_a_results$calibration$mean_prediction$std_error,
  bmi_a_results$calibration$differential_prediction$estimate, bmi_a_results$calibration$differential_prediction$std_error,
  format_p(bmi_a_results$calibration$differential_prediction$p_value_two_sided),
  # Specification (b)
  bmi_b_results$ate_details$estimate, bmi_b_results$ate_details$std_error,
  bmi_b_results$ate_details$ci_95_lower, bmi_b_results$ate_details$ci_95_upper,
  format_p(bmi_b_results$ate_details$p_value),
  bmi_b_results$calibration$mean_prediction$estimate, bmi_b_results$calibration$mean_prediction$std_error,
  bmi_b_results$calibration$differential_prediction$estimate, bmi_b_results$calibration$differential_prediction$std_error,
  format_p(bmi_b_results$calibration$differential_prediction$p_value_two_sided),
  # Notes top 5
  format_top5(bmi_a_results$top_variables),
  format_top5(bmi_b_results$top_variables)
)

writeLines(s14_md, table_s14_path)
cat("Table S14 written successfully to:", table_s14_path, "\n")

cat("\n==============================================================================\n")
cat("Analysis Measured Child BMI Sensitivity Re-fitting Complete.\n")
cat("==============================================================================\n")
