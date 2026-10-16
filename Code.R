library(grf)


############################## @ GRF analysis ##############################
# choose dataset
# four models: 1) environmental covariates only, 2) environmental and genetic covariates, 3) environmental and brain covariates, 4) environmental, genetic, and brain covariates
# two samples: European and multiethnic
data <- data_multi_env_gen

# choose variables that do not have to be z-scaled
not_cont <- c('subjectkey','parent_adhd','site_id_l',
              'famhx_ss_momdad_dprs_p','famhx_ss_momdad_dg_p','allergies','head_injury',
              'ethnic_g_AFR','ethnic_g_AMR','ethnic_g_EAS','ethnic_g_EUR',
              'preeclampsia','preterm_birth','planned_preg','alcohol','tobacco','marijuana','coc_crack','sex',
              'demo_prnt_marital_v2_1','demo_prnt_marital_v2_2',
              'demo_prnt_marital_v2_3','demo_prnt_marital_v2_4',
              'demo_prnt_marital_v2_5','demo_prnt_marital_v2_6')

# choose variables that have to be z-scaled and scale them
cont_cols <- names(data)[!(names(data) %in% not_cont)]
for(i in 2:ncol(data)){
  if(names(data)[i] %in% cont_cols){
    data[,i] <- as.vector(scale(data[,i]))
  }
}

# save computation time
start_time <- Sys.time()

# create an empty list for the results
results <- list()

# choose a GRF setting
Y <- data$child_adhd  #outcome
W <- data$parent_adhd #treatment
X <- subset(data, select=-c(child_adhd, parent_adhd, site_id_l))  #covariate matrix
site_id <- data$site_id_l  #clustering variable

# set a number of seeds in the seed ensemble
num_seeds <- 5

# generate random seeds and save them
seeds <- sample(1:10000, size=num_seeds, replace=FALSE)
results[[paste0('seeds')]] <- seeds

# create an empty list for GRF models
forest_list <- list()

# for each seed, build a GRF model with 5,000 trees
for(seed_ind in 1:length(seeds)){
  
  print(paste('seed_ind', seed_ind, sep=' = '))  #print the seed number to track the process
  
  forest <- causal_forest(X, Y, W, 
                          Y.hat=NULL, W.hat=NULL,  #Y.hat and W.hat are estimated using separate regression forests 
                          compute.oob.predictions = FALSE, #out-of-bag predictions are made based on the merged GRF model below
                          clusters = site_id,  #apply clustering
                          num.trees = 5000, 
                          tune.parameters = 'all',  #tune all tunable parameters by cross-validation
                          num.threads = 16,
                          seed = seeds[seed_ind])
  forest_list[[paste0('seed', seed_ind)]] <- forest  #save the models to the list
}

# merge the resulting models to one big forest
forest_merged <- merge_forests(forest_list, compute.oob.predictions = TRUE)

# estimate average treatment effect and save the results
ate <- average_treatment_effect(forest_merged, method='AIPW', target.sample='all')
ate_t <- ate[1] / ate[2]
ate_p <- 1.96 * (1 - pnorm(abs(ate_t)))
ate_ci_low <- ate[1] - 1 * qnorm(0.975) * ate[2] 
ate_ci_high <- ate[1] + 1 * qnorm(0.975) * ate[2] 
ate_results <- data.frame(ate[1], ate[2], ate_t, ate_p, ate_ci_low, ate_ci_high)
rownames(ate_results) <- 'ATE'
colnames(ate_results)<-c("estimate", "std.error", "t-score", "P-value", "[95% ", " CI]")
print(ate_results)
results[[paste0('ate')]] <- ate_results

# compute and save test calibration results
test_cal <- test_calibration(forest_merged)
print(test_cal)
results[[paste0('test_cal')]] <- test_cal

# compute and save variable importance
var_imp <- c(variable_importance(forest_merged))
names(var_imp) <- colnames(X)
sorted_var_imp <- sort(var_imp, decreasing=TRUE)
sorted_var_imp <- data.frame(importance=sorted_var_imp)
print(sorted_var_imp)
results[[paste0('var_imp')]] <- sorted_var_imp

# name the list with the results
multi_env_gen <- results

# save computation time
end_time <- Sys.time()
end_time - start_time

# remove forest models to empty memory space
rm(forest, forest_list, forest_merged)



############################## @ Variable importance significance test (based on Gong et al. (2021)) ##############################
# save computation time
start_time <- Sys.time()

# create an empty list for variable importance
var_imp_sign <- list()

# set a number of iterations
num_iter <- 500

# for each iteration, compute variable importance
for(iter in 1:num_iter){
  print(paste('iter', iter, sep=' = '))  #print the iteration to track the process
  
  W_perm <- permute(W)  #permute the treatment variable
  
  forest <- causal_forest(X, Y, W_perm,  #build a forest model with 5,000 trees
                          Y.hat=NULL, W.hat=NULL, 
                          compute.oob.predictions = TRUE, 
                          clusters = site_id,
                          num.trees = 5000, 
                          tune.parameters = 'all',
                          num.threads = 16,
                          seed = sample(1:10000, size=1))
  
  var_imp_iter <- c(variable_importance(forest))  #compute variable importance
  var_imp_sign[[paste0('iter=', iter)]] <- var_imp_iter  #save variable importance
  
}

# create a vector with variable importance scores from all iterations
var_imp_vec <- unlist(var_imp_sign)

# calculate a 95% threshold and save it
threshold <- quantile(var_imp_vec, probs=0.95)
results[[paste0('threshold')]] <- threshold

# find and save variables with importance scores higher than the threshold
imp_names <- rownames(sorted_var_imp)[sorted_var_imp$importance > threshold]
results[[paste0('imp_names')]] <- imp_names

# name the list with the results
multi_env_gen <- results

# save computation time
end_time <- Sys.time()
end_time - start_time

# remove forest models to empty memory space
rm(forest, forest_merged, forest_list)



############################## @ Partial dependence plots ##############################
# create a function to identify a trend in a numeric vector
is.sorted <- function(num.vec){
  if(length(num.vec) > 1 & all(sort(num.vec, decreasing=F) == num.vec)){
    return('asc')
  } else if(length(num.vec) > 1 & all(sort(num.vec, decreasing=T) == num.vec)){
    return('des')
  } else {
    return('none')
  }
}

# choose the data and results to work with
data <- data_multi_env_gen
results <- multi_env_gen

# choose variables with significant importance scores
features <- results$imp_names

# create a data frame for the results
univ.sim <- data.frame(features)
names(univ.sim) <- "feature"
univ.sim[, "tau.hat.Q0"] <- 0
univ.sim[, "tau.hat.Q1"] <- 0
univ.sim[, "tau.hat.Q2"] <- 0
univ.sim[, "tau.hat.Q3"] <- 0
univ.sim[, "tau.hat.Q4"] <- 0
univ.sim[, "trend"] <- 0

# choose variables that do not have to be z-scaled
not_cont <- c('subjectkey','parent_adhd','site_id_l',
                 'famhx_ss_momdad_dprs_p','famhx_ss_momdad_dg_p','allergies','head_injury',
                 'ethnic_g_EUR','ethnic_g_AMR','ethnic_g_AFR','ethnic_g_EAS',
                 'preeclampsia','preterm_birth','planned_preg','alcohol','tobacco','marijuana','coc_crack','sex',
                 'demo_prnt_marital_v2_1','demo_prnt_marital_v2_2',
                 'demo_prnt_marital_v2_3','demo_prnt_marital_v2_4',
                 'demo_prnt_marital_v2_5','demo_prnt_marital_v2_6')

# choose variables that have to be z-scaled and scale them
cont_cols <- names(data)[!(names(data) %in% not_cont)]
for(i in 2:ncol(data)){
  if(names(data)[i] %in% cont_cols){
    data[,i] <- as.vector(scale(data[,i]))
  }
}

# calculate the data for partial dependence plots
for(i in 1:nrow(univ.sim)){
  
  # select one important covariate at a time
  selected.covariate <- univ.sim$feature[i]
  X <- subset(data, select=-c(subjectkey, child_adhd, parent_adhd, site_id_l))
  covariates <- names(X)
  other.covariates <- covariates[which(covariates != selected.covariate)]
  
  # calculate percentiles for the selected covariate
  percentiles <- c(0, .25, .5, .75, 1)
  grid.size <- length(percentiles)
  covariate.grid <- quantile(data[,selected.covariate], probs=percentiles)
  
  # take median of other covariates 
  medians <- apply(data[, other.covariates, F], 2, median)
  
  # construct a dataset
  data.grid <- data.frame(sapply(medians, function(x) rep(x, grid.size)), covariate.grid)
  colnames(data.grid) <- c(other.covariates, selected.covariate)
  
  # expand the data
  fmla <- formula(paste0("~ 0 + ", paste0(covariates, collapse="+")))
  X.grid <- model.matrix(fmla, data.grid)
  
  # calculate point predictions of the CATE and standard errors, when other covariates are fixed at the median 
  forest.pred <- predict(results$forest, newdata = X.grid, estimate.variance=TRUE)
  tau.hat <- forest.pred$predictions
  tau.hat.se <- sqrt(forest.pred$variance.estimates)
  tau.hat.tstat <- tau.hat / tau.hat.se
  tau.hat.p <- 1.96 * (1 - pnorm(abs(tau.hat.tstat)))
  tau.hat.p.fdr <- p.adjust(tau.hat.p)
  quantiles <- c('Q0','Q1','Q2','Q3','Q4')
  univ.sim[i,c(2:6)] <- tau.hat
  
  # identify the trend in CATE estimates
  univ.sim$trend[i] <- is.sorted(tau.hat)
}
