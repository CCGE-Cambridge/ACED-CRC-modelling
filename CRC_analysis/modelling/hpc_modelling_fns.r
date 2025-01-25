# =====================================================================================
# Script Title: Survival Analysis Utility Functions
# Author: Sam Ip
# Description:
# Provides functions for fitting Cox PH models (no selection, bidirectional selection, naive LASSO, and group LASSO),
# scaling data, and utility methods for extracting and processing model coefficients.
# * Cox PH model (no selection) -- fit_1a_coxph_vanilla_model
# * Cox PH model (bidirectional selection) -- fit_1b_coxph_bidirselect_model
# * Cox PH model (naive LASSO selection) -- fit_2_lasso_model
# * Cox PH model (group LASSO selection) -- fit_2_grouplasso_model, select_lambda_coeffs, get_selected_features_groupedlevels
# * Scaling data -- calculate_scaling_params, scale_with_traindata_params
# * Optimising memory usage (strip Cox model) -- strip_coxph
# =====================================================================================

# Fit Cox PH model (no selection)
# Inputs:
#   - data: Dataset containing survival information and predictors.
#   - time_col: Name of the time-to-event column.
#   - event_col: Name of the event occurrence column.
#   - id_col: ID column for clustering.
# Output: Fitted Cox model.
fit_1a_coxph_vanilla_model <- function(data=train_data_coalition, time_col="futime", event_col="status", id_col="eid") {
  # data=train_data_coalition; time_col="futime"; event_col="status"; id_col="eid"; feature_names=setdiff(names(train_data_coalition), c("futime", "status", "eid"))
  feature_names <- setdiff(names(data), c("futime", "status", "eid"))

  # Prepare full formula
  full_fml <- ifelse(length(feature_names) == 0,
                  paste0("Surv(time = ", time_col, ", event = ", event_col, ") ~ 1"),
                  paste0("Surv(time = ", time_col, ", event = ", event_col, ") ~ ",
                     paste(setdiff(feature_names, id_col), collapse = "+"), "+", glue::glue("cluster({id_col})")))

  cat("\nfull_fml......", full_fml, "\n")

  # Fit final model
  final_model <- coxph(as.formula(full_fml), data = data, x = TRUE)
  
  return(final_model)
}

# Fit Cox PH model (bidirectional selection)
# Inputs:
#   - data: Dataset containing survival information and predictors.
#   - time_col: Name of the time-to-event column.
#   - event_col: Name of the event occurrence column.
#   - id_col: ID column for clustering.
# Output: Fitted Cox model.
fit_1b_coxph_bidirselect_model <- function(data=train_data_coalition, time_col="futime", event_col="status", id_col="eid") {
  # data=train_data_coalition; time_col="futime"; event_col="status"; id_col="eid"; feature_names=setdiff(names(train_data_coalition), c("futime", "status", "eid"))
  feature_names <- setdiff(names(data), c("futime", "status", "eid"))
  # Prepare full formula
  full_fml <- ifelse(length(feature_names) == 0,
                  paste0("Surv(time = ", time_col, ", event = ", event_col, ") ~ 1"),
                  paste0("Surv(time = ", time_col, ", event = ", event_col, ") ~ ",
                     paste(setdiff(feature_names, id_col), collapse = "+")))

  cat("\nfull_fml......", full_fml, "\n")
  # Prepare minimal formula
  min_cols <- feature_names[grepl("^gen_|^core_age|^core_sex", feature_names)]
  min_fml <- ifelse(length(min_cols) == 0,
                    paste0("Surv(time = ", time_col, ", event = ", event_col, ") ~ 1"),
                    paste0("Surv(time = ", time_col, ", event = ", event_col, ") ~ ",
                           paste(min_cols, collapse = "+")))

  cat("\nmin_fml......", min_fml, "\n")

  # Fit full and minimal models
  full_fit <- coxph(as.formula(full_fml), data = data)
  min_fit <- coxph(as.formula(min_fml), data = data)

  # Stepwise bidirectional selection
  fit_stepwise_bidirectional <- step(
    min_fit,
    scope = list(lower = min_fit, upper = full_fit),
    data = data, steps = 500, direction = "both"
  )
  
  # Extract selected covariates
  covars_kept <- attributes(fit_stepwise_bidirectional$terms)$term.labels %>% unique()
  cat("covars_kept (inc forced)......\n")
  print(covars_kept)
  
  # Prepare reduced formula
  red_fml <- if (length(covars_kept) == 0) {
    paste0("Surv(time = ", time_col, ", event = ", event_col, ") ~ 1")
  } else {
    paste0("Surv(time = ", time_col, ", event = ", event_col, ") ~ ",
           paste(c(covars_kept, glue::glue("cluster({id_col})")), collapse = "+"))
  }
  
  cat("red_fml......\n")
  print(red_fml)
  
  # Fit final model
  final_model <- coxph(as.formula(red_fml), data = data, x = TRUE)
  
  return(final_model)
}

# Identifies the lambda value based on the specified selection criterion 
# ("1se") from cross-validation results and extracts the non-zero 
# coefficients associated with the chosen lambda.
# Inputs:
#   - cvfit: Fitted cross-validation object from grpreg.
#   - lambda_min_or_1se: Criterion for lambda selection:
#       - "1se": Selects the most regularised lambda within one standard error 
#                of the minimum cross-validation error (for a simpler model).
#       - "min" (not implemented): Selects the lambda that minimises the cross-validation error.
# Outputs:
#   - A list containing:
#       - `lambda`: The chosen lambda value based on the selection criterion.
#       - `selected_features`: Features with non-zero coefficients for the chosen lambda.
#       - `binned_features`: Features with zero coefficients for the chosen lambda.
select_lambda_coeffs <- function(cvfit, lambda_min_or_1se = "1se") { 
  # Identify the minimum cross-validation error, its index, and standard error
  cve.min <- min(cvfit$cve)
  index.min <- which.min(cvfit$cve); print(cvfit$lambda[cvfit$cve==cve.min])
  cvse.min <- cvfit$cvse[index.min]

  # Identify the largest lambda within 1se threshold
  lambda.1se <- max(cvfit$lambda[cvfit$cve <= (cve.min + cvse.min)])

  # Extract coefficients for the selected lambda (1se)
  coefficients_1se <- coef(cvfit$fit, lambda = lambda.1se)

  # Identify non-zero and zero coefficients
  selected_features <- names(coefficients_1se)[coefficients_1se != 0]
  binned_features <- names(coefficients_1se)[coefficients_1se == 0]; cat("\nbinned_features......", binned_features, "\n")

  # Return 1se results as a list
  return(list(
    lambda = lambda.1se,
    selected_features=selected_features, binned_features=binned_features))
}

# Function: Map selected feature levels (resp. features) to their corresponding parent features (resp. group)
# Description:
#   Matches selected features to their corresponding grouped levels (longest matching prefix).
# Inputs:
#   - selected_features: Vector of selected feature names.
#   - feature_names: Vector of all feature names.
# Outputs:
#   - Vector of unique feature names that correspond to grouped levels.
get_selected_features_groupedlevels <- function(selected_features, feature_names) {
  # Initialise a vector to store matched feature names
  matched_features <- character(length(selected_features))
  
  for (i in seq_along(selected_features)) {
    selected_feature <- selected_features[i]
    
    # Find matching feature names based on longest prefix
    matches <- sapply(feature_names, function(fn) {
      if (str_detect(selected_feature, paste0("^", fn))) {
        return(nchar(fn))
      } else {
        return(0)
      }
    })
    if (any(matches > 0)) {
      longest_match <- feature_names[which.max(matches)]
      matched_features[i] <- longest_match
    } else {
      matched_features[i] <- NA
    }
  }
  
  # Remove NA values and get unique feature names
  unique_matched_features <- unique(matched_features[!is.na(matched_features)])
  
  return(unique_matched_features)
}

# Fit Cox PH model (naive LASSO selection), where internal CV is ID-aware
fit_2_lasso_model <- function(df=train_data_coalition,lambda_min_or_1se="1se") {
  # Identify predictor feature names (exclude survival and ID columns)
  feature_names <- setdiff(names(df), c("futime", "status", "eid"))

  # Create fold IDs for internal cross-validation (ID-aware)
  nfolds<-10
  df_foldid <- data.frame(
    eid = unique(df$eid),
    foldid = sample(rep_len(1:nfolds, dplyr::n_distinct(df$eid)))
  )
  foldid <- df %>% dplyr::left_join(df_foldid, by = "eid") %>% dplyr::pull(foldid)
  
  # Prepare the design matrix
  df_onlyfeats <- df %>% dplyr::select(feature_names)
  covariate_matrix <- model.matrix(~ . , data = df_onlyfeats)[, -1]
    
  # Set penalties: force certain features (e.g., PCs, age, sex) to remain in the model
  idx_forcekeep <- grepl("^gen_pc|^core_age|^core_sex", colnames(covariate_matrix))
  penalty <- rep(1, ncol(covariate_matrix)); penalty[idx_forcekeep] <- 0

  # Fit LASSO model with internal cross-validation
  fit_lasso <- glmnet::cv.glmnet(
    x = covariate_matrix, y = Surv(df$futime, df$status),
    type.measure = "C", family = "cox", foldid = foldid,
    penalty.factor=penalty
    )
  cat("fit_lasso......\n"); print(fit_lasso)

  # Extract coefficients for the chosen lambda (1se or min)
  idx_lambda.1se <- which(fit_lasso$lambda %in% fit_lasso$lambda.1se)
  res_lasso <- fit_lasso$glmnet.fit$beta[,idx_lambda.1se] %>% as.data.frame()
  # idx_lambda.min <- which(fit_lasso$lambda %in% fit_lasso$lambda.min)
  # res_lasso <- fit_lasso$glmnet.fit$beta[,idx_lambda.min] %>% as.data.frame() 
  cat("res_lasso......\n"); print(res_lasso)

  # Filter for non-zero coefficients and retrieve variable names
  colnames(res_lasso) <- "estimate"
  res_lasso <- res_lasso %>% dplyr::mutate(term=rownames(res_lasso)) %>% dplyr::filter(estimate!=0L)
  covars_kept <- res_lasso %>% dplyr::pull(term)

  # Map feature names to readable labels (if applicable)
  if ("core_sex_genetic1" %in% covars_kept) {
    covars_kept <- sapply(covars_kept, function(covar) {
    best_match <- readable_names$covar_code2[which.max(stringdist::stringsim(covar, readable_names$covar_code2))]
    return(best_match)
  })} 

  # Fit the final Cox PH model with the selected features
  final_model <- fit_1a_coxph_vanilla_model(data=df %>% dplyr::select("eid", "futime", "status", all_of(covars_kept)), time_col="futime", event_col="status", id_col="eid")
  
  return(final_model)
}

# Fit Cox PH model (group LASSO selection), where internal CV is ID-aware
fit_2_grouplasso_model <- function(df=train_data_coalition, grp_featlvls, lambda_min_or_1se="1se", nlambda=100, time_col="futime", event_col="status", id_col="eid") {
  # Identify predictor feature names (exclude survival and ID columns)
  feature_names <- setdiff(names(df), c("futime", "status", "eid"))
  
  # Prepare design matrix and response
  x <- model.matrix(as.formula(paste("~", paste(feature_names, collapse = "+"), "- 1")), data = df)
  y <- Surv(df[["futime"]], df[["status"]])
  
  # Create eid-aware fold IDs for internal cross-validation
  nfolds <- 10
  df_foldid <- data.frame(
      eid = unique(df$eid),
      foldid = sample(rep_len(1:nfolds, dplyr::n_distinct(df$eid)))
    )
  foldid <- df %>% dplyr::left_join(df_foldid, by = "eid") %>% dplyr::pull(foldid)

  # Fit group LASSO model with internal cross-validation
  cvfit <- grpreg::cv.grpsurv(x, y, group=grp_featlvls, nlambda=nlambda, nfolds=10, fold=foldid) 
  res_lambda1se <- select_lambda_coeffs(cvfit, lambda_min_or_1se = "1se")

  # Map selected feature levels to their corresponding parent features
  selected_features_groupedlevels <- get_selected_features_groupedlevels(res_lambda1se$selected_features, feature_names); 
  cat("\nselected_features_groupedlevels......", glue::glue("{length(selected_features_groupedlevels)}/{length(feature_names)}"), "\n", selected_features_groupedlevels, "\n")
  cat("\nbinned grouped levels......", selected_features_groupedlevels, "\n")

  # Get indices of features (e.g., PCs, age, sex) forced to remain in the model
  idx_forcekeep <- grepl("^gen_pc|^core_age|^core_sex", feature_names); print(feature_names[idx_forcekeep])

  # Fit final Cox PH model with selected features
  final_model <- fit_1a_coxph_vanilla_model(data=df %>% dplyr::select("eid", "futime", "status", all_of(unique(c(feature_names[idx_forcekeep], selected_features_groupedlevels)))), time_col="futime", event_col="status", id_col="eid")
  
  return(final_model)
}


# Computes mean and standard deviation for scaling
calculate_scaling_params <- function(data) {
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  numeric_cols <- setdiff(numeric_cols, c("status", "eid", "futime")) 
  scaling_params <- data.table(
    column = numeric_cols,
    mean = sapply(data[, ..numeric_cols], mean, na.rm = TRUE),
    sd = sapply(data[, ..numeric_cols], sd, na.rm = TRUE)
  )
  return(scaling_params)
}

# Scales numeric predictors using given means and standard deviations.
scale_with_traindata_params <- function(data, scaling_params) {
    if (nrow(scaling_params) == 0) {
    return(data)
  }
  
  scaled_data <- copy(data)
  for (i in 1:nrow(scaling_params)) {
    col <- scaling_params$column[i]
    col_mean <- scaling_params$mean[i]
    col_sd <- scaling_params$sd[i]
    scaled_data[, (col) := (get(col) - col_mean) / col_sd]
  }
  return(scaled_data)
}

# Strip unnecessary components from a fitted Cox model to optimise memory usage
strip_coxph <- function(cox_model){
  cox_model_stripped <- cox_model
  cox_model_stripped$model <- NULL           # Remove the model frame (actual data)
  cox_model_stripped$y <- NULL               # Remove the response variable
  cox_model_stripped$x <- NULL               # Remove the model matrix
  cox_model_stripped$call$data <- NULL       # Remove the data reference in the call
  cox_model_stripped$residuals <- NULL       # Remove residuals if not needed
  cox_model_stripped$linear.predictors <- NULL # Remove linear predictors if not needed
  cox_model_stripped$fitted.values <- NULL   # Remove fitted values if not needed
  cox_model_stripped$weights <- NULL         # Remove weights if not needed
  cox_model_stripped$offset <- NULL          # Remove offset if not needed
  return(cox_model_stripped)
}