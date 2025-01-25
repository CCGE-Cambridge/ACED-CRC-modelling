#!/usr/bin/env Rscript
# =====================================================================================
# Script Title: Bootstrap resampling for Cox PH model with no or group-LASSO seelction
# Author: Sam Ip
# Description:
# This script performs the following steps:
# 1. Bootstrap-resamples the full cohort (ID-aware, s.t. all records of an individual are either fully included or fully excluded -- sampling by ID, not by record).
# 2. Splits data into train-test sets (ID-aware s.t. each ID is exclusive to either the train or test set).
# 3. Fits models for each predictor coalition on the training set.
# 4. Tests models on the corresponding test set.
# =====================================================================================


# =============================================================================
# Setup ----
# =============================================================================
rm(list=ls())
options(future.globals.maxSize = 5 * 1024^3); cat("future.globals.maxSize is set to: ", getOption("future.globals.maxSize"), "\n")

# Initialise----
args = commandArgs(trailingOnly=TRUE)
# args <- c(5,"_OnlySymptomaticPop", 1, 200, "", "", "_2_lasso")
n_cores <- as.numeric(args[1])
assign("tag", as.character(args[2]), envir = .GlobalEnv)  #"", _OnlySymptomaticPop, _OnlySymptomaticPop_SA_AllSympSansFatigue
idx_bootsmp_start <- as.numeric(args[3])
idx_bootsmp_end <- as.numeric(args[4])
assign("tag_Vision", as.character(args[5]), envir = .GlobalEnv) # _sansVision
assign("tag_LDpred", as.character(args[6]), envir = .GlobalEnv) # _LDpred
cat("\n\nn_cores: ", n_cores, "\ntag: ", tag, ", tag_Vision: ", tag_Vision, ", tag_LDpred: ", tag_LDpred, "\nidx_bootsmp_start-idx_bootsmp_end:  ", idx_bootsmp_start,"-", idx_bootsmp_end,  "\n\n")

set.seed(137) 

# Load external scripts ----
dir_ACED <- "~/ACED-CRC-modelling"
scripts_modelling_dir <- file.path(dir_ACED, "CRC_analysis", "modelling")
source(file.path(scripts_modelling_dir, "hpc_init_getcoxdat.r"))
source(file.path(scripts_modelling_dir, "hpc_modelling_fns.r"))
source(file.path(codelist_dir, "defn_data_sources.r"))

# Specify predictor selection method ----
mdl <- as.character(args[7]) # Overwrite default bidirectional selection from in hpc_config.r -- choose from _1a_coxph_vanilla, _2_lasso

setDT(coxdat_stacked_allages_num)

# Define directories ----
res_plots_dir <- "~/rds/hpc-work/BMJOnc_reviewerresp/tmp_boot_res" #HERE 
dpath_tmp_boot_res <- file.path(res_plots_dir, glue::glue("bootsmps_combos{tag}"))
if (exists("tag_LDpred", envir = .GlobalEnv)) {dpath_tmp_boot_res <- glue::glue("{dpath_tmp_boot_res}{tag_LDpred}")} 
if (exists("tag_Vision", envir = .GlobalEnv)) {dpath_tmp_boot_res <- glue::glue("{dpath_tmp_boot_res}{tag_Vision}")} 
dir.create(dpath_tmp_boot_res, recursive = TRUE)
cat("\nResults will be saved to ......", dpath_tmp_boot_res, "\n\n")

# =============================================================================
# Define Predictor Coalitions ----
# =============================================================================
data_sources <- c("core", "gen", "medhist", "symp", "biom", "lifestyle")
combos <- paste(sort(data_sources), collapse = "_")

# Identify remaining combinations to process (fmt: {bootsmpnum}_{coalition}) ----
csv_files <- list.files(dpath_tmp_boot_res, pattern = "^cindex_(.*).csv$")
csv_models <- sub("^cindex_(.*)\\.csv$", "\\1", csv_files)
required_combos <- expand.grid(idx_bootsmp_start:idx_bootsmp_end, combos)
required_combos <- paste(required_combos$Var1, required_combos$Var2, sep = "_")
bootsmpnum_coalition_to_process <- setdiff(required_combos, csv_models)

cat("\n\nIDENTIFIED bootsmpnum * combos to process: ", length(bootsmpnum_coalition_to_process), "......\n\n")


# =============================================================================
# Functions ----
# =============================================================================
# Extract bootstrap samples ----
get_bootstrap_sample_i <- function(IDs_traintest_bootstrap_sample, coxdat_stacked_allages_num) {
  train_data <- coxdat_stacked_allages_num[IDs_traintest_bootstrap_sample$train, on = .(eid), nomatch = 0]
  test_data <- coxdat_stacked_allages_num[IDs_traintest_bootstrap_sample$test, on = .(eid), nomatch = 0]
  cat("\nDONE Extract bootstrap samples ......\n")
  return(list(train_data = train_data, test_data = test_data))
}
is_null_model <- function(model_formula) {
  as.character(model_formula)[3] == "1"
}

# Get predictor levels groupings for group-LASSO selection option (same as in hpc_2_lasso.r)
generate_group_vector <- function(data) {
  feature_names <- setdiff(names(data), c("futime", "status", "eid"))
  model_matrix <- model.matrix(as.formula(paste("~", paste(feature_names, collapse = "+"), "- 1")), data = data)
  group_vector <- integer(ncol(model_matrix))
  
  current_group <- 1
  
  # Group gen_pc predictors
  gen_pc_columns <- grep("^gen_pc\\d+$", colnames(model_matrix))
  if (length(gen_pc_columns) > 0) {
    group_vector[gen_pc_columns] <- current_group
    current_group <- current_group + 1
  }
  
  # Group other predictors
  for (base_name in feature_names) {
    # Skip gen_pc predictors as they are already grouped
    if (grepl("^gen_pc\\d+$", base_name)) next
    
    if (is.factor(data[[base_name]]) || is.character(data[[base_name]])) {
      matched_columns <- grep(paste0("^", base_name), colnames(model_matrix))
      group_vector[matched_columns] <- current_group
      current_group <- current_group + 1
    } else {
      group_vector[which(colnames(model_matrix) == base_name)] <- current_group
      current_group <- current_group + 1
    }
  }
  
  return(group_vector)
}

# Fit models ----
process_bootstrap_sample <- function(bootsmpnum_coalition, coxdat_stacked_allages_num, tag, mdl) {
  # bootsmpnum_coalition <- bootsmpnum_coalition_to_process[1]
  cat("\nSTART process_bootstrap_sample ......", bootsmpnum_coalition, "\n")

  # Parse bootsmpnum_coalition
  parts <- strsplit(bootsmpnum_coalition, "_")[[1]]
  i <- as.numeric(parts[1])
  combo <- paste(parts[-1], collapse = "_")

  # Load bootstrap sample train & test IDs
  cat("\n IDs: ", file.path(dpath_IDs_200bootstrapsamples_indiv, glue::glue("bootstrap_sample_{i}.rds")), "\n") # Generated by hpc_bootstrap_samples.r
  bootsamp_i <-  readRDS(file.path(dpath_IDs_200bootstrapsamples_indiv, glue::glue("bootstrap_sample_{i}.rds")))

  # Get train & test datasets of the bootstrap sample -- use DT for speed
  IDs_traintest_bootstrap_sample <- list(
    train = data.table(eid = as.character(bootsamp_i$train)),
    test = data.table(eid = as.character(bootsamp_i$test))
  )

  ls_bootstrap_train_test <- get_bootstrap_sample_i(IDs_traintest_bootstrap_sample, coxdat_stacked_allages_num)
  train_data <- ls_bootstrap_train_test$train_data; 
  test_data <- ls_bootstrap_train_test$test_data
  cat("\nDONE Get train & test datasets of the bootstrapped sample ......\n")


  cat("\nSTART combo ......", combo, "\n")
  fpath_i_combo <- file.path(dpath_tmp_boot_res, glue::glue("cindex_{bootsmpnum_coalition}.csv"))
  
  if(file.exists(fpath_i_combo)){return(NULL)}

  # Get coalition-specific train & test datasets
  predictor_prefixes <- strsplit(combo, "_")[[1]]
  train_data_coalition <- train_data %>%
    dplyr::select(matches(glue::glue("^({paste(c('eid', 'futime', 'status', predictor_prefixes), collapse='|')})"))) # keep only predictors in specified coalition 
  test_data_coalition <- test_data %>%
    dplyr::select(matches(glue::glue("^({paste(c('eid', 'futime', 'status', predictor_prefixes), collapse='|')})"))) 

  # Scale data
  scaling_params <- calculate_scaling_params(train_data_coalition)
  train_data_coalition <- scale_with_traindata_params(train_data_coalition, setDT(scaling_params))
  test_data_coalition <- scale_with_traindata_params(test_data_coalition, setDT(scaling_params))

  # *** TRAIN ***
  source(file.path(scripts_modelling_dir, "hpc_modelling_fns.r"), local=TRUE)

  if (mdl=="_1a_coxph_vanilla"){
    cat("\nFitting model ......", mdl, "\n\n")
    mdl_fit <- fit_1a_coxph_vanilla_model(data=train_data_coalition)
  } else if (mdl=="_2_lasso"){
    cat("\nFitting model ......", mdl, "\n\n")
    library(grpreg)
    grp_featlvls <- generate_group_vector(train_data_coalition)
    mdl_fit <- fit_2_grouplasso_model(df=train_data_coalition, grp_featlvls, lambda_min_or_1se="1se", nlambda=100, time_col="futime", event_col="status", id_col="eid")
  }
  cat("\nDONE Train model on train_data_coalition ......\n")
  
  # *** TEST ***
  c_index <- survival::concordance(mdl_fit, newdata = test_data_coalition)$concordance
  cat("\nC-index ......\n"); print(c_index)
  cat("\nDONE Predict on test_data_coalition ......\n")
  
  # Save model and metrics ----
  
  # Format model
  if (!is_null_model(mdl_fit$formula)){
    # Process the model using broom::tidy
    tidy_model <- broom::tidy(mdl_fit, conf.int = TRUE, exponentiate = TRUE)
    tidy_model <- tidy_model %>%
        dplyr::select(term, estimate, conf.low, conf.high, p.value)

    # Add the bootsmpnum_coalition column
    tidy_model <- tidy_model %>%
        mutate(bootsmpnum_coalition = bootsmpnum_coalition)%>%
        dplyr::select(bootsmpnum_coalition, everything())  
  } else (tidy_model <- data.frame())

  # Format metric
  cindex_result <- data.frame(
    bootsmpnum_coalition = bootsmpnum_coalition,
    cindex = c_index
  )

  # Save
  temp_cindex_file <- file.path(dpath_tmp_boot_res, paste0("cindex_", bootsmpnum_coalition, ".csv"))
  temp_tidy_model_file <- file.path(dpath_tmp_boot_res, paste0("model_", bootsmpnum_coalition, ".csv"))
  write.table(cindex_result, file = temp_cindex_file, sep = ",", col.names = TRUE, row.names = FALSE)
  write.table(tidy_model, file = temp_tidy_model_file, sep = ",", col.names = TRUE, row.names = FALSE)

  cat(glue::glue("\nDONE {bootsmpnum_coalition} ......\n\n"))

}

# =============================================================================
# Main Execution ----
# =============================================================================
future::plan(future::multisession, workers = n_cores)
future.apply::future_lapply(bootsmpnum_coalition_to_process, function(bootsmpnum_coalition) {
    cat(glue::glue("\nSTART {bootsmpnum_coalition} ......\n"))
    process_bootstrap_sample(bootsmpnum_coalition, coxdat_stacked_allages_num, tag, mdl)
    cat(glue::glue("\nDONE {bootsmpnum_coalition} ......\n"))
})

