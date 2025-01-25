#!/usr/bin/env Rscript
# =====================================================================================
# Script Title: Overall Cox PH model with group-lasso selection
# Author: Sam Ip
# Description: Fits an overall Cox PH model with group-lasso selection 
#              (grouping predictor levels during selection) to the full cohort without resampling.
# =====================================================================================

# Clear workspace and set global options ----
rm(list=ls())
# options(future.globals.maxSize = 5 * 1024^3); cat("future.globals.maxSize is set to: ", getOption("future.globals.maxSize"), "\n")

# Initialise ----
args = commandArgs(trailingOnly=TRUE)
# args <- c("","","")
assign("tag", as.character(args[1]), envir = .GlobalEnv)  #"", _OnlySymptomaticPop, _OnlySymptomaticPop_SA_AllSympSansFatigue
assign("tag_Vision", as.character(args[2]), envir = .GlobalEnv) # _sansVision
assign("tag_LDpred", as.character(args[3]), envir = .GlobalEnv) # _LDpred

cat("\n\ntag, tag_Vision, tag_LDpred ......", tag, tag_Vision, tag_LDpred, "\n\n")

# File paths & initialisation ----
library(grpreg)
dir_ACED <- "~/ACED-CRC-modelling"
scripts_modelling_dir <- file.path(dir_ACED, "CRC_analysis", "modelling")
source(file.path(scripts_modelling_dir, "hpc_init_getcoxdat.r"))
source(file.path(scripts_modelling_dir, "hpc_modelling_fns.r"))
setDT(coxdat_stacked_allages_num) # Prepare data 
mdl <- "_2_lasso" # Model of interest -- overwrites default 1b_coxph_stepwise in hpc_config.r
dir.create(file.path(res_plots_dir_tagfree, glue::glue("overallmdl")), recursive = TRUE)
cat("\nResults will be saved to ......", file.path(res_plots_dir_tagfree, glue::glue("overallmdl")), "\n\n")


# Fit model ----
# Scale 
scaling_params <- calculate_scaling_params(coxdat_stacked_allages_num)
coxdat_stacked_allages_num <- scale_with_traindata_params(coxdat_stacked_allages_num, setDT(scaling_params))

# Get predictor levels groupings
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
grp_featlvls <- generate_group_vector(coxdat_stacked_allages_num)

# Fit
mdl_lasso <- fit_2_grouplasso_model(df=coxdat_stacked_allages_num, grp_featlvls, lambda_min_or_1se="1se", nlambda=100, time_col="futime", event_col="status", id_col="eid")
cat("\nDONE Scale and Train model ......\n")
print(mdl_lasso)

# Save model ----
mdl_lasso <- strip_coxph(mdl_lasso) # Reduces file size
saveRDS(mdl_lasso, file = file.path(res_plots_dir_tagfree, glue::glue("overallmdl"), glue::glue("res_overall{mdl}{tag}{tag_LDpred}{tag_Vision}.rds")))

cat("\nDONE Saved overall model for tag=", tag, " ......", file.path(res_plots_dir_tagfree, glue::glue("overallmdl"), glue::glue("res_overall{mdl}{tag}{tag_LDpred}{tag_Vision}.rds")), "\n")
