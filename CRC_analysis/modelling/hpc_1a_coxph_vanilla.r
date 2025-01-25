#!/usr/bin/env Rscript
# =====================================================================================
# Script Title: Overall Cox PH model without selection
# Author: Sam Ip
# Description: Fits an overall Cox PH model without predictor selection to the full cohort without resampling.
# =====================================================================================

# Clear workspace ----
rm(list=ls())

# Initialise ----
args = commandArgs(trailingOnly=TRUE)
# args <- c("","","")
assign("tag", as.character(args[1]), envir = .GlobalEnv)  #"", _OnlySymptomaticPop, _OnlySymptomaticPop_SA_AllSympSansFatigue
assign("tag_Vision", as.character(args[2]), envir = .GlobalEnv) # _sansVision
assign("tag_LDpred", as.character(args[3]), envir = .GlobalEnv) # _LDpred

cat("\n\ntag, tag_Vision, tag_LDpred ......", tag, tag_Vision, tag_LDpred, "\n\n")

# File paths & initialisation ----
dir_ACED <- "~/ACED-CRC-modelling"
scripts_modelling_dir <- file.path(dir_ACED, "CRC_analysis", "modelling")
source(file.path(scripts_modelling_dir, "hpc_init_getcoxdat.r"))
source(file.path(scripts_modelling_dir, "hpc_modelling_fns.r"))
setDT(coxdat_stacked_allages_num) # Prepare data 
mdl <- "_1a_coxph_vanilla" # Model of interest -- overwrites default 1b_coxph_stepwise in hpc_config.r
dir.create(file.path(res_plots_dir_tagfree, glue::glue("overallmdl")), recursive = TRUE) # Set up result directory
cat("\nResults will be saved to ......", file.path(res_plots_dir_tagfree, glue::glue("overallmdl")), "\n\n")


# Fit model ----
# Scale 
scaling_params <- calculate_scaling_params(coxdat_stacked_allages_num)
coxdat_stacked_allages_num <- scale_with_traindata_params(coxdat_stacked_allages_num, setDT(scaling_params))
# Fit
mdl_vanilla <- fit_1a_coxph_vanilla_model(data=coxdat_stacked_allages_num)
cat("\nDONE Scale and Train model ......\n")
print(mdl_vanilla)

# Save model ----
mdl_vanilla <- strip_coxph(mdl_vanilla) # Reduces file size
saveRDS(mdl_vanilla, file = file.path(res_plots_dir_tagfree, glue::glue("overallmdl"), glue::glue("res_overall{mdl}{tag}{tag_LDpred}{tag_Vision}.rds")))

cat("\nDONE Saved overall model for tag=", tag, " ......", file.path(res_plots_dir_tagfree, glue::glue("overallmdl"), glue::glue("res_overall{mdl}{tag}{tag_LDpred}{tag_Vision}.rds")), "\n")
