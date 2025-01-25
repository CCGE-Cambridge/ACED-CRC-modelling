#!/usr/bin/env Rscript
# =====================================================================================
# Script Title: Define "Symptomatic"
# Author: Sam Ip
# Description: Select symptoms for defining "symptomatic" using stepwise Cox regression.
# =====================================================================================

# Clear workspace and set global options ----
rm(list=ls())
options(future.globals.maxSize = 5 * 1024^3); cat("future.globals.maxSize is set to: ", getOption("future.globals.maxSize"), "\n")

# Initialise ----
args = commandArgs(trailingOnly=TRUE)
# args <- c("")
assign("tag", as.character(args[1]), envir = .GlobalEnv)  #"", _OnlySymptomaticPop, _OnlySymptomaticPop_SA_AllSympSansFatigue
cat("\n\ntag ......", tag, "\n\n")

# File paths & initialisation ----
dir_ACED <- "~/ACED-CRC-modelling"
scripts_modelling_dir <- file.path(dir_ACED, "CRC_analysis", "modelling")
source(file.path(scripts_modelling_dir, "hpc_init_getcoxdat.r"))
source(file.path(scripts_modelling_dir, "hpc_modelling_fns.r"))
setDT(coxdat_stacked_allages_num) # Prepare data 
dir.create(file.path(res_plots_dir, glue::glue("define_symptomatic")), recursive = TRUE) # Set up result directory
cat("\nDONE Set results dir ......", file.path(res_plots_dir, glue::glue("define_symptomatic")), "\n\n")

# Restrict predictors to age (age2), sex, and symptoms only ----
vars <- paste0("^(", paste(c("eid", "futime", "status", "symp", "core_age", "core_sex"), collapse="|"), ")")
coxdat_stacked_allages_num <- dplyr::select(coxdat_stacked_allages_num, matches(vars)) # keep only predictors in specified coalition 
str(coxdat_stacked_allages_num)

# Fit model ----
# Scale 
scaling_params <- calculate_scaling_params(coxdat_stacked_allages_num)
coxdat_stacked_allages_num <- scale_with_traindata_params(coxdat_stacked_allages_num, setDT(scaling_params))
# Fit 
mdl_bidirselect <- fit_1b_coxph_bidirselect_model(data=coxdat_stacked_allages_num)
cat("\nDONE Scale and Train model ......\n")
print(mdl_bidirselect)

# Save model ----
mdl_bidirselect <- strip_coxph(mdl_bidirselect) # Reduces file size
saveRDS(mdl_bidirselect, file = file.path(res_plots_dir, glue::glue("define_symptomatic"), glue::glue("res_define_symptomatic{tag}.rds")))

cat("\nDONE Saved overall model for tag=", tag, " ......", file.path(res_plots_dir, glue::glue("define_symptomatic"), glue::glue("res_define_symptomatic{tag}.rds")), "\n")


# Extract selected symptoms ----
selected_symptoms <- grep("^symp_", names(mdl_bidirselect$coefficients), value = TRUE); 
cat("\n\nselected_symptoms ......", selected_symptoms, "\n")
selected_symptoms <- sub("1$", "", selected_symptoms) # Clean symptom names
cat("\n\n CLEANED selected_symptoms ......", selected_symptoms, "\n")

# Save selected symptoms ----
saveRDS(selected_symptoms, file.path(res_plots_dir, glue::glue("define_symptomatic"), glue::glue("symps_define_symptomatic{tag}.rds")))
