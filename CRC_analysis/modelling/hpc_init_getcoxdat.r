#!/usr/bin/env Rscript

# =====================================================================================
# Script Title: Fitting an overall Cox PH model
# Author: Sam Ip
# Description: Fits an overall Cox PH model to the full cohort without resampling.
# =====================================================================================

suppressPackageStartupMessages(library(dplyr)) 
gc()

# Define Directories ----
dir_ACED <- "~/ACED-CRC-modelling"
setwd(file.path(dir_ACED, "CRC_analysis"))
scripts_modelling_dir <- file.path(dir_ACED, "CRC_analysis", "modelling")

# Load Configurations ----
source(file.path(scripts_modelling_dir, "hpc_config.r"))

# Load Dataset (coxdat_stacked_allages_num) ----
fpath_coxdat <- file.path(shared_dir_hpc, paste0("coxdat_stacked_allages_num_", date_str_tagfree, ".rds")); cat("\nfpath_coxdat......", fpath_coxdat, "\n")
coxdat_stacked_allages_num <- readRDS(fpath_coxdat) %>% as.data.frame() %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    dplyr::select(-any_of(c("nelaa")))

# Remove Vision people if applicable ----
if (grepl("sansVision", date_str)){
  cat("\nSTART rm Vision eids ......", dim(coxdat_stacked_allages_num), "\n")
  eids_Vision <- fread(file.path(shared_dir, "eids_Vision.csv"))$eid
  coxdat_stacked_allages_num <- coxdat_stacked_allages_num %>% dplyr::filter(! eid %in% eids_Vision)
  cat("\nDONE rm Vision eids ......", dim(coxdat_stacked_allages_num), "\n")
} 

# Swap PRS (to LDPred) if applicable ----
if (grepl("LDpred", date_str)){
  coxdat_stacked_allages_num <- coxdat_stacked_allages_num %>% dplyr::select(-gen_prs)
  df_genetic <- readRDS(file.path(root_privatedata_dir, glue::glue("df_genetic{tag_LDpred}.rds")))
  df_genetic <- df_genetic %>%
    rename_with(~ paste0("gen_", .), -eid) %>%
    mutate(eid = as.character(eid))
  coxdat_stacked_allages_num <- coxdat_stacked_allages_num %>%
    left_join(df_genetic, by = "eid")
  cat("\nSwapped PRS to LDPred.\n")
}


# Rename predictors to indicate type----
source(file.path(codelist_dir, "defn_data_sources.r"))

colnames_df <- names(coxdat_stacked_allages_num)
data.table::setnames(
  coxdat_stacked_allages_num,
  old = c(
    c("age", "age2"),
    names_core[!names_core %in% c("age", "age2")],
    names_medhist,
    colnames_df[grep( "^mrk_" , colnames_df)],
    names_lifestyle[!names_lifestyle %in% c("v0_tot_METs_wk")]
  ), 
  new = c(
    paste0("core_", c("age", "age2")),
    # paste0("age_", c("age", "age2")),
    sub("^[^_]*_", "core_", names_core[!names_core %in% c("age", "age2")]),
    sub("^[^_]*_", "medhist_", names_medhist),
    sub("^[^_]*_", "symp_", colnames_df[grep( "^mrk_" , colnames_df)]),
    sub("^[^_]*_", "lifestyle_", names_lifestyle[!names_lifestyle %in% c("v0_tot_METs_wk")])
    ))


# Filter symptomatic populations if applicable ----
if (grepl("OnlySymptomaticPop", tag)){
  coxdat_stacked_allages_num <- filter_coxdat_OnlySymptomaticPop_hpc(coxdat_stacked_allages_num, tag) #requires selected symptoms -- graph_1b_coxph_stepwise_bidirectional_selectsymptoms_agesex.rds
}
if (grepl("OnlySymptomaticPop_SA_AllSympSansFatigue", tag)){
    coxdat_stacked_allages_num <- filter_coxdat_SA_AllSympSansFatigue(coxdat_stacked_allages_num, tag)
}

# Collapse predictor levels due to small counts if applicable ----
coxdat_stacked_allages_num <- coxdat_stacked_allages_num %>%
  dplyr::select(-one_of(c(
    "medhist_DIB_T1D_ever", "symp_ABDO_LUMP_in_lb", "symp_RECTAL_MASS_in_lb",
    "symp_PELVIC_PAIN_freq_in_lb", "symp_WEIGHT_LOSS_new_onset", "symp_JAUNDICE_new_onset"
    )))  %>% 
  dplyr::mutate(
    core_ethn_selfrep = recode_factor(
        core_ethn_selfrep,
        "Chinese" = "Other/Unknown",
        "Missing" = "Other/Unknown",
        "Other" = "Other/Unknown"
        ),
    symp_ABDO_BLOAT_freq_in_lb = ifelse(symp_ABDO_BLOAT_freq_in_lb > 0, "1+", symp_ABDO_BLOAT_freq_in_lb),
    symp_ABDO_PAIN_freq_in_lb = ifelse(symp_ABDO_PAIN_freq_in_lb > 0, "1+", symp_ABDO_PAIN_freq_in_lb),
    symp_STOM_DIS_freq_in_lb =ifelse(symp_STOM_DIS_freq_in_lb > 0, "1+", symp_STOM_DIS_freq_in_lb)
    )
  cat("\nSmall predictor levels collapsed.\n")

# Finalise dataset ----

# Identify binary columns in the dataset
binarycols <- coxdat_stacked_allages_num %>% dplyr::select_if(
    function(col) {
        (n_distinct(col) == 2)
    })  %>% dplyr::select(!status) %>% names()

# Columns to turn into factor
factor_cols <- c(binarycols, "lifestyle_edu_quals") %>% unique()
coxdat_stacked_allages_num <- coxdat_stacked_allages_num %>% dplyr::mutate(across(all_of(factor_cols), function(col) {as.factor(col)}))
str(coxdat_stacked_allages_num)

# Load a dictionary of covariate names ----
# The dictionary file contains mappings for covariates, including:
# - "covariate_readable": Human-readable description of the covariate
# - "covariate_code": Original column name in the dataset
# - "covar_code2": Modified or standardized column name used in analyses
readable_names <- fread(file.path(codelist_dir, "dict_covariate_names.csv"))

# Add a new column "data_source" to classify each covariate based on its origin
readable_names <- readable_names %>% dplyr::mutate(data_source = case_when(
    covariate_code %in% names_core ~ "core",
    covariate_code %in% names_medhist ~ "medhist",
    covariate_code %in% names_lifestyle ~ "lifestyle",
    TRUE ~ sub("\\_.*", "", covariate_code)
))


# Identify eventless risk factor levels ----
# The eventless_riskfactor_levels function identifies categorical variables (factors)
# that have levels with no cases (where `status == 1`) in the dataset.
# These "eventless" levels are not informative for event-based models (e.g., survival analysis)
# and can be flagged for removal or adjustment.
eventless_cols <- eventless_riskfactor_levels(ind_symptomatic_only=ifelse(grepl("OnlySymptomaticPop", tag), TRUE, FALSE))


