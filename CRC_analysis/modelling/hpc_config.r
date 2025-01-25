# ========================================================================
# Script Title: Configuration and Setup 
# Author: Sam Ip
# ========================================================================

# Toggles for dataset ----
tag <-  get("tag", envir = .GlobalEnv) #"", _OnlySymptomaticPop, _OnlySymptomaticPop_SA_AllSympSansFatigue
if (exists("tag_Vision", envir = .GlobalEnv)) {tag_Vision <- get("tag_Vision", envir = .GlobalEnv)} # _sansVision
if (exists("tag_LDpred", envir = .GlobalEnv)) {tag_LDpred <- get("tag_LDpred", envir = .GlobalEnv)} # _LDpred
date_str <- glue::glue("20240616_EditedFerritin"); date_str_tagfree <- date_str #20230726_PRScsx_CRC_OkSkinNM
if (exists("tag_LDpred", envir = .GlobalEnv)) {date_str <- glue::glue("{date_str}{tag_LDpred}")} 
if (exists("tag_Vision", envir = .GlobalEnv)) {date_str <- glue::glue("{date_str}{tag_Vision}")} 
cat("date_str......", date_str, "\n")

mdl <- "1b_coxph_stepwise"
ctrl2case <- NA #5, NA

# Hyperparameters ----
nfolds_rsmp <- 10 #3, 5 -- inner and outer resampling folds

# Directory paths ----
dir_ACED <- "~/ACED-CRC-modelling"
root_privatedata_dir <- "~/rds/rds-ccge1-hdMXhK21vco/biobank/si_aced_crc"
cat("dir_ACED:", dir_ACED, "\nroot_privatedata_dir:", root_privatedata_dir, "\n")

# libraries ----
vec_libraries <- c("arrow", "bbotk", "checkmate", "data.table", "DiceKriging", "dplyr", "forestmodel",  
    "GGally", "ggplot2", "ggpubr", "gtools", "lubridate", "paradox", "ranger",  "readxl", "reshape2", 
    "stringi", "survival", "tidyverse", "zoo")
needed_libraries <- vec_libraries[!(vec_libraries %in% installed.packages()[,"Package"])]
cat("\n\nmissing but needed libraries......\n")
needed_libraries %>% print()

ls_library <- as.list(vec_libraries[!vec_libraries %in% needed_libraries])
lapply(ls_library, function(lib) {
    suppressPackageStartupMessages(library(lib, character.only = TRUE))
})

cat("\n\nPackages loaded ......\n")

# Study configuration ----
study_end_date <- as.Date("2020-01-01", format = "%Y-%m-%d")
num_yrs_back <- 2
num_yrs_fwd <- 2
num_yrs_comorb_lookback <- 6
lma_spacing <- 1
lmas <- seq(40,75, lma_spacing) 

# Directories ----
res_dir <- file.path(root_privatedata_dir, glue::glue("res/{date_str}_nback{num_yrs_back}_nfwd{num_yrs_fwd}"))
codelist_dir <- file.path(dir_ACED, "CRC_analysis", "codelists")
scripts_modelling_dir <- file.path(dir_ACED, "CRC_analysis", "modelling")
scripts_dataprep_dir <- file.path(dir_ACED, "CRC_analysis", "data_preparation")
shared_dir <- file.path("~/rds/rds-ccge1-hdMXhK21vco/biobank", "ACED-RREDD-EHR_colorectal")
shared_dir_hpc <- file.path("~/rds/rds-ccge1-hdMXhK21vco/biobank", "si_aced_crc")

res_plots_dir <- file.path(res_dir, "results")
res_plots_dir_tagfree <- sub(date_str, date_str_tagfree, res_plots_dir)

ls_dirs <- list(
    res_dir,
    file.path(res_dir, "results"),
    file.path(res_dir, "survobj"),
    file.path(res_dir, "biomarkers"),
    file.path(res_dir, "coxdat"),
    file.path(res_dir, "plots"),
    file.path(res_dir, "multimorb"),
    file.path(res_dir, "other_events"), 
    file.path(res_dir, "modifier_events"), 
    file.path(res_dir, "marker_events")
)

lapply(ls_dirs, dir.create, recursive = TRUE)


# Codelists ----
fpath_codelist_cancers <- file.path(codelist_dir, "codelist_cancers.csv")
fpath_codelist_ukb_baseline <- file.path(codelist_dir, "ukb_baseline_core.csv")

# Data filepaths ----
dpath_hpc_healthrecords <- "~/rds/rds-ccge1-hdMXhK21vco/biobank/health_records"
fpath_gp <- file.path(dpath_hpc_healthrecords, "gp_clinical.txt")
fpath_gp_reg <- file.path(dpath_hpc_healthrecords, "gp_registrations.txt")
fpath_baseline <- file.path(dpath_hpc_healthrecords, "ukb50372.csv") 
fpath_genetic_pc_ethn <- file.path("~/dfs", "jonathan_pc_ethn.parquet") 
fpath_genetic_prs <- file.path( 
    dir_ACED, "CRC_analysis", "codelists/genetic", 
    case_when(
        date_str %in% c("20220718", "20221212") ~ "crc_LDPred",
        grepl("LDpred", date_str) ~ "crc_LDPred",
        date_str %in% c("20230423_PRScsx_CRC", "20240616_EditedFerritin") ~ "crc_PRScsx_CRC",
        date_str == "20230423_best_COADREAD" ~ "crc_best_COADREAD",
        date_str == "20230423_PRS_CRC" ~ "crc_PRS_CRC",
        TRUE ~ "crc_huyghe"
    ), 
    "id_prs.txt") #id_prsivw
fpath_first_deg_relatives <- file.path(normalizePath(file.path(dpath_hpc_healthrecords, "..")), "phenotypes/ukb_samples_post_genotyping_qc.csv")
fpath_ukb_dob <- file.path(res_dir, "dob_from_baseline.parquet")
fpath_key_fam <- file.path(dpath_hpc_healthrecords, "cambridge_ucl_id_lookup.txt")

# Util scripts ----
dpath_hpc_CRCvariables <- file.path(dir_ACED, "CRC_variables/dependencies")
wrangling_functions <- file.path(scripts_modelling_dir, "0_wrangling_functions.r")
get_ukb_dob <- file.path(scripts_dataprep_dir, "0_get_ukb_dob.r")
biomarker_manage <- file.path(dpath_hpc_CRCvariables, "biomarker_manage_march2022.R")
gp_data_funcs <- file.path(dpath_hpc_CRCvariables, "gp_data_funcs.R")
gen_functions <- file.path(dpath_hpc_CRCvariables, "gen_functions.R")
base_vars_funcs <- file.path(dpath_hpc_CRCvariables, "base_vars_functions.R")
biomarker_clean_funcs <- file.path(dpath_hpc_CRCvariables, "biomarker_cleaning.R")
CRC_biomarker_funcs <- file.path(dpath_hpc_CRCvariables, "colorectal_biomarker_covariates.R")
CRC_threshold_defs <- file.path(dpath_hpc_CRCvariables, "CRC_thresholds.R")
mmscore_defs_and_funcs <- file.path(dpath_hpc_CRCvariables, "multimorb_score_info_and_funcs_v2.R")
wide_baseline_funcs <- file.path(dpath_hpc_CRCvariables, "wide_baseline_funcs.R")


lapply(
    list(
        wrangling_functions,
        get_ukb_dob,
        biomarker_manage,
        gp_data_funcs,
        gen_functions,
        base_vars_funcs,
        biomarker_clean_funcs,
        CRC_biomarker_funcs,
        CRC_threshold_defs,
        mmscore_defs_and_funcs,
        wide_baseline_funcs
    ),
    function(file_path) {
        cat("Sourcing:", file_path, "\n")
        source(file_path)}
)

# Display options ----
options(max.print = 1000)
# options(error = recover)

cat("\nDONE hpc_config ......\n")