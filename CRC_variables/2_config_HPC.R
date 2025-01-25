##File name: 0_config.R
##Author: Hannah Harrison
##Last Edit: 09/04/20204
##Description:using as config file to update biomarker variables on HPC (June 2023), install libraries, sets paths and loads functions, sets global variables
##configured for use on HPC


##LIBRARIES
#library(tidyverse) #didn't get working
#library(readxl) #didn't get working
library(tidyr)
library(dplyr)
library(data.table)
library(lubridate) #required local "timechange" install?
#library(remotes)
library(parallel)
library(survival)
library(stringr) #/lib64/libstdc++.so.6: version `CXXABI_1.3.8' not found, required local stringi install

##GLOBAL VARIABLES
cancer_reg_censor_date <- as.Date("2022-11-30") 
gp_study_end_date <- as.Date("2018-08-01", format = "%Y-%m-%d")#?
study_end_date <- as.Date("2018-08-01", format = "%Y-%m-%d")#?

gp_censor_date_p1 <- as.Date("2017-05-01")#england vision
gp_censor_date_p2 <- as.Date("2017-04-01")#scotland
gp_censor_date_p3 <- as.Date("2016-07-01")#england TPP
gp_censor_date_p4 <- as.Date("2017-08-01")#wales

gp_oldest_permitted <- as.Date("1985-01-01")
earliest_baseline <- as.Date("2006-01-01")

#PATHS TO DATA (uploaded by Joe)
#fpath_baseline_old_has_TDS <- "/store/biobank/health_records/ukb50372.csv" #before update, most recent doesn't include TDS?
fpath_gp_reg <- file.path(UCL_raw, "gp_registrations.txt")
fpath_gp <- file.path(UCL_raw, "gp_clinical.txt") #need to query externally
fpath_gp_scripts <- file.path(UCL_raw, "gp_scripts.txt") #need to query externally
fpath_first_deg_relatives <- file.path(CCGE_raw, "ukb_samples_post_genotyping_qc.csv")
fpath_key_fam <- file.path(UCL_raw,"cambridge_ucl_id_lookup.txt") #links CCGE and UCL biobank applications
fpath_withdrawn <- file.path(CCGE_raw, "withdraw_id", "exclude_id_20230420.csv")


##PATHS TO DATA (pre-processed uploaded by Hannah)
fpath_baseline <- file.path(HH_storage, "pheno_small.csv") #after update in April 2023 - UCL data
fpath_baseline_CCGE <- file.path(HH_storage, "pheno_small_CCGE.csv")  #after update in April 2023

##PATHS TO CODELISTS
fpath_codelists <- file.path(ACED_git_folder, "codelists")
fpath_ICD10_codes <- file.path(fpath_codelists,  "ICD10_cancer.csv")
fpath_ICD9_codes <- file.path(fpath_codelists,  "ICD9_cancer.csv")

##PATHS TO OTHER INPUTS
fpath_baseline_vars <- file.path(fpath_codelists, "baseline_vars_inc_HES_list.csv")

##LOAD FUNCTION FILES
fpath_deps <- file.path(ACED_git_folder, project_name, "dependencies")
source(file.path(fpath_deps, "gen_functions.R"))
source(file.path(fpath_deps, "wide_baseline_funcs.R"))
source(file.path(fpath_deps, "base_vars_functions.R"))
source(file.path(fpath_deps, "gp_data_funcs.R"))
source(file.path(fpath_deps, "biomarker_cleaning.R"))
source(file.path(fpath_deps, "colorectal_biomarker_covariates.R"))
source(file.path(fpath_deps, "biomarker_manage_march2022.R"))


##load data setting files
source(file.path(fpath_deps, "CRC_thresholds.R"))

##Configure key
df_key <- fread(fpath_key_fam) %>% rename(eid=UCL_ID)

##PATHS TO SAVE data for this project
fpath_save <- HH_data_temp_store

##super-landmark framework settings
lma_spacing <- 1
lmas <- seq(40,75, lma_spacing)

yrs_back <- 2
yrs_fwd <-  2



