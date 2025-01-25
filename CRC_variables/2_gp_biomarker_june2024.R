##File name: 2_gp_biomarker_june2024.R
##Author: Hannah Harrison
##Last Edit: 11/06/2024
##Description: extracts biomarker data and converts to variables for use in colorectal cancer model, using lma framework, following bug fix in June 2023, reran on HPC 
#(this should be run as a job - I split into two sections 1. query the clinical record with high memory/cpu usage, then 2. use multiple cores to run the lma loop)
#rm(list = ls())

##top level paths (to be removed before code release)
UCL_raw <- "~/rds/rds-ccge1-hdMXhK21vco/biobank/health_records"
CCGE_raw <- "~/rds/rds-ccge1-hdMXhK21vco/biobank/phenotypes"
HH_data_temp_store <- "~/rds/rds-ccge1-hdMXhK21vco/biobank/HH_biomarker_data"
SI_data_temp_store <- "~/rds/rds-ccge1-hdMXhK21vco/biobank/ACED-RREDD-EHR_colorectal"

HH_storage <- "~/rds/hpc-work"
ACED_git_folder <- "~/ACED_project/ACED-RREDD-EHR" 
project_name <- "CRC_variables"

#temp code (also in config)
#source("~/ACED-RREDD-EHR/sam/scripts/0_config.r")
source(file.path(ACED_git_folder, project_name, "2_config_HPC.R"))

#codelist
codelist_biomarkers <- fread(file.path(fpath_codelists, "bowel_cancer_specific", "CRC_codelist_biomarkers_formatted.csv"), colClasses = "character")
codelist_biomarkers <- rename(codelist_biomarkers, biomarker = my_group)

#get gp data (run as job)
df_gp_clinical <- fread(fpath_gp, colClasses = "character")
df_gp_clinical <- df_gp_clinical %>% mutate(read_2 = na_if(read_2, ""), read_3 = na_if(read_3, ""))

#find all biomarker events (run as job)
df_gp_biom <- gp_clinical_query(df_gp_clinical, codelist_biomarkers)
df_gp_biom <- value_extract(df_gp_biom)
save(df_gp_biom, file = file.path(HH_data_temp_store, "df_gp_all_biom_values_RAW.Rdata")) # done as a job to here


#load(file = file.path(HH_data_temp_store, "df_gp_all_biom_values_RAW.Rdata"))#loads df_gp_biom
df_gp_biom <- biomarker_basic_clean_CRC(df_gp_biom) ##cleans (trims and rescales) the biomarkers required for CRC analyses
 
##load baseline data (we need both age and sex), merge with biomarker and calculate dob and age at event
df_baseline <- fread(fpath_baseline, select = c("eid", "31-0.0", "34-0.0", "52-0.0"), colClasses = "character") # extracts baseline data (eid, sex, year of bith, month of birth)
names(df_baseline) <- make.names(c("eid", "sex", "yob", "mob"))
df_baseline$dob <- gen_birth_dates(df_baseline$yob, df_baseline$mob)
df_gp_biom_plus <- merge(df_gp_biom, df_baseline, by = "eid")
df_gp_biom_plus$event_age <- gen_event_age(df_gp_biom_plus$dob, df_gp_biom_plus$event_dt)
save(df_gp_biom_plus, file = file.path(HH_data_temp_store, "df_gp_all_biom_values_plus_age_and_sex.Rdata"))

# 
# ###check against gp registrations to get EHR dataset only (need to check this process with Sam)
load(file.path(SI_data_temp_store, "df_gp_reg.Rdata")) #loads dataframe gp_reg (hh only)
#gp_reg <- read_parquet(file.path(res_dir, "df_gp_reg.parquet"))
gp_reg$eid <- as.character(gp_reg$eid)
df_EHR <- right_join(df_baseline, gp_reg, by = "eid")
save(df_EHR, file = file.path(HH_data_temp_store, "df_EHR_biomarker.Rdata"))

# ##quick load (ran as a job to this point)
base::load(file.path(HH_data_temp_store, "df_gp_all_biom_values_plus_age_and_sex.Rdata")) #loads pre-extracted CRC biomarkers, with merged age and sex info
base::load(file.path(HH_data_temp_store, "df_EHR_biomarker.Rdata")) #loads pre-configure EHR eligible eids, with sex and age data

get_df_gp_biom_vars <- function(lma) {
  #limit EHR eids to individuals with included data in lma window
  # limit biomarkers events to dates within window ----
  df_gp_biom_plus$window_min <-df_gp_biom_plus$dob + years(lma) - years(num_years_back)
  df_gp_biom_plus$window_max <- df_gp_biom_plus$dob + years(lma) # lma date
  df_gp_biom_window <-df_gp_biom_plus[(df_gp_biom_plus$event_dt >= df_gp_biom_plus$window_min) & (df_gp_biom_plus$event_dt <= df_gp_biom_plus$window_max)]
  ##we should also drop more than 2 years pre-baseline assessment? - need to check this

  #identify number of measurements and abnormal results for each biomarker per eid
  df_gp_biom_meta <- biomarker_meta_CRC(df_gp_biom_window) #by eid, number of tests and number of abnormal results (both per test type)

  #model variables for iron deficiency
  df_gp_biom_meta <- iron_def_vars(df_gp_biom_meta)
  #model variables for inflammation
  df_gp_biom_meta <- inflam_vars(df_gp_biom_meta)
  df_gp_biom_meta <- subset(df_gp_biom_meta,  select = -sex)

  #find all individuals with at least 6 months primary care data in the correct timeperiod??
  df_EHR$cutoff_min <- df_EHR$dob + years(lma) - years(num_years_back) # lma date (cut-off)
  df_EHR$cutoff_max <- df_EHR$dob + years(lma) # lma date (cut-off)
  #drop from df_EHR if not registered at a gp at lms age (check with Sam!) - they need at least 6months right?
  df_EHR <- df_EHR[(df_EHR$cutoff_max >= df_EHR$gp_start_date + months(6)) & (df_EHR$cutoff_min <= df_EHR$gp_end_date - months(6))]
  df_EHR_biom <- right_join(df_gp_biom_meta, df_EHR, by = "eid")
  df_EHR_biom[is.na(df_EHR_biom)] <- 0
  df_EHR_biom <- df_EHR_biom %>% select("eid", "iron_def", "iron_def_measured", "inflam_all", "inflam_measured")
  df_EHR_biom$lma <- lma
  return(df_EHR_biom)
  #save(df_EHR_biom, file = file.path(res_dir, paste0("biomarkers/df_gp_biom_lma_", lma, "_nback", ".Rdata")))
}

num_years_back <- yrs_back #this should be fixed at 2 years for biomarkers
df_lma_EHR_biom <- mclapply(as.list(lmas), mc.cores=4, function(lma) get_df_gp_biom_vars(lma)) %>% rbindlist(fill = TRUE)
save(df_lma_EHR_biom, file = file.path(HH_data_temp_store, "df_biomarker_for_lma.Rdata"))

load(file = file.path(HH_data_temp_store, "df_biomarker_for_lma.Rdata"))
load(file = file.path(SI_data_temp_store, "coxdat_stacked_allages_num_20230726_PRScsx_CRC_OkSkinNM.Rdata"))
df_lma_EHR_biom <- df_lma_EHR_biom %>% rename(age = lma)
df_joined <- left_join(coxdat_stacked_allages_num, df_lma_EHR_biom, by = c("eid", "age"))

wid <- as.data.frame(fread(fpath_withdrawn))
df_key <- fread(fpath_key_fam)
wid_plus_key <- left_join(wid, df_key, by = c("x" = "Cambridge_ID"))
df_joined_drop_wid <- df_joined %>% filter(!eid %in% as.character(wid_plus_key$UCL_ID))


