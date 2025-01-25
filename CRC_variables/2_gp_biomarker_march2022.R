##File name: 2_gp_biomarkers_march2022.r
##Author: Hannah Harrison
##Last Edit: 30/03/2022
##Description: extracts biomarker data and converts to variables for use in colorectal cancer model, using lma framework, updated variables in march 2022

#temp code (also in config)
source("~/ACED-RREDD-EHR/sam/scripts/0_config.r")

# fpath_baseline <- "/store/biobank/health_records/ukb44651.csv"
fpath_gp <- "/store/biobank/health_records/gp_clinical.txt" #GP clinical data
# res_dir <- "~/df"
res_dir_tmp <- file.path(res_dir, "tmp")
dir.create(res_dir_tmp, recursive = TRUE)

#define functions
#source("~/ACED-RREDD-EHR/CRC_variables/dependencies/biomarker_manage_march2022.r")
#source("~/ACED-RREDD-EHR/CRC_variables/dependencies/gp_data_funcs.r")
#source("~/ACED-RREDD-EHR/CRC_variables/dependencies/biomarker_cleaning.r")

#codelist
fpath_codelists <- "~/ACED-RREDD-EHR/codelists" 
codelist_biomarkers <- fread(file.path(fpath_codelists, "bowel_cancer_specific", "CRC_codelist_biomarkers_formatted.csv"), colClasses = "character")
codelist_biomarkers <- rename(codelist_biomarkers, biomarker = my_group)

#get gp data
df_gp_clinical <- fread(fpath_gp, colClasses = "character")
df_gp_clinical <- df_gp_clinical %>% mutate(read_2 = na_if(read_2, ""), read_3 = na_if(read_3, ""))

#find all biomarker events
df_gp_biom <- gp_clinical_query(df_gp_clinical, codelist_biomarkers)
df_gp_biom <- value_extract(df_gp_biom)
df_gp_biom <- biomarker_basic_clean_CRC(df_gp_biom) ##cleans (trims and rescales) the biomarkers required for CRC analyses

##load baseline data (we need both age and sex), merge with biomarker and calculate dob and age at event
df_baseline <- fread(fpath_baseline, select = c("eid", "31-0.0", "34-0.0", "52-0.0"), colClasses = "character") # extracts baseline data (eid, sex, year of bith, month of birth)
names(df_baseline) <- make.names(c("eid", "sex", "yob", "mob"))
df_baseline$dob <- gen_birth_dates(df_baseline$yob, df_baseline$mob)
df_gp_biom_plus <- merge(df_gp_biom, df_baseline, by = "eid")
df_gp_biom_plus$event_age <- gen_event_age(df_gp_biom_plus$dob, df_gp_biom_plus$event_dt)
save(df_gp_biom_plus, file = file.path(res_dir_tmp, "df_gp_biom_values_plus_age_and_sex.Rdata"))

###check against gp registrations to get EHR dataset only (need to check this process with Sam)
# load(file.path(res_dir, "df_gp_reg.Rdata")) #loads dataframe gp_reg (hh only)
gp_reg <- read_parquet(file.path(res_dir, "df_gp_reg.parquet"))
gp_reg$eid <- as.character(gp_reg$eid)
df_EHR <- right_join(df_baseline, gp_reg, by = "eid")
save(df_EHR, file = file.path(res_dir_tmp, "df_EHR_biomarker.Rdata"))

# ##quick load
# base::load(file.path(res_dir_tmp, "df_gp_biom_values_plus_age_and_sex.Rdata")) #loads pre-extracted CRC biomarkers, with merged age and sex info
# base::load(file.path(res_dir_tmp, "df_EHR_biomarker.Rdata")) #loads pre-configure EHR eligible eids, with sex and age data
num_years_back <- 2 #this should be fixed at 2 years for biomarkers

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
  ##model variables for rasied CA125
  df_gp_biom_meta <- CA125_raised_vars(df_gp_biom_meta)
  df_gp_biom_meta <- subset(df_gp_biom_meta,  select = -sex)

  #find all individuals with at least 6 months primary care data in the correct timeperiod??
  df_EHR$cutoff_min <- df_EHR$dob + years(lma) - years(num_years_back) # lma date (cut-off)
  df_EHR$cutoff_max <- df_EHR$dob + years(lma) # lma date (cut-off)
  #drop from df_EHR if not registered at a gp at lms age (check with Sam!) - they need at least 6months right?
  df_EHR <- df_EHR[(df_EHR$cutoff_max >= df_EHR$gp_start_date + months(6)) & (df_EHR$cutoff_min <= df_EHR$gp_end_date - months(6))]
  df_EHR_biom <- right_join(df_gp_biom_meta, df_EHR, by = "eid")
  df_EHR_biom[is.na(df_EHR_biom)] <- 0
  df_EHR_biom <- df_EHR_biom %>% select("eid", "iron_def", "iron_def_measured", "inflam_all", "inflam_measured")
  save(df_EHR_biom, file = file.path(res_dir, paste0("biomarkers/df_gp_biom_lma_", lma, "_nback", ".Rdata")))
}

parallel::mclapply(as.list(lmas), mc.cores=16, function(lma) get_df_gp_biom_vars(lma))


#get_df_gp_biom_vars(60)
#load(file = file.path(res_dir, paste0("biomarkers/df_gp_biom_lma_60_nback", ".Rdata")))