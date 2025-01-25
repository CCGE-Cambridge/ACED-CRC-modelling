##File name: 2_gp_multimorb_HH_May22.r
##Authors: Sam Ip and Hannah Harrison
##Last Edit: 15/11/2022
##Description: identifies events from EHR (gp clinical)  within lma framework
##these events must have no time restriction other than happened before date of interest 
##and no additional formatting or conditions are applied (e.g. event happened ever)
root_dir <- "~/ACED-RREDD-EHR/sam"
script_dir <- file.path(root_dir, "scripts")
source(file.path(script_dir, "0_config.r"))


# =================================  PARAMS  ==================================
# lma <- 40
# ===============================  read in ====================================
# source("~/ACED-RREDD-EHR/CRC_variables/dependencies/gp_data_funcs.r")
# source("~/ACED-RREDD-EHR/CRC_variables/dependencies/gen_functions.r")
# source("~/ACED-RREDD-EHR/CRC_variables/dependencies/multimorb_score_info_and_funcs_v2.r")
my_conditions <- "~/ACED-RREDD-EHR/codelists/mm_score/cam_comorb_score.csv"

#required codelists
df_comorbs <- fread(my_conditions) #list of included comorbidities
fpath_codelists <- "~/ACED-RREDD-EHR/codelists"
df_codes_clinical <- fread(file.path(fpath_codelists, "mm_score", "cam_comorb_clinical.csv"), colClasses = "character")
df_codes_script <- fread(file.path(fpath_codelists, "mm_score", "cam_comorb_meds.csv"), colClasses = "character")

#load each dataset and query with respective codelists
df_gp_clinical <- fread(fpath_gp, colClasses = "character")
df_gp_clinical <- df_gp_clinical %>% mutate(read_2 = na_if(read_2, ""), read_3 = na_if(read_3, ""))
df_clinical_diags <- gp_clinical_query(df_gp_clinical, df_codes_clinical)

df_gp_scripts <- fread(fpath_gp_med, colClasses = "character")
df_gp_scripts <- df_gp_scripts %>% mutate(read_2 = na_if(read_2, ""),
                                          bnf_code = na_if(bnf_code, ""),
                                          dmd_code = na_if(dmd_code, "NA"))
df_script_events <- gp_script_query(df_gp_scripts, df_codes_script)
df_script_events <- df_script_events %>% select(eid, event_dt, med_type, Type)

#load a version of baseline data which included pre-calculated dob
base::load(file.path(res_dir, "df_dob.Rdata"))
# load("~/df/df_dob.Rdata") #loads pre-calculated dob (hh only)

#check against gp registrations to get EHR dataset only 
##?sam do you do this later on when you combine in lma structure...
#fpath_gp_reg <- "/store/biobank/health_records/gp_registrations.txt"
#source("~/ACED-RREDD-EHR/hannah/sam_scripts_mod/1_gp_reg_v0_nov0.r")
# load("~/df/df_gp_reg.Rdata") #loads dataframe gp_reg (hh only)
gp_reg <- read_parquet(file.path(res_dir, "df_gp_reg.parquet"))
gp_reg$eid <- as.character(gp_reg$eid)
df_EHR <- right_join(df_dob, gp_reg, by = "eid")

#remove data provider column (if present) to avoid nasty merge bug later
df_EHR <- df_EHR %>% select(-any_of(c("data_provider"))) 


get_df_gp_multimorb <- function(lma) {
    #lma<- 60
  #use same condition for gp records (at least 6 months within lb) as for other variables - not a requirement for mm_score in general, but reasonable for model development
  df_EHR$cutoff_min <- df_EHR$dob + years(lma) - years(2) # lma date (cut-off)
  df_EHR$cutoff_max <- df_EHR$dob + years(lma) # lma date (cut-off)
  df_EHR <- df_EHR[(df_EHR$cutoff_max >= df_EHR$gp_start_date + months(6)) & (df_EHR$cutoff_min <= df_EHR$gp_end_date - months(6))]

  # limit events to dates before current lma age ----
  df_clinical_diags_before <- events_relative_to_date(df_clinical_diags, df_EHR, "cutoff_max")
  df_script_events_before <- events_relative_to_date(df_script_events, df_EHR, "cutoff_max")

  ##extract clinical events and pivot to wide format (see "multimorb_score_info_and_funcs.r" for some of these variable definitions)
  df_ckd_only <- define_CKD(df_clinical_diags_before) #apply specific conditions for kidney disease
  df_can_only <- define_CAN(df_clinical_diags_before) #apply specific conditions for cancer diagnoses
  df_clinical_others <- df_clinical_diags_before %>% 
                            filter((Comorb %in% comorb_ever) | 
                                  ((Comorb %in% comorb_in_last_year) & (event_dt >= my_date_minus_1_year))) %>%
                                distinct(eid, Comorb, .keep_all = FALSE) 
  df_clinical_all <- rbind(df_clinical_others, df_ckd_only, df_can_only) %>%
                        mutate(dummy = 1) %>% 
                            pivot_wider(names_from = "Comorb", values_from = "dummy",  names_prefix = "clinic_", values_fill = 0)

  #extract prescription events and pivot to wide format (see "multimorb_score_info_and_funcs.r" for some of these variable definitions)
  df_script_all <- df_script_events_before %>% 
                    filter(event_dt >= my_date_minus_1_year | med_type == "SCZ") %>%
                        add_count(eid, med_type, name = "num_scripts") %>%
                            distinct(eid, med_type, .keep_all = TRUE) %>% 
                                select(eid, med_type, num_scripts)
  df_script_all <- df_script_all %>% 
                    pivot_wider(names_from = "med_type", values_from = "num_scripts",  names_prefix = "script_", values_fill = 0)

   ##combine clinical and script with the full EHR database events
   df_EHR_events <- right_join(df_clinical_all, df_EHR, by = "eid")
   df_EHR_events <- right_join(df_script_all, df_EHR_events, by = "eid")
   df_EHR_events[is.na(df_EHR_events)] <- 0

   ##calculate mm age categories (note I don't include the under 40 age categories)
    df_EHR_events$age <- gen_event_age(df_EHR_events$dob, df_EHR_events$assessment_date)
    df_EHR_events <- gen_mm_age_categories(df_EHR_events)

    ##define 37 comorbidites needed for the model
    df_mm_covars <- define_37_mm_covariates(df_EHR_events)

    ##create and save consulation model version (37 conditions) 
    #df_mm_consultation_long <- mm_score_HR(df_mm_covars, df_cam_multi_morb, "cons_NB_37_RR") %>% select("eid", "residual_fit") %>% rename(mm_score_res_fit = residual_fit)
    #save(df_mm_consultation_long, file = file.path(res_dir, paste0("multimorb/df_gp_multimorb_lma" , lma, ".Rdata"))) #test save
    df_mm_hosp_long <- mm_score_HR(df_mm_covars, df_cam_multi_morb, "hospital_37_HR") %>% select("eid", "prognostic_index") %>% rename(mm_score_PI = prognostic_index)
    save(df_mm_hosp_long, file = file.path(res_dir, paste0("multimorb/df_gp_multimorb_lma" , lma, ".Rdata"))) #test save
}

# get_df_gp_multimorb(60) #test
# df_mm_consultation_long

# cat("get_df_gp_multimorb(40)......\n")
# df_mm_consultation_long get_df_gp_multimorb(40) #test

# base::load(file.path(res_dir, "multimorb", paste0("df_gp_multimorb_lma" , 40, ".Rdata"))) %>% get() %>% head() %>% print()
parallel::mclapply(as.list(lmas), mc.cores=16, function(lma) {
    get_df_gp_multimorb(lma)
    cat("Done ", lma, "\n")
    })

