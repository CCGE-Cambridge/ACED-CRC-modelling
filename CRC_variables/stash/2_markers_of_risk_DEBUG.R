##File name: 2_markers_of_risk.r
##Authors: Hannah Harrison
##Last Edit: 28/02/2023
##Description: identifies markers of risk from EHR (gp_clincial and gp_scripts)  within lma framework
##events in lookback period, with other conditions added in some cases
source("~/ACED-RREDD-EHR/sam/scripts/0_config.r")

#load dependencies
#source("~/ACED-RREDD-EHR/CRC_variables/dependencies/gp_data_funcs.r")

##get codelists 
fpath_codelists <- "~/ACED-RREDD-EHR/codelists"
codelist_markers <- fread(file.path(fpath_codelists, "bowel_cancer_specific", "CRC_codelist_markers_clinical_formatted.csv"))
symptoms_list <- unique(codelist_markers %>% select(symptom))

codelist_medDIARR <-  fread(file.path(fpath_codelists, "medication/diarrhoea_final.csv"))
codelist_medHAEMR <-  fread(file.path(fpath_codelists, "medication/haemorrhoid_final.csv"))
codelist_medCONS <-  fread(file.path(fpath_codelists, "medication/laxatives_final.csv"))

###CLINICAL RECORDS
##read in clinical data
df_gp <- fread(fpath_gp, colClasses = "character")
df_gp <- df_gp %>% mutate(read_2 = na_if(read_2, ""), read_3 = na_if(read_3, ""))
df_gp$event_dt <- as.Date(df_gp$event_dt, format = "%d/%m/%Y")

##query clinical data (only do this once!)
df_clinical_events <- gp_clinical_query(df_gp, codelist_markers) %>% rename(event_type = symptom)
save(df_clinical_events, file = file.path(res_dir, "df_clinical_events.Rdata"))
base::load(file.path(res_dir, "df_clinical_events.Rdata"))

any_in_lb_events = c("ABDO_LUMP", "RECTAL_MASS", "CHANGE_BH", "RECTAL_BLEED") # alarm symptoms
freq_in_lb_events = c("ABDO_BLOAT", "ABDO_PAIN", "PELVIC_PAIN", "STOM_DIS") #non-alarm symptoms
new_onset_events = c("CONS", "DIARR", "WEIGHT_LOSS", "JAUNDICE", "FATIGUE", "DIV", "IBS", "HAEMR") #new onset sysmptoms or possible mis-diagnoses

###PRESCRIPTION RECORDS

##read in prescription data
df_scripts <- fread(fpath_gp_med, colClasses = "character")
df_scripts$issue_date <- as.Date(df_scripts$issue_date, format = "%d/%m/%Y")

##query prescription data
df_scripts_DIARR <- gp_script_query(df_scripts, codelist_medDIARR) %>% mutate(script_type = "DIARR")
df_scripts_HAEMR <- gp_script_query(df_scripts, codelist_medHAEMR) %>% mutate(script_type = "HAEMR")
df_scripts_CONS <- gp_script_query(df_scripts, codelist_medCONS) %>% mutate(script_type = "CONS") %>% select(-c(med_type))
df_script_events <- rbind(df_scripts_DIARR, df_scripts_HAEMR, df_scripts_CONS)
#save(df_script_events, file = file.path(res_dir, "df_script_events_markers.Rdata"))
#base::load(file.path(res_dir, "df_script_events_markers.Rdata"))
new_onset_scripts <- c("DIARR", "HAEMR", "CONS")

###GP REGISTRATION RECORDS
#get gp registration and date of birth info
base::load(file.path(res_dir, "df_dob.Rdata")) #loads pre-calculated dob
#gp_reg <- read_parquet(file.path(res_dir, "df_gp_reg.parquet"))
base::load(file.path(res_dir, "df_gp_reg.Rdata")) #loads dataframe gp_reg (hh only)
gp_reg$eid <- as.character(gp_reg$eid)
df_EHR <- right_join(df_dob, gp_reg, by = "eid")
df_EHR <- df_EHR %>% select(-any_of(c("data_provider")))  #remove data provider column (if present) to avoid nasty merge bug later
df_EHR <- df_EHR_savepoint 


get_df_gp_marker_events <- function(lma) {
  #look back window - between age of 18 (records of adults only) and lma date
  df_EHR$cutoff_min <- df_EHR$dob + years(18) 
  df_EHR$cutoff_max <- df_EHR$dob + years(lma) # lma date (cut-off)

  df_EHR$two_yr_lb <- df_EHR$dob + years(lma) - years(2)
  df_EHR$five_yr_lb <- df_EHR$dob + years(lma) - years(5)

  ##CLINICAL EVENTS
  # limit GP events to dates before lma age
  df_events_plus <- merge(df_clinical_events, df_EHR, all.x = TRUE)
  df_events_plus <- df_events_plus[(df_events_plus$event_dt >= df_events_plus$cutoff_min) & 
                                    (df_events_plus$event_dt <= df_events_plus$cutoff_max)]
  
  #group data by person
  df_events_by_eid <- df_events_plus %>%  group_by(eid) 
  
  #alarm symptoms (any in lb)
  for (i in any_in_lb_events) {
    df_events_by_eid <- event_in_lb(df_events_by_eid, i)
  }
  #non-alarm symtpoms (freq in lb)
  for (i in freq_in_lb_events) {
    df_events_by_eid <- event_freq_in_lb(df_events_by_eid, i)
  }
  #new onset symptoms (or misdiagnoses)
  for (i in new_onset_events) {
    df_events_by_eid <- event_new_onset(df_events_by_eid, i)
  }
  
  #get df with one row per eid
  df_events_by_eid_ungroup <- df_events_by_eid %>% 
                                filter(row_number()==1) %>% 
                                  select("eid",  contains(c("_ever", "_in_lb", "new_onset"))) %>% 
                                    ungroup    
  
  #merge in rest of primary care pop
  df_EHR <- df_EHR[(df_EHR$cutoff_max >= df_EHR$gp_start_date + months(6)) & (df_EHR$two_yr_lb <= df_EHR$gp_end_date - months(6))]
  df_gp_marker_events <- right_join(df_events_by_eid_ungroup, df_EHR, by = "eid")
  df_gp_marker_events[is.na(df_gp_marker_events)] <- 0
  df_gp_marker_events <- df_gp_marker_events %>% select(eid, assessment_date, cutoff_max, ends_with("in_lb"), ends_with("new_onset"))

  ##PRESCRIPTION EVENTS
  # limit prescription events to dates before lma age
  df_scripts_plus <- merge(df_script_events, df_EHR, all.x = TRUE)
  df_scripts_plus <- df_scripts_plus[(df_scripts_plus$event_dt >= df_scripts_plus$cutoff_min) & 
                                    (df_scripts_plus$event_dt <= df_scripts_plus$cutoff_max)]
  #group data by person
  df_scripts_by_eid <- df_scripts_plus %>%  group_by(eid) 

  #find new onset prescriptions
  for (i in new_onset_scripts) {
    df_scripts_by_eid <- script_new_onset(df_scripts_by_eid, i)
  }
  
  #get df by eid
  df_scripts_by_eid_ungroup <- df_scripts_by_eid %>% 
                                filter(row_number()==1) %>% 
                                  select("eid",  contains(c("new_onset"))) %>% 
                                    ungroup 

  #merge in with whole EHR pop for lma
  df_gp_marker_scripts <- right_join(df_scripts_by_eid_ungroup, df_EHR, by = "eid")
  df_gp_marker_scripts[is.na(df_gp_marker_scripts)] <- 0
  df_gp_marker_scripts <- df_gp_marker_scripts %>% select(eid,  ends_with("new_onset"))

  #COMBINE CLINICAL AND PRESCRIPTION INFO
  df_gp_marker_all <- full_join(df_gp_marker_events, df_gp_marker_scripts, by = "eid")
  df_gp_marker_all <- df_gp_marker_all %>% mutate(all_CONS_new_onset = case_when((CONS_new_onset == 1 | MED_CONS_new_onset) ~ 1, TRUE ~ 0), 
                                            all_DIARR_new_onset = case_when((DIARR_new_onset == 1 | MED_DIARR_new_onset) ~ 1, TRUE ~ 0), 
                                            all_HAEMR_new_onset = case_when((HAEMR_new_onset == 1 | MED_HAEMR_new_onset) ~ 1, TRUE ~ 0))
  df_gp_marker_all <- df_gp_marker_all %>% select(-c(CONS_new_onset, MED_CONS_new_onset, DIARR_new_onset, MED_DIARR_new_onset, HAEMR_new_onset, MED_HAEMR_new_onset))
    
  #save final dataframe of marker variables
  save(df_gp_marker_all, file = file.path(res_dir, paste0("marker_events/df_gp_marker_events_lma" , lma, ".Rdata"))) #test save
  #save(df_other_events, file.path(res_dir, paste0("other_events/df_gp_other_events_lma", lma, ".Rdata")))
}

#lapply(as.list(lmas), function(lma) {
 # get_df_gp_marker_events(lma)
 # cat("Done ", lma, "\n")
  #})

##debug
res_dir = "~/df"
get_df_gp_marker_events(50)
load(file = file.path(res_dir, paste0("marker_events/df_gp_marker_events_lma50.Rdata")))
df_temp_50 <- df_gp_marker_all
df_temp_50_eligible <- df_temp_50 %>% filter(assessment_date <= cutoff_max)



#get_df_gp_marker_events(60)
get_df_gp_marker_events(60)
load(file = file.path(res_dir, paste0("marker_events/df_gp_marker_events_lma60.Rdata")))
df_temp_60 <- df_gp_marker_all
df_temp_60_eligible <- df_temp_60 %>% filter(assessment_date <= cutoff_max)

get_df_gp_marker_events(70)
load(file = file.path(res_dir, paste0("marker_events/df_gp_marker_events_lma70.Rdata")))
df_temp_70 <- df_gp_marker_all
df_temp_70_eligible <- df_temp_70 %>% filter(assessment_date <= cutoff_max)


df_temp_50 %>% summarise(across(where(is.numeric), list(sum = ~sum(.x), num = ~n()))) %>%
                  pivot_longer(everything(), names_to = c("variable",".value"), names_pattern = "(.+)_(.+)") %>% 
                    mutate(percent = (sum/num)*100) %>% print()

df_temp_60 %>% summarise(across(where(is.numeric), list(sum = ~sum(.x), num = ~n()))) %>%
                  pivot_longer(everything(), names_to = c("variable",".value"), names_pattern = "(.+)_(.+)") %>% 
                    mutate(percent = (sum/num)*100) %>% print()

df_temp_70 %>% summarise(across(where(is.numeric), list(sum = ~sum(.x), num = ~n()))) %>%
            pivot_longer(everything(), names_to = c("variable",".value"), names_pattern = "(.+)_(.+)") %>% 
              mutate(percent = (sum/num)*100) %>% print()

df_temp_50_eligible %>% summarise(across(where(is.numeric), list(sum = ~sum(.x), num = ~n()))) %>%
          pivot_longer(everything(), names_to = c("variable",".value"), names_pattern = "(.+)_(.+)") %>% 
            mutate(percent = (sum/num)*100) %>% print()

df_temp_60_eligible %>% summarise(across(where(is.numeric), list(sum = ~sum(.x), num = ~n()))) %>%
          pivot_longer(everything(), names_to = c("variable",".value"), names_pattern = "(.+)_(.+)") %>% 
            mutate(percent = (sum/num)*100) %>% print()

df_temp_70_eligible %>% summarise(across(where(is.numeric), list(sum = ~sum(.x), num = ~n()))) %>%
          pivot_longer(everything(), names_to = c("variable",".value"), names_pattern = "(.+)_(.+)") %>% 
            mutate(percent = (sum/num)*100) %>% print()