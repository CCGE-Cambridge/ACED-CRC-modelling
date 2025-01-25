##File name: 2_modifiers_of_risk.r
##Authors: Hannah Harrison
##Last Edit: 26/05/2022
##Description: identifies modifiers of risk from EHR (gp_clincial)  within lma framework
##these events must have no time restriction other than happened before date of interest 
##NB: additional filtering applied to diabetes codes - to mange the separation into T1D and T2D, while still using the information from unspecified diabetes type codes
##NB: in terms of analysis should be grouped with variables from 2_gp_meds.r (NSAIDs and aspirin) - could add here for simplicity?
source("~/ACED-RREDD-EHR/sam/scripts/0_config.r")

##get codelists
#source("~/ACED-RREDD-EHR/CRC_variables/dependencies/gp_data_funcs.r")
fpath_codelists <- "~/ACED-RREDD-EHR/codelists"
codelist_modifiers <- fread(file.path(fpath_codelists, "bowel_cancer_specific", "CRC_codelist_modifiers_clinical_formatted.csv"))
codelist_NSAIDs <- fread(file.path(fpath_codelists, "medication/NSAIDs_final.csv"))
codelist_aspirin <- fread(file.path(fpath_codelists, "medication/aspirin_final.csv"))

###CLINICAL RECORDS
##read in clinical data
df_gp <- fread(fpath_gp, colClasses = "character")
df_gp <- df_gp %>% mutate(read_2 = na_if(read_2, ""), read_3 = na_if(read_3, ""))
df_gp$event_dt <- as.Date(df_gp$event_dt, format = "%d/%m/%Y")

##query clinical data (only do this once!)
df_clinical_events <- gp_clinical_query(df_gp, codelist_modifiers) %>% rename(event_type = condition)
save(df_clinical_events, file = file.path(res_dir, "df_clinical_events_modifiers.Rdata"))
base::load(file.path(res_dir, "df_clinical_events_modifiers.Rdata"))
ever_events = c("IBD", "DIB_T1D", "DIB_T2D", "GALL") # long term conditions

###PRESCRIPTION RECORDS
##read in prescription data
df_scripts <- fread(fpath_gp_med, colClasses = "character")
df_scripts$issue_date <- as.Date(df_scripts$issue_date, format = "%d/%m/%Y")

##query prescription data
df_scripts_NSAIDs <- gp_script_query(df_scripts, codelist_NSAIDs) %>% mutate(script_type = "NSAIDs")
df_scripts_aspirin <- gp_script_query(df_scripts, codelist_aspirin) %>% mutate(script_type = "ASPIRIN")
df_script_events <- rbind(df_scripts_NSAIDs, df_scripts_aspirin)
save(df_script_events, file = file.path(res_dir, "df_script_events.Rdata"))
base::load(file.path(res_dir, "df_script_events.Rdata"))
regular_scripts <- c("NSAIDs", "ASPIRIN")

#get gp registration and date of birth info
base::load(file.path(res_dir, "df_dob.Rdata")) #loads pre-calculated dob 
gp_reg <- read_parquet(file.path(res_dir, "df_gp_reg.parquet"))
# base::load(file.path(res_dir, "df_gp_reg.Rdata")) #loads dataframe gp_reg (hh only)
gp_reg$eid <- as.character(gp_reg$eid)
df_EHR <- right_join(df_dob, gp_reg, by = "eid")
df_EHR <- df_EHR %>% select(-any_of(c("data_provider")))  #remove data provider column (if present) to avoid nasty merge bug later

# for test run only

get_df_gp_modifier_events <- function(lma) {
  #look back window - between age of 18 (records of adults only) and lma date
  df_EHR$cutoff_min <- df_EHR$dob + years(18) #ever occured as an adult 
  df_EHR$cutoff_max <- df_EHR$dob + years(lma) # lma date
  df_EHR$two_yr_lb <- df_EHR$dob + years(lma) - years(2) # lma date
  ##CLINICAL EVENTS
  # limit GP events to: occured when adult and before lma age
  df_events_plus <- merge(df_clinical_events, df_EHR, all.x = TRUE)
  df_events_plus <- df_events_plus[(df_events_plus$event_dt >= df_events_plus$cutoff_min) & 
                                    (df_events_plus$event_dt <= df_events_plus$cutoff_max)]

  #group data by eid
  df_events_by_eid <- df_events_plus %>%  group_by(eid) 

  #long term conditions (event ever)
  for (i in ever_events) {
    df_events_by_eid <- event_ever(df_events_by_eid, i)
  }

  #get df with one row per eid
  df_events_by_eid_ungroup <- df_events_by_eid %>% 
                                filter(row_number()==1) %>% 
                                  select("eid",  contains(c("_ever"))) %>% 
                                    ungroup    
  
  #merge in rest of primary care pop
  df_EHR <- df_EHR[(df_EHR$cutoff_max >= df_EHR$gp_start_date + months(6)) & (df_EHR$two_yr_lb <= df_EHR$gp_end_date - months(6))]

  df_gp_modifier_events <- right_join(df_events_by_eid_ungroup, df_EHR, by = "eid")
  df_gp_modifier_events[is.na(df_gp_modifier_events)] <- 0
  df_gp_modifier_events <- df_gp_modifier_events %>% select(eid, ends_with("ever"))

  ##PRESCRIPTION EVENTS
  # limit prescription events to dates before lma age
  df_scripts_plus <- merge(df_script_events, df_EHR, all.x = TRUE)
  df_scripts_plus <- df_scripts_plus[(df_scripts_plus$event_dt >= df_scripts_plus$cutoff_min) & 
                                    (df_scripts_plus$event_dt <= df_scripts_plus$cutoff_max)]
  #group data by person
  df_scripts_by_eid <- df_scripts_plus %>%  group_by(eid) 

  #find regular prescriptions
  for (i in regular_scripts) {
    df_scripts_by_eid <- script_reg_ever(df_scripts_by_eid, i, 4, 1)
  }
  
  #get df by eid
  df_scripts_by_eid_ungroup <- df_scripts_by_eid %>% 
                                filter(row_number()==1) %>% 
                                  select("eid",  contains(c("reg_ever"))) %>% 
                                    ungroup 

  #merge in with whole EHR pop for lma
  df_gp_modifier_scripts <- right_join(df_scripts_by_eid_ungroup, df_EHR, by = "eid")
  df_gp_modifier_scripts[is.na(df_gp_modifier_scripts)] <- 0
  df_gp_modifier_scripts <- df_gp_modifier_scripts %>% select(eid,  ends_with("reg_ever"))

      
  #COMBINE CLINICAL AND PRESCRIPTION INFO
  df_gp_modifier_all <- full_join(df_gp_modifier_events, df_gp_modifier_scripts, by = "eid")

  # save ----
  save(df_gp_modifier_all, file = file.path(res_dir, paste0("modifier_events/df_gp_modifier_events_lma" , lma, ".Rdata"))) #test save
  #save(df_other_events, file.path(res_dir, paste0("other_events/df_gp_other_events_lma", lma, ".Rdata")))
}


# lmas <- seq(66,75, 2)

parallel::mclapply(as.list(lmas), mc.cores = 16, function(lma) {
  get_df_gp_modifier_events(lma)
  cat("Done ", lma, "\n")
  })

# get_df_gp_modifier_events(50) #test
# get_df_gp_modifier_events(60) #test
#load(file = file.path(res_dir, paste0("modifier_events/df_gp_modifier_events_lma60.Rdata")))
# get_df_gp_modifier_events(70) #test
# get_df_gp_modifier_events(80) #test

