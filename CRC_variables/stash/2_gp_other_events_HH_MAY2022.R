##File name: 2_gp_other_events.r
##Authors: Sam Ip and Hannah Harrison
##Last Edit: 27/05/2021
##Description: identifies "other" events from EHR (gp clinical)  within lma framework
##colonoscopy in records and eligibility for bc screening
source("~/ACED-RREDD-EHR/sam/scripts/0_config.r")

##get codelists
#source("~/ACED-RREDD-EHR/CRC_variables/dependencies/gp_data_funcs.r")
fpath_codelists <- "~/ACED-RREDD-EHR/codelists"
codelist_colonoscopy <- fread(file.path(fpath_codelists, "bowel_cancer_specific", "colonoscopy_codes_final.csv"))

##CLINICAL RECORDS 
#read in clinical data 
df_gp <- fread(fpath_gp, colClasses = "character")
df_gp <- df_gp %>% mutate(read_2 = na_if(read_2, ""), read_3 = na_if(read_3, ""))
df_gp$event_dt <- as.Date(df_gp$event_dt, format = "%d/%m/%Y")
##query clinical data (only do this once!)
df_clinical_events <- gp_clinical_query(df_gp, codelist_colonoscopy) %>% mutate(event_type = "colonoscopy")
save(df_clinical_events, file = file.path(res_dir, "df_clinical_events.Rdata"))
base::load(file.path(res_dir, "df_clinical_events.Rdata"))
any_in_10yr_lb_events = c("colonoscopy") 


##get gp registration and date of birth info
#df_dob <- read_parquet(fpath_ukb_dob)
base::load(file.path(res_dir, "df_dob.Rdata")) #loads pre-calculated dob (hh only)
gp_reg <- read_parquet(file.path(res_dir, "df_gp_reg.parquet"))
# base::load(file.path(res_dir, "df_gp_reg.Rdata")) #loads dataframe gp_reg (hh only)
gp_reg$eid <- as.character(gp_reg$eid)
df_EHR <- right_join(df_dob, gp_reg, by = "eid")
df_EHR <- df_EHR %>% select(-any_of(c("data_provider")))  #remove data provider column (if present) to avoid nasty merge bug later


get_df_gp_other_events <- function(lma) {
    # lma<- 50
  #restrict to dates between age 18 and lma age. 
  df_EHR$cutoff_min <- df_EHR$dob + years(18) #ever occured as an adult 
  df_EHR$cutoff_max <- df_EHR$dob + years(lma) # lma date
  df_EHR$two_yr_lb <- df_EHR$dob + years(lma) - years(2) # lma date

  ##CLINICAL EVENTS
  # limit GP events to dates before lma age
  df_events_plus <- merge(df_clinical_events, df_EHR, all.x = TRUE)
  df_events_plus <- df_events_plus[(df_events_plus$event_dt >= df_events_plus$cutoff_min) & 
                                    (df_events_plus$event_dt <= df_events_plus$cutoff_max)]
  
  #group data by person
  df_events_by_eid <- df_events_plus %>%  group_by(eid) 
  
  #any colonoscopy in last 10 years
  for (i in any_in_10yr_lb_events) {
    df_events_by_eid <- event_in_var_lb(df_events_by_eid, i, 10) #use 10year lookback
  }
  
  #get df with one row per eid
  df_events_by_eid_ungroup <- df_events_by_eid %>% 
                                filter(row_number()==1) %>% 
                                  select("eid",  contains(c("colonoscopy"))) %>% 
                                    ungroup    
  
  #merge in rest of primary care pop
  df_EHR <- df_EHR[(df_EHR$cutoff_max >= df_EHR$gp_start_date + months(6)) & (df_EHR$two_yr_lb <= df_EHR$gp_end_date - months(6))]
  df_gp_other_events <- right_join(df_events_by_eid_ungroup, df_EHR, by = "eid")
  df_gp_other_events[is.na(df_gp_other_events)] <- 0

  ##SCREENING ELIGIBILITY (was individual eligible 2 years before lma date)
  #criteria: aged between 60-74 years 2 years before lma and it's after 2008 (mid-point of roll-out)
  ###do we have data provider at this point?
  df_gp_other_events <- df_gp_other_events %>% 
            mutate(bcscreen_eligible = case_when((lma >= 62 & lma <= 76 & two_yr_lb >= as.Date("2008-01-01")) ~ 1, TRUE ~0))  

  df_gp_other_events <- df_gp_other_events %>% select(eid, bcscreen_eligible, colonoscopy_in_10_yr_lb)
  

  # save ----
  save(df_gp_other_events, file = file.path(res_dir, paste0("other_events/df_gp_other_events_lma" , lma, ".Rdata"))) #test save

  #save(df_other_events, file.path(res_dir, paste0("other_events/df_gp_other_events_lma", lma, ".Rdata")))
}

lapply(as.list(lmas), function(lma) {
  get_df_gp_other_events(lma)
  cat("Done ", lma, "\n")
  })

# get_df_gp_other_events(60)
#load(file = file.path(res_dir, paste0("other_events/df_gp_other_events_lma60.Rdata")))
# get_df_gp_other_events(80)
#load(file = file.path(res_dir, paste0("other_events/df_gp_other_events_lma80.Rdata")))