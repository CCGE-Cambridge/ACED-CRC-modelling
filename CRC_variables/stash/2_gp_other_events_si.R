##File name: 2_gp_other_events_HH_NOV2022.R
##Authors: Sam Ip and Hannah Harrison
##Last Edit: 11/11/2022
##Description: identifies events from EHR (gp clinical)  within lma framework
##these events must have no time restriction other than happened before date of interest 
##and no additional formatting or conditions are applied (e.g. event happened ever)

# =================================  PARAMS  ==================================
# lma <- 40
# ===============================  read in ====================================
root_dir <- "~/ACED-RREDD-EHR/sam"
script_dir <- file.path(root_dir, "scripts")
source(file.path(script_dir, "0_config.r"))
library(parallel)

fpath_codelists <- "~/ACED-RREDD-EHR/codelists"
codelist_colonoscopy <- fread(file.path(fpath_codelists, "bowel_cancer_specific", "colonoscopy_codes_final.csv"))
HES_colonoscopy_codes <- fread(file.path(fpath_codelists, "bowel_cancer_specific", "colonoscopy_HES.csv"))
patt_colonoscopy_codes <- (HES_colonoscopy_codes %>% summarise(OPCS_Code = paste(OPCS_Code, collapse = "|")))$OPCS_Code[1]
patt_colonoscopy_codes <- gsub(x = patt_colonoscopy_codes, pattern = "\\.", replacement = "")

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

##get HES events info and dataframe
#baseline variable lists
fpath_baseline_vars_inc_HES <- "~/ACED-RREDD-EHR/codelists/baseline/baseline_vars_inc_HES_list.csv"
all_baseline_fields <- fread(fpath_baseline_vars_inc_HES)
HES_oper_fields <- all_baseline_fields %>% filter(type == "HES_oper") 
df_baseline_oper <- fread(fpath_baseline, select = c("eid", HES_oper_fields$field))
df_baseline_oper$eid <- as.character(df_baseline_oper$eid)

df_baseline_oper[sapply(df_baseline_oper, `%in%`, c("", "NA"))] <- NA
colnames(df_baseline_oper) <- gsub(x = names(df_baseline_oper), pattern = "-|\\.", replacement = "_") 
colnames(df_baseline_oper) <- gsub(x = names(df_baseline_oper), pattern = "^([0-9])", replacement = "v\\1") 

df_EHR <- df_EHR %>% full_join(df_baseline_oper, by="eid") #make sure extra fields are avaliable

patt_oper_var <- "v41272_0_"
patt_operdate_var <- "v41282_0_"

get_df_gp_other_events <- function(lma) {
  # lma<- 60
  #restrict to dates between age 18 and lma age. 
  df_EHR$cutoff_min <- df_EHR$dob + years(18) #ever occured as an adult 
  df_EHR$cutoff_max <- df_EHR$dob + years(lma) # lma date
  df_EHR$two_yr_lb <- df_EHR$dob + years(lma) - years(2) # lma date

  ##CLINICAL EVENTS
  # limit GP events to dates before lma age
  df_EHR_temp <- df_EHR%>%  select(!starts_with("v"))
  df_events_plus <- merge(df_clinical_events, df_EHR_temp, all.x = TRUE)
  df_events_plus <- df_events_plus[(df_events_plus$event_dt >= df_events_plus$cutoff_min) & 
                                    (df_events_plus$event_dt <= df_events_plus$cutoff_max)]
  
  #group data by person
  df_events_by_eid <- df_events_plus %>%  group_by(eid) 
  
  #any colonoscopy in last 10 years (gp records)
  for (i in any_in_10yr_lb_events) {
    df_events_by_eid <- event_in_var_lb(df_events_by_eid, i, 10) #use 10year lookback
  }
  
  #get dataframe with one row per eid
  df_events_by_eid_ungroup <- df_events_by_eid %>% 
                                filter(row_number()==1) %>% 
                                  select("eid",  contains(c("colonoscopy"))) %>% 
                                    ungroup    
  
  #get the HES colonoscopy records
  df_EHR <- UKB_wide_event_in_lb(df_EHR, patt_oper_var, patt_operdate_var, 123, "cutoff_max", 10, patt_colonoscopy_codes, "HES_colonoscopy")
  df_EHR <- df_EHR %>%  select(!starts_with("v"))

  #merge in rest of primary care pop
  df_EHR <- df_EHR[(df_EHR$cutoff_max >= df_EHR$gp_start_date + months(6)) & (df_EHR$two_yr_lb <= df_EHR$gp_end_date - months(6))]
  df_gp_other_events <- right_join(df_events_by_eid_ungroup, df_EHR, by = "eid")
  df_gp_other_events[is.na(df_gp_other_events)] <- 0

  #combine Primary care and HES records for colonsocopy
  df_gp_other_events <- df_gp_other_events %>% 
                mutate(colonoscopy_ALL_in_10_yr_lb = case_when(colonoscopy_in_10_yr_lb == 1 |  HES_colonoscopy_in_lb == 1 ~ 1, TRUE ~ 0))
  ##SCREENING ELIGIBILITY (was individual eligible 2 years before lma date)
  #criteria: aged between 60-74 years 2 years before lma and it's after 2008 (mid-point of roll-out)
  ###do we have data provider at this point?
  df_gp_other_events <- df_gp_other_events %>% 
            mutate(bcscreen_eligible = case_when((lma >= 62 & lma <= 76 & two_yr_lb >= as.Date("2008-01-01")) ~ 1, TRUE ~0))  

  #df_gp_other_events <- df_gp_other_events %>% select(eid, bcscreen_eligible, colonoscopy_in_10_yr_lb)
  df_gp_other_events <- df_gp_other_events %>% select(eid, bcscreen_eligible, colonoscopy_ALL_in_10_yr_lb)

  # save ----
  save(df_gp_other_events, file = file.path(res_dir, paste0("other_events/df_gp_other_events_lma" , lma, ".Rdata"))) #test save

  #save(df_other_events, file.path(res_dir, paste0("other_events/df_gp_other_events_lma", lma, ".Rdata")))
}

# get_df_gp_other_events(60) #test
# get_df_gp_other_events(70) #test
# get_df_gp_other_events(80) #test
# base::load(file.path(res_dir, "colonoscopy", paste0("df_gp_other_events_lma" , 60, ".Rdata")))
parallel::mclapply(as.list(lmas), mc.cores=16, function(lma) {
  get_df_gp_other_events(lma)
  cat("Done ", lma, "\n")
  })
