# ==============================================================================
# Script Title: Process GP Registrations Data and Create Continuous GP Record Periods [gp_ scripts 1 of 3]
# Author: Sam Ip
# - Imposes provider-specific censoring dates for `deduct_date`
# - Removes invalid registration periods (e.g., `reg_date` > `deduct_date`)
# - Removes periods that end before the earliest relevant date
# - Filters out discontinuous registration periods (>90-day gaps)
# - Keeps only continuous GP records that overlap with UK Biobank assessment
# - Applys further eligibility filters (e.g., removing prevalent cancer cases)
# ==============================================================================
rm(list=ls())

# Define Directories ----
dir_ACED <- "~/ACED-CRC-modelling"
setwd(file.path(dir_ACED, "CRC_analysis"))
scripts_modelling_dir <- file.path(dir_ACED, "CRC_analysis", "modelling")

# Load Configurations ----
source(file.path(scripts_modelling_dir, "hpc_config.r"))

# Load Data ----
gp_reg <- fread(fpath_gp_reg) # GP registration data
df_baseline <- read_parquet( # UKB baseline data
    file.path(res_dir, "df_baseline.parquet"), 
    col_select = c("eid", "date_attendassessmentcentre", "first_cancer_date", "first_cancer_age", "first_type_crc", "date_firstcancer_or_death")) 

# Merge GP registration data with baseline data ----
gp_reg <-  merge(x = gp_reg, y = df_baseline, by = "eid") # inner join 

## Apply Provider-Specific `deduct_date` Censoring ----
gp_reg <- gp_reg %>% mutate_at(vars(grep("date", names(gp_reg), value = TRUE)), funs(as.Date(.,"%d/%m/%Y"))) # Convert date columns to Date format
gp_reg$deduct_date <- fifelse((is.na(gp_reg$deduct_date) | (gp_reg$deduct_date > as.Date("2017-05-01"))) & gp_reg$data_provider==1, 
    as.Date("2017-05-01", format="%Y-%m-%d"),
    gp_reg$deduct_date
    )
gp_reg$deduct_date <- fifelse((is.na(gp_reg$deduct_date) | (gp_reg$deduct_date > as.Date("2017-04-01"))) & gp_reg$data_provider==2, 
    as.Date("2017-04-01", format="%Y-%m-%d"),
    gp_reg$deduct_date
    )

gp_reg$deduct_date <- fifelse((is.na(gp_reg$deduct_date) | (gp_reg$deduct_date > as.Date("2016-07-01"))) & gp_reg$data_provider==3, 
    as.Date("2016-07-01", format="%Y-%m-%d"),
    gp_reg$deduct_date
    )
gp_reg$deduct_date <- fifelse((is.na(gp_reg$deduct_date) | (gp_reg$deduct_date > as.Date("2017-08-01"))) & gp_reg$data_provider==4, 
    as.Date("2017-08-01", format="%Y-%m-%d"),
    gp_reg$deduct_date
    )

# Filter Invalid Registration Periods ----
# Remove entries with `reg_date` after `deduct_date`
gp_reg <- gp_reg %>% filter(reg_date<deduct_date)

# Remove nested registration periods
gp_reg <- gp_reg %>% group_by(eid) %>% 
    arrange(reg_date, .by_group = TRUE) %>% 
    mutate(next_reg_date = lead(reg_date), next_deduct_date = lead(deduct_date)) %>% ungroup() 
nested <- gp_reg %>% filter((next_reg_date >= reg_date) & (next_deduct_date <= deduct_date))
gp_reg <-  gp_reg %>% dplyr::anti_join(nested) 

# Flag discontinuous registration periods (>90-day gaps)
gp_reg <- gp_reg %>% group_by(eid) %>% 
    arrange(reg_date, .by_group = TRUE) %>% 
    mutate(discont = ifelse(!is.na(next_reg_date), 1*((next_reg_date-deduct_date)>90), 0)) %>% ungroup()

# Keep the last continuous period for each individual
deduct_date_lastdiscont <- gp_reg %>% filter(discont == 1)%>% 
    group_by(eid) %>% summarise(deduct_date_lastdiscont = max(deduct_date[discont==1])) 
gp_reg <- merge(x = gp_reg, y = deduct_date_lastdiscont, by = "eid", all.x = TRUE) %>% 
    dplyr::filter((reg_date > deduct_date_lastdiscont) | is.na(deduct_date_lastdiscont))


# Summarise Continuous GP Periods ----
gp_reg <- gp_reg %>% group_by(eid) %>% summarise(gp_start_date = min(reg_date), gp_end_date = max(deduct_date)) 

# Apply Eligibility Filters ----
# Remove entries where GP records end before UKB assessment date
gp_reg <- gp_reg %>% left_join(df_baseline, by="eid") %>% filter(gp_end_date > date_attendassessmentcentre)

# Remove individuals with prevalent cancers or deaths
gp_reg <- gp_reg %>% filter(is.na(date_firstcancer_or_death) | (date_firstcancer_or_death > gp_start_date))

# Keep GP records of at least 6 months
gp_reg <- as.data.frame(gp_reg) %>% filter(as.numeric(gp_end_date-gp_start_date) >= 365.25*0.5) 

# Save ----
gp_reg <- gp_reg %>% dplyr::select(eid, gp_start_date, gp_end_date) %>% distinct()
write_parquet(gp_reg, file.path(res_dir, "df_gp_reg.parquet"))
base::save(gp_reg, file = file.path(shared_dir, "df_gp_reg.Rdata"))
