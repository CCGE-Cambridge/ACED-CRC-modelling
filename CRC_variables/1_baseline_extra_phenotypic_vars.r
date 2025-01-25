##File name: extra_phenotypic_vars.r
##Author: Hannah Harrison
##Last Edit: 04/03/2022
##Description: creates composite variables for colorectal cancer using phenotypic UKB dataset

#libraries if needed here (most now added in config)
#library(tidyverse) 
#library(readxl) 
#library(dplyr)
#library(data.table)
#library(lubridate)
library(remotes)
# fpath_baseline <- "/store/biobank/health_records/ukb44651.csv" #load from config -- updated data!
# path_to_local_data <- "~/df"
source("~/ACED-RREDD-EHR/sam/scripts/0_config.r")

path_to_local_data <- res_dir
##functions
###adaptations to rowsums, handles weighting and manages missing values (if all(is.na(.)) ~ NA, otherwise treat NA = 0)
#source("~/ACED-RREDD-EHR/CRC_variables/dependencies/gen_functions.r") 
#source("~/ACED-RREDD-EHR/CRC_variables/dependencies/base_vars_functions.r")

#baseline variable lists
fpath_baseline_vars <- "~/ACED-RREDD-EHR/codelists/baseline/baseline_vars_hannah.csv"
all_baseline_fields <- fread(fpath_baseline_vars)

##Education
edu_fields  <- all_baseline_fields %>% filter(type == "education")
df_baseline_edu <- fread(fpath_baseline, select = c("eid", edu_fields$field))
names(df_baseline_edu) <- c("eid", edu_fields$name)
df_baseline_edu <- gen_edu_length(df_baseline_edu)
df_baseline_edu <- df_baseline_edu %>% select(eid, edu_quals)

##Physical Activity
####generates mets_per_day (continuous), IPAQ_cat (ordered categorical) and recommonded_pa (binary)
pa_fields <- all_baseline_fields %>% filter(type == "pa" & name != "")
df_baseline_pa <- fread(fpath_baseline, select = c("eid", pa_fields$field))
names(df_baseline_pa) <- c("eid", pa_fields$name)
df_baseline_pa <- pa_do_all(df_baseline_pa)
df_baseline_pa <- df_baseline_pa %>% select(eid, tot_METs_wk)

##Alcohol consumption)
alcohol_fields <- all_baseline_fields %>% filter(type == "alcohol")
df_baseline_alcohol <- fread(fpath_baseline, select = c("eid", alcohol_fields$field))
names(df_baseline_alcohol) <- c("eid", alcohol_fields$name)
df_baseline_alcohol <- alc_daily_units(df_baseline_alcohol)
#df_baseline_alcohol <- alc_exceed(df_baseline_alcohol) not using
df_baseline_alcohol <- df_baseline_alcohol %>% select(eid, alcohol_units_daily)

##Dietary Variables
dietary_fields <- all_baseline_fields %>% filter(type == "dietary")
df_baseline_dietary <- fread(fpath_baseline, select = c("eid", dietary_fields$field))
names(df_baseline_dietary) <- c("eid", dietary_fields$name)

df_baseline_dietary <- gen_tot_red_meat(df_baseline_dietary)
df_baseline_dietary <- gen_tot_proc_meat(df_baseline_dietary)
#df_baseline_dietary <- gen_tot_all_meat(df_baseline_dietary) # don't need for colorectal cancer
df_baseline_dietary <- gen_partial_fibre_score(df_baseline_dietary)
#df_baseline_dietary <- gen_ca_score(df_baseline_dietary) #not using
df_baseline_dietary <- df_baseline_dietary %>% 
  select(eid, tot_process_meat_wk, tot_red_meat_wk, partial_fibre_score)

###sam's family history code

###sam's ethnicity code

##sam's code for BMI (link from other dataset)
                    

##merge into "complete dataset"
#put all data frames into list
df_list <- list(df_baseline_edu, df_baseline_pa, df_baseline_alcohol, df_baseline_dietary)
#merge all data frames in list
df_all <- df_list %>% reduce(full_join, by='eid')

##save database here 
save(df_all, file = file.path(path_to_local_data, "df_baseline_extra_phenotypic_vars.Rdata"))


# load(file = "~/df/df_baseline_colorectal.Rdata")




