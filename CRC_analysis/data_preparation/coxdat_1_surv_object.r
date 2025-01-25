# =============================================================================
# Script Title: Generate Survival Object [coxdat_ scripts 1 of 3]
# Author: Sam Ip
# - Combines baseline, GP registration, genetic, and DOB data into an outcome dataset.
# - Filters and prepares follow-up data based on specified follow-up windows.
# - Creates survival objects (tstart, tstop, status) for specified LMA (landmark age) values.
# - Outputs survival objects as parquet files for downstream analysis.
# =============================================================================

rm(list=ls())

# Define Directories ----
dir_ACED <- "~/ACED-CRC-modelling"
setwd(file.path(dir_ACED, "CRC_analysis"))
scripts_modelling_dir <- file.path(dir_ACED, "CRC_analysis", "modelling")

# Load Configurations ----
source(file.path(scripts_modelling_dir, "hpc_config.r"))

# Read Data ----
df_baseline <- read_parquet(file.path(res_dir, "df_baseline.parquet"), col_select = c(eid, death_date, first_cancer_date, first_type_crc, date_attendassessmentcentre))%>% mutate(across(eid, as.character))
df_gp_reg <- read_parquet(file.path(res_dir, "df_gp_reg.parquet"))%>% mutate(across(eid, as.character))
df_genetic <- read_parquet(file.path(res_dir, "df_genetic.parquet"), col_select = c(eid))%>% mutate(across(eid, as.character))
df_dob <- read_parquet(file.path(res_dir, "dob_from_baseline.parquet"), col_select = c(eid, dob)) %>% mutate(across(eid, as.character))

# Combine Data into df_outcome ----
df_outcome <- df_gp_reg %>% inner_join(df_baseline, by="eid") %>% inner_join(df_genetic, by="eid") %>% left_join(df_dob, by="eid")

# Filter out individuals with gp_start_date after imposed study_end_date (based on distribution of cases) ----
df_outcome <- df_outcome %>% filter(gp_start_date < study_end_date)

# Create overall follow-up start and stop times for each individual (not LMA-specific) ----
df_outcome <- df_outcome %>% mutate(
  tstart_baseline_gpstart = pmax(gp_start_date, date_attendassessmentcentre),
  tstop_death_firstcancer = pmin(death_date, first_cancer_date, study_end_date, na.rm=TRUE)
  )
  

# Create survival objects for a specific LMA
get_survobj <- function(lma){
  # Define LMA-specific follow-up window
  df_outcome$lma_followupwindow_start <- df_outcome$dob + years(lma) # Date of LMA
  df_outcome$lma_followupwindow_stop <- df_outcome$dob + years(lma) + years(num_yrs_fwd) # Date of follup-up end (for LMA)

  # Follow up (excl lookback) must start after UKB baseline
  df_outcome <- df_outcome %>% filter((lma_followupwindow_start > tstart_baseline_gpstart) & (lma_followupwindow_start < tstop_death_firstcancer))

  # Define time-to-event variables
  df_outcome <- df_outcome %>% dplyr::mutate(
    tstart = lma,
    tstop = gen_event_age(dob, pmin(lma_followupwindow_stop, tstop_death_firstcancer)),
    lookback_cutoff_min = lma_followupwindow_start- years(num_yrs_back), # Date of LMA-lookback time = earliest date lookback starts
    lookback_cutoff_max = lma_followupwindow_start # Date of LMA
    )

  # Filter records based on GP registration and follow-up criteria
  df_outcome <- df_outcome %>% filter(
    # Must have at least 6 months of lookback
    (lookback_cutoff_max >= gp_start_date + months(6)) &  
    (lookback_cutoff_min <= gp_end_date - months(6)) & 
    # Must have (post LMA) followup time
    (lookback_cutoff_max <= pmin(study_end_date - years(num_yrs_fwd), tstop_death_firstcancer))
  )

  # Event status ----
  df_outcome <- df_outcome %>% dplyr::mutate(first_cancer_age = gen_event_age(dob, first_cancer_date))
  df_outcome <- df_outcome %>% dplyr::mutate(status = ifelse((first_type_crc==1) & (first_cancer_age >= tstart) & (first_cancer_age <= tstop), 1, 0))

  # Save survival object as parquet ---- 
  # Select relevant columns for survival analysis
  df_outcome <- df_outcome %>% dplyr::select(eid, tstart, tstop, status)
  write_parquet(df_outcome, file.path(res_dir, paste0("survobj/survobj_lma", lma, "_nfwd", num_yrs_fwd, ".parquet")))
  
  return(df_outcome)
}

# MAIN: Apply Function for All LMAs ----
ls_surv <- lapply(as.list(lmas), function(lma) get_survobj(lma))
