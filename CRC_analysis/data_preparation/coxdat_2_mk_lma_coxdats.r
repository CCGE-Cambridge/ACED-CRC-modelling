# ==============================================================================
# Script Title: Create Landmarking Cox Regression Dataset [coxdat_ scripts 2 of 3]
# Author: Sam Ip
# - Preceded by coxdat_1_surv_object.r
# Processes multiple data sources (baseline, biomarkers, genetic, etc.) 
# for a specified landmark age (LMA) and combines them into a single dataset 
# (`coxdat`). The dataset includes:
# - Relevant predictors (baseline, genetic, biomarker, multimorbidity, etc.)
# - Survival information (e.g., time to event, event status)
# - Derived variables for survival modeling
# ==============================================================================

library(parallel)

# Load Baseline Data ----
df_baseline <- read_parquet(
  file.path(res_dir, "df_baseline.parquet"),
  col_select=c("eid", "alcohol_units_daily", "birth_year", 
  "bmi", "edu_quals", "ethn_selfrep", "famhist_bowel_cancer", 
  "famhist_breast_cancer", "famhist_lung_cancer", "partial_fibre_score", "sex_genetic",
  "smoking_status", "tot_METs_wk", "tot_process_meat_wk", "tot_red_meat_wk", "townsend")
)

# Prefix variable names with `v0_` to indicate baseline values ----
names(df_baseline)[names(df_baseline) != "eid"] <- paste0("v0_", names(df_baseline)[names(df_baseline) != "eid"])


## Process Data for Each Landmark Age (LMA) ----
# For each specified landmark age, process the datasets and create a combined Cox dataset.
mclapply(as.list(lmas), mc.cores=36, function(lma) {
# lapply(as.list(lmas), function(lma) {  
  cat("Processing data for landmark age (LMA):", lma, "\n")
  
  # Load datasets specific to the landmark age
  df_biom <- base::get(base::load(file.path(res_dir, paste0("biomarkers/df_gp_biom_lma_", lma, "_nback", ".Rdata"))))
  names(df_biom)[names(df_biom) != "eid"] <- paste0("biom_", names(df_biom)[names(df_biom) != "eid"])
  ### "biom_iron_def", "biom_iron_def_measured", "biom_inflam_all", "biom_inflam_measured"

  df_multimorb <- base::get(base::load(file.path(res_dir, paste0("multimorb/df_gp_multimorb_lma" , lma, ".Rdata"))))
  names(df_multimorb)[names(df_multimorb) != "eid"] <- paste0("multimorb_", names(df_multimorb)[names(df_multimorb) != "eid"])
  ### multimorb_mm_score_PI -- after regression on age, sex -- multimorb_mm_score_res_fit

  df_other_events <- base::get(base::load(file.path(res_dir, paste0("other_events/df_gp_other_events_lma" , lma, ".Rdata"))))
  names(df_other_events)[names(df_other_events) != "eid"] <- paste0("other_", names(df_other_events)[names(df_other_events) != "eid"])
  ### other_bcscreen_eligible, other_colonoscopy_in_10_yr_lb

  df_modifier <- base::get(base::load(file.path(res_dir, paste0("modifier_events/df_gp_modifier_events_lma" , lma, ".Rdata"))))
  df_modifier  <- df_modifier %>% dplyr::select(eid, IBD_ever, DIB_T1D_ever, DIB_T2D_ever, GALL_ever, MED_NSAIDs_reg_ever, MED_ASPIRIN_reg_ever)
  names(df_modifier)[names(df_modifier) != "eid"] <- paste0("mod_", names(df_modifier)[names(df_modifier) != "eid"])
  ### mod_IBD_ever, mod_DIB_T1D_ever, mod_DIB_T2D_ever, mod_GALL_ever, mod_MED_NSAIDs_reg_ever, mod_MED_ASPIRIN_reg_ever

  df_marker <- base::get(base::load(file.path(res_dir, paste0("marker_events/df_gp_marker_events_lma" , lma, ".Rdata"))))
  names(df_marker)[names(df_marker) != "eid"] <- paste0("mrk_", names(df_marker)[names(df_marker) != "eid"])
  ### mrk_ABDO_LUMP_in_lb, mrk_RECTAL_MASS_in_lb, mrk_CHANGE_BH_in_lb, mrk_RECTAL_BLEED_in_lb, 
  ### mrk_ABDO_BLOAT_freq_in_lb, mrk_ABDO_PAIN_freq_in_lb, mrk_PELVIC_PAIN_freq_in_lb, 
  ### mrk_STOM_DIS_freq_in_lb, mrk_WEIGHT_LOSS_new_onset, mrk_JAUNDICE_new_onset, mrk_FATIGUE_new_onset, 
  ### mrk_DIV_new_onset, mrk_IBS_new_onset, mrk_all_CONS_new_onset, mrk_all_DIARR_new_onset, mrk_all_HAEMR_new_onset
  
  df_genetic <- read_parquet(file.path(res_dir, "df_genetic.parquet"), col_select=c("eid", "prs", paste0("pc", 1:10)))
  names(df_genetic)[names(df_genetic) != "eid"] <- paste0("gen_", names(df_genetic)[names(df_genetic) != "eid"])
  ### gen_prs, gen_pc1-10

  survobj <- read_parquet(file.path(res_dir, paste0("survobj/survobj_lma", lma, "_nfwd", num_yrs_fwd, ".parquet")))
  ### tstart, tstop, status, lookback


  # eid as.character
  ls_df_names <- list(survobj, df_baseline, df_genetic, df_biom, df_multimorb, df_other_events, df_modifier, df_marker)
  ls_dfs <- lapply(ls_df_names, function(df){dplyr::mutate_at(df, vars(matches("^eid")), as.character)})

  # Combine datasets
  coxdat <- ls_dfs[1:3] %>% reduce(inner_join, by = "eid") # Inner join core datasets (cohort spine) -- only include people with entries in all of these dfs
  coxdat <- append(list(coxdat),  ls_dfs[4:length(ls_df_names)]) %>% reduce(left_join, by = "eid") # Left join additional data
  sort(names(coxdat))

  # Compute derived variables for Cox regression (tstart2, tstop != tstart) ----
  coxdat <- coxdat %>%
    dplyr::mutate(
      age = tstart,
      age2 = age^2,
      futime = ifelse(tstop == tstart, 0.001, tstop - tstart)  # Ensure non-zero follow-up time
    )

  # Select survival info & covariates
  covars <- names(coxdat)[grepl("^v0_|^gen_|^biom_|^mod_|^mrk_|^multimorb_|^other_", names(coxdat))]
  coxdat <- coxdat %>% dplyr::select(eid, futime, status, age, age2, all_of(covars))

  # Save LMA dataset as parquet ----
  write_parquet(coxdat, file.path(res_dir, paste0("coxdat/coxdat_lma", lma, "_back", num_yrs_back, "_fwd", num_yrs_fwd, ".parquet")))
})
