# ==============================================================================
# Script Title: Creates and cleans superlandmark dataset (stacked LMA Cox datasets) [coxdat_ scripts 3 of 3]
# Author: Sam Ip
# - Preceded by coxdat_2_mk_lma_coxdats.r
# - Loading age-specific datasets and combining them into a single stacked dataset
# - Regressing out age/sex from the multimorbidity score
# - Identifying and handling missing data
# - Removing withdrawn participants and ensuring complete cases
# ==============================================================================

# Read and combine Cox datasets for all landmark ages ----
coxdat_stacked_allages_num <- lapply(as.list(lmas), function(lma) {read_parquet(file.path(res_dir, paste0("coxdat/coxdat_lma", lma, "_back", num_yrs_back, "_fwd", num_yrs_fwd, ".parquet")))}) %>% rbindlist()

# Adjust multimorbidity score by regressing out the effects of age and sex ----
coxdat_stacked_allages_num <-  mm_score_take_res(coxdat_stacked_allages_num, "multimorb_mm_score_PI") %>%
    dplyr::rename(multimorb_mm_score_res_fit=residual_fit) %>% dplyr::select(-one_of(c("multimorb_mm_score_PI", "age_sex_fit")))

# (Added) Remove physical activity column (v0_tot_METs_wk) and apply complete case filtering ----
coxdat_stacked_allages_num <- coxdat_stacked_allages_num %>% 
    dplyr::select(!v0_tot_METs_wk) %>%
    dplyr::filter(complete.cases(.))

# Check and remove participants who have subsequently withdrawn from the UK Biobank ----
ids_withdrawn <- fread("/store/biobank/phenotypes/withdraw_id/w28126_20220222.csv")
if (any(coxdat_stacked_allages_num$eid %in% ids_withdrawn$V1)) stop("Contains withdrawn participants\n")

# Save the superlandmark dataset as RData ----
save(coxdat_stacked_allages_num, file=file.path(shared_dir, paste0("coxdat_stacked_allages_num_", date_str, ".Rdata")))
