# =====================================================================================
# Script Title: Shapley Value Analysis for Data Sources with Bootstrap Sampling
# Author: Sam Ip
# Description: 
# - Calculates (classic) Shapley values for each predictor group of each bootstrap sample
# - Calculates mean Shapley values and confidence intervals
# =====================================================================================

rm(list=ls())

# Load external scripts ----
dir_ACED <- "~/ACED-CRC-modelling"
scripts_modelling_dir <- file.path(dir_ACED, "CRC_analysis", "modelling")
source(file.path(scripts_modelling_dir, "hpc_init_getcoxdat.r"))

# Define directories ----
res_plots_dir <- "~/rds/hpc-work/BMJOnc_reviewerresp/tmp_boot_res" 
dpath_tmp_boot_res_tag <- file.path(res_plots_dir, glue::glue("bootsmps_combos{tag}"))
if (exists("tag_LDpred", envir = .GlobalEnv)) {dpath_tmp_boot_res <- glue::glue("{dpath_tmp_boot_res}{tag_LDpred}")} 
if (exists("tag_Vision", envir = .GlobalEnv)) {dpath_tmp_boot_res <- glue::glue("{dpath_tmp_boot_res}{tag_Vision}")} 

# Load bootstrap results and format ----
# Bootstrap results generated by hpc_mlr3_contributions_shapley_spreadovernodes.r
fpaths_cindices_all_combos_200btsmp <- list.files(dpath_tmp_boot_res_tag, pattern = "^cindex_", full.names = TRUE) 
cindices_all_combos_200btsmp <- parallel::mclapply(fpaths_cindices_all_combos_200btsmp, fread, mc.cores = 20) %>%
  bind_rows()
cindices_all_combos_200btsmp <- cindices_all_combos_200btsmp %>%
  dplyr::mutate(
    bootsmp_num = sub("_.*", "", bootsmpnum_coalition), # Extract bootstrap sample number
    coalition = sub("^[^_]*_", "", bootsmpnum_coalition) # Extract coalition
  )

# Functions ----

# Derive the base coalition name (e.g., "core") from the input coalition (e.g., "biom_core")
# by removing the predictor group of interest (data_source, e.g., "biom").
extract_base_coalition <- function(coalition, data_source) {
  coalition <- str_replace_all(coalition, paste0("(^|_)(", data_source, ")(_)?"), "\\1")
  coalition <- str_replace_all(coalition, "__", "_")  # Replace double underscores with a single underscore
  str_replace_all(coalition, "^_|_$", "")  # Remove leading or trailing underscores
}

# Calculate Shapley value for a predictor group (one bootstrap sample)
calc_shap_of_datasource <- function(data_source, cindices_all_combos_200btsmp) {
  # Create a new column (`base_coalition`) by removing the predictor group (`data_source`) 
  # from the coalition. This represents the coalition without the predictor group of interest.
  df_with_base <- cindices_all_combos_200btsmp %>%
    dplyr::mutate(base_coalition = extract_base_coalition(coalition, data_source))
  
  # Determine the total number of predictor groups (data sources).
  N_all <- df_source_dict$data_source %>% length
  
  # Identify coalition pairs:
  # - One coalition contains the predictor group of interest (e.g., "biom_core").
  # - The other is the same coalition without the predictor group (e.g., "core").
  paired_df <- df_with_base %>%
    inner_join(df_with_base, by = c("base_coalition" = "base_coalition", "bootsmp_num" = "bootsmp_num"), suffix = c(".x", ".y")) %>% 
    filter(coalition.x != coalition.y & str_detect(coalition.x, data_source)) %>%
    dplyr::select(
      bootsmp_num,
      coalition_with_ds = coalition.x, cindex_with_ds = cindex.x,
      coalition_without_ds = coalition.y, cindex_without_ds = cindex.y
      ) %>% 
    # Calculate the difference in c-index between coalitions, compute the Shapley weight, and derive the Shapley term.
    dplyr::mutate(
      cindex_diff = cindex_with_ds - cindex_without_ds, # Incremental contribution of the predictor group.
      N_subset_S = str_count(coalition_without_ds, "_") + 1, # Size of the coalition without the predictor group.
      weight = factorial(N_subset_S) * factorial(N_all - N_subset_S - 1) / factorial(N_all), # Compute the Shapley weight based on coalition size.
      shap_term = weight * cindex_diff # Weighted contribution (Shapley term).
      )
  
  # Aggregate the Shapley terms for the predictor group across bootstrap samples.
  shap_of_datasource <- paired_df %>%
    group_by(bootsmp_num) %>% # Group by the same bootstrap sample
    summarize(shap_of_datasource = sum(shap_term, na.rm = TRUE)) %>% 
    dplyr::mutate(data_source=data_source)

  return(shap_of_datasource)
}

# Calculate Shapley values for all predictor groups
df_shap_of_datasources <- lapply(as.list(df_source_dict$data_source), function(data_source) calc_shap_of_datasource(data_source, cindices_all_combos_200btsmp)) %>% 
  rbindlist()

# Calculate the 95% confidence intervals
calculate_ci <- function(x) {
  alpha <- 0.05
  lower <- quantile(x, alpha / 2)
  upper <- quantile(x, 1 - alpha / 2)
  return(c(lower, upper))
}


# Main ----
# Express Shapley contributions within each bootstrap sample as percentages
df_shap_of_datasources <- df_shap_of_datasources %>%
  group_by(bootsmp_num) %>%
  mutate(perc = shap_of_datasource / sum(shap_of_datasource)*100)

# Calculate mean Shapley value percentages and confidence intervals for each predictor group.
df_shap_summary <- df_shap_of_datasources %>%
  group_by(data_source) %>%
  summarize(
    mean_shap_perc = mean(perc, na.rm = TRUE),
    ci_lower_shap_perc = calculate_ci(perc)[1],
    ci_upper_shap_perc = calculate_ci(perc)[2],
    mean_CI = glue::glue("{scales::percent(mean_shap_perc / 100, accuracy = 0.1)} ({scales::percent(ci_lower_shap_perc / 100, accuracy = 0.1)},{scales::percent(ci_upper_shap_perc / 100, accuracy = 0.1)})")
  ) %>% 
  # Add colour and readable labels for predictor groups
  dplyr::mutate(
    color_datasource = dplyr::case_when(
        data_source == "core" ~ "black", 
        data_source == "symp" ~ "#C27C1A",
        data_source == "gen" ~ "#64734B",
        data_source == "lifestyle" ~ "#216974",
        data_source == "biom" ~ "#983B3B",
        data_source == "medhist" ~ "#6D2952",
        TRUE ~ NA_character_
    ))  %>%
  dplyr::mutate(
    data_source = case_when(
        data_source == "core" ~ "Core demographics", 
        data_source == "symp" ~ "Symptoms",
        data_source == "gen" ~ "Polygenic score",
        data_source == "lifestyle" ~ "Lifestyle",
        data_source == "biom" ~ "PC blood tests",
        data_source == "medhist" ~ "Medical history"
    )
  )%>%
    # Reorder data sources by mean Shapley value percentage in descending order
    dplyr::mutate(data_source = fct_reorder(data_source, mean_shap_perc, .desc = TRUE))

# Save the Shapley summary table (mean and CI) to a CSV file. ----
df_shap_summary <- df_shap_summary %>% 
  dplyr::arrange(data_source) %>% 
  dplyr::select(data_source, mean_CI)
fwrite(df_shap_summary, file.path(res_plots_dir, glue::glue("tbl_shap_meanCI{tag}{tag_LDpred}{tag_Vision}.csv")), row.names=FALSE)
cat("Saved df_shap_summary ", file.path(res_plots_dir, glue::glue("tbl_shap_meanCI{tag}{tag_LDpred}{tag_Vision}.csv")), "......\n")
