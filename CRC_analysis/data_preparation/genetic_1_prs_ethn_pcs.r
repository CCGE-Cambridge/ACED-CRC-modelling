# =====================================================================================
# Script Title: Genetic Data Compilation (PRS, Ethnicity, Top 10 PCs)
# Author: Sam Ip
## - Combines data from `jonathan_pc_ethn.parquet` and `id_prs.txt`, which contain top 10
##   Principal Components (PCs), genetic ethnicity (based on the top 4 PCs), and PRS for each individual.
## - Requires a .fam file for joining with datasets.
# =====================================================================================

## Notes
## - Combines data from `jonathan_pc_ethn.parquet` and `id_prs.txt`, which contain 
##   Principal Components (PCs), genetic ethnicity (based on the top 4 PCs), and PRS for each individual.
## - Requires a .fam file for joining with datasets.
## =============================================================================

# Clear workspace ----
rm(list=ls())

# Initialise ----
args <- c("","","_LDpred")
assign("tag", as.character(args[1]), envir = .GlobalEnv)  # "", _OnlySymptomaticPop, _OnlySymptomaticPop_SA_AllSympSansFatigue
assign("tag_Vision", as.character(args[2]), envir = .GlobalEnv) # _sansVision
assign("tag_LDpred", as.character(args[3]), envir = .GlobalEnv) # _LDpred

cat("\n\ntag, tag_Vision, tag_LDpred ......", tag, tag_Vision, tag_LDpred, "\n\n")

# File paths & initialisation ----
dir_ACED <- "~/ACED-CRC-modelling"
scripts_modelling_dir <- file.path(dir_ACED, "CRC_analysis", "modelling")
source(file.path(scripts_modelling_dir, "hpc_init_getcoxdat.r"))
source(file.path(scripts_modelling_dir, "hpc_modelling_fns.r"))

# Link genetic data with main data using eid mapping key ----
# Read in files
key <- fread(fpath_key_fam) %>% dplyr::rename(ID_v0=Cambridge_ID, eid=UCL_ID)
prs <- fread(fpath_genetic_prs)
pcs <- read_parquet(fpath_genetic_pc_ethn)

# Unify eids
df_genetic <- merge(x = prs, y = key, by = "ID_v0", all.x = TRUE) %>% 
  dplyr::select(!"ID_v0") 

# Remove rows with missing data
cat("\n\nRemoving rows with missing data.\n")
df_genetic <- df_genetic[rowSums(is.na(df_genetic))==0,]
df_genetic <- df_genetic %>% dplyr::mutate(eid = as.character(eid))

# Join genetic data with main data
distinct_coxdat_eids <- coxdat_stacked_allages_num %>%
  dplyr::select(eid) %>%
  distinct()

eids_in_coxdat_genetic <- distinct_coxdat_eids %>% inner_join(df_genetic, by = "eid")

# Save ----
saveRDS(eids_in_coxdat_genetic, file.path(root_privatedata_dir, glue::glue("df_genetic{tag_LDpred}.rds")))
