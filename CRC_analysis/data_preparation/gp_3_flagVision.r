# ==============================================================================
# Script Title: Extract Data Provider "Vision" Patients' EIDs [gp_ scripts 3 of 3]
# Author: Sam Ip
# - Not dependent on preceding gp_2_primarycarecovariates.r
# - Extracts EIDs for Vision patients from GP registration data
#   Ref: https://biobank.ndph.ox.ac.uk/ukb/refer.cgi?id=591 Section 6 P.16
#   data_provider 1= England (Vision), 2= Scotland, 3 = England (TPP), 4 = Wales
# - Creates a new results directory for analyses excluding Vision patients
# ==============================================================================
rm(list=ls())

# Define Directories ----
dir_ACED <- "~/ACED-CRC-modelling"
setwd(file.path(dir_ACED, "CRC_analysis"))
scripts_modelling_dir <- file.path(dir_ACED, "CRC_analysis", "modelling")

# Load Configurations ----
source(file.path(scripts_modelling_dir, "hpc_config.r"))

# Load GP registration data ----
gp_reg <- fread(fpath_gp_reg)

# Extract EIDs of Vision patients (data_provider == 1) ----
eids_Vision <- gp_reg %>% dplyr::filter(data_provider=="1") %>% dplyr::select(eid)

# Save Vision patients' EIDs as a CSV file ----
fwrite(eids_Vision, file.path(shared_dir, "eids_Vision.csv"), row.names=FALSE)
