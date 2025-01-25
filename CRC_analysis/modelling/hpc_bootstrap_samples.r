# =====================================================================================
# Script Title: Generating Bootstrap Samples and Train-Test Splits
# Author: Sam Ip
# Description: Generates bootstrap samples and splits data into training and testing sets, 
#              saving the results for downstream analysis.
# =====================================================================================
# Define Directories ----
dir_ACED <- "~/ACED-CRC-modelling"
setwd(file.path(dir_ACED, "CRC_analysis"))
scripts_modelling_dir <- file.path(dir_ACED, "CRC_analysis", "modelling")

# Load dependencies ----
source(file.path(scripts_modelling_dir, "hpc_init_getcoxdat.r"))

# Configuration ----
n_bootstrap_samples <- 200 # Number of bootstrap samples
tag <-  "" #"", _OnlySymptomaticPop, _OnlySymptomaticPop_SA_AllSympSansFatigue
tag_Vision <- "" # _sansVision
tag_LDpred <- "" # _LDpred

# Generate bootstrap samples and save eids ----
set.seed(137)
bootstrap_samples_eid <- list()
train_test_split_ratio <- 0.7  # Define the train-test split ratio

for (i in 1:n_bootstrap_samples) {
  # Split unique IDs into train and test sets, then sample with replacement in reach set
  unique_ids <- unique(coxdat_stacked_allages_num$eid)
  n_train <- floor(length(unique_ids) * train_test_split_ratio)
  train_ids <- sample(unique_ids, n_train, replace = TRUE)
  test_ids <- sample(unique_ids[!(unique_ids %in% train_ids)], length(unique_ids) - n_train, replace = TRUE)

  # Save the train and test IDs
  bootstrap_samples_eid[[i]] <- list(
    train = train_ids,
    test = test_ids
  )}

saveRDS(bootstrap_samples_eid, file = file.path(res_dir, "bootstrap_samples_eid.rds"))
