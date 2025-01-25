#!/usr/bin/env Rscript
# Free up memory; improve space efficiency for bootstrap res

rm(list=ls())
options(future.globals.maxSize = 20 * 1024^3); cat("future.globals.maxSize is set to: ", getOption("future.globals.maxSize"), "\n")
# Initialise----
args = commandArgs(trailingOnly=TRUE)
# args <- c(5, "bootsmps_combos")
n_cores <- as.numeric(args[1])
assign("res_dir_basename", as.character(args[2]), envir = .GlobalEnv) #"bootsmps_combos_OnlySymptomaticPop"


library(broom)
library(dplyr)
library(parallel)
library(survival)

# Directory containing the RDS files
directory_path <- file.path("~/rds/rds-ccge1-hdMXhK21vco/biobank/si_aced_crc/res/20240616_EditedFerritin_nback2_nfwd2/results", res_dir_basename)
dpath_tmp_boot_res <- file.path("~/rds/hpc-work/BMJOnc_reviewerresp/tmp_boot_res", res_dir_basename)

# List all RDS files in the directory
rds_files <- list.files(directory_path, pattern = "\\.rds$", full.names = TRUE)
# rds_files <- rds_files[1:3] #HERE small for debugging

# Restrict to unprocessed files
csv_files <- list.files(dpath_tmp_boot_res, pattern = "^cindex_(.*).csv$")
csv_models <- sub("^cindex_(.*)\\.csv$", "\\1", csv_files); #head(csv_models)
rds_models <- sub("^res_bootsmps(.*)\\.rds$", "\\1", basename(rds_files)); #head(rds_models)
rds_files <- rds_files[!rds_models %in% csv_models]

cat("\nDONE List all RDS files ...... left to process:", length(rds_files), "\n")#; print(head(rds_files))

# Function to process each RDS file
is_null_model <- function(model_formula) {
  as.character(model_formula)[3] == "1"
}

process_rds_file <- function(fname) {
  # fname <- rds_files[2]; fname
  # Extract the bootsmpnum_coalition value
  bootsmpnum_coalition <- gsub("^res_bootsmps|\\.rds$", "", basename(fname))
  
  # Read the RDS file
  check <- readRDS(fname)
  
  # if (!grepl("^\\d+_$", bootsmpnum_coalition)){
  if (!is_null_model(check$model$formula)){
    # Process the model using broom::tidy
    tidy_model <- broom::tidy(check$model, conf.int = TRUE, exponentiate = TRUE)
    tidy_model <- tidy_model %>%
        dplyr::select(term, estimate, conf.low, conf.high, p.value)

    # Add the bootsmpnum_coalition column
    tidy_model <- tidy_model %>%
        mutate(bootsmpnum_coalition = bootsmpnum_coalition)%>%
         dplyr::select(bootsmpnum_coalition, everything())
    # cat("\nDONE tidy_model ......\n")#; print(tidy_model)
  } else (tidy_model <- data.frame())

  # Process the cindex value
  cindex_result <- data.frame(
    bootsmpnum_coalition = bootsmpnum_coalition,
    cindex = check$cindex
  )
  # cat("\nDONE cindex_result ......\n")#; print(cindex_result)

  # Save ----
  temp_cindex_file <- file.path(dpath_tmp_boot_res, paste0("cindex_", bootsmpnum_coalition, ".csv"))
  temp_tidy_model_file <- file.path(dpath_tmp_boot_res, paste0("model_", bootsmpnum_coalition, ".csv"))
  write.table(cindex_result, file = temp_cindex_file, sep = ",", col.names = TRUE, row.names = FALSE)
  write.table(tidy_model, file = temp_tidy_model_file, sep = ",", col.names = TRUE, row.names = FALSE)
  cat("\nDONE process_rds_file ......", bootsmpnum_coalition, " ......\n")
}

# Process each file
parallel::mclapply(rds_files, function(file_path) {
    tryCatch({
      process_rds_file(file_path)
    }, error = function(e) {
      cat("Error processing file:", file_path, "\n", conditionMessage(e), "\n")
    })
  }, mc.cores = n_cores)

cat("\nDONE mclapply process_rds_file ......\n")

# # Debug ----
# for (file_path in rds_files) {
#   cat("\n", basename(file_path), "\n")
#   tryCatch({
#       process_rds_file(file_path)
#     }, error = function(e) {
#       cat("Error processing file:", file_path, "\n", conditionMessage(e), "\n")
#     })
#   }

# # Merge ----
# merge_csv_files <- function(temp_dir, pattern, output_file) {
#   temp_files <- list.files(dpath_tmp_boot_res, pattern = pattern, full.names = TRUE)
#   combined_data <- do.call(rbind, lapply(temp_files, read.csv, stringsAsFactors = FALSE))
#   write.csv(combined_data, output_file, row.names = FALSE)
# }
# merge_csv_files(dpath_tmp_boot_res, "cindex_*.csv", glue::glue("all_cindex_results_{res_dir_basename}.csv"))
# merge_csv_files(dpath_tmp_boot_res, "model_*.csv", glue::glue("all_model_results_{res_dir_basename}.csv"))
