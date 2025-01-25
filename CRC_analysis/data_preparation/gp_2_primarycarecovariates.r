# ==============================================================================
# Script Title: Extraction of Primary Care Variables [gp_ scripts 2 of 3]
# Author: Sam Ip
# - Not dependent on preceding gp_1_reg.r
# - Executes multiple R scripts each extracting a specific type of variable from primary care records
# ==============================================================================
rm(list=ls())

# Define Directories ----
dir_ACED <- "~/ACED-CRC-modelling"
setwd(file.path(dir_ACED, "CRC_analysis"))
scripts_modelling_dir <- file.path(dir_ACED, "CRC_analysis", "modelling")

# Load Configurations ----
source(file.path(scripts_modelling_dir, "hpc_config.r"))

# Paths to R scripts for preprocessing tasks ----
scripts <- list(
  biomarkers = "~/ACED-RREDD-EHR/CRC_variables/2_gp_biomarker_march2022.R",       # Biomarker data
  symptoms = "~/ACED-RREDD-EHR/CRC_variables/2_markers_of_risk.R",                # Symptoms data
  medhist_ever = "~/ACED-RREDD-EHR/CRC_variables/2_modifiers_of_risk.R",          # Medical history (ever)
  multimorb = "~/ACED-RREDD-EHR/CRC_variables/2_gp_multimorb_HH_May22.R",         # Multimorbidity data (Prior to age/sex regression)
  other_events = "~/ACED-RREDD-EHR/CRC_variables/2_gp_other_events_HH_NOV2022.R"  # Other events (e.g., colonoscopy)
)

# Execute Scripts in Parallel ----
# Directory for nohup output logs
nohup_dir <- "~/ACED-RREDD-EHR/sam/nohups/"
dir.create(nohup_dir, showWarnings = FALSE, recursive = TRUE)  # Create directory if it doesn't exist

# Execute Scripts (Each parallelised within itself) ----
lapply(scripts, function(fpath) {
  # Extract script name without extension for log file naming
  script_name <- tools::file_path_sans_ext(basename(fpath))
  
  # Construct nohup command
  nohup_command <- paste0(
    "nohup Rscript ", fpath, 
    " > ", nohup_dir, script_name, ".out 2>&1 &"
  )
  
  # Execute command
  system(nohup_command)
  
  # Log execution
  cat("Started:", fpath, "\nLog file:", nohup_dir, script_name, ".out\n\n")
})

