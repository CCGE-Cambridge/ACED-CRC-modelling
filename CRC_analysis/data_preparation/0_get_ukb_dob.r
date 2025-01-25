# =====================================================================================
# Script Title: Extract and process UKB baseline data to calculate the date of birth (DoB)
# Author: Sam Ip
# =====================================================================================

# Define Directories ----
dir_ACED <- "~/ACED-CRC-modelling"
setwd(file.path(dir_ACED, "CRC_analysis"))
scripts_modelling_dir <- file.path(dir_ACED, "CRC_analysis", "modelling")

# Load Configurations ----
source(file.path(scripts_modelling_dir, "hpc_config.r"))

# Extract DoB from Baseline Data ----
# Read UKB baseline data with relevant columns: eid, year of birth (34-0.0), month of birth (52-0.0)
df_baseline <- fread(fpath_baseline, select = c("eid", "34-0.0", "52-0.0"), colClasses = "character")

# Rename columns for clarity and consistency
names(df_baseline) <- make.names(c("eid", "yob", "mob"))

# Calculate DoB using year and month of birth (assuming the first day of the month)
df_baseline$dob <- as.Date(ISOdate(df_baseline$yob, df_baseline$mob, 1), "%Y-%m-%d")

# Save as a parquet file ----
write_parquet(df_baseline, file.path(res_dir, "dob_from_baseline.parquet"))


# Create a version of HH's df_dob (using gen_birth_dates) for covariate scripts ----
source(file.path(dpath_hpc_CRCvariables, "gen_functions.R")) 

# Read baseline data with additional columns: sex (31-0.0), assessment date (53-0.0)
df_dob <- fread(fpath_baseline, select = c("eid", "31-0.0", "34-0.0", "52-0.0", "53-0.0"), colClasses = "character")

names(df_dob) <- make.names(c("eid", "sex", "yob", "mob", "assessment_date"))
df_dob$dob <- gen_birth_dates(df_dob$yob, df_dob$mob)
df_dob$assessment_date <- as.Date(df_dob$assessment_date, format = "%Y-%m-%d")
df_dob <- df_dob %>% select(eid, sex, assessment_date, dob)

# Save as a Rdata file ----
save(df_dob, file = file.path(res_dir, "df_dob.Rdata"))