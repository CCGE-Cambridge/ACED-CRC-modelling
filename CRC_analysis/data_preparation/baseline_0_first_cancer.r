# ========================================================================
# Script Title: Creates a dataset of the earliest recorded cancer (type, age, date) 
# and death (age, date) [baseline_ scripts 1 of 3]
# Author: Sam Ip
# - Processes raw UK Biobank (UKB) data to extract first cancer and death records.
# - Flags whether the first cancer is colorectal cancer (CRC).
# ========================================================================

rm(list=ls())

# Define Directories ----
dir_ACED <- "~/ACED-CRC-modelling"
setwd(file.path(dir_ACED, "CRC_analysis"))
scripts_modelling_dir <- file.path(dir_ACED, "CRC_analysis", "modelling")

# Load Configurations ----
source(file.path(scripts_modelling_dir, "hpc_config.r"))

# Load cancer code lists for CRC and Skin Non-Melanoma
medcodes_crc <- fread(fpath_codelist_cancers, select = "CRC") %>% filter(!is.na(CRC))
medcodes_crc <- as.character(medcodes_crc$CRC)
medcodes_skinNM <- fread(fpath_codelist_cancers, select = "Skin_NM") %>% filter(!is.na(Skin_NM))
medcodes_skinNM <- as.character(medcodes_skinNM$Skin_NM)


# Extract and process baseline data related to cancer and death incidences ----
# "eid", "first_cancer_date", "first_cancer_age", "death_date", "death_age", "first_type_crc"
baseline_fields <- fread(fpath_codelist_ukb_baseline) %>% 
  filter(grepl("death|cancer|date_attendassessmentcentre", desc)) # get fields with death/cancer in description
df_cancer <- get_baselinefields(fpath_baseline, baseline_fields) # eid | death_date-0.0/1.0 | cancer_date-0.0/1.0/.../17.0 | cancer_type--0.0/1.0/.../16.0 | resp. age fields
df_cancer <- as.data.frame(df_cancer)

# Remove redundant cancer_date-17.0 (no corresponding type field) ----
# cancer_type-<0.0-16.0> vs cancer_date-<0.0-17.0> -- date 17 w/o corresp type field 
df_cancer %>% filter(!is.na(`cancer_date-17.0`)) # only ONE person has a non-NA cancer_date-17.0, which is in 2019 and lacks a corresponding cancer-type
df_cancer <- df_cancer %>% dplyr::select(!`cancer_date-17.0`) 


# make sure cancer_types are as.character ----
df_cancer[grep("^cancer_type", names(df_cancer))] <- sapply(df_cancer[grep("^cancer_type", names(df_cancer))], as.character)
sapply(df_cancer, class) # check resulting class

# Coalesce ICD9 & ICD10 cancer types (they share cancer date and age fields) ----
ls_col_combined <- lapply(as.list(0:16), function(idx) {
  col9 <- paste0("cancer_type9-", idx, ".0")
  col10 <- paste0("cancer_type10-", idx, ".0")
  if (idx %in% c(0:14)) {
    col_combined <- as.data.frame(coalesce(df_cancer[, col9], df_cancer[, col10]))
  } else {
    col_combined <- as.data.frame(df_cancer[, col10])
  }
  names(col_combined) <- paste0("cancer_type-", idx, ".0")
  return(col_combined)
})

cols_cancer_type <- ls_col_combined %>% reduce(cbind)
str(cols_cancer_type)

df_cancer <- df_cancer %>% dplyr::select(!c(
  names(.)[grep("^cancer_type9", names(.))],
  names(.)[grep("^cancer_type10", names(.))]
))
df_cancer <- df_cancer %>% cbind(., cols_cancer_type)


# Cleaning incomplete/irrelevant cancer records ----
df_cancer[sapply(df_cancer, `%in%`, c("", "NA", medcodes_skinNM))] <- NA
cols_cancerdate <-  grep("^cancer_date", names(df_cancer), value=TRUE)
cols_corresp_cancertype <-  grep("^cancer_type", names(df_cancer), value=TRUE)
for (idx in 1:length(cols_cancerdate)){
  df_cancer[[cols_cancerdate[idx]]][is.na(df_cancer[[cols_corresp_cancertype[idx]]])] <- NA
}

# Check that no cancer type NA and cancer date not NA
ls_na_type_notna_date <- map2(as.list(cols_corresp_cancertype), as.list(cols_cancerdate), function(type, date){
  df_cancer %>% filter(is.na(type) & !is.na(date)) %>% dplyr::select({{type}}, {{date}})
})


# Get earliest cancer and death date recorded (& corresp age) ----
df_cancer <- collapse_cancer_records_1stcancerdeathdates(df_cancer) %>% as.data.frame()


# Extract the earliest cancer type and date for each individual
first_cancer_dates_typecolname <- apply(df_cancer[cols_cancerdate], 1, function(row) {
  first_date <- min(row, na.rm = TRUE) # get 0.0/1.0/.../16.0 for earliest cancer date
  col_first_type <- names(row)[which(row == first_date)] %>% gsub("cancer_date-", "cancer_type-", .) # could be more than 1 entry
  return(list(first_date, col_first_type))
})



# Flag if first incident cancer is CRC ---- 
library(future)
library(purrr)
library(furrr)

ls_first_typecolname <- sapply(first_cancer_dates_typecolname, "[[", 2) 

extract_first_type_crc <- function(row, type_col) {
  row_nonNAs <- row %>% dplyr::select(where(~!all(is.na(.))))
  first_type_crc <- 1 * any(as.character(unlist(row_nonNAs[, type_col])) %in% medcodes_crc)
  return(first_type_crc)
}

plan(multisession, workers = 15) # Parallelised processing with batching
  batch_size <- 1000
  result <- bind_rows(
    future_map_dfr(
      split(df_cancer, ceiling(seq_len(nrow(df_cancer)) / batch_size)), 
      function(batch) future_pmap_dfr(
        list(split(batch, seq(nrow(batch))), ls_first_typecolname[seq_len(nrow(batch))]), 
        extract_first_type_crc
      )
    )
  )
plan(sequential) # Reset the plan to sequential processing after completing the parallel tasks

df_cancer$first_type_crc <- unlist(result)
sum(df_cancer$first_type_crc) #7341
cat("DONE first_type_crc......", format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"), "\n")

# Final Dataset and Censoring ----
df_first_cancer <- df_cancer %>% dplyr::select(eid, first_cancer_date, first_cancer_age, death_date, death_age, first_type_crc)

# Censor at min of first cancer or DoD (later to add GP deduct date when linked)
df_first_cancer$censor_date <- pmin(df_first_cancer$first_cancer_date, df_first_cancer$death_date, na.rm = TRUE)

# Save as a Parquet file ----
write_parquet(df_first_cancer, file.path(res_dir, "df_first_cancer.parquet"))

cat("DONE baseline_0_first_cancer.r......\n")
print(summary(df_first_cancer))