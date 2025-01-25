# ========================================================================
# Script Title: Extract & wrangle UKB baseline assessment data [baseline_ scripts 3 of 3]
# Author: Sam Ip
# - Preceded by baseline_1_loaddata.r
# - Further extracts and processes UKB baseline assessment predictors
# - Outputs df_baseline.parquet
# ========================================================================

# Source UKB baseline variables recoding dictionaries ----
source(file.path(codelist_dir, "ukb_baseline_dicts.r"))

# Functions  ----

# Ensure factor levels are ordered and fill NAs with "missing"
mk_factor_orderlevels <- function(df, colname){
  df <- df %>% dplyr::mutate(
    !!sym(colname) := factor(!!sym(colname), levels = str_sort(unique(df[[colname]]), numeric = TRUE)))
  df[[colname]] <- na.fill(df[[colname]], "missing") 
  return(df)
}

# Recodes categorical variables and sets levels, reference, and fills missing values
categorical_recode_relevel_fillna <- function(df, colname, dict, ref_level, na_level, order_levels=NA){
  df <- df %>% mutate(
    !!sym(colname) := dplyr::recode(!!sym(colname), !!!setNames(as.character(dict$new), dict$original)))
  df[[colname]] <- na.fill(df[[colname]], na_level) 
  if (length(order_levels) > 1){
    df[[colname]] <- factor(df[[colname]], levels = order_levels, ordered=TRUE)
  } else {
    df[[colname]] <- as.factor(df[[colname]])
    if(!is.na(ref_level)){
      df[[colname]] <- relevel(df[[colname]], ref = ref_level)
    }
  }
  print(levels(df[[colname]]))
  return(df)
}

# Wrangle predictors ----

### Smoking status ----
df_baseline <- categorical_recode_relevel_fillna(
  df_baseline, "smoking_status-0.0", dict_smoking_status, ref_level = "never", na_level="missing")

### Ethnicity (self-reported) ----
df_baseline <- categorical_recode_relevel_fillna(df_baseline, "ethn_selfrep-0.0", dict_ethn_selfrep, ref_level = "White", na_level="missing")

### Townsend score ----
df_baseline$townsend <- na.fill(df_baseline$`townsend-0.0`, median(df_baseline$`townsend-0.0`, na.rm=TRUE)) 
any(is.na(df_baseline$townsend))
df_baseline <- df_baseline %>% dplyr::select(!`townsend-0.0`)

### Family history ----
cols_famhist <- df_baseline %>% dplyr::select(grep_colnames("^famhist", df_baseline)) 
colnames_famhist_og <- grep_colnames("^famhist", cols_famhist)
for (colname in names(cols_famhist)){ # Recode family history using dictionary
  cols_famhist <- cols_famhist %>% mutate(
  !!sym(colname) := dplyr::recode(!!sym(colname), !!!setNames(as.character(dict_famhist$new), dict_famhist$original)))
}

# Filter and create binary indicators for relevant conditions 
relevant_famhist <- c("lung_cancer", "bowel_cancer", "breast_cancer")
irrelevant_famhist <- unique(dict_famhist$new )[!unique(dict_famhist$new ) %in% relevant_famhist]

for (feat in irrelevant_famhist){cols_famhist[cols_famhist==feat] <- NA}

for (feat in relevant_famhist){
  cols_famhist <- cols_famhist %>% 
    dplyr::mutate(!!sym(paste0("famhist_", feat)) := ifelse(apply(cols_famhist == feat, 1, any), 1, 0))
}
cols_famhist <- cols_famhist %>% dplyr::select(! colnames_famhist_og)
cols_famhist[is.na(cols_famhist)] <- 0
df_baseline <- df_baseline %>% dplyr::select(! colnames_famhist_og)
df_baseline <- cbind(df_baseline, cols_famhist)

# Formatting ----
### Remove -0.0 from column names ----
names(df_baseline) <- sub("\\-0.0", "", names(df_baseline))


# Add phenotypic variables ----
# Generated by CRC_variables/1_baseline_extra_phenotypic_vars.r
df_extra_pheno <- get(base::load(file.path(res_dir, "df_baseline_extra_phenotypic_vars.Rdata"))) 
df_extra_pheno <- df_extra_pheno %>% dplyr::select(
  eid, edu_quals, tot_METs_wk, alcohol_units_daily, tot_process_meat_wk, tot_red_meat_wk, partial_fibre_score)
df_baseline <- df_baseline %>% left_join(df_extra_pheno, by="eid")

# Save as Parquet file ----
write_parquet(df_baseline, file.path(res_dir, "df_baseline.parquet"))
