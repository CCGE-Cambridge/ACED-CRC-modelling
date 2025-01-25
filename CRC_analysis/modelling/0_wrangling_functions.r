# =====================================================================================
# Script Title: Data Preprocessing Utilities
# Author: Sam Ip
# Description:
# This script provides a set of utility functions for preprocessing survival and phenotypic data.
# =====================================================================================


# Read baseline assessment data file, extracts fields based on a mapping table, and 
# returns a cleaned dataframe with renamed columns.
# - baseline_fields -- dictionary with columns c("field","desc") 
#                      corresponding values c(34,"birth_year"), for example
get_baselinefields <- function(fpath_baseline, baseline_fields) {
    # Read phenotype numbers ("field") from the baseline file
    pheno_nums <- names(fread(fpath_baseline, nrows=1))

    # Generate a list of numeric and descriptive field names
    ls_field_nums_and_names <- lapply(
        split(baseline_fields, seq(nrow(baseline_fields))), 
        function(baseline_field) {
            # Extract field numbers from `pheno_nums` that match the given field pattern.
            # E.g.: For `field = 34`, matches columns like "34-0.0", "34-1.0", etc.
            field_nums <- pheno_nums[grepl(paste0("^", baseline_field$field, "-"), pheno_nums)]

            # Replace the numeric `field` with the descriptive `desc` in the matched field names.
            # E.g.: Replace "34-0.0" with "birth_year-0.0".
            field_names <- gsub(baseline_field$field, baseline_field$desc, field_nums)
            
            # Return a df containing both field numbers and descriptive field names.
            return(data.frame(field_nums = field_nums, field_names = field_names))     
             }
    )

    fields_df <- do.call(rbind, ls_field_nums_and_names)
    field_nums <- fields_df$field_nums
    field_names <- fields_df$field_names

    # Read the selected columns from the baseline file
    df_baseline <- fread(fpath_baseline, select=c("eid", field_nums))
    names(df_baseline) <- c("eid", field_names) # Rename with descriptive field names

    return(df_baseline)
}

# Collapse multiple cancer and death date fields into single fields for
# the earliest recorded cancer and death.
collapse_cancer_records_1stcancerdeathdates <- function(df_cancer) {
    # Convert date columns to Date type
    df_cancer <- df_cancer %>% mutate_at(
        vars(grep("date", names(df_cancer), value = TRUE)), funs(as.Date(.,"%Y-%m-%d")))

    # Calculate the first cancer and death dates/ages
    df_cancer <- df_cancer %>% mutate(
        first_cancer_date = pmin(!!!rlang::syms(grep_colnames("cancer_date-", df_cancer)), na.rm=TRUE), 
        first_cancer_age = pmin(!!!rlang::syms(grep_colnames("cancer_age-", df_cancer)), na.rm=TRUE),
        death_date = pmin(!!!rlang::syms(grep_colnames("death_date-", df_cancer)), na.rm=TRUE),
        death_age = pmin(!!!rlang::syms(grep_colnames("death_age-", df_cancer)), na.rm=TRUE)
        ) 
    return(df_cancer)
}

# Helper function to search for column names containing a specific substring.
grep_colnames <- function(col_string, df) {
    return(names(df)[grepl(col_string, names(df))])
}

# Generates the powerset of a given set (all possible non-empty subsets).
get_powerset <- function(set) {
  n <- length(set) # Number of elements in the input set
  masks <- 2^(1:n-1) # Binary masks for each element in the set
  lapply( 1:2^n-1, function(u) set[ bitwAnd(u, masks) != 0 ] ) # Generate all subsets except the empty set
}

# Identifies columns in a dataframe that contain missing values.
which_cols_hv_NAs <- function(df){
    return(names(which(colSums(is.na(df))>0)))
}


# Filters data to include individuals with any ‘symptoms’ predictors, excluding fatigue
filter_coxdat_SA_AllSympSansFatigue <- function(coxdat_stacked_allages_num, tag){
    # str(coxdat_stacked_allages_num)
    if (tag=="_OnlySymptomaticPop_SA_AllSympSansFatigue"){
        # Select relevant symptom columns
        symp_selected <- grep("^symp_", colnames(coxdat_stacked_allages_num), value=TRUE) # All symptoms
        symp_selected <- symp_selected[!symp_selected %in% c("symp_FATIGUE_new_onset")] # Except fatigue
        cat("\n\nsymp_selected ......", symp_selected, "\n")

        # Identify symptomatic individuals
        ind_symptomatic <- coxdat_stacked_allages_num %>% dplyr::select(all_of(c("eid", "core_age", symp_selected))) %>%
            dplyr::mutate(across(all_of(symp_selected), function(x) {as.numeric(as.character(x))}))  %>%
            dplyr::mutate(any_mrk = select(., all_of(symp_selected)) %>% rowSums(na.rm = TRUE), 
                symptomatic=factor(1*(any_mrk >  0), levels=c("0", "1"))) 

        # Filter data for symptomatic individuals
        coxdat_stacked_allages_num <- coxdat_stacked_allages_num %>% 
            left_join(ind_symptomatic %>% dplyr::select(eid,core_age,symptomatic), by = c("eid", "core_age")) %>%
            dplyr::filter(symptomatic=="1") %>%
            dplyr::select(!symptomatic)

        # Summary
        cat("N_symp: ", coxdat_stacked_allages_num$eid %>% unique() %>% length(), "; N_symp_cases: ",
        (coxdat_stacked_allages_num %>% filter(status==1))$eid %>% unique() %>% length(), "; dim: ", dim(coxdat_stacked_allages_num),"\n")
    }

    return(coxdat_stacked_allages_num)
}


# Filters the dataset to include only symptomatic individuals based on a predefined 
# list of symptoms loaded from an external file. (On HPC)
filter_coxdat_OnlySymptomaticPop_hpc <- function(coxdat_stacked_allages_num, tag){
    if (tag=="_OnlySymptomaticPop"){
        # Load predefined symptom list from external file
        symp_selected <- readRDS(file=file.path(res_plots_dir_tagfree, glue::glue("define_symptomatic"), glue::glue("symps_define_symptomatic.rds"))) #HERE
        cat("\n\nsymp_selected ......", symp_selected, "\n")

        # Identify individuals with symptoms
        ind_symptomatic <- coxdat_stacked_allages_num %>% dplyr::select(all_of(c("eid", "core_age", symp_selected))) %>%
            dplyr::mutate(across(symp_selected, function(x) {as.numeric(as.character(x))}))  %>%
            dplyr::mutate(any_mrk = select(., all_of(symp_selected)) %>% rowSums(na.rm = TRUE), 
                symptomatic=factor(1*(any_mrk >  0), levels=c("0", "1"))) 

        # Filter the dataset to retain only symptomatic individuals
        coxdat_stacked_allages_num <- coxdat_stacked_allages_num %>% 
            left_join(ind_symptomatic %>% dplyr::select(eid,core_age,symptomatic), by = c("eid", "core_age")) %>%
            dplyr::filter(symptomatic=="1") %>%
            dplyr::select(!symptomatic)

        # Summary
        cat("N_symp: ", coxdat_stacked_allages_num$eid %>% unique() %>% length(), "; N_symp_cases: ",
        (coxdat_stacked_allages_num %>% filter(status==1))$eid %>% unique() %>% length(), "; dim: ", dim(coxdat_stacked_allages_num),"\n")
    }

    return(coxdat_stacked_allages_num)
}

# Rounds a numeric value to the specified number of decimal places.
my_round <- function(x, dp){format(round(x,dp), nsmall = dp) }

# Computes the 95% confidence interval for a numeric vector, rounds the 
# results to the specified number of decimal places, and formats the output.
CI_my_round <- function(x, dp){
  x <- Rmisc::CI(x,ci=0.95)
  format(round(x,dp), nsmall = dp) }

# Replaces numeric education codes with corresponding human-readable education levels.
education_dict <- c(
    "1" = "College/University", "2" = "A/AS levels", "3" = "O levels/GCSEs", 
    "4" = "CSEs", "5" = "NVQ/HND/HNC", "6" = "Other", "7" =  "None")

replace_education <- function(df, column_name, ind_abbrev) {
    # Update the dictionary for abbreviated labels if applicable
    if (ind_abbrev){
        education_dict <- c("1" = "Col/Univ", "2" = "A/AS", "3" = "O/GCSE", "4" = "CSEs", "5" = "NVQ/HND/HNC", "6" = "Other", "7" =  "None")
    }

    # Replace numeric codes with corresponding education levels
    replace_education_in_column <- function(x) {
        is_education <- str_detect(x, "^Education\\.")
        num <- str_extract(x[is_education], "\\d+$")
        new_value <- education_dict[num]
        x[is_education] <- ifelse(is.na(new_value), x[is_education], str_replace(x[is_education], paste0("\\.", num), paste0(".", new_value)))
        return(x)
    }
    df[[column_name]] <- as.character(df[[column_name]])
    df[[column_name]] <- replace_education_in_column(df[[column_name]])

    # Convert education column to a factor with updated levels
    df <- df %>% dplyr::mutate(covariate_readable=factor(covariate_readable, levels=unique(df$covariate_readable)))
    
    return(df)
}

# Formats a dataframe containing model results by cleaning term names, 
# assigning readable names, mapping data sources, and setting color schemes 
# for visualisation.
fmt_df_mdl_datasrc <- function(df_mdl, mdl){
    # Add significance indicators and sort by significance and estimate
    df_mdl <- df_mdl %>% dplyr::mutate(
        sig = 1*(p.value <0.05),
        estimate_sign_sig = as.factor(sign(estimate)*sig)
    ) %>% dplyr::arrange(estimate_sign_sig, estimate)

    # Add readable names and data sources
    readable_names <- readable_names %>% dplyr::mutate(data_source=sub("\\_.*", "", covar_code2))
    if ("core_sex_genetic1" %in% df_mdl$term) { # Correct term names if necessary
        df_mdl <- df_mdl %>% dplyr::mutate(term = sapply(term, correct_term, covar_codes = readable_names$covar_code2))
    }
    df_mdl <- df_mdl %>%
        dplyr::mutate(term_root = sub("\\..*", "", term)) %>%
        dplyr::left_join(readable_names %>% select(covar_code2, data_source), by = c("term_root" = "covar_code2")) %>%
        dplyr::left_join(readable_names %>% select(covar_code2, covariate_readable), by = c("term_root" = "covar_code2")) %>%
        dplyr::mutate(
            covariate_readable = coalesce(covariate_readable, term),
            covariate_readable = stringr::str_replace_all(
                covariate_readable, 
                setNames(readable_names$covariate_readable, readable_names$covariate_code)
            ),
            data_source = coalesce(data_source, sub("\\_.*", "", term))
        )

    # Filter out unwanted terms and adjust readable names
    df_mdl <- df_mdl %>%
        dplyr::filter(!grepl("^gen_pc", term_root)) %>%
        dplyr::mutate(
            covariate_readable_root = covariate_readable,
            covariate_readable = paste0(covariate_readable, trimws(term, "left", "\\w")),
            covariate_readable = ifelse(term_root %in% binarycols, sub("\\..*", "", covariate_readable), covariate_readable),
            covariate_readable = ifelse(term_root == "core_sex_genetic", "Male (genetic)", covariate_readable)
        ) #https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=22001

    # Add colours based on predictor groups and adjust factors
    df_mdl <- df_mdl %>%
        dplyr::mutate(
            color_datasource = dplyr::case_when(
                data_source == "core" ~ "black", 
                data_source == "symp" ~ "#C27C1A",
                data_source == "gen" ~ "#64734B",
                data_source == "lifestyle" ~ "#216974",
                data_source == "biom" ~ "#983B3B",
                data_source == "medhist" ~ "#6D2952",
                TRUE ~ NA_character_
            )) %>%
        dplyr::mutate(
            data_source = ifelse(data_source == "mrk", "symp", data_source),
            covariate_readable = factor(covariate_readable, levels = df_mdl$covariate_readable),
            data_source = factor(data_source, levels = unique(df_mdl$data_source))
        ) %>%
        dplyr::arrange(data_source, as.character(covariate_readable))
    
    df_mdl <- df_mdl %>%
        dplyr::mutate(color_datasource = factor(color_datasource, levels = unique(df_mdl$color_datasource)))


    cat("\ndf_mdl......\n"); print(df_mdl)
    cat("\nData sources......", levels(df_mdl$data_source), "...... \n")
    
    return(df_mdl)
}

# Identifies predictor levels with no associated CRC events.
eventless_riskfactor_levels <- function(ind_symptomatic_only){
    # Convert character columns to factors
    coxdat_stacked_allages_num <- coxdat_stacked_allages_num %>% 
        dplyr::mutate_if(is.character, as.factor) %>% 
        dplyr::select(-any_of(c("nelaa")))

    # Define relevant factor columns and ensure they are factors
    factor_cols <- c(binarycols, "core_smoking_status", "core_ethn_selfrep", "lifestyle_edu_quals", "symp_ABDO_BLOAT_freq_in_lb", "symp_ABDO_PAIN_freq_in_lb", "symp_PELVIC_PAIN_freq_in_lb", "symp_STOM_DIS_freq_in_lb") %>% intersect(.,colnames(coxdat_stacked_allages_num))
    coxdat_stacked_allages_num <- coxdat_stacked_allages_num %>% dplyr::mutate(across(all_of(factor_cols), function(col) {as.factor(col)}))
    str(coxdat_stacked_allages_num)

    # Count cases for each factor level
    cases_onehot_count <- colSums(model.matrix(~., coxdat_stacked_allages_num %>% dplyr::filter(status==1) %>% dplyr::select(factor_cols)))
    
    # Identify factor levels with zero cases
    df_eventless_cols <- data.frame(eventless_cols=cases_onehot_count[cases_onehot_count==0] %>% names())
    
    # Correct column names for readability if necessary
    if ("core_sex_genetic1" %in% names(cases_onehot_count)) {
        df_eventless_cols <- df_eventless_cols %>% dplyr::mutate(
        eventless_cols=stringr::str_replace_all(df_eventless_cols$eventless_cols, setNames(glue::glue("{readable_names$covar_code2}."), 
            readable_names$covar_code2)),
        eventless_cols=gsub("[.]$", "", eventless_cols))
        }
    
    return(df_eventless_cols$eventless_cols)
}

# Capitalizes the first letter of a string while trimming any leading 
# or trailing whitespace.
my_cap <- function(x) {
  x <- trimws(x) # Remove leading and trailing whitespace
  glue::glue("{toupper(substring(x, 1,1))}{substring(x, 2)}") # Capitalise the first letter
}
