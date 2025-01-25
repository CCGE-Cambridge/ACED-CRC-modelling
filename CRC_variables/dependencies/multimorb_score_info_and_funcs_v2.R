
##File name:multimorb_score_info_and_funcs.r
##Author: Hannah Harrison
##Last Edit: 28/01/2022
##Description: imports cam_multi_morb score info including coeffs and 
##defines functions for generating multimorb score

##load libraries (now loaded in config)
#library(tidyverse) 
#library(readxl) 
#library(dplyr)
#library(data.table)

##load any functions files needed (now added in config file)
#source("~/ACED-RREDD-EHR/CRC_variables/dependencies/biomarker_cleaning.r") #defines functions that clean biomarker variables

#general description of cam multi-morbidity models
dir_hh <- file.path(dir_ACED, "hannah") 
cam_multi_morb_file <- file.path(dir_hh, "comorbidity/cam_multi_morb_score_coeffs.csv")

df_cam_multi_morb <- fread(cam_multi_morb_file, header = TRUE)

#definition of events
##from clinical data
comorb_ever <- c("ALC", "ANO", "AST", "ATR", "BLI", "BRO", "CLD", "SIN", "COP", "CHD", "DEM", "DIB", 
                    "DIV", "EPI", "HEL", "HEF", "HYP", "IBD", "IBS", "LEA", "MSC", "PRK", "PEP", "PVD", 
                    "PRO", "PSO", "PSM", "RHE", "SCZ", "STR", "THY")
comorb_in_last_year <- c("ANX", "DEP")
comorb_in_last_5_years <- c("CAN") #keep only first cancer diagnosis
comorb_from_measures <- c("CKD") #chronic kideny disease looks at filtration rate

##from prescription data
script_ever <- c("SCZ")
script_1_in_last_year <- c("AST", "EPI")
script_4_in_last_year <- c("ANX",  "DEP", "CON", "IBS", "MIG", "PNC", "PSO")


events_relative_to_date <- function(df_events, df_EHR, date_column){
    ##for a particular date, removes events after date and generates date variables 1 year and 5 years in past
    df_events_pre_cutoff <- left_join(df_events, df_EHR, by = "eid") %>% 
                                filter(event_dt < .data[[date_column]]) %>% 
                                    mutate(my_date_minus_1_year = as.Date(.data[[date_column]]) - 365) %>%
                                        mutate(my_date_minus_5_years = as.Date(.data[[date_column]]) - 5*365)

    return(df_events_pre_cutoff)
}

define_CKD <- function(df_clinical_pre_cutoff) {
    ##defines chronic kidney disease for mm score
    ##CKD: event if highest of last two GFR tests < 60 mL/min
    ##throws "NAs introduced by coercion warning" 
    df_ckd_only <- df_clinical_pre_cutoff %>% 
                        filter(Comorb == "CKD") %>% 
                            mutate(biomarker = Comorb) 

    #extracts and manages any rescaling of GFR test results
    df_ckd_only <- value_extract(df_ckd_only) 
    df_ckd_only <- trim_dist_and_rm(df_ckd_only, "CKD", c(0.1, 150)) #removes extreme GFR measurements

    #drops any events that doen't meet CKD criteria
    df_ckd_only <- df_ckd_only %>% 
                    group_by(eid) %>% 
                        arrange(desc(event_dt), by_group = TRUE) %>% 
                            filter(n() != 1) %>%
                                filter(row_number() %in% c(1, 2)) %>%
                                    filter(value == max(value)) %>%
                                        filter(value < 60) %>% 
                                            ungroup %>%
                                                distinct(eid, Comorb, .keep_all = FALSE) 
    return(df_ckd_only)
}

define_CAN <- function(df_clinical_pre_cutoff){
    ##defines "cancer" event for mm score
    ##CAN: event if first cancer diagnosis was in last 5 years
    df_can_only <- df_clinical_pre_cutoff %>% 
                    filter(Comorb == "CAN") %>% 
                        group_by(eid) %>%
                            arrange(event_dt, by_group = TRUE) %>% 
                                filter(n() == 1) %>%
                                    filter(event_dt > my_date_minus_5_years) %>% 
                                        ungroup %>%
                                            select(eid, Comorb) 

    return(df_can_only)
}

gen_mm_age_categories <- function(df_inc_age) {
    #creates categorical age variables that are needed for the mm socre
    df_inc_age_cats <- df_inc_age %>% 
                mutate(age_31_40 = ifelse(age <= 40, 1, 0),
                        age_41_50 = ifelse(age > 40 & age <= 50, 1, 0),
                        #age_51_60 = ifelse(age > 50 & age <= 60, 1, 0),
                        age_61_70 = ifelse(age > 60 & age <= 70, 1, 0),
                        age_71_80 = ifelse(age > 70 & age <= 80, 1, 0),
                        age_81_max = ifelse(age > 80 , 1, 0))
    return(df_inc_age_cats)
}

define_37_mm_covariates <- function(df_clinical_and_script_events) {
    #takes both the clinical (gp) and prescription events and combines them into the 37 morbidity variables needed in mm_score
    list_of_columns <- c("clinic_ALC", "clinic_ANO", "clinic_AST", "script_AST", "script_DEP", "script_ANX",
     "clinic_DEP",  "clinic_ANX","clinic_ATR", "clinic_BLI", "clinic_BRO", "clinic_CAN", "clinic_CKD", 
     "clinic_CLD", "clinic_SIN", "script_CON", "clinic_COP", "clinic_CHD", "clinic_DEM", "clinic_DIB", 
     "clinic_DIV", "clinic_EPI", "script_EPI", "clinic_HEL", "clinic_HEF", "clinic_HYP", "clinic_IBD", 
      "clinic_IBS", "script_IBS", "clinic_LEA", "script_MIG", "clinic_MSC", "script_PNC", "clinic_PRK", 
     "clinic_PEP", "clinic_PVD", "clinic_PRO", "clinic_PSM", "clinic_RHE", "clinic_SCZ", "script_SCZ", 
     "clinic_STR", "clinic_THY")
    
    missing_cols <- list_of_columns[!list_of_columns %in% colnames(df_clinical_and_script_events)]
    for(name in missing_cols) df_clinical_and_script_events <- df_clinical_and_script_events %>% mutate(!! name := 0)

    df_37_mm_covars <- df_clinical_and_script_events %>%
        mutate(event_ALC = clinic_ALC,
                event_ANO = clinic_ANO,
                event_AST = case_when(clinic_AST == 1 & script_AST >= 1 ~ 1, TRUE ~ 0), 
                event_ANX_DEP = case_when((script_DEP >= 4 | script_ANX >= 4 | clinic_DEP == 1 | clinic_ANX == 1) ~ 1, TRUE ~ 0),
                event_ATR = clinic_ATR,
                event_BLI = clinic_BLI,
                event_BRO = clinic_BRO,
                event_CAN = clinic_CAN,
                event_CKD = clinic_CKD,
                event_CLD = clinic_CLD,
                event_SIN = clinic_SIN,
                event_CON = case_when(script_CON>= 4 ~ 1, TRUE ~0),
                event_COP = clinic_COP,
                event_CHD = clinic_CHD,
                event_DEM = clinic_DEM,
                event_DIB = clinic_DIB,
                event_DIV = clinic_DIV,
                event_EPI = case_when(clinic_EPI == 1 & script_EPI >=1 ~ 1, TRUE ~ 0),
                event_HEL = clinic_HEL,
                event_HEF = clinic_HEF,
                event_HYP = clinic_HYP,
                event_IBD = clinic_IBD, 
                event_IBS = case_when(clinic_IBS == 1 | script_IBS >= 4 ~ 1, TRUE ~ 0),
                event_LEA = clinic_LEA, 
                event_MIG = case_when(script_MIG >= 4 ~ 1, TRUE ~ 0),
                event_MSC = clinic_MSC,
                event_PNC = case_when(script_PNC >= 4 | (script_EPI >= 4 & clinic_EPI == 0) ~ 1, TRUE ~ 0), #incorrect number of epi scripts
                event_PRK = clinic_PRK,
                event_PEP = clinic_PEP, 
                event_PVD = clinic_PVD, 
                event_PRO = clinic_PRO, 
                event_PSO = case_when(clinic_PSO == 1 & script_PSO >= 4 ~ 1, TRUE ~ 0),
                event_PSM = clinic_PSM,
                event_RHE = clinic_RHE, 
                event_SCZ = case_when(clinic_SCZ == 1 | script_SCZ > 0 ~ 1, TRUE ~ 0),
                event_STR = clinic_STR,
                event_THY = clinic_THY) %>%
        select(starts_with(c("eid", "event", "age", "sex")))
    return(df_37_mm_covars)
}

mm_score_HR <- function(mm_covars, score_info, model_col) {
#applies the coefficients of the specified type of mm_score to the pre-calculated morbidity variables
#returns all 37 event variables, mm score and residual score (removing effects of age and sex)    
    #  mm_covars<-df_mm_covars; score_info<-df_cam_multi_morb; model_col<-"hospital_37_HR"
    
    my_model <- score_info %>% dplyr::select(codes, !!model_col)
    colnames(my_model) <- c("codes", "HR")
    my_model <- my_model %>% mutate(coeff = log(HR)) %>% mutate(coeff = ifelse(is.na(coeff), 0, coeff))
    comorb_col_idx <- grep(pattern="^event", x=colnames(mm_covars))
    mm_covars <- mm_covars %>% mutate(sum_temp = 0)
    

    for (i in comorb_col_idx) {
        event_name <- colnames(mm_covars[, i])
        comorb_code <- gsub("event_", "", event_name)
        coeff_temp <- my_model %>% filter(codes == comorb_code) %>% dplyr::select(coeff)
        mm_covars <- mm_covars %>% mutate(sum_temp = sum_temp + (!!as.name(event_name) *  unlist(coeff_temp)))
    }


    age_col_idx <- grep(pattern="^age_", x=colnames(mm_covars))
    #  mm_covars[,age_col_idx] %>% View()
    for (i in age_col_idx) {
        # i<- 40
        age_type <- colnames(mm_covars[, i])
        coeff_temp <- my_model %>% filter(codes == age_type) %>% dplyr::select(coeff)
        mm_covars <- mm_covars %>% mutate(sum_temp = sum_temp + (!!as.name(age_type) * unlist(coeff_temp)))
    }

    coeff_temp <- my_model %>% filter(codes == "sex") %>% dplyr::select(coeff)
    mm_covars <- mm_covars %>% mutate(sum_temp = sum_temp + (as.numeric(sex) * unlist(coeff_temp)))

    coeff_temp <- my_model %>% filter(codes == "constant") %>% dplyr::select(coeff)
    mm_covars <- mm_covars %>% mutate(sum_temp = sum_temp + unlist(coeff_temp))
    
    #fitting by age and sex to get residuals
    my_fit <- lm(sum_temp ~ factor(as.numeric(sex)) + age, data = mm_covars)
    my_coeffs <- summary(my_fit)$coefficients 
    mm_covars <- mm_covars %>% mutate(age_sex_fit = my_coeffs[3, 1]*age + my_coeffs[2,1]*as.numeric(sex) + my_coeffs[1,1])
    mm_covars <- mm_covars %>% mutate(residual_fit = sum_temp - age_sex_fit)

    #rename column
    mm_covars <- mm_covars %>% mutate(prognostic_index = sum_temp)
    return(mm_covars)
}

mm_score_take_res <- function(df_inc_mm_var, mm_score_col_name) {
#requires df_inc_mm_var to include columns named "age" and "sex" as well as mm_score column
    df_res_temp <- as.data.frame(df_inc_mm_var)[c(mm_score_col_name, "v0_sex_genetic", "age")]
    my_fit <- lm(df_res_temp)
    my_coeffs <- summary(my_fit)$coefficients 
    df_inc_mm_var <- df_inc_mm_var %>% mutate(age_sex_fit = my_coeffs[3, 1]*age + my_coeffs[2,1]*as.numeric(v0_sex_genetic) + my_coeffs[1,1])
    df_inc_mm_var <- df_inc_mm_var %>% mutate(residual_fit = !!as.name(mm_score_col_name) - age_sex_fit)
    return(df_inc_mm_var)
}
