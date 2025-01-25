##File name: colorectal_biomarker_covariates.r
##Author: Hannah Harrison
##Last Edit: 28/03/2022
##Description: defines biomakers (from gp data) used for colorectal cancer models
#libraries if needed here


##simple functions to identify abnormal measurements of biomarkers separately

##below threshold - simple
flag_abnormal_albumin <- function(df, albumin_thresh){
    return(df %>% 
        mutate(albumin_flag = case_when((biomarker == "Albumin" & !is.na(value) & value < albumin_thresh) ~ 1, TRUE ~ 0)))
} 

flag_abnormal_Ferritin_low <- function(df, Ferritin_thresh_low){
    return(df %>% 
        mutate(Ferritin_low_flag = case_when((biomarker == "Ferritin" & !is.na(value) & value < Ferritin_thresh_low) ~ 1, TRUE ~ 0)))
} 

flag_abnormal_MCV <- function(df, MCV_thresh){
    return(df %>% 
        mutate(MCV_flag = case_when((biomarker == "MCV" & !is.na(value) &  value < MCV_thresh) ~ 1, TRUE ~ 0)))
}

##below threshold - with sex condition
flag_abnormal_Ferritin_high <- function(df, Ferritin_thresh_high){
    return(df %>% 
        mutate(Ferritin_high_flag = case_when((biomarker == "Ferritin" & !is.na(value) & 
                                    ((value > Ferritin_thresh_high[1] & sex == 0 ) | (value > Ferritin_thresh_high[2] & sex == 1))) ~ 1, TRUE ~ 0)))
} 

flag_abnormal_Hgb <- function(df, Hgb_thresh){
    return(df %>% 
        mutate(Hgb_flag = case_when((biomarker == "Haemoglobinconc" & !is.na(value) & 
                                    ((value < Hb_thresh[1] & sex == 0 ) | (value < Hb_thresh[2] & sex == 1))) ~ 1, TRUE ~ 0)))
} 

flag_abnormal_Hmc <- function(df, Hmc_thresh){
    return(df %>% 
        mutate(Hmc_flag = case_when((biomarker == "Haematocritperc" & !is.na(value) &
                                 ((value < Hmc_thresh[1] & sex == 0 ) | (value < Hmc_thresh[2] & sex == 1))) ~ 1, TRUE ~ 0)))
} 

##above threshold - simple
flag_abnormal_CRP <- function(df, CRP_thresh){
    return(df %>% 
        mutate(CRP_flag = case_when((biomarker == "CRP" & !is.na(value) & value > CRP_thresh ) ~ 1, TRUE ~ 0)))
} 

flag_abnormal_PV <- function(df, PV_thresh){
    return(df %>% 
        mutate(PV_flag = case_when((biomarker == "PV" & !is.na(value) &  value > PV_thresh) ~ 1, TRUE ~ 0)))
}

flag_abnormal_Plate <- function(df, Plate_thresh){
    return(df %>% 
        mutate(Plate_flag = case_when((biomarker == "Platelets" & !is.na(value) &  value > Plate_thresh) ~ 1, TRUE ~ 0)))
}


flag_abnormal_CA125 <- function(df, CA125_thresh){
    return(df %>% 
        mutate(CA125_flag = case_when((biomarker == "CA125" & !is.na(value) &  value > CA125_thresh) ~ 1, TRUE ~ 0)))
}

##above threshold - with age and sex condition
flag_abnormal_ESR <- function(df, df_ESR){
    df <- df %>% mutate(ESR_age_groups = cut(event_age,
                                               breaks = c(0, 40, 50, 60, 70, 80, 110),
                                               labels = c(1, 2, 3, 4, 5, 6))) #age groups for merge
    ####merge in ESR thresholds by age group and sex
    df <- merge(df, df_ESR, by = c("sex", "ESR_age_groups"))

    return(df %>% mutate(ESR_flag = case_when((biomarker == "ESR" & !is.na(value) & value > ESR_thresh) ~ 1, TRUE ~ 0)))
}


###creates model variables (for df with one row per eid, runs after have summarise test meta-data)
####IRON DEFICIENCY 
iron_def_vars <- function(df){
    df<- df %>% mutate(iron_def_direct = case_when((Hgb_tests > 0 & Ferritin_tests > 0) & #both Hbg and Ferritin have been measured
                                                    (Hgb_abnormal > 0 & Ferritin_low_abnormal > 0) ~ 1, TRUE ~ 0)) # and both have at least one abnormally low measurment 

    df<- df %>% mutate(iron_def_proxy1 = case_when((Ferritin_tests == 0 & Hgb_tests > 0 & MCV_tests > 0) & #Ferritin hasn't been measured
                                                    (Hgb_abnormal > 0 & MCV_abnormal > 0) ~ 1, TRUE ~ 0)) # but hbg is abnormally low

    df<- df %>% mutate(iron_def_proxy2 = case_when((Ferritin_tests == 0 & MCV_tests == 0 & Hgb_tests > 0) & #Ferritin hasn't been measured
                                                    (Hgb_abnormal > 0) ~ 1, TRUE ~ 0)) # but hbg is abnormally low

    #df<- df %>% mutate(iron_def_proxy3 = case_when((Ferritin_tests == 0 & Hgb_tests == 0 & Hmc_tests > 0) & #ferritin and hgb haven't been measured
                                       #(Hmc_abnormal > 0) ~ 1, TRUE ~ 0)) # but haemocrit is abnormally low (have removed haemocrit)
                                   

    df <- df %>% mutate(iron_def = case_when((iron_def_direct == 1 | iron_def_proxy1 == 1 | iron_def_proxy2 == 2) ~ 1, TRUE ~0))
    df <- df %>% mutate(iron_def_measured = case_when(Hgb_tests > 0 ~ 1, TRUE ~0))

    return(df)
}

####INFLAMMATION
inflam_vars <- function(df){
    df <- df %>% mutate(inflam_direct = case_when((CRP_tests > 0 | ESR_tests > 0 | PV_tests > 0) & # one of the inflammation tests has been measured
                                                    (CRP_abnormal > 0 | ESR_abnormal > 0 | PV_abnormal > 0) ~ 1, TRUE ~ 0)) # and an abnormal result measured

     df <- df %>% mutate(inflam_proxy1 = case_when((CRP_tests == 0 & ESR_tests == 0 & PV_tests == 0 & Plate_tests > 0) & #no inflammation tests
                                                    (Plate_abnormal > 0) ~ 1, TRUE ~ 0))                                #but platelets are high    

    df <- df %>%  mutate(inflam_proxy2 = case_when((CRP_tests == 0 & ESR_tests == 0 & PV_tests == 0 & albumin_tests > 0) & #no inflammation tests
                                                    (albumin_abnormal > 0) ~ 1, TRUE ~ 0))                                #but albumin is low    

    df <- df %>%  mutate(inflam_proxy3 = case_when((CRP_tests == 0 & ESR_tests == 0 & PV_tests == 0 & Ferritin_tests > 0) & #no inflammation tests, but ferritin high                                                        #but Ferritin is high
                                                    (Ferritin_high_abnormal > 0) ~ 1, TRUE ~ 0))         

    df <- df %>% mutate(inflam_all = case_when((inflam_direct == 1 | inflam_proxy1 == 1 | inflam_proxy2 == 1 | inflam_proxy3 == 1) ~ 1, TRUE ~0))
    df <- df %>% mutate(inflam_measured = case_when((CRP_tests > 0 | ESR_tests > 0 | PV_tests > 0 | albumin_tests > 0 | Ferritin_tests > 0) ~ 1, TRUE ~0))

    return(df) 
}

CA125_raised_vars <- function(df){
    df <- df %>% mutate(CA125_raised_men = case_when((CA125_tests > 0 & CA125_abnormal > 0 & sex == 1) ~ 1, TRUE ~ 0)) # data check for CA125 measurements in men
    df <- df %>% mutate(CA125_raised_women = case_when((CA125_tests > 0 & CA125_abnormal > 0 & sex == 0) ~ 1, TRUE ~ 0)) # expected for some women
    df <- df %>% mutate(CA125_raised_all = case_when((CA125_raised_men == 1 | CA125_raised_women == 1) ~ 1, TRUE ~ 0))
    return(df)
}
##archived old functions for combined biomarker results
##Inflammatory Markers Function
archive_inflam_var <- function(my_df, CRP_thresh, PV_thresh, df_ESR_thresh) {
    ##checks if either CRP, ESR or PV above threshold
    ####formats main dataframe for merge
    my_df <- my_df %>% mutate(ESR_age_groups = cut(event_age,
                                               breaks = c(0, 40, 50, 60, 70, 80, 110),
                                               labels = c(1, 2, 3, 4, 5, 6))) #age groups for merge
    ####merge in ESR thresholds by age group and sex
    my_df <- merge(my_df, df_ESR_thresh, by = c("sex", "ESR_age_groups"))

    ####defines new variable for "any inflammatory marker is high"
    my_df <- my_df %>% 
                mutate(gp_inflam = case_when(
                    ((biomarker == "CRP" & !is.na(value) & value > CRP_thresh)  |
                    (biomarker == "ESR" & !is.na(value) & value > ESR_thresh) |
                    (biomarker == "PV" & !is.na(value) &  value > PV_thresh)) ~ 1, TRUE ~ 0))
    ####tidy up and return
    return(my_df)

}

archive_anaemia_var <- function(my_df, Hb_thresh, Hmc_thresh, Ferritin_thresh) {
    ##Iron deficiency/anemia
    ####based on haemoglobin, haematocrit or ferritin thresholds
    return(my_df %>% 
        mutate(gp_anaemia = case_when(
                                ((biomarker == "Haemoglobinconc" & !is.na(value) &
                                ((value < Hb_thresh[1] & sex == 0 ) | (value < Hb_thresh[2] & sex == 1))) |
                                (biomarker == "Haemocritperc" & !is.na(value) &
                                ((value < Hmc_thresh[1] & sex == 0 ) | (value < Hmc_thresh[2] & sex == 1))) |
                                (biomarker == "Ferritin" & !is.na(value) &
                                (value < Ferritin_thresh))) ~ 1, TRUE ~ 0 )))
}

archive_thrombo_var <- function(my_df, Plate_thresh) {
    ##Thrombocytosis
    ####based on high platlets
    return(my_df %>% 
        mutate(gp_thrombo = case_when(biomarker == "Platelets" & !is.na(value) &
                                    value > Plate_thresh  ~ 1, TRUE ~ 0)))
}

