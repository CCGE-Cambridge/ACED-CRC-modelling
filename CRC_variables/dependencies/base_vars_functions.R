##File name: base_vars_functions.r
##Author: Hannah Harrison
##Last Edit: 04/03/2022
##Description: functions to define baseline variables

##Smoking
smoking_stat_clean <- function(df_smoke){
    df_smoke <- df_smoke %>%
        mutate(smoking_status= na_if(smoking_status, -3))
    return(df_smoke)
}
###Education 
clean_edu_quals <- function(df_edu) {
    df_edu <- df_edu %>% 
        mutate(edu_quals= na_if(edu_quals, -3)) %>% #missing if prefer not to answer
        mutate(edu_quals = as.numeric(edu_quals)) %>%
        mutate(edu_quals = recode(edu_quals, `-7` = 7)) # category 7 for none of the above
    return(df_edu)
}

gen_edu_length <- function(df_edu){
    #handle missing values for education end age 
    df_edu <- df_edu %>% 
        mutate(edu_end_age = na_if(edu_end_age, -1)) %>%
        mutate(edu_end_age= na_if(edu_end_age, -3)) %>%
        mutate(edu_end_age = as.numeric(edu_end_age)) %>%
        mutate(edu_end_age = recode(edu_end_age, `-2` = 5 )) #recode never went to school to left school at 5
    
    #calculate education length based on education leaving age
    df_edu <- df_edu %>% 
        mutate(edu_length = edu_end_age - 5) #assume everyone starts school at 5
    
    ##clean edu quals as needed
    df_edu <- clean_edu_quals(df_edu)

    #fill in missing length data using education type
    df_edu <- df_edu %>%
        mutate(edu_length_from_quals = recode(edu_quals , `7` = 10, `1` = 16, `2` = 13, `3` = 11, `4` = 11, `5` = 14, `6` = 16 )) %>%
        mutate(edu_length = case_when(is.na(edu_length) ~ edu_length_from_quals, TRUE ~ edu_length)) 
    
    return(df_edu)
}

###Dietary
fibre_from_fv <- function(df_diet) {
    #information needed
    fv_fields <- c("cook_veg", "raw_veg", "fresh_fruit", "dried_fruit")
    fv_fibre <- c(1, 1, 2, 0.5) #fibre content per portion in grams
    df_fv_fibre <- data.frame(fv_fields, fv_fibre)

    ##handle missing data
    df_diet <- df_diet %>% 
        mutate(across(all_of(df_fv_fibre$fv_fields), ~na_if(. ,-1))) %>%
        mutate(across(all_of(df_fv_fibre$fv_fields), ~na_if(. ,-3))) %>%
        mutate(across(all_of(df_fv_fibre$fv_fields), as.numeric)) %>%
        mutate(across(all_of(df_fv_fibre$fv_fields), ~recode(. , `-10` = 0.5)))   

    ##sum fv fields weighted by fibre content        
     df_diet <- df_diet %>% 
        mutate(fruit_and_veg_fibre_daily =  wtd_and_na_handle.rowSums(select(., df_fv_fibre$fv_fields), wts = df_fv_fibre$fv_fibre)) 

    return(df_diet)
}

fibre_from_bread <- function(df_diet){
    #information needed
    bread_cats <- c("white", "brown", "wholemeal", "other") # ALERT: MISSING IS CODED TO 1.25 FOR BREAD
    bread_fibre <- c(0.68, 1.26, 1.80, 1.25)
    df_bread_fibre <- data.frame(bread_cats,  bread_fibre)

    #handle missing data
    df_diet <- df_diet %>% 
        mutate(bread_intake = na_if(bread_intake, -1)) %>%
        mutate(bread_intake = na_if(bread_intake, -1)) %>%
        mutate(bread_type = as.numeric(bread_type)) %>%
        mutate(bread_type = recode(bread_type, `-3` = 4, `-1` = 4)) 
    
    #extract fibre content (including conversion from weekly to daily)
    df_diet <- df_diet %>% 
        mutate(bread_fibre_per_portion = recode(bread_type, !!!df_bread_fibre$bread_fibre)) %>%
        mutate(bread_fibre_daily = (bread_intake*bread_fibre_per_portion)/7)
    
    return(df_diet)
}

fibre_from_cereal <- function(df_diet){
    ##information on cereal variables
    cereal_cats <- c("bran", "biscuit", "oat", "muesli", "other", "unknown") # ALERT: MISSING IS CODED TO 3.34 FOR CEREAL
    cereal_fibre <- c(7.16, 2.92, 1.92, 4.18, 0.54, 3.34)
    df_cereal_fibre <- data.frame(cereal_cats, cereal_fibre)

    ##handle missing data
    df_diet <- df_diet %>% 
        mutate(cereal_intake = na_if(cereal_intake, -1)) %>%
        mutate(cereal_intake = na_if(cereal_intake, -3)) %>%
        mutate(cereal_intake = as.numeric(cereal_intake)) %>%
        mutate(cereal_intake= recode(cereal_intake, `-10` = 0.25)) %>%  #<1/wk is coded as 0.25/wk
        mutate(cereal_type = as.numeric(cereal_type)) %>%
        mutate(cereal_type = recode(cereal_type, `-3` = 6, `-1` = 6))  #missing recoded to unknown
        
    ##extract fibre content (and convert from weekly to daily)
    df_diet <- df_diet %>%     
        mutate(cereal_fibre_per_portion = recode(cereal_type, !!!df_cereal_fibre$cereal_fibre)) %>%
        mutate(cereal_fibre_daily = (cereal_intake*cereal_fibre_per_portion)/7)

    return(df_diet)
}
gen_partial_fibre_score <- function(df_diet){
#calculates fibre from fruit/veg, bread and cereal, then combines into fibre score.
    df_diet <- fibre_from_fv(df_diet)
    df_diet <- fibre_from_bread(df_diet)
    df_diet <- fibre_from_cereal(df_diet)
    fibre_fields <- c("fruit_and_veg_fibre_daily", "bread_fibre_daily", "cereal_fibre_daily")
    df_diet <- df_diet %>%  
        mutate(partial_fibre_score = na_handle.rowSums(select(., all_of(fibre_fields))))
    return(df_diet)
}

gen_tot_red_meat <- function(df_diet){
    #calculates weekly red meat consumption
    ##fields of interest
    red_meat_fields <- c("beef_meat", "lamb_meat", "pork_meat") # all use data coding 100377

    ##handle missing data
    df_diet <- df_diet %>%  
        mutate(across(all_of(red_meat_fields), ~na_if(. ,-1))) %>%
        mutate(across(all_of(red_meat_fields), ~na_if(. ,-3))) 
    
    ##recode frequency categories to number of portions
    df_diet <- df_diet %>%  
        mutate(across(all_of(red_meat_fields), as.numeric)) %>%
        mutate(across(all_of(red_meat_fields), ~recode(. , `1` = 0.5, `2` = 1, `4` = 5.5, `5` = 7), .names = "{col}_portions")) 

    #sum all red meat portion fields
    red_meat_portions <-  paste(red_meat_fields, "portions", sep = "_")
    df_diet <- df_diet %>%  
        mutate(tot_red_meat_wk = na_handle.rowSums(select(., all_of(red_meat_portions))))

    return(df_diet)
}
gen_tot_proc_meat <- function(df_diet){
    ##missing data
    df_diet <- df_diet %>%  
        mutate(process_meat = na_if(process_meat, -1)) %>%
        mutate(process_meat = na_if(process_meat, -3)) 

    ##recode frequency categories to portions     
    df_diet <- df_diet %>%  
        mutate(process_meat = as.numeric(process_meat)) %>%
        mutate(tot_process_meat_wk = recode(process_meat , `1` = 0.5, `2` = 1, `4` = 5.5, `5` = 7))
                           
    return(df_diet)
}

gen_tot_all_meat <- function(df_diet){
    #calculates weekly meat consumption
    ##fields of interest
    all_meat_fields <- c("beef_meat", "lamb_meat", "pork_meat", "process_meat", "poultry_meat")

    ##handle missing data
    df_diet <- df_diet %>%  
        mutate(across(all_of(all_meat_fields), ~na_if(. ,-1))) %>%
        mutate(across(all_of(all_meat_fields), ~na_if(. ,-3))) 
    
    ##recode frequency categories to number of portions
    df_diet <- df_diet %>%      
        mutate(across(all_of(all_meat_fields), as.numeric)) %>%
        mutate(across(all_of(all_meat_fields), ~recode(. , `0` = 0, `1` = 0.5, `2` = 1,  `4` = 5.5, `5` = 7), .names = "{col}_portions")) 
    
    ##sum across all meat fields
    all_meat_portions <- month_fields = paste(all_meat_fields, "portions", sep = "_")
    df_diet <- df_diet %>%  
        mutate(tot_all_meat_wk = na_handle.rowSums(select(., all_of(all_meat_portions))))

    return(df_diet)
}

ca_from_cereal <- function(df_diet){
##calcium from milk in cereal (bradbury et al. 2020)
    calcium_per_cereal <- 120 #in mg assuming 100ml of milk per bowl

    ##handle missing data
    df_diet <- df_diet %>% 
        mutate(cereal_intake = na_if(cereal_intake, -1)) %>%
        mutate(cereal_intake = na_if(cereal_intake, -3)) 

    ##recode special values to bowls of cereal               
    df_diet <- df_diet %>% 
    mutate(cereal_intake = as.numeric(cereal_intake)) %>%
    mutate(cereal_intake= recode(cereal_intake, `-10` = 0.25))  ##<1 to 0.25
    
    ##estimate of calcium from milk in cereal
    df_diet <- df_diet %>% 
        mutate(calcium_from_cereal_milk_wk = cereal_intake*calcium_per_cereal)

    return(df_diet)
}

ca_from_tea <- function(df_diet){
##estimate of calcium intake from milk in tea (bradbury et al. 2020)
    calcium_per_tea <- 42 #in mg, assuming 35ml of milk per tea

    #handle missing data
    df_diet <- df_diet %>% 
        mutate(tea_intake = na_if(tea_intake, -1)) %>%
        mutate(tea_intake = na_if(tea_intake, -3)) 
    
    #recode special values to cups of tea
    df_diet <- df_diet %>%
        mutate(tea_intake = as.numeric(tea_intake)) %>%
        mutate(tea_intake= recode(tea_intake, `-10` = 0.25))  ##<1 to 0.25
    
    #calculate calcium per cup of tea
    df_diet <- df_diet %>% 
        mutate(calcium_from_tea_milk_wk = (tea_intake*calcium_per_tea)*7)

    return(df_diet)
}

ca_from_cheese <- function(df_diet){
#estimate of calcium intake from cheese (bradbury et al. 2020)
    calcium_per_cheese <- 150 #739mg per 100g, assume that a serving is approx 20g

    #handle missing data
    df_diet <- df_diet %>% 
        mutate(cheese_intake = na_if(cheese_intake, -1)) %>%
        mutate(cheese_intake = na_if(cheese_intake, -3))
    
    #recode frequency categories to portions consumed
    df_diet <- df_diet %>% 
        mutate(cheese_intake = as.numeric(cheese_intake)) %>%
        mutate(cheese_intake_portion= recode(cheese_intake, `1` = 0.25, `2` = 1, `4` = 5.5, `5` = 7)) ##<1 to 0.25
    
    #calculate calcium from cheese
    df_diet <- df_diet %>% 
        mutate(calcium_from_cheese_wk = cheese_intake_portion*calcium_per_cheese)

    return(df_diet)
}

ca_from_supp <- function(df_diet){
    #calculate calcium from supplement, if taking estimated daily dose of 250mg
    calcium_in_supp <- 250
    df_diet <- df_diet %>% 
        mutate(calcium_from_supp_wk = case_when(min_other_supps == 3 ~ calcium_in_supp*7, TRUE ~ 0 ))

    return(df_diet)
}

gen_ca_score <-function(df_diet) {
    df_diet <- ca_from_cereal(df_diet)
    df_diet <- ca_from_tea(df_diet)
    df_diet <- ca_from_cheese(df_diet)
    df_diet <- ca_from_supp(df_diet)
    weekly_diet_calcium_fields <- c("calcium_from_cereal_milk_wk", "calcium_from_tea_milk_wk", "calcium_from_cheese_wk")
    df_diet <- df_diet %>% 
        mutate(tot_calcium_wk = na_handle.rowSums(select(., all_of(weekly_diet_calcium_fields)))) %>%
        mutate(tot_calcium_wk = tot_calcium_wk +  calcium_from_supp_wk)

    return(df_diet)
}


###Alcohol
alc_clean <- function(df_alc, df_fields){
    df_alc <- df_alc %>% 
        mutate(across(all_of(c(df_fields$month_fields, df_fields$week_fields)), ~na_if(. ,-1))) %>%
        mutate(across(all_of(c(df_fields$month_fields, df_fields$week_fields)), ~na_if(. ,-3))) %>%        
    return(df_alc)
}

alc_info <- function() {
    alcohol_types <- c("rwine", "wwine", "fwine", "beer", "spirits", "other_alcohol")
    units <- c(2.1, 2.1, 2.4, 2.5, 1.0, 1.5) # units per portion 
    df_alcohol <- data.frame(alcohol_types, units) %>%
        mutate(month_fields = paste(alcohol_types, "mnth", sep = "_"), 
               week_fields = paste(alcohol_types, "wk", sep = "_"))
    rownames(df_alcohol) <- df_alcohol$alcohol_types
    return(df_alcohol)
}

alc_daily_units <- function(df_alc) {
#n.b. ukb collects data differently for alcohol consumption depending on how an earlier question was answers
#so, units weekly and units monthly only have values for half of the cohort each

    df_fields <- alc_info()
    df_alc <- alc_clean(df_alc, df_fields)
    df_alc <- df_alc %>%
        mutate(units_monthly = wtd_and_na_handle.rowSums(select(., df_fields$month_fields), wts = df_fields$units)) %>%
        mutate(units_weekly = wtd_and_na_handle.rowSums(select(., df_fields$week_fields), wts = df_fields$units)) %>% 
        mutate(alcohol_units_daily = wtd_and_na_handle.rowSums(select(., units_monthly, units_weekly), wts = c(1/30.44, 1/7))) %>%
        mutate(alcohol_units_daily = case_when(alcohol_freq == 6 ~ 0,
                                       (is.na(alcohol_units_daily) & alcohol_freq == 5) ~ 0,
                                       (is.na(alcohol_units_daily) & alcohol_freq == 4) ~ 0.15, 
                                        TRUE ~ alcohol_units_daily)) 
    return(df_alc)
}

alc_exceed <- function(df_alc) {
    #flags if consuming more than 14 units daily - expects that you have already run alc_daily_units()
    df_alc <- df_alc %>%
        mutate(exceed_alcohol_rec = case_when((alcohol_units_daily >=  0 &  alcohol_units_daily*7 <= 14) ~ 0, alcohol_units_daily*7 > 14 ~ 1,  TRUE ~ NA_real_))
    return(df_alc)
}

##Physical Activity
pa_clean <- function(df_pa, num_days, length_in_mins) {
#handle missing values, recode special values and truncate excessively high length of exercise
    df_pa <- df_pa %>% 
        mutate(across(all_of(c(num_days, length_in_mins)), ~na_if(. ,-1))) %>%
        mutate(across(all_of(c(num_days, length_in_mins)), ~na_if(. ,-3))) %>%
        mutate(across(all_of(c(num_days, length_in_mins)), as.numeric)) %>%
        mutate(cannot_walk = case_when(walk_pa_freq == -2 ~ 1, TRUE ~ 0)) %>%
        mutate(across(all_of(num_days), ~recode(. , `-2` = 0))) %>% #recode unable to walk as 0 days
        mutate(across(all_of(length_in_mins), ~ case_when(. > 180 ~ 180, . <= 180 ~ .))) ##truncate exercise lengths of more than 3 hours in a day
    
    return(df_pa)
}

pa_mets <- function(df_pa){
#caluclate extimate of mets/week, expects that yuou will have already done pa_clean
    walk_mets_per_min <- 3.3
    mod_mets_per_min <- 4.0 
    vig_mets_per_min <- 8.0

    df_pa <- df_pa %>% 
        mutate(walk_week_mets = walk_mets_per_min*walk_pa_freq*walk_pa_length,
                mod_week_mets = mod_mets_per_min*mod_pa_freq*mod_pa_length, 
                vig_week_mets = vig_mets_per_min*vig_pa_freq*walk_pa_length) %>%
        mutate(tot_METs_wk = na_handle.rowSums(select(. , all_of(c("walk_week_mets", "mod_week_mets" , "vig_week_mets"))))) %>%
        mutate(tot_METs_wk = case_when((is.na(tot_METs_wk) & cannot_walk == 1) ~ 0, TRUE ~ tot_METs_wk))
    
    return(df_pa)
}

pa_IPAQ <- function(df_pa, num_days_pa_vars){
###IPAQ score for activity (1 = inactive, 2 = moderate, 3 = active), expects that you have already done pa_clean() and pa_mets()
    df_pa <- df_pa %>% 
        mutate(active_days = rowSums(select(. , all_of(num_days_pa_vars)), na.rm = TRUE)) %>%
        mutate(IPAQ_score = case_when(is.na(tot_METs_wk) ~ NA_real_ ,
            ((active_days >= 7 & tot_METs_wk >= 3000) | (vig_pa_freq >= 7 & tot_METs_wk >= 1500)) ~ 3, 
            ((vig_pa_freq >= 3 & vig_pa_length >= 20) | (mod_pa_freq >= 5 & mod_pa_length >= 30) | (active_days >= 5 & tot_METs_wk >= 600)) ~ 2 ,
        TRUE ~ 1))
    return(df_pa)
}

pa_do_all <- function(df_pa){
    num_days_pa_vars <- c("walk_pa_freq", "mod_pa_freq", "vig_pa_freq")
    length_in_mins_pa_vars <- c("walk_pa_length", "mod_pa_length", "vig_pa_length")
    
    df_pa <- pa_clean(df_pa, num_days_pa_vars, length_in_mins_pa_vars)
    df_pa <- pa_mets(df_pa)
    df_pa <- pa_IPAQ(df_pa, num_days_pa_vars)
    return(df_pa)
}