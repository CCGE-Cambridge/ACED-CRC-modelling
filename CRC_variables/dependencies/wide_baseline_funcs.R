##File name: wide_baseline_funs.r
##Author: Hannah Harrison
##Last Edit: 27/07/2022
##Description:functions to manage "wide" data in UKB baseline

#prep codelists for search

#function to check if there is any record of any event
UKB_wide_event_any <- function(df, patt_event_var, num_vars, type_name) {

    #create new column to hold result
    new_var_name <- paste0(type_name, "_any_ever")
    df[[new_var_name]] <- 0


    #iterate through event variables
    for (i in 0:num_vars) {
        curr_event_var <- paste0(patt_event_var, i)
        df <- df %>% mutate(!!new_var_name := 
                        case_when((!is.na(!!as.name(curr_event_var)) | !!as.name(new_var_name) == 1) ~ 1, TRUE ~ 0))
    }
    return(df)
}


UKB_wide_event_any_except <- function(df, patt_event_var, num_vars, patt_ignore, type_name) {
 #for example, ignore nm skin cancer

    #create new column to hold result
    new_var_name <- paste0(type_name, "_any_except")
    df[[new_var_name]] <- 0

    #iterate through event variables
    for (i in 0:num_vars) {
        curr_event_var <- paste0(patt_event_var, i)
        df <- df %>% mutate(!!new_var_name := 
                        case_when(((!is.na(!!as.name(curr_event_var)) & 
                                    !str_detect(!!as.name(curr_event_var), patt_ignore)) |
                                        !!as.name(new_var_name) == 1) ~ 1, TRUE ~ 0))
    }
    return(df)
}

UKB_wide_event_any_except_after <- function(df, patt_event_var, patt_date_var, num_vars, patt_ignore, type_name, col_min_date, study_end_date) {
 #for example, ignore nm skin cancer

    #create new column to hold result
    new_var_name <- paste0(type_name, "_any_except_after")
    df[[new_var_name]] <- 0

    #iterate through event variables
    for (i in 0:num_vars) {
        curr_event_var <- paste0(patt_event_var, i)
        curr_date_var <- paste0(patt_date_var, i)

        df <- df %>% mutate(!!new_var_name := 
                        case_when(((!is.na(!!as.name(curr_event_var)) & 
                                 (as.Date(!!as.name(curr_date_var)) > !!as.name(col_min_date) & 
                                 as.Date(!!as.name(curr_date_var)) <= study_end_date) &   
                                    !str_detect(!!as.name(curr_event_var), patt_ignore)) |
                                        !!as.name(new_var_name) == 1) ~ 1, TRUE ~ 0))
    }
    return(df)
}

UKB_wide_event_any_except_in_range <- function(df, patt_event_var, patt_date_var, num_vars, patt_ignore, type_name, col_min_date, years_follow_up) {
 #for example, ignore nm skin cancer

    #create new column to hold result
    new_var_name <- paste0(type_name, "_any_except_", years_follow_up, "yr")
    df[[new_var_name]] <- 0

    #iterate through event variables
    for (i in 0:num_vars) {
        curr_event_var <- paste0(patt_event_var, i)
        curr_date_var <- paste0(patt_date_var, i)

        df <- df %>% mutate(!!new_var_name := 
                        case_when(((!is.na(!!as.name(curr_event_var)) & 
                                 (as.Date(!!as.name(curr_date_var)) > !!as.name(col_min_date) & 
                                 as.Date(!!as.name(curr_date_var)) <= !!as.name(col_min_date) + years(years_follow_up)) &   
                                    !str_detect(!!as.name(curr_event_var), patt_ignore)) |
                                        !!as.name(new_var_name) == 1) ~ 1, TRUE ~ 0))
    }
    return(df)
}


#UKB_wide_event_fetch_first #get event code and date for first event type

#UKB_wide_event_fetch_first_except #for example, ignore nm skin cancer
UKB_wide_event_fetch_first_except <- function(df, patt_event_var, patt_date_var, num_vars, patt_ignore, type_name, starter = 0) {
 #for example, ignore nm skin cancer

    #create new column to hold result
    new_var_event <- paste0(type_name, "_code_first")
    new_var_date <- paste0(type_name, "_date_first")
    df[[new_var_event]] <- NA_character_
    df[[new_var_date]] <- NA_Date_

    #iterate through event variables 
    for (i in starter:num_vars) {
        curr_event_var <- paste0(patt_event_var, i)
        curr_date_var <- paste0(patt_date_var, i)

        df <- df %>% mutate(last_event = !!as.name(new_var_event)) 
        
        df <- df %>% mutate(!!new_var_event := 
                        case_when(((!is.na(!!as.name(curr_event_var)) & 
                                    !str_detect(!!as.name(curr_event_var), patt_ignore)) &
                                    (is.na(!!as.name(new_var_event)) | 
                                    as.Date(!!as.name(curr_date_var)) < !!as.name(new_var_date))) ~ as.character(!!as.name(curr_event_var)), 
                                    TRUE ~ !!as.name(new_var_event)))
        
       
        df <- df %>% mutate(!!new_var_date := 
                        case_when(((!is.na(!!as.name(curr_event_var)) & 
                                    !str_detect(!!as.name(curr_event_var), patt_ignore)) &
                                    (is.na(last_event) | 
                                    as.Date(!!as.name(curr_date_var)) < !!as.name(new_var_date))) ~ as.Date(!!as.name(curr_date_var)),
                                    TRUE ~ !!as.name(new_var_date)))

    }

    df <- df %>% select(!last_event)
    return(df)
}

UKB_wide_event_fetch_first_except_after <- function(df, patt_event_var, patt_date_var, num_vars, patt_ignore, type_name, col_min_date, study_end_date, starter = 0) {
 #for example, ignore nm skin cancer

    #create new column to hold result
    new_var_event <- paste0(type_name, "_code_first_after")
    new_var_date <- paste0(type_name, "_date_first_after")

    df[[new_var_event]] <- NA_character_
    df[[new_var_date]] <- NA_Date_

    #iterate through event variables 
    for (i in starter:num_vars) {
        curr_event_var <- paste0(patt_event_var, i)
        curr_date_var <- paste0(patt_date_var, i)

        df <- df %>% mutate(last_event = !!as.name(new_var_event)) 
        
        df <- df %>% mutate(!!new_var_event := 
                        case_when(((!is.na(!!as.name(curr_event_var)) & 
                                (as.Date(!!as.name(curr_date_var)) > !!as.name(col_min_date) & 
                                 as.Date(!!as.name(curr_date_var)) <= study_end_date) &
                                    !str_detect(!!as.name(curr_event_var), patt_ignore)) &
                                    (is.na(!!as.name(new_var_event)) | 
                                    as.Date(!!as.name(curr_date_var)) < !!as.name(new_var_date))) ~ as.character(!!as.name(curr_event_var)), 
                                    TRUE ~ !!as.name(new_var_event)))
        
       
        df <- df %>% mutate(!!new_var_date := 
                        case_when(((!is.na(!!as.name(curr_event_var)) & 
                                    (as.Date(!!as.name(curr_date_var)) > !!as.name(col_min_date) & 
                                 as.Date(!!as.name(curr_date_var)) <= study_end_date) &    
                                    !str_detect(!!as.name(curr_event_var), patt_ignore)) &
                                    (is.na(last_event) | 
                                    as.Date(!!as.name(curr_date_var)) < !!as.name(new_var_date))) ~ as.Date(!!as.name(curr_date_var)),
                                    TRUE ~ !!as.name(new_var_date)))

    }

    df <- df %>% select(!last_event)
    return(df)
}

UKB_wide_event_fetch_first_except_in_range <- function(df, patt_event_var, patt_date_var, num_vars, patt_ignore, type_name, col_min_date, years_follow_up, starter = 0) {
 #for example, ignore nm skin cancer

    #create new column to hold result
    new_var_event <- paste0(type_name, "_code_first_", years_follow_up, "yr")
    new_var_date <- paste0(type_name, "_date_first_", years_follow_up, "yr")

    df[[new_var_event]] <- NA_character_
    df[[new_var_date]] <- NA_Date_

    #iterate through event variables 
    for (i in starter:num_vars) {
        curr_event_var <- paste0(patt_event_var, i)
        curr_date_var <- paste0(patt_date_var, i)

        df <- df %>% mutate(last_event = !!as.name(new_var_event)) 
        
        df <- df %>% mutate(!!new_var_event := 
                        case_when(((!is.na(!!as.name(curr_event_var)) & 
                                (as.Date(!!as.name(curr_date_var)) > !!as.name(col_min_date) & 
                                 as.Date(!!as.name(curr_date_var)) <= !!as.name(col_min_date) + years(years_follow_up)) &
                                    !str_detect(!!as.name(curr_event_var), patt_ignore)) &
                                    (is.na(!!as.name(new_var_event)) | 
                                    as.Date(!!as.name(curr_date_var)) < !!as.name(new_var_date))) ~ as.character(!!as.name(curr_event_var)), 
                                    TRUE ~ !!as.name(new_var_event)))
        
       
        df <- df %>% mutate(!!new_var_date := 
                        case_when(((!is.na(!!as.name(curr_event_var)) & 
                                    (as.Date(!!as.name(curr_date_var)) > !!as.name(col_min_date) & 
                                 as.Date(!!as.name(curr_date_var)) <= !!as.name(col_min_date) + years(years_follow_up)) &    
                                    !str_detect(!!as.name(curr_event_var), patt_ignore)) &
                                    (is.na(last_event) | 
                                    as.Date(!!as.name(curr_date_var)) < !!as.name(new_var_date))) ~ as.Date(!!as.name(curr_date_var)),
                                    TRUE ~ !!as.name(new_var_date)))

    }

    df <- df %>% select(!last_event)
    return(df)
}



#function to check if there is any record of the specified event
UKB_wide_event_ever <- function(df, patt_event_var, num_vars, event_codes, event_name, starter = 0) {

    #create new column to hold result
    new_var_name <- paste0(event_name, "_ever")
    df[[new_var_name]] <- 0

    #iterate through event variables
    for (i in starter:num_vars) {
        curr_event_var <- paste0(patt_event_var, i)
        df <- df %>% mutate(!!new_var_name := 
                        case_when(str_detect(!!as.name(curr_event_var), event_codes) ~ 1, 
                                    TRUE ~ !!as.name(new_var_name)))
    }
    return(df)
}

#function to check how many times an event happened
UKB_wide_event_count <- function(df, patt_event_var, num_vars, event_codes, event_name) {

    #create new column to hold result
    new_var_name <- paste0(event_name, "_count")
    df[[new_var_name]] <- 0

    #iterate through event variables
    for (i in 0:num_vars) {
        curr_event_var <- paste0(patt_event_var, i)
        df <- df %>% mutate(!!new_var_name := 
                        case_when(str_detect(!!as.name(curr_event_var), event_codes) ~ !!as.name(new_var_name) +1, 
                                TRUE ~ !!as.name(new_var_name)))
    }
    return(df)
}
#function to get the date of the first event
UKB_wide_event_first_date <- function(df, patt_event_var, patt_date_var,  num_vars, event_codes, event_name) {

    #create new column to hold result
    new_var_name <- paste0(event_name, "_first_date")
    df[[new_var_name]] <- NA_Date_

    #iterate through event variables
    for (i in 0:num_vars) {
        curr_event_var <- paste0(patt_event_var, i)
        curr_date_var <- paste0(patt_date_var, i)

        df <- df %>% mutate(!!new_var_name := 
                        case_when(str_detect(!!as.name(curr_event_var), event_codes) &
                                (is.na(!!as.name(new_var_name)) |  
                                as.Date(!!as.name(curr_date_var)) < !!as.name(new_var_name)) ~ as.Date(!!as.name(curr_date_var)),
                                        TRUE ~ !!as.name(new_var_name)))
    }
    return(df)
}


#function to get the date of the last event
UKB_wide_event_last_date <- function(df, patt_event_var, patt_date_var,  num_vars, event_codes, event_name) {

    #create new column to hold result
    new_var_name <- paste0(event_name, "_last_date")
    df[[new_var_name]] <- NA_Date_

    #iterate through event variables
    for (i in 0:num_vars) {
        curr_event_var <- paste0(patt_event_var, i)
        curr_date_var <- paste0(patt_date_var, i)

        df <- df %>% mutate(!!new_var_name := 
                        case_when(str_detect(!!as.name(curr_event_var), event_codes) &
                                (is.na(!!as.name(new_var_name)) |  
                                as.Date(!!as.name(curr_date_var)) > !!as.name(new_var_name)) ~ as.Date(!!as.name(curr_date_var)),
                                        TRUE ~ !!as.name(new_var_name)))
    }
    return(df)
}


#function to check if event happened in a lookback period
UKB_wide_event_in_lb <- function(df, patt_event_var, patt_date_var,  num_vars, date_col_name, lb_in_years, event_codes, event_name) {

    #create new column to hold result
    new_var_name <- paste0(event_name, "_in_lb")
    df[[new_var_name]] <- 0

    #iterate through event variables
    for (i in 0:num_vars) {
        curr_event_var <- paste0(patt_event_var, i)
        curr_date_var <- paste0(patt_date_var, i)

        df <- df %>% mutate(!!new_var_name := 
                        case_when((str_detect(!!as.name(curr_event_var), event_codes) &
                                as.Date(!!as.name(curr_date_var)) < as.Date(!!as.name(date_col_name)) &
                                as.Date(!!as.name(curr_date_var)) >= as.Date(!!as.name(date_col_name)) - years(lb_in_years)) ~ 1,
                                        TRUE ~ !!as.name(new_var_name)))
    }
    return(df)
}

#function to check if event happened after a specific date
UKB_wide_event_post <- function(df, patt_event_var, patt_date_var,  num_vars, date_col_name, event_codes, event_name) {

    #create new column to hold result
    new_var_name <- paste0(event_name, "_post")
    df[[new_var_name]] <- 0

    #iterate through event variables
    for (i in 0:num_vars) {
        curr_event_var <- paste0(patt_event_var, i)
        curr_date_var <- paste0(patt_date_var, i)

        df <- df %>% mutate(!!new_var_name := 
                        case_when((str_detect(!!as.name(curr_event_var), event_codes) &
                                as.Date(!!as.name(curr_date_var)) > as.Date(!!as.name(date_col_name))) ~ 1,
                                        TRUE ~ !!as.name(new_var_name)))
    }
    return(df)
}

#function to check if event happend before (or on) a specific date
UKB_wide_event_pre <- function(df, patt_event_var, patt_date_var,  num_vars, date_col_name, event_codes, event_name) {

    #create new column to hold result
    new_var_name <- paste0(event_name, "_pre")
    df[[new_var_name]] <- 0

    #iterate through event variables
    for (i in 0:num_vars) {
        curr_event_var <- paste0(patt_event_var, i)
        curr_date_var <- paste0(patt_date_var, i)

        df <- df %>% mutate(!!new_var_name := 
                        case_when((str_detect(!!as.name(curr_event_var), event_codes) &
                                as.Date(!!as.name(curr_date_var)) <= as.Date(!!as.name(date_col_name))) ~ 1,
                                        TRUE ~ !!as.name(new_var_name)))
    }
    return(df)
}