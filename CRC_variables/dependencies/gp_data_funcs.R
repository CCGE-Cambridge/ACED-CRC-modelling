##File name: gp_data_code_funcs.r
##Author: Hannah Harrison
##Last Edit: 25/11/2021
##Description: functions to for cleaning and querying EHR (basic and generalisable)
##(handling for both gp_clinical and gp_script)

gp_date_censor <- function(gp_dataset) { 
    ##start date: readv2 and read ctv3 codes in use since 1985
    gp_dataset <- gp_dataset %>% filter(event_dt > as.Date("1985-01-01"))

    ##end date varies by data provider, most conservative date (from Table 6, UKB primary care doc) used
    ##English Vision data [25th May 2017]
    ##Scottish data (EMIS/Vision) [19th April 2017]
    ##English TTP data [14th Jun 2016]
    ##welsh data (EMIS/Vision) [18th September 2017]
    gp_dataset <- gp_dataset %>% 
                    filter((data_provider == 1 & event_dt < as.Date("2017-05-25")) |
                           (data_provider == 2 & event_dt < as.Date("2017-04-19")) |
                           (data_provider == 3 & event_dt < as.Date("2016-06-14")) |
                           (data_provider == 4 & event_dt < as.Date("2017-09-18")))
    return(gp_dataset)
}

gp_script_query <- function(df_gp_scripts,  df_codelist) {
    #takes gp_script and pre-formatted codelist
    #searches by data provider and recombines final list of events
    #formats issue date to event date
    #call function to censor by data provider 
    df_script_events_scot <- inner_join(df_gp_scripts[data_provider == 2], df_codelist[Type == "bnf_6_char"], by = c("bnf_code" = "Code"))
    df_script_events_wales <- inner_join(df_gp_scripts[data_provider == 4], df_codelist[Type == "readv2"], by = c("read_2" = "Code"))
    df_script_events_engvis <- inner_join(df_gp_scripts[data_provider == 1], df_codelist[Type == "dmd"], by = c("dmd_code" = "Code"))
    df_script_events_engttp <- inner_join(df_gp_scripts[data_provider == 3], df_codelist[Type == "bnf_chap_id"], by = c("bnf_code" = "Code"))
    df_script_events <- rbind(df_script_events_scot, df_script_events_wales, df_script_events_engvis, df_script_events_engttp)

    df_script_events$event_dt <- as.Date(df_script_events$issue_date, format = "%d/%m/%Y")
    df_script_events <- gp_date_censor(df_script_events)
    return(df_script_events)
}

gp_clinical_query <- function(df_gp_clinical,  df_codelist) {
    #takes gp_script and pre-formatted codelist
    #searches readcode type and recombines final list of events
    #formats event date
    #call function to censor by data provider 
    df_events_rv2 <- inner_join(df_gp_clinical, df_codelist[Type == "read_v2"], by = c("read_2" = "Code"))
    df_events_rctv3 <- inner_join(df_gp_clinical, df_codelist[Type == "read_ctv3"], by = c("read_3" = "Code"))
    df_events <- rbind(df_events_rv2, df_events_rctv3)

    df_events$event_dt <- as.Date(df_events$event_dt, format = "%d/%m/%Y")
    df_events <- gp_date_censor(df_events)
    return(df_events)
}

##The following functions
##Accept: for a dataframe grouped by eid, with a list of events, date of event included. 
##Generate: for each person, a column indicating if their set of events meet specific criteria

event_ever <- function(df, event_name) {#ever had clinical event of interest
  varname <- paste(event_name, "ever", sep="_")
  df <- df %>% mutate(!!varname := case_when(any(event_type == event_name) ~ 1, TRUE ~0))
  return(df)
}

event_in_lb <- function(df, event_name) {#had clinical event of interest in fixed 2 year lookback period
  varname <- paste(event_name, "in_lb", sep="_")
  df <- df %>% mutate(temp_2_yr_lb = case_when((event_type == event_name & event_dt >= two_yr_lb & event_dt <= cutoff_max) ~ 1, TRUE ~ 0))
  df <- df %>% mutate(!!varname := case_when(any(temp_2_yr_lb == 1) ~ 1, TRUE ~0))
  return(df)
  }

event_in_var_lb <- function(df, event_name, lb) {#had clinical event of interest in dynamic lookback period (given in years)
  varname <- paste(event_name, "in", lb, "yr_lb", sep="_")
  df <- df %>% mutate(temp_x_yr_lb = case_when((event_type == event_name & event_dt >= (cutoff_max - years(lb)) & event_dt <= cutoff_max) ~ 1, TRUE ~ 0))
  df <- df %>% mutate(!!varname := case_when(any(temp_x_yr_lb == 1) ~ 1, TRUE ~0))
  return(df)
  }

event_freq_in_lb <- function(df, event_name) {#frequency of clinical event in lookback period
  varname <- paste(event_name, "freq_in_lb", sep="_")
  df <- df %>% mutate(temp_2_yr_lb = case_when((event_type == event_name & event_dt >= two_yr_lb & event_dt <= cutoff_max) ~ 1, TRUE ~ 0))
  df <- df %>% mutate(temp = case_when(any(temp_2_yr_lb == 1) ~ 1, TRUE ~0))
  df <- df %>% mutate(temp_sum = ifelse(temp == 1, sum(temp_2_yr_lb), 0))
  df <- df %>% mutate(!!varname := ifelse(temp_sum > 4, 4, temp_sum)) # frequency is truncated to 4. 
  return(df)
}

event_new_onset <- function(df, event_name) {#reported instance of clincal event in lb, but no reported instances in preceeding 3 years
  varname <- paste(event_name, "new_onset", sep="_")
  df <- df %>% mutate(temp_2_yr_lb = case_when((event_type == event_name & event_dt >= two_yr_lb & event_dt <= cutoff_max) ~ 1, TRUE ~ 0))
  df <- df %>% mutate(temp_pre_lb = case_when((event_type == event_name & event_dt >= five_yr_lb & event_dt < two_yr_lb)  ~ 1, TRUE ~ 0))
  df <- df %>% mutate(!!varname := case_when((any(temp_2_yr_lb == 1) & !any(temp_pre_lb == 1)) ~ 1, TRUE ~ 0))
  return(df)
}

event_freq_new_onset <- function(df, event_name) {#frequency of new onset clinical event in lookback period
  varname <- paste(event_name, "freq_new_onset", sep="_")
  df <- df %>% mutate(temp_2_yr_lb = case_when((event_type== event_name & event_dt >= two_yr_lb & event_dt <= cutoff_max) ~ 1, TRUE ~ 0))
  df <- df %>% mutate(temp_pre_lb = case_when((event_type == event_name & event_dt >= five_yr_lb & event_dt < two_yr_lb)  ~ 1, TRUE ~ 0))
  df <- df %>% mutate(temp = case_when((any(temp_2_yr_lb == 1) & !any(temp_pre_lb == 1)) ~ 1, TRUE ~ 0))                                  
  df <- df %>% mutate(!!varname := ifelse(temp == 1, sum(temp_2_yr_lb), 0))
  return(df)
}

script_ever <- function(df, event_name) {#ever has recorded instance of a prescription for event 
  varname <- paste("MED", event_name, "ever", sep="_")
  df <- df %>% mutate(!!varname := case_when(any(script_type == event_name) ~ 1, TRUE ~0))
  return(df)
}
script_in_lb <- function(df, event_name) {#had a relevant prescription in the lookback period (2 years)
  varname <- paste("MED", event_name, "in_lb", sep="_")
  df <- df %>% mutate(temp_2_yr_lb = case_when((script_type == event_name & issue_date >= two_yr_lb & issue_date <= cutoff_max) ~ 1, TRUE ~ 0))
  df <- df %>% mutate(!!varname := case_when(any(temp_2_yr_lb == 1) ~ 1, TRUE ~0))
  return(df)
}

script_freq_in_lb <- function(df, event_name) {#number of relevant prescriptions in the llokback period
  varname <- paste("MED", event_name, "freq_in_lb", sep="_")
  df <- df %>% mutate(temp_2_yr_lb = case_when((script_type == event_name & issue_date >= two_yr_lb & issue_date <= cutoff_max) ~ 1, TRUE ~ 0))
  df <- df %>% mutate(temp = case_when(any(temp_2_yr_lb == 1) ~ 1, TRUE ~0))
  df <- df %>% mutate(!!varname := ifelse(temp == 1, sum(temp_2_yr_lb), 0))
  return(df)
}

script_reg_ever <- function(df, event_name, script_num, period_in_years) {#more than X prescriptions in any X year period - caution this takes a while
  varname <- paste("MED", event_name, "reg_ever", sep="_")
  #df_subset <- df %>% select(script_type == event_name)
  df <- df %>% arrange(script_type, issue_date, .by_group=TRUE)
  df <- df %>% mutate(lag_date = lag(issue_date, n = script_num - 1))
  df <- df %>% mutate(comp_type = lag(script_type, n = script_num - 1))
  df <- df %>% mutate(time_diff = as.numeric((issue_date-lag_date)/365.25)) #time difference in years 

  df <- df %>% mutate(temp = case_when((script_type == event_name & comp_type == event_name & time_diff <= period_in_years) ~ 1, TRUE ~ 0))
  df <- df %>% mutate(!!varname := case_when(any(temp == 1) ~ 1, TRUE ~0))

  return(df)
}

script_new_onset <- function(df, event_name) {#had a new prescription in the lookback period (none in preceeding 3 years)
  varname <- paste("MED", event_name, "new_onset", sep="_")
  df <- df %>% mutate(temp_2_yr_lb = case_when((script_type == event_name & issue_date >= two_yr_lb & issue_date <= cutoff_max) ~ 1, TRUE ~ 0))
  df <- df %>% mutate(temp_pre_lb = case_when((script_type == event_name & issue_date >= five_yr_lb & issue_date < two_yr_lb)  ~ 1, TRUE ~ 0))
  df <- df %>% mutate(!!varname := case_when((any(temp_2_yr_lb == 1) & !any(temp_pre_lb == 1)) ~ 1, TRUE ~ 0))
  return(df)
}