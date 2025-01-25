##File name: biomarker_cleaning.r
##Author: Hannah Harrison
##Last Edit: 29/10/2021
##Description: cleans specified variables (e.g. biomakers) from gp data 
##Possible to expand to clean symptom data too?

#libraries if needed here

rescale_dist <- function(my_df, my_var, my_range, my_factor) {
    return(my_df %>% 
        mutate(value = ifelse(biomarker == my_var & !is.na(value) &
                                value > my_range[1] & value < my_range[2],
                                value * my_factor, value)))

}

trim_dist <- function(my_df, my_var, my_range) {
    return(my_df %>%
         mutate(value = ifelse(biomarker == my_var & (value < my_range[1] | value > my_range[2]), NA, value)))
}

trim_dist_and_rm <- function(my_df, my_var, my_range) {
    return(my_df %>%
         mutate(value = ifelse(biomarker == my_var & (value < my_range[1] | value > my_range[2]), NA_real_, value)) %>%
            filter(!is.na(value)))
}

value_extract <- function(my_df) { 
#only pass dataframe for which you want values, other types of readcodes are discarded
#may generate coerced NA warning
    return(my_df %>% 
            filter(!is.na(value1) | !is.na(value2) | !is.na(value3)) %>%
                mutate(value = ifelse(data_provider == 2, value2, value1)) %>%
                    mutate(value = as.numeric(na_if(value, ""))) %>%
                        filter(!is.na(value)))
}