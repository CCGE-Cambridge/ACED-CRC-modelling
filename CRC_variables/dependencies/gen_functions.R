##File name: gen_functions.r
##Author: Hannah Harrison
##Last Edit: 03/03/2022
##Description: functions for performing simple operations on UKB data

###General date and age functions
gen_birth_dates <- function(year, month) {
    #gives everyone in UKB a birthdate for the 1st of the month they were born in
    return(as.Date(ISOdate(year, month, 1), "%Y-%m-%d"))
}

gen_event_age <- function(dob, event_dt) {
    #calculates age at a particular event using date of birth and event date
    return(as.numeric(event_dt - dob) / 365.25)
}


##custom rowwise operations
###rowwise sum with weighting (e.g. ax1 + bx2 + cx3 ....)
wtd.rowSums <- function(x, wts=1, na.rm=TRUE) {
  rowSums(t(t(x) * wts), na.rm=na.rm) }

###rowwise sum with missing data handling, NA if all are missing, otherwise ignore NA
na_handle.rowSums <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  z <- base::rowSums(x, na.rm = TRUE)
  z[!base::rowSums(!is.na(x))] <- NA
  z}

###combined above functions, weighted sum and missing data handling
wtd_and_na_handle.rowSums <- function(x, wts=1) {
  na_handle.rowSums(t(t(x) * wts)) }