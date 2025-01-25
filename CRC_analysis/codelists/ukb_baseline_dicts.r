# ==============================================================================
# Script Title: Dictionaries and Level Order for Recoding Categorical Variables
# Author: Sam Ip
# ==============================================================================

# smoke ----
dict_smoking_status <- data.frame("original" =  c("0", "1","2", "-3"), 
                                  "new" = c("never", "previous", "current", "missing"))

# alcohol_freq ----
dict_alcohol_freq <- data.frame("original" =  c("1","2", "3", "4", "5", "6", "-3"), 
                                "new" = c("daily", "three_four_week", "once_twice_week", "monthly_or_special", "monthly_or_special", "never", "missing"))

# edu_quals ----
dict_edu_quals <- data.frame("original" =  c("1","2", "3", "4", "5", "6", "-7", "-3"), 
                             "new" = c("university_college", 
                                       "Alevels", 
                                       "GCSE", "GCSE", 
                                       "technical_college", 
                                       "GCSE", "none", "missing"))
order_levels_edu_qual <- c("university_college", 
                           "technical_college",
                           "Alevels", 
                           "GCSE", 
                           "none", 
                           "missing")

# ibs ----
dict_ibs <- data.frame("original" =  c("0", "1","-121", "-818"), 
                       "new" = c("0", "1", "0", "0"))

# nsaid ----
dict_nsaid <- data.frame("original" =  c("1", "2", "3", "4", "5", "6", "-7", "-1", "-3"), 
                         "new" = c("pain_med", "pain_med", "pain_med", 
                                   "heartburn_med", "heartburn_med",
                                   "constipation_med",
                                   "none",
                                   "none", "none"))
# abs ----
dict_abs_nonfreq <- data.frame("original" =  c("0", "1", "-818"), 
                         "new" = c("0", "1", "0"))
dict_abs_freq <- data.frame("original" =  c("-818"), 
                         "new" = c("0"))

# famhist ----
dict_famhist <- data.frame("original" =  c(1:6, 8:13, -11, -13, -17, -21, -23, -27), 
                           "new" = c(
                             "heart",
                             "heart",
                             "lung_cancer",
                             "bowel_cancer",
                             "breast_cancer",
                             "chron_bronchi_emphy",
                             "heart",
                             "diabetes",
                             "alz_demen_park",
                             "alz_demen_park",
                             "severe_depression",
                             "prostate_cancer",
                             "missing1",
                             "missing1",
                             "none1",
                             "missing2",
                             "missing2",
                             "none2"
                           ))

# meds_hrt_bpchol_insulin ----
dict_meds_hrt_bpchol_insulin <- data.frame("original" =  c(1:5, -7, -1, -3), 
                               "new" = c(
                                 "med_chol",
                                 "med_bp",
                                 "med_insulin",
                                 "med_hrt",
                                 "med_cocp",
                                 "none",
                                 "missing",
                                 "missing"
                               ))

# ethn_selfrep ----
dict_ethn_selfrep <- data.frame(
  "original" =  c(
    "1", "1001", "1002", "1003", # White
    "2", "2001", "2002", "2003", "2004", # Mixed
    "3", "3001", "3002", "3003", "3004", # Asian
    "4", "4001", "4002", "4003", # Black
    "5", # Chinese
    "6", # Other
    "-1", "-3" # Missing
  ), 
  "new" = c(
    rep("White", 4), rep("Mixed", 5), rep("Asian", 5), rep("Black", 4), "Chinese", "Other", rep("missing", 2))
)