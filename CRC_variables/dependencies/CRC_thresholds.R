##File name: CRC_thresholds.r
##Author: Hannah Harrison
##Last Edit: 28/03/2022
##Description: Thresholds used for CRC biomakers

CRP_thresh <- 6.8 #mg/L
PV_thresh <- 1.72 #mPa.s

ESR_age_groups <- factor(c(1, 2, 3, 4, 5, 6))# "<40", "40-49", "50-59", "60-69", "70-79", ">80"
men <- c(11, 12, 14, 14, 20, 20)#mm/hr
women <- c(14, 15, 17, 18, 22, 23)#mm/hr
df_ESR <- data.frame(ESR_age_groups, men, women)
df_ESR <- df_ESR %>% 
            pivot_longer(c("men", "women"), names_to = "sex", values_to = "ESR_thresh") %>% 
                mutate(sex =recode(sex, "men" = "1", "women" = "0")) #unique value by age group and sex

Hb_thresh <- c(120, 130) #g/L - women then men
Hmc_thresh <- c(37, 42) #  as percentage, women then men
Ferritin_thresh_low <- 30 #mg/L
Ferritin_thresh_high <- c(200, 300) #ug/L  
Plate_thresh <- 450 #(10^9)/L

CA125_thresh <- 35 # U/ml, equal to or more than
albumin_thresh <- 35 # g/L (less than)
MCV_thresh <- 85 #fl (taken from Hopkins 2022 et al.) - microcytosis is less than this

#NB: using low haemotcrit but ref is bad
#NB: no good threshold for CRP/albumin ratio, leave as continuous??