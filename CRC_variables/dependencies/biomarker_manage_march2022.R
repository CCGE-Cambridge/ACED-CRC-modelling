##File name: biomarker_manager.r
##Author: Hannah Harrison (Sam edited for merge with pipeline)
##Last Edit: 28/03/2022
##Description:Preps EHR biomaker dataframe then defines variables for colorectal cancer

####Load any libraries (most now added in config)
#library(tidyverse) 
#library(readxl)  
#library(tidyr)
#library(dplyr)
#library(readr)
#library(data.table)

##Define any functions or variables needed (now added in config file)
#source(file.path("~/ACED-RREDD-EHR/CRC_variables/dependencies/colorectal_biomarker_covariates.r")) #defines functions that create summary biomarker vars
#source(file.path("~/ACED-RREDD-EHR/CRC_variables/dependencies/biomarker_cleaning.r")) #defines functions that clean biomarker variables
#source(file.path("~/ACED-RREDD-EHR/CRC_variables/dependencies/CRC_thresholds.r")) #defines thresholds and ranges for biomarker variables
#source(file.path("~/ACED-RREDD-EHR//CRC_variables/dependencies/gen_functions.r"))

biomarker_basic_clean_CRC <- function(df_gp_biom){
    ####start with gp biom data (raw biomarkers already defined), remove bad values 
    df_gp_biom <- df_gp_biom %>% filter(
        (!is.na(event_dt)) & 
        (event_dt > as.Date("1960-01-01")) &
        (event_dt < as.Date("2017-09-18")) &  #last date allowed by UKB
        (!is.na(value))
        ) #each row must have a date and a measurement

    ##All biomarkers, remove zeros, NAs
    df_gp_biom <- df_gp_biom %>% filter(value != 0) 
     
    #remove measurements for which zero was recorded
    ##Iron deficiency variables
    ####Hb - two distinct distributions, rescaled to g/L
    df_gp_biom <- rescale_dist(df_gp_biom, "Haemoglobinconc", c(0, 30), 10^1)#moves measurements in g/dL to g/L
    df_gp_biom <- trim_dist(df_gp_biom, "Haemoglobinconc", c(30, 300))#removes extremes

    ####Haemocrit - two distinct ditributions, rescaled decimal values (<1) to percentage values 
    df_gp_biom <- rescale_dist(df_gp_biom, "Haematocritperc", c(0, 1), 10^2)#moves measurements in g/dL to g/L
    df_gp_biom <- trim_dist(df_gp_biom, "Haematocritperc", c(0, 100)) #removes extremes

    ####Ferritin - unclear if there is a secondary distribution at low values(?), value3 units absolutely no help
    df_gp_biom <- trim_dist(df_gp_biom, "Ferritin", c(0, 800)) #removes extremes

    ##Thrombocytosis variables
    ####Platelets
    df_gp_biom <- rescale_dist(df_gp_biom, "Platelets", c(10^6, 10^10), 10^(-6))#moves measurements in /l to 10^9/l
    df_gp_biom <- trim_dist(df_gp_biom, "Platelets", c(0, 1500)) #removes extremes
    
    ##Inflammatory Markers
    ####C-reactive Protein (CRP) - seems like all data is in mg/L - but long tail so hard to tell
    df_gp_biom <- trim_dist(df_gp_biom, "CRP", c(0, 200)) #removes extremes

    ####Erythrocyte Sedimentation Rate (ESR) - seems like all data is in mm/hr - but long tail so hard to tell
    df_gp_biom <- trim_dist(df_gp_biom, "ESR", c(0, 250)) #removes extremes

    ####Plasma Viscosity (PV) - found three distributions, rescaled to mPa.s
    df_gp_biom <- rescale_dist(df_gp_biom, "PV", c(10, 100), 10^(-1))
    df_gp_biom <- rescale_dist(df_gp_biom, "PV", c(99.9, 300), 10^(-2))
    df_gp_biom <- trim_dist(df_gp_biom, "PV", c(0, 3)) #removes extremes

    ###Albumin - found some measurements <10, probably in g/dL (instead of g/L)
    df_gp_biom <- rescale_dist(df_gp_biom, "Albumin", c(0, 10), 10^1)
    df_gp_biom <- trim_dist(df_gp_biom, "Albumin", c(10, 100)) #removes extremes

    ###CA125 - unclear what is a spuriously high measure, seems to vary by lab (300 U/ml is right order of magnitude)
    df_gp_biom <- trim_dist(df_gp_biom, "CA125", c(1, 300)) #removes extremes

    ##MCV, typically use units of fL (expected range is 80-100fL), these is a small extra distribution around 30 - can't work out what the units would be - 
    ##doesn't look quite like hmc entered in error, will remove 
    df_gp_biom <- trim_dist(df_gp_biom, "MCV", c(50, 200)) #removes extremes, range taken from hopkins et al. 2020


    return(df_gp_biom)
}

biomarker_meta_CRC <- function(df_gp_biom){
    ###load threshold information
    #source(file.path("~/ACED-RREDD-EHR/hannah/colorectal_biomarkers/CRC_thresholds.r"))
    #has now moved and is loaded in config

    ##identify "abnormal" measurements for each biomarker
    df_gp_biom_plus <- flag_abnormal_albumin(df_gp_biom_plus, albumin_thresh) 
    df_gp_biom_plus <- flag_abnormal_CRP(df_gp_biom_plus, CRP_thresh) 
    df_gp_biom_plus <- flag_abnormal_ESR(df_gp_biom_plus, df_ESR) 
    df_gp_biom_plus <- flag_abnormal_PV(df_gp_biom_plus, PV_thresh) 
    df_gp_biom_plus <- flag_abnormal_Hgb(df_gp_biom_plus, Hgb_thresh) 
    df_gp_biom_plus <- flag_abnormal_Hmc(df_gp_biom_plus, Hmc_thresh) 
    df_gp_biom_plus <- flag_abnormal_Ferritin_low(df_gp_biom_plus, Ferritin_thresh_low) 
    df_gp_biom_plus <- flag_abnormal_Ferritin_high(df_gp_biom_plus, Ferritin_thresh_high) 
    df_gp_biom_plus <- flag_abnormal_Plate(df_gp_biom_plus, Plate_thresh) 
    df_gp_biom_plus <- flag_abnormal_MCV(df_gp_biom_plus, MCV_thresh) 
    df_gp_biom_plus <- flag_abnormal_CA125(df_gp_biom_plus, CA125_thresh) 

    df_gp_biom_plus_group <- df_gp_biom_plus %>% group_by(eid)
    df_gp_biom_plus_group <- df_gp_biom_plus_group %>% 
            mutate(albumin_tests = sum(biomarker == "Albumin"), 
                  CRP_tests = sum(biomarker == "CRP"), 
                  ESR_tests = sum(biomarker == "ESR"), 
                  Hgb_tests =sum(biomarker == "Haemoglobinconc"), 
                  Hmc_tests = sum(biomarker =="Haematocritperc"),
                  Ferritin_tests = sum(biomarker =="Ferritin"),
                  PV_tests = sum(biomarker == "PV"),
                  Plate_tests = sum(biomarker == "Platelets"),
                  MCV_tests = sum(biomarker == "MCV"),  
                  CA125_tests = sum(biomarker == "CA125"))

    df_gp_biom_plus_group <- df_gp_biom_plus_group %>% 
            mutate(albumin_abnormal = sum(albumin_flag), 
                  CRP_abnormal = sum(CRP_flag), 
                  ESR_abnormal = sum(ESR_flag), 
                  Hgb_abnormal =sum(Hgb_flag), 
                  Hmc_abnormal = sum(Hmc_flag),
                  Ferritin_low_abnormal = sum(Ferritin_low_flag),
                  Ferritin_high_abnormal = sum(Ferritin_high_flag),
                  PV_abnormal = sum(PV_flag),
                  Plate_abnormal = sum(Plate_flag),
                  MCV_abnormal = sum(MCV_flag),  
                  CA125_abnormal = sum(CA125_flag))

    ##keep one row per eid with all summary data    
    df_gp_biom_by_eid <- df_gp_biom_plus_group %>% filter(row_number()==1) %>% select("eid", "sex", contains(c("_tests", "_abnormal"))) %>% ungroup
    
}
    
    ##for each person and each biomarker - did measurement take place, was result "abnormal"

    ####creates column gp_anaemia (0 or 1, where 1 indicates that this measurement shows iron deficiency)
    #df_gp_biom_plus <- anaemia_var(df_gp_biom_plus, Hb_thresh, Hmc_thresh, Ferritin_thresh) 
    ####creates column gp_thrombo (0 or 1, where 1 indicates that this measurement shows thrombocytosis (high platelet count)))
    #df_gp_biom_plus <- thrombo_var(df_gp_biom_plus, Platelet_thresh)
    ####creates column gp_inflam (0 or 1 , where 1 indicates a raised inflammatory marker)
    #df_gp_biom_plus <- inflam_var(df_gp_biom_plus, CRP_thresh, PV_thresh, df_ESR)

    ##for each person.... 
    ### indiviual specific flags ----
    #flags <- c("gp_anaemia", "gp_thrombo", "gp_inflam")
    #df_gp_biom_plus <- df_gp_biom_plus %>% dplyr::select( c("eid", flags))
    #df_gp_biom_indiv <- df_gp_biom_plus %>% group_by(eid) %>% summarise_at(flags, sum)
    #names(df_gp_biom_indiv) <- c("eid", "flag_anaemia", "flag_thrombo", "flag_inflam")


##example histogram to visualise distribution (MCV)

# source("~/ACED-RREDD-EHR/hannah/wd/graph_theme.r")
# my_dir <- "~/ACED-RREDD-EHR/hannah/colorectal_biomarkers"
# temp_MCV <- df_gp_biom %>% filter(biomarker == "MCV") %>% filter(value < 500)
# p1 <- ggplot(temp_MCV) + 
#         geom_histogram(aes(x = value), bins = 55, color = "black", fill = "white")  + 
#         coord_cartesian(xlim = c(5, 150)) + 
#         xlab("MCV (fL)") + 
#         ylab("Number of Measurements") + 
#         theme_hh()
# ggsave(file.path(my_dir, "MCV_graphs/graph_1.png"),  plot = p1, dpi = 300, type = "cairo")

# #test low dist

# temp_MCV_low <- temp_MCV %>% filter(value > 10 & value <50)
# p2 <- ggplot(temp_MCV_low) + 
#         geom_histogram(aes(x = value), bins = 55, color = "black", fill = "white")  + 
#         xlab("MCV (fL)") + 
#         ylab("Number of Measurements") + 
#         theme_hh()
# ggsave(file.path(my_dir, "MCV_graphs/graph_2.png"),  plot = p2, dpi = 300, type = "cairo")

# #test haemocrit
# temp_hmc <- df_gp_biom %>% filter(biomarker == "Haematocritperc") %>% filter(value > 10 & value < 50)
# p3 <- ggplot(temp_hmc) + 
#         geom_histogram(aes(x = value), bins = 55, color = "black", fill = "white")  + 
#         xlab("Haematocrit") + 
#         ylab("Number of Measurements") + 
#         theme_hh()
# ggsave(file.path(my_dir, "MCV_graphs/graph_3.png"),  plot = p3, dpi = 300, type = "cairo")