


biomarker_manage <- "~/ACED-RREDD-EHR/CRC_variables/dependencies/biomarker_manage_march2022.r"
gp_data_funcs <- "~/ACED-RREDD-EHR/CRC_variables/dependencies/gp_data_funcs.r"
gen_functions <- "~/ACED-RREDD-EHR/CRC_variables/dependencies/gen_functions.r"
base_vars_funcs <- "~/ACED-RREDD-EHR/CRC_variables/dependencies/base_vars_functions.r"
biomarker_clean_funcs <- "~/ACED-RREDD-EHR/CRC_variables/dependencies/biomarker_cleaning.r"
CRC_biomarker_funcs <- "~/ACED-RREDD-EHR/CRC_variables/dependencies/colorectal_biomarker_covariates.r"
CRC_threshold_defs <- "~/ACED-RREDD-EHR/CRC_variables/dependencies/CRC_thresholds.r"
mmscore_defs_and_funcs <- "~/ACED-RREDD-EHR/CRC_variables/dependencies/multimorb_score_info_and_funcs_v2.r"
wide_baseline_funcs <- "~/ACED-RREDD-EHR/CRC_variables/dependencies/wide_baseline_funcs.R"


lapply(
    list(
        #wrangling_functions,
        #get_ukb_dob,
        biomarker_manage,
        gp_data_funcs,
        gen_functions,
        base_vars_funcs,
        biomarker_clean_funcs, 
        CRC_biomarker_funcs,
        CRC_threshold_defs,
        mmscore_defs_and_funcs,
        wide_baseline_funcs
    ),
    source
)
