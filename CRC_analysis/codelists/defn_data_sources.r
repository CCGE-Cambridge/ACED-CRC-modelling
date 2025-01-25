# =============================================================================
# Script Title: Regroup variables from legacy categories
# Author: Sam Ip
# =============================================================================



# core
names_core <- c(
    "age", "age2",          # Age and its squared term
    "v0_birth_year",        # Year of birth
    "v0_bmi",               # Body Mass Index (BMI) at baseline assessment
    "v0_ethn_selfrep",      # Self-reported ethnicity at baseline assessment
    "v0_sex_genetic",       # Genetically determined sex
    "v0_smoking_status",    # Smoking status at baseline assessment
    "v0_townsend"           # Townsend deprivation index at baseline assessment
    )

# Genetics Variables ----
# All variables related to genetics are prefixed with "gen_".

# Medical History Variables ----
names_medhist <- c(
  "v0_famhist_bowel_cancer",      # Family history of bowel cancer
  "v0_famhist_breast_cancer",     # Family history of breast cancer
  "v0_famhist_lung_cancer",       # Family history of lung cancer
  "multimorb_mm_score_res_fit",   # Multimorbidity score
  "other_bcscreen_eligible",      # Eligibility for bowel cancer screening
  "other_colonoscopy_ALL_in_10_yr_lb", # Prior colonoscopy in the last 10 years
  "mod_IBD_ever",                 # History of inflammatory bowel disease
  "mod_DIB_T1D_ever",             # History of Type 1 diabetes
  "mod_DIB_T2D_ever",             # History of Type 2 diabetes
  "mod_GALL_ever",                # History of gallstones
  "mod_MED_NSAIDs_reg_ever",      # Regular NSAID use (ever)
  "mod_MED_ASPIRIN_reg_ever"      # Regular aspirin use (ever)
)

# Symptoms Variables ----
# All variables related to symptoms are prefixed with "mrk_".


# Biomarkers Variables ----
# All variables related to biomarkers are prefixed with "biom_".

# Lifestyle Variables ----
names_lifestyle <- c(
  "v0_alcohol_units_daily",  # Daily alcohol consumption
  "v0_edu_quals",            # Highest education qualification
  "v0_partial_fibre_score",  # Dietary fiber intake score
  "v0_tot_METs_wk",          # Weekly physical activity (METs)
  "v0_tot_process_meat_wk",  # Weekly processed meat intake
  "v0_tot_red_meat_wk"       # Weekly red meat intake
)