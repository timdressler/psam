# tid_psam_erp_analysis.R
#
# Performs analysis of behavioural data.
#
# Tim Dressler, 18.04.25

#-------------------------------------Set up------------------------------------
#load packages
library(tidyr) 
library(afex)
library(emmeans)
library(readxl)
library(car) 
library(corrplot) 
library(dplyr)
library(stringr)
library(ez)
library(ggplot2)
library(ggstatsplot)
library(DescTools)
library(ggpubr)
library(ggeffects)
library(cowplot) 
library(tidyverse)
library(psych)
library(RVAideMemoire)
library(rstatix)
library(lme4)

rm(list=ls())

options(scipen = 999)

# Set costum colors
colors <- list()
colors$main_blue <- "#004F9F"
colors$main_red <- "#D53D0E"
colors$main_green <- "#00786B"
colors$light_blue <- "#5BC5F2"
colors$main_yellow <- "#FDC300"
colors$UI <- "grey"

# Set up paths
if (interactive()) {
  # Running inside RStudio
  library(rstudioapi)
  SCRIPTPATH <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  # Running via Rscript (called from tid_psam_run_pipeline.py)
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
  SCRIPTPATH <- dirname(script_path)
}
if (grepl("psam/analysis_script/R", SCRIPTPATH)) {
  cat("Path OK") 
} else {
  stop("Path not OK")  
}

MAINPATH <- gsub("/analysis_script/R", "", SCRIPTPATH)
INPATH <- file.path(MAINPATH, "data", "analysis_data", "erp_analysis")
OUTPATH <- file.path(MAINPATH, "data", "analysis_data", "stats_erp_analysis")

FUNPATH <- file.path(MAINPATH, "functions")
source(file.path(FUNPATH, "tid_psam_check_folder_TD.R"))
source(file.path(FUNPATH, "tid_psam_clean_up_folder_TD.R"))

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

setwd(OUTPATH)

#-----------------------------------Load data-----------------------------------

# Load data
df_erp <- read_excel(file.path(INPATH, "all_subj_cor_erp_data.xlsx"))
df_erp$task_instruction <-as.factor(df_erp$task_instruction)
df_erp$probe_onset_cat <-as.factor(df_erp$probe_onset_cat)
df_erp$probe_type <-as.factor(df_erp$probe_type)

#-------------------------------Create needed dfs-------------------------------


# PSAM effect for probe-onset and probe_type for each subject
df_erp <- df_erp %>%
  mutate(temp_matching_condition = str_remove(condition_full, "^(act|pas)_")) # Create a temporary matching variable
temp_df_active <- df_erp %>%
  filter(task_instruction == "Active") %>%
  select(subj, temp_matching_condition, probe_onset_cat, probe_type,
         erp_amp, erp_lat, erp_from, erp_till, condition_full) %>%
  rename_with(~ paste0("active_", .), -c(subj, temp_matching_condition, probe_onset_cat, probe_type)) # Create a temporary df which includes only active conditions

temp_df_passive <- df_erp %>%
  filter(task_instruction == "Passive") %>%
  select(subj, temp_matching_condition, probe_onset_cat, probe_type,
         erp_amp, erp_lat, erp_from, erp_till, condition_full) %>%
  rename_with(~ paste0("passive_", .), -c(subj, temp_matching_condition, probe_onset_cat, probe_type)) # Create a temporary df which includes only passive conditions

df_psam <- left_join(temp_df_active, temp_df_passive,
                     by = c("subj", "temp_matching_condition", "probe_onset_cat", "probe_type")) %>%
  mutate(
    psam_amp = active_erp_amp - passive_erp_amp,
    psam_lat = active_erp_lat - passive_erp_lat
  ) %>%
  transmute(
    subj,
    erp_from = active_erp_from,
    erp_till = active_erp_till,
    psam_amp,
    psam_lat,
    condition_full = active_condition_full,
    probe_onset_cat,
    probe_type) # Get PSAM value and merge dfs

# Clean-Up
df_erp <- df_erp %>% select(-temp_matching_condition) # Delete the temporary matching variable
rm(temp_df_active, temp_df_passive) # Delete the temporary dfs


#------------------------------------Analysis-----------------------------------

# MAIN_ERPBASE
# Linear Mixed Model (Random Intercepts, Fixed Slopes) (DV = N1 Amplitude)
# Analysis MAIN_ERPBASE is included as a baseline model.
MAIN_ERPBASE <- lme4::lmer(erp_amp ~ 1 + (1|subj), data = df_erp)
summary(MAIN_ERPBASE)
performance::icc(MAIN_ERPBASE)

# Assumptions
performance::check_convergence(MAIN_ERPBASE)
performance::check_model(MAIN_ERPBASE)

#------------------------------------------------------------------------------#
#
#

# Assumptions
# (- Convergence: )
# - Normally distributed random effects:
# - Normally distributed residuals:  
# - Homogeneity of variance of the residuals: 
#------------------------------------------------------------------------------#

# MAIN_ERP1
# Linear Mixed Model (Random Intercepts, Fixed Slopes) (DV = N1 Amplitude, within = Task (active, passive), Probe-type (unaltered, altered), Probe-onset (early, late))
# Analysis MAIN_ERP1 concerns how N1 ERP amplitudes are influenced by probe-type, probe-onset and task.
MAIN_ERP1 <- lme4::lmer(erp_amp ~ task_instruction*probe_onset_cat*probe_type + (1|subj), data = df_erp)
summary(MAIN_ERP1)

# Plot: N1 ERP amplitude by probe-onset and probe-type
ezPlot(
  data = df_erp 
  , dv = erp_amp 
  , wid = subj  
  , within= .(task_instruction,probe_type, probe_onset_cat)
  , x = .(probe_type)
  , split   = .(task_instruction)
  , row = .(probe_onset_cat)
)

# Descriptive statistics
psych::describeBy(
  df_erp$erp_amp,
  list(df_erp$task_instruction,df_erp$probe_type, df_erp$probe_onset_cat)
)

# Follow-Up Tests
MAIN_ERP1_FUT <- emmeans(MAIN_ERP1, ~ task_instruction)
contrast(MAIN_ERP1_FUT, method = "pairwise", adjust = "bonferroni")

# Assumptions
performance::check_convergence(MAIN_ERP1)
performance::check_model(MAIN_ERP1)

#------------------------------------------------------------------------------#
#
#

# Assumptions
# (- Convergence: )
# - Normally distributed random effects:
# - Normally distributed residuals:  
# - Homogeneity of variance of the residuals: 
#------------------------------------------------------------------------------#

# MAIN_ERP2
# Linear Mixed Model (Random Intercepts, Fixed Slopes) (DV = PSAM effect, within = Probe-type (unaltered, altered), Probe-onset (early, late))
# Analysis MAIN_ERP2 concerns how the PSAM effect is influenced by probe-type, probe-onset. 
MAIN_ERP2 <- lme4::lmer(psam_amp ~ probe_onset_cat*probe_type + (1|subj), data = df_psam)
summary(MAIN_ERP2)

# Plot: PSAM effect by probe-onset and probe-type
ezPlot(
  data = df_psam 
  , dv = psam_amp 
  , wid = subj  
  , within= .(probe_type, probe_onset_cat)
  , x = .(probe_type)
  , split   = .(probe_onset_cat)
)

# Descriptive statistics
psych::describeBy(
  df_psam$psam_amp,
  list(df_psam$probe_type, df_psam$probe_onset_cat)
)

# Follow-Up Tests
MAIN_ERP2_FUT <- emmeans(MAIN_ERP2, ~ probe_type | probe_onset_cat)
contrast(MAIN_ERP2_FUT, method = "pairwise", adjust = "bonferroni")



# Assumptions
performance::check_convergence(MAIN_ERP1)
performance::check_model(MAIN_ERP1)

#------------------------------------------------------------------------------#
#
#

# Assumptions
# (- Convergence: )
# - Normally distributed random effects:
# - Normally distributed residuals:  
# - Homogeneity of variance of the residuals: 
#------------------------------------------------------------------------------#

# MAIN_ERP3
# Linear Mixed Model (Random Intercepts, Fixed Slopes) (DV = N1 Latency, within = Task (active, passive), Probe-type (unaltered, altered), Probe-onset (early, late))
# Analysis MAIN_ERP3 concerns how the N1 ERP latency is influenced by probe-type, probe-onset and task. 
MAIN_ERP3 <- lme4::lmer(erp_lat ~ task_instruction*probe_onset_cat*probe_type + (1|subj), data = df_erp)
summary(MAIN_ERP3)

# Plot: N1 ERP latency by probe-onset and probe-type
ezPlot(
  data = df_erp 
  , dv = erp_lat 
  , wid = subj  
  , within= .(task_instruction,probe_type)
  , within_full= .(task_instruction,probe_type, probe_onset_cat)
  , x = .(probe_type)
  , split   = .(task_instruction)
)

# Descriptive statistics# Descriptive statiprobe_onset_catstics
psych::describeBy(
  df_erp$erp_lat,
  list(df_erp$task_instruction,df_erp$probe_type)
)

# Follow-Up Tests
MAIN_ERP3_FUT <- emmeans(MAIN_ERP3, ~ probe_type | task_instruction)
contrast(MAIN_ERP3_FUT, method = "pairwise", adjust = "bonferroni")


# Assumptions
performance::check_convergence(MAIN_ERP3)
performance::check_model(MAIN_ERP3)

#------------------------------------------------------------------------------#
#
#

# Assumptions
# (- Convergence: )
# - Normally distributed random effects:
# - Normally distributed residuals:  
# - Homogeneity of variance of the residuals: 
#------------------------------------------------------------------------------#

















