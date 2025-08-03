# tid_psam_erp_analysis.R
#
# Performs analysis of ERP data.
#
# Tim Dressler, 18.04.25

#-------------------------------------Set up------------------------------------
# Load packages
library(Rmisc)
library(tidyr) 
library(afex)
library(emmeans)
library(readxl)
library(car) 
library(corrplot) 
library(dplyr)
library(stringr)
library(ez)
library(lmerTest)
library(ggplot2)
library(ggstatsplot)
library(DescTools)
library(ggpubr)
library(ggeffects)
library(cowplot) 
library(tidyverse)
library(psych)
library(rstatix)
library(lme4)
library(RVAideMemoire)

rm(list=ls())
set.seed(123)
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

# Normalize path separators for consistent matching
normalized_path <- gsub("[\\]", "/", SCRIPTPATH)

# Check if path is correct
if (grepl("psam/analysis_script/R", normalized_path)) {
  cat("Path OK\n") 
} else {
  stop("Path not OK")  
}

MAINPATH <- sub("/analysis_script/R$", "", normalized_path)
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

# Within SEs for plots
df_erp_withinSE <- summarySEwithin(
  df_erp,
  measurevar = "erp_amp",
  withinvars = c("task_instruction", "probe_type", "probe_onset_cat"),
  idvar = "subj"
)

df_psam_withinSE <- summarySEwithin(
  df_psam,
  measurevar = "psam_amp",
  withinvars = c("probe_type", "probe_onset_cat"),
  idvar = "subj"
)

df_erp_lat_withinSE <- summarySEwithin(
  df_erp,
  measurevar = "erp_lat",
  withinvars = c("task_instruction", "probe_type", "probe_onset_cat"),
  idvar = "subj"
)

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
# Assumptions
# (- Convergence: OK)
# - Normally distributed random effects: OK
# - Normally distributed residuals: OK
# - Homogeneity of variance of the residuals: OK
#------------------------------------------------------------------------------#

# MAIN_ERP1
# Linear Mixed Model (Random Intercepts, Fixed Slopes) (DV = N1 Amplitude, within = Task (active, passive), Probe-type (unaltered, altered), Probe-onset (early, late))
# Analysis MAIN_ERP1 concerns how N1 ERP amplitudes are influenced by probe-type, probe-onset and task.
MAIN_ERP1 <- lmer(erp_amp ~ task_instruction*probe_onset_cat*probe_type + (1|subj), data = df_erp)
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
# task_instruction sig.; all other n.s.

# Assumptions
# (- Convergence: OK)
# - Normally distributed random effects: OK
# - Normally distributed residuals: OK
# - Homogeneity of variance of the residuals: OK 
#------------------------------------------------------------------------------#

# MAIN_ERP2
# Linear Mixed Model (Random Intercepts, Fixed Slopes) (DV = PSAM effect, within = Probe-type (unaltered, altered), Probe-onset (early, late))
# Analysis MAIN_ERP2 concerns how the PSAM effect is influenced by probe-type, probe-onset. 
MAIN_ERP2 <- lmer(psam_amp ~ probe_onset_cat*probe_type + (1|subj), data = df_psam)
summary(MAIN_ERP2)

# Plot: PSAM effect by probe-onset and probe-type
ezPlot(
  data = df_psam 
  , dv = psam_amp 
  , wid = subj  
  , within= .(probe_type, probe_onset_cat)
  , x = .(probe_onset_cat)
  , split   = .(probe_type)
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
performance::check_convergence(MAIN_ERP2)
performance::check_model(MAIN_ERP2)

#------------------------------------------------------------------------------#
# All n.s.

# Assumptions
# (- Convergence: OK)
# - Normally distributed random effects: OK
# - Normally distributed residuals: OK
# - Homogeneity of variance of the residuals: OK
#------------------------------------------------------------------------------#

# MAIN_ERP3
# Linear Mixed Model (Random Intercepts, Fixed Slopes) (DV = N1 Latency, within = Task (active, passive), Probe-type (unaltered, altered), Probe-onset (early, late))
# Analysis MAIN_ERP3 concerns how the N1 ERP latency is influenced by probe-type, probe-onset and task. 
MAIN_ERP3 <- lmer(erp_lat ~ task_instruction*probe_onset_cat*probe_type + (1|subj), data = df_erp)
summary(MAIN_ERP3)

# Plot: N1 ERP latency by probe-onset and probe-type
ezPlot(
  data = df_erp 
  , dv = erp_lat 
  , wid = subj  
  , within= .(task_instruction,probe_type, probe_onset_cat)
  , x = .(probe_type)
  , split   = .(task_instruction)
  , row = .(probe_onset_cat)
)

# Descriptive statistics
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
# all n.s.

# Assumptions
# (- Convergence: OK)
# - Normally distributed random effects: OK
# - Normally distributed residuals: OK
# - Homogeneity of variance of the residuals: OK 
#------------------------------------------------------------------------------#


#-------------------------------------Plots-------------------------------------

# Plot 1: ERP Amplitudes by Probe Type, Task Instruction, and Probe Onset
P1 <- ggplot(df_erp_withinSE, aes(x = task_instruction, y = erp_amp, color = probe_type, group = probe_type)) +
  geom_line(aes(y = erp_amp), size = 1.2) +
  geom_errorbar(aes(ymin = erp_amp - se, ymax = erp_amp + se), width = 0.2, size = 0.8) +
  geom_point(size = 3) +
  facet_wrap(~ probe_onset_cat, nrow = 1,
             labeller = as_labeller(c("Early" = "Early Probe Onset",
                                      "Late" = "Late Probe Onset"))) +
  scale_color_manual(values = c(
    "Altered" = colors$main_red,
    "Unaltered" = colors$main_blue
  )) +
  labs(
    #title = "N1 ERP Amplitudes by Probe Type, Task Instruction, and Probe Onset",
    y = "N1 ERP Amplitude [µV]",
    x = "Task Condition",
    color = "Probe Type"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "top"
  ) + 
  scale_x_discrete(labels = c("Active" = "Active",
                              "Passive" = "Passive"))
P1

# Save plot
ggsave(
  filename = "tid_psam_erp_amp.png", 
  plot = P1,
  width = 8,      
  height = 6,     
  dpi = 900,
  bg = "white"
)

# Plot 2: PSAM Effect Amplitudes by Probe Type and Probe Onset
P2 <- ggplot(df_psam_withinSE, aes(x = probe_onset_cat, y = psam_amp, color = probe_type, group = probe_type)) +
  geom_line(aes(y = psam_amp), size = 1.2) +
  geom_errorbar(aes(ymin = psam_amp - se, ymax = psam_amp + se), width = 0.2, size = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Altered" = colors$main_red,
    "Unaltered" = colors$main_blue
  )) +
  labs(
    #title = "PSAM Effect Amplitudes by Probe Type and Probe Onset",
    y = "PSAM Effect Amplitude [µV]",
    x = "Probe Onset",
    color = "Probe Type"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "top"
  ) + 
  scale_x_discrete(labels = c("Early" = "Early",
                              "Late" = "Late"))
P2

# Save plot
ggsave(
  filename = "tid_psam_psam_effect.png", 
  plot = P2,
  width = 8,      
  height = 6,     
  dpi = 900,
  bg = "white"
)

# Plot 3: PSAM Effect Amplitudes by Probe Type and Probe Onset (spreaded-out)
# Preparation
df_psam$probe_combination <- interaction(df_psam$probe_onset_cat, df_psam$probe_type, sep = "_")
df_psam$probe_combination <- factor(df_psam$probe_combination,
                                    levels = c("Early_Altered", "Early_Unaltered", "Late_Altered", "Late_Unaltered"),
                                    labels = c("Early Onset\nAltered", "Early Onset\nUnaltered", "Late Onset\nAltered", "Late Onset\nUnaltered"))

max_amp <- max(df_psam$psam_amp, na.rm = TRUE)
y_pos_level1 <- max_amp + 6 
y_pos_level2 <- max_amp + 7.5 
y_pos_level3 <- max_amp + 9 

# Plot
P3 <- ggplot(df_psam, aes(x = probe_combination, y = psam_amp, fill = probe_type)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  # Comparison 1: Early Onset, Altered vs Unaltered
  geom_signif(comparisons = list(c("Early Onset\nAltered", "Early Onset\nUnaltered")),
              annotations = "n.s.",
              y_position = y_pos_level1,
              tip_length = 0.02,
              textsize = 2.75,
              vjust = 0.1) +
  # Comparison 2: Late Onset, Altered vs Unaltered
  geom_signif(comparisons = list(c("Late Onset\nAltered", "Late Onset\nUnaltered")),
              annotations = "n.s.",
              y_position = y_pos_level1, 
              tip_length = 0.02,
              textsize = 2.75,
              vjust = 0.1) +
  scale_fill_manual(values = c("Altered" = colors$main_red, "Unaltered" = colors$main_blue)) +
  labs(
    #title = "PSAM Effect Amplitudes by Probe Type and Probe Onset Combinations",
    y = "PSAM Effect Amplitude [µV]",
    x = NULL,
    fill = "Probe Type"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 0, hjust = 1)
  )
P3

# Save plot
ggsave(
  filename = "tid_psam_psam_effect_violin.png", 
  plot = P3,
  width = 8,      
  height = 6,     
  dpi = 900,
  bg = "white"
)

# Plot 4: ERP Latencies by Probe Type, Task Instruction, and Probe Onset
P4 <- ggplot(df_erp_lat_withinSE, aes(x = task_instruction, y = erp_lat, color = probe_type, group = probe_type)) +
  geom_line(aes(y = erp_lat), size = 1.2) +
  geom_errorbar(aes(ymin = erp_lat - se, ymax = erp_lat + se), width = 0.2, size = 0.8) +
  geom_point(size = 3) +
  facet_wrap(~ probe_onset_cat, nrow = 1,
             labeller = as_labeller(c("Early" = "Early Probe Onset",
                                      "Late" = "Late Probe Onset"))) +
  scale_color_manual(values = c(
    "Altered" = colors$main_red,
    "Unaltered" = colors$main_blue
  )) +
  labs(
    #title = "N1 ERP Latencies by Probe Type, Task Instruction, and Probe Onset",
    y = "N1 ERP Latency [ms]",
    x = "Task Condition",
    color = "Probe Type"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "top"
  ) + 
  scale_x_discrete(labels = c("Active" = "Active",
                              "Passive" = "Passive"))
P4

# Save plot
ggsave(
  filename = "tid_psam_erp_lat.png", 
  plot = P4,
  width = 8,      
  height = 6,     
  dpi = 900,
  bg = "white"
)

# Combine Plot 2 and Plot 3
P2_P3 <- ggdraw() +
  draw_plot(P2, x = 0, y = 0, width = .5, height = 1) +
  draw_plot(P3, x = .5, y = 0, width = .5, height = 1) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0.5), y = c(1, 1))
P2_P3

ggsave(
  filename = "tid_psam_psam_effect_combined.png", 
  plot = P2_P3,
  width = 8,      
  height = 6,     
  dpi = 900,
  bg = "white"
)


