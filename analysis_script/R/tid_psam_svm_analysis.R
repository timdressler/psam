# tid_psam_beh_analysis.R
#
# Performs analysis of SVM data.
#
# Tim Dressler, 18.04.25

#-------------------------------------Set up------------------------------------
# Load packages
library(tidyr) 
library(afex)
library(emmeans)
library(readxl)
library(car) 
library(corrplot) 
library(dplyr) 
library(ez)
library(ggplot2)
library(ggstatsplot)
library(DescTools)
library(ggpubr)
library(ggeffects)
library(cowplot) 
library(tidyverse)
library(psych)
library(rstatix)
library(RVAideMemoire)
library(effectsize)
library(smplot2)

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
if (grepl("psam/analysis_script/R", SCRIPTPATH)) {
  cat("Path OK") 
} else {
  stop("Path not OK")  
}

MAINPATH <- gsub("/analysis_script/R", "", SCRIPTPATH)
INPATH <- file.path(MAINPATH, "data", "analysis_data", "svm_analysis")
OUTPATH <- file.path(MAINPATH, "data", "analysis_data", "stats_svm_analysis")

FUNPATH <- file.path(MAINPATH, "functions")
source(file.path(FUNPATH, "tid_psam_check_folder_TD.R"))
source(file.path(FUNPATH, "tid_psam_clean_up_folder_TD.R"))

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

setwd(OUTPATH)

#-----------------------------------Load data-----------------------------------

# Load data
df_svm_wide <- read_excel(file.path(INPATH, "all_subj_accuracy_proportions.xlsx"))

df_svm <- df_svm_wide %>% # Same df in long format 
  pivot_longer(cols = starts_with("prop"), names_to = "feature_extraction_window", values_to = "percent_above_chance")
df_svm$feature_extraction_window <- as.factor(df_svm$feature_extraction_window)

#------------------------------------Analysis-----------------------------------

# MAIN_SVM1
# Paired T-Test (DV = Percentage of above chance level performing hyperparameter pairs, within = Temporal window (early, late))
# Analysis MAIN_SVM1 concerns how the percentage of above chance level performing hyperparameter pairs is influenced by the time window used for extracting the features. 
MAIN_SVM1 <- t.test(
  df_svm_wide$prop_sig_early,
  df_svm_wide$prop_sig_late,
  paired = TRUE
)
MAIN_SVM1

MAIN_SVM1_ES <- as.data.frame(df_svm) %>% 
  rstatix::cohens_d(percent_above_chance ~ feature_extraction_window, paired = TRUE)
MAIN_SVM1_ES

# Plot: Probe F0 by probe-type
ggplot(df_svm, aes(x = feature_extraction_window, y = percent_above_chance, fill = feature_extraction_window)) +
  geom_boxplot(show.legend = T) +  
  theme_ggstatsplot() 

# Descriptive statistics
psych::describeBy(df_svm$percent_above_chance,
                  group = df_svm$feature_extraction_window)

# Assumptions 
# Normal distibution
df_svm %>%
  ggplot(aes(x = percent_above_chance)) +
  geom_histogram(bins = 50) +
  facet_wrap(df_svm$feature_extraction_window) +
  theme_ggstatsplot()

byf.shapiro(percent_above_chance ~ feature_extraction_window, 
            data = df_svm)

#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
#------------------------------------------------------------------------------#

# MAIN_SVM1_ALT
# Wilcoxon Signed-Rank Test (DV = Percentage of above chance level performing hyperparameter pairs, within = Temporal window (early, late))
# Analysis MAIN_SVM1_ALT concerns how the percentage of above chance level performing hyperparameter pairs is influenced by the time window used for extracting the features. 
# Since the assumptions were not met for MAIN_SVM1, a non-parametric alternative is used.
MAIN_SVM1_ALT <- wilcox.test(df_svm_wide$prop_sig_early, df_svm_wide$prop_sig_late, paired = TRUE)
MAIN_SVM1_ALT

MAIN_SVM1_ES_ALT <- effectsize::rank_biserial(df_svm_wide$prop_sig_early, df_svm_wide$prop_sig_late, paired = TRUE)
MAIN_SVM1_ES_ALT

# Plot: Probe F0 by probe-type
ggplot(df_svm, aes(x = feature_extraction_window, y = percent_above_chance, fill = feature_extraction_window)) +
  geom_boxplot(show.legend = T) +  
  theme_ggstatsplot() 

# Descriptive statistics
psych::describeBy(df_svm$percent_above_chance,
                  group = df_svm$feature_extraction_window)

#------------------------------------------------------------------------------#
#
#
#------------------------------------------------------------------------------#

# SVM2 # SVM3
# Descriptive only
# Analysis SVM2 concerns the chance-level thresholds for each subject.
# Analysis SVM3 concerns the number of included trials per subject.

# Descriptive statistics
describe(df_svm_wide)

#------------------------------------------------------------------------------#
#
#
#------------------------------------------------------------------------------#

#-------------------------------------Plots-------------------------------------

# Plot 1: 
P1 <- df_svm %>%
  ggplot(aes(x = feature_extraction_window, y = percent_above_chance*100, fill = feature_extraction_window)) + # *100 to convert to percent
  sm_raincloud(boxplot.params = list(fill = "white", outlier.shape = NA), violin.params = list(alpha = .6)
               , point.params = list(), legends = F) +
  scale_fill_manual(values = c(colors$main_yellow, colors$main_blue)) +
  scale_x_discrete(labels = c('prop_sig_early' = 'Early Window', 'prop_sig_late' = 'Late Window')) +
  scale_y_continuous(n.breaks = 3) +
  labs(x = NULL, y = "Percentage of possible Hyperparameter-Pairs leading to an above-chance Classification") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.position="none") +
  geom_signif(comparisons=list(c("prop_sig_early", "prop_sig_late")), annotations="n.s.",
              y_position = 100, tip_length = 0.02,  vjust=0.4) 
P1

# Save plot
ggsave(
  filename = "tid_psam_hyperparamters_violin.png", 
  plot = P1,
  width = 8,      
  height = 6,     
  dpi = 300,
  bg = "white"
)

