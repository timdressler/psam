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
library(rstudioapi)
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
SCRIPTPATH <- dirname(rstudioapi::getSourceEditorContext()$path)
if (grepl("psam/analysis_script/R", SCRIPTPATH)) {
  cat("Path OK\n") 
} else {
  stop("Path not OK")  
}

MAINPATH <- gsub("/analysis_script/R", "", SCRIPTPATH)
INPATH <- file.path(MAINPATH, "data", "analysis_data", "erp_analysis")
OUTPATH <- file.path(MAINPATH, "data", "analysis_data", "stats_erp_analysis")

if (!dir.exists(OUTPATH)) {
  dir.create(OUTPATH, recursive = TRUE)
}

setwd(OUTPATH)

#-----------------------------------Load data-----------------------------------

# Load data
df_erp <- read_excel(file.path(INPATH, "all_subj_erp_data.xlsx"))


# MAIN_ERP1
# Linear Mixed Model (Random Intercepts, Fixed Slopes) (DV = N1 Amplitude, within = Task (active, passive), Probe-type (unaltered, altered), Probe-onset (early, late))
# Analysis MAIN_ERP1 concerns how N1 ERP amplitudes are influenced by probe-type, probe-onset and task.
MAIN_ERP1 <- lme4::lmer(erp_amp ~ task_instruction + (1|subj), data = df_erp)
summary(MAIN_ERP1)

# Assumptions 


#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
#------------------------------------------------------------------------------#

















