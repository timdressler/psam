# tid_psam_beh_analysis.R
#
# Performs analysis of behavioural data.
#
# Tim Dressler, 18.04.25

#-------------------------------------Set up------------------------------------
#load packages
library(tidyr) 
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
INPATH <- file.path(MAINPATH, "data", "processed_data", "beh_preprocessed_clean")
OUTPATH <- file.path(MAINPATH, "data", "analysis_data", "stats_beh_analysis")

if (!dir.exists(OUTPATH)) {
  dir.create(OUTPATH, recursive = TRUE)
}

setwd(OUTPATH)

#-----------------------------------Load data-----------------------------------

# Load data
df_beh <- read_excel(file.path(INPATH, "all_subj_beh_preprocessed_clean.xlsx"))


#-------------------------------Create needed dfs-------------------------------

# Z-standardized probe F0 values for unaltered and altered probes for each subject
df_probe_properties <- df_beh %>%
  distinct(subj, probe_f0_unaltered_z, probe_f0_altered_z)
df_probe_properties <- df_probe_properties %>%
  pivot_longer(
    cols = c("probe_f0_unaltered_z", "probe_f0_altered_z"),
    names_to = "probe_type",
    values_to = "f0_z"
  )

# Z-standardized vocal responses F0 values for probe and no-probe trials for each subject
df_probe_f0_z <- df_beh %>% 
  group_by(probe, subj) %>%
  summarise(mean(recording_f0_z, na.rm = T))  
colnames(df_probe_f0_z) <- c("probe", "subj", "recording_f0_z")

# Z-standardized vocal responses F0 values for probe and no-probe trials for each subject
df_probe_type_onset_f0_z <- filter(df_beh, probe == "Yes") %>% 
  group_by(probe_type, probe_onset_cat, subj) %>%
  summarise(mean(recording_f0_z, na.rm = T))  
colnames(df_probe_type_onset_f0_z) <- c("probe_type", "probe_onset_cat", "subj", "recording_f0_z")

# Z-standardized VOTs for probe and no-probe trials for each subject
df_probe_vot_z <- df_beh %>% 
  group_by(probe, subj) %>%
  summarise(mean(recording_vot_z, na.rm = T))  
colnames(df_probe_vot_z) <- c("probe", "subj", "recording_vot_z")

# Z-standardized VOTs for probe and no-probe trials for each subject
df_probe_type_onset_vot_z <- filter(df_beh, probe == "Yes") %>% 
  group_by(probe_type, probe_onset_cat, subj) %>%
  summarise(mean(recording_vot_z, na.rm = T))  
colnames(df_probe_type_onset_vot_z) <- c("probe_type", "probe_onset_cat", "subj", "recording_vot_z")

# Z-standardized vocal responses F0 values for each block and subject
df_block_f0_z <- df_beh %>% 
  group_by(block, subj) %>%
  summarise(mean(recording_f0_z, na.rm = T))  
colnames(df_block_f0_z) <- c("block", "subj", "recording_f0_z")

# Z-standardized VOTs for each block and subject
df_block_vot_z <- df_beh %>% 
  group_by(block, subj) %>%
  summarise(mean(recording_vot_z, na.rm = T))  
colnames(df_block_vot_z) <- c("block", "subj", "recording_vot_z")

#------------------------------------Analysis-----------------------------------

# BEH1

# BEH2

# BEH3
# rmANOVA (DV = Z-transformed F0 value, Within = Probe-type (unaltered, altered), Probe-onset (early, late)) 
# Analysis BEH3 concerns how the F0 of the vocal responses during the experiment is influenced by probe-type and probe-onset.

ezANOVA(df_probe_type_onset_f0_z,dv=.(recording_f0_z),within=.(probe_type, probe_onset_cat),wid=.(subj), type = 3) 

#PLOT: Hitrate by cond group 
ezPlot(
  data = df_probe_type_onset_f0_z 
  , dv = recording_f0_z 
  , wid = subj  
  , within= .(probe_type, probe_onset_cat)
  , x = .(probe_type)
)

#Follow-Up T-Tests
pairwise.t.test(df_cat$Hitrate, df_cat$cond_group, p.adjust.method = "bonf", paired = TRUE) 

pwc_PREB1 <- df_cat %>% 
  pairwise_t_test(Hitrate ~ cond_group, paired = T, p.adjust.method = "bonferroni") 
pwc_PREB1

#PLOT: Hitrate by cond group 
pwc_PREB1 <- pwc_PREB1 %>%
  add_xy_position(x = "cond_group")
ggboxplot(df_cat, x = "cond_group", y = "Hitrate", add = "point") +
  stat_pvalue_manual(pwc_PREB1)


#Assupmtions 
#Sphericity: checked in rmANOVA 

#Normal distibution:
psych::describeBy(df_probe_type_onset_f0_z$recording_f0_z, df_probe_type_onset_f0_z$probe_onset_cat)

psych::describeBy(
  df_probe_type_onset_f0_z$recording_f0_z,
  list(df_probe_type_onset_f0_z$probe_type, df_probe_type_onset_f0_z$probe_onset_cat)
)

df_probe_type_onset_f0_z %>%
  ggplot(aes(x=Hitrate)) +
  geom_histogram(bins = 50) +
  facet_wrap(df_probe_type_onset_f0_z$cond_group) +
  theme_ggstatsplot() 

byf.shapiro(df_probe_type_onset_f0_z$recording_f0_z~df_probe_type_onset_f0_z$probe_onset_cat, data = df_probe_type_onset_f0_z) 

#Balance of the design
ezDesign(df_probe_type_onset_f0_z, x = cond_group, y = subj) 

#------------------------------------------------------------------------------#
#significant differences between categories on Hitrate
#all comparisons sig.

#assumptions
# - normal distribution: not OK 
# - Sphericity: not OK, corrected with HF (GGe > .75)
# - Balance of the design: OK
#------------------------------------------------------------------------------#



