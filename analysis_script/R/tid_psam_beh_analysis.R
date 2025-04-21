# tid_psam_beh_analysis.R
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
df_probe_properties$probe_type <- as.factor(df_probe_properties$probe_type)

# Z-standardized vocal responses F0 values for probe and no-probe trials for each subject
df_probe_f0_z <- df_beh %>% 
  group_by(probe, subj) %>%
  summarise(mean(recording_f0_z, na.rm = T))  
colnames(df_probe_f0_z) <- c("probe", "subj", "recording_f0_z")
df_probe_f0_z$probe <- as.factor(df_probe_f0_z$probe)

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
df_probe_vot_z$probe <- as.factor(df_probe_vot_z$probe)

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
df_block_f0_z$block <- as.factor(df_block_f0_z$block)

# Z-standardized VOTs for each block and subject
df_block_vot_z <- df_beh %>% 
  group_by(block, subj) %>%
  summarise(mean(recording_vot_z, na.rm = T))  
colnames(df_block_vot_z) <- c("block", "subj", "recording_vot_z")
df_block_vot_z$block <- as.factor(df_block_vot_z$block)

#------------------------------------Analysis-----------------------------------

# BEH1
# Paired T-Test: (DV = Z-transformed probe F0 value, within = Probe-type)
# Analysis BEH1 concerns the F0 of the auditory probes (unaltered and altered) relative to the distribution of the F0 of the vocal responses during the experiment.
BEH1 <- t.test(data = df_probe_properties, f0_z ~ probe_type, paired = TRUE) 
BEH1

BEH1_ES <- as.data.frame(df_probe_properties) %>% 
  cohens_d(f0_z ~ probe_type, paired = TRUE)
BEH1_ES

# Plot: Probe F0 by probe-type
ggplot(df_probe_properties, aes(x = probe_type, y = f0_z, fill = probe_type)) +
  geom_boxplot(show.legend = T) +  
  theme_ggstatsplot() 

# Descriptive statistics
psych::describeBy(df_probe_properties$f0_z,
                  group = df_probe_properties$probe_type)

# Assupmtions 
# Normal distibution
df_probe_properties %>%
  ggplot(aes(x = f0_z)) +
  geom_histogram(bins = 50) +
  facet_wrap(df_probe_properties$probe_type) +
  theme_ggstatsplot()

byf.shapiro(f0_z ~ probe_type, 
            data = df_probe_properties)

#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
#------------------------------------------------------------------------------#

# BEH2
# Paired T-Test: (DV = Z-transformed F0 value, within = Probe)
# Analysis BEH2 concerns how the F0 of the vocal responses during the experiment is influenced by a probe being presented.
BEH2 <- t.test(data = df_probe_f0_z, recording_f0_z ~ probe, paired = TRUE) 
BEH2

BEH2_ES <- as.data.frame(df_probe_f0_z) %>% 
  cohens_d(recording_f0_z ~ probe, paired = TRUE) 
BEH2_ES

# Plot: F0 by probe
ggplot(df_probe_f0_z, aes(x = probe, y = recording_f0_z, fill = probe)) +
  geom_boxplot(show.legend = T) +  
  theme_ggstatsplot() 

# Descriptive statistics
psych::describeBy(df_probe_f0_z$recording_f0_z,
                  group = df_probe_f0_z$probe)

# Assupmtions 
# Normal distibution
df_probe_f0_z %>%
  ggplot(aes(x = recording_f0_z)) +
  geom_histogram(bins = 50) +
  facet_wrap(df_probe_f0_z$probe) +
  theme_ggstatsplot()

byf.shapiro(recording_f0_z ~ probe, 
            data = df_probe_f0_z)

#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
#------------------------------------------------------------------------------#

# BEH3
# rmANOVA (DV = Z-transformed F0 value, Within = Probe-type (unaltered, altered), Probe-onset (early, late)) 
# Analysis BEH3 concerns how the F0 of the vocal responses during the experiment is influenced by probe-type and probe-onset.
BEH3 <- aov_ez(id = "subj",
                dv = "recording_f0_z",
                data = df_probe_type_onset_f0_z,
                within = c("probe_type", "probe_onset_cat"),
                type = 3)
summary(BEH3)
##anova_test(data = df_probe_type_onset_f0_z, dv = recording_f0_z, wid = subj, within = c(probe_type,probe_onset_cat), type = 3, effect.size = "ges") # To get corrected dfs

# Plot: F0 by probe-onset and probe-type
ezPlot(
  data = df_probe_type_onset_f0_z 
  , dv = recording_f0_z 
  , wid = subj  
  , within= .(probe_type, probe_onset_cat)
  , x = .(probe_onset_cat)
  , split   = .(probe_type)
)

# Descriptive statistics
psych::describeBy(
  df_probe_type_onset_f0_z$recording_f0_z,
  list(df_probe_type_onset_f0_z$probe_type, df_probe_type_onset_f0_z$probe_onset_cat)
)

# Follow-Up T-Tests
BEH3_EM <- emmeans(BEH3, ~ probe_type | probe_onset_cat)
BEH3_PWC <- pairs(BEH3_EM, adjust = "bonferroni")
BEH3_PWC

# Assupmtions 
# Sphericity
# Not applicable due to only 2 factor levels per factor

# Normal distibution
df_probe_type_onset_f0_z %>%
  ggplot(aes(x = recording_f0_z)) +
  geom_histogram(bins = 50) +
  facet_grid(probe_type ~ probe_onset_cat) +
  theme_ggstatsplot()

byf.shapiro(recording_f0_z ~ probe_type * probe_onset_cat, 
            data = df_probe_type_onset_f0_z)

#Balance of the design
ezDesign(df_probe_type_onset_f0_z, x = probe_type, y = subj, row = probe_onset_cat) 

#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
# - Sphericity: Not applicable due to only 2 factor levels per factor
# - Balance of the design: OK
#------------------------------------------------------------------------------#

# BEH4
# Paired T-Test: (DV = Z-transformed vocal onset time, within = Probe)
# Analysis BEH4 concerns how the vocal onset time is influenced by a probe being presented.
BEH4 <- t.test(data = df_probe_vot_z, recording_vot_z ~ probe, paired = TRUE) 
BEH4

BEH4_ES <- as.data.frame(df_probe_vot_z) %>% 
  cohens_d(recording_vot_z ~ probe, paired = TRUE) 
BEH4_ES

# Plot: vot by probe
ggplot(df_probe_vot_z, aes(x = probe, y = recording_vot_z, fill = probe)) +
  geom_boxplot(show.legend = T) +  
  theme_ggstatsplot() 

# Descriptive statistics
psych::describeBy(df_probe_vot_z$recording_vot_z,
                  group = df_probe_vot_z$probe)

# Assupmtions 
# Normal distibution
df_probe_vot_z %>%
  ggplot(aes(x = recording_vot_z)) +
  geom_histogram(bins = 50) +
  facet_wrap(df_probe_vot_z$probe) +
  theme_ggstatsplot()

byf.shapiro(recording_vot_z ~ probe, 
            data = df_probe_vot_z)

#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
#------------------------------------------------------------------------------#

# BEH5
# rmANOVA (DV = Z-transformed vocal onset time, Within = Probe-type (unaltered, altered), Probe-onset (early, late)) 
# Analysis BEH5 concerns how the vocal onset time is influenced by probe-type and probe-onset.
BEH5 <- aov_ez(id = "subj",
               dv = "recording_vot_z",
               data = df_probe_type_onset_vot_z,
               within = c("probe_type", "probe_onset_cat"),
               type = 3)
summary(BEH5)


# Plot: Vocal onset time by probe-onset and probe-type
ezPlot(
  data = df_probe_type_onset_vot_z 
  , dv = recording_vot_z 
  , wid = subj  
  , within= .(probe_type, probe_onset_cat)
  , x = .(probe_onset_cat)
  , split   = .(probe_type)
)

# Descriptive statistics
psych::describeBy(
  df_probe_type_onset_f0_z$recording_f0_z,
  list(df_probe_type_onset_vot_z$probe_type, df_probe_type_onset_vot_z$probe_onset_cat)
)

# Follow-Up T-Tests
BEH5_EM <- emmeans(BEH5, ~ probe_onset_cat| probe_type)
BEH5_PWC <- pairs(BEH5_EM, adjust = "bonferroni")
BEH5_PWC

# Assupmtions 
# Sphericity
# Not applicable due to only 2 factor levels per factor

# Normal distibution
df_probe_type_onset_vot_z %>%
  ggplot(aes(x = recording_vot_z)) +
  geom_histogram(bins = 50) +
  facet_grid(probe_type ~ probe_onset_cat) +
  theme_ggstatsplot()

byf.shapiro(recording_vot_z ~ probe_type * probe_onset_cat, 
            data = df_probe_type_onset_vot_z)

#Balance of the design
ezDesign(df_probe_type_onset_vot_z, x = probe_type, y = subj, row = probe_onset_cat) 

#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
# - Sphericity: Not applicable due to only 2 factor levels per factor
# - Balance of the design: OK
#------------------------------------------------------------------------------#

# BEH6
# rmANOVA (DV = Z-transformed F0 value, Within = Block (1,2,3,4,5,6,7,8)) 
# Analysis BEH6 concerns how the F0 of the vocal responses during the experiment is influenced by the block of the experiment.
BEH6 <- aov_ez(id = "subj",
               dv = "recording_f0_z",
               data = df_block_f0_z,
               within = c("block"),
               type = 3)
summary(BEH6)


# Plot: Vocal onset time by probe-onset and probe-type
ezPlot(
  data = df_block_f0_z 
  , dv = recording_f0_z 
  , wid = subj  
  , within= .(block)
  , x = .(block)
  , type = 1
)

# Descriptive statistics
psych::describeBy(
  df_block_f0_z$recording_f0_z,
  list(df_block_f0_z$block)
)

# Follow-Up T-Tests
BEH6_EM <- emmeans(BEH6, ~ block)
BEH6_PWC <- pairs(BEH6_EM, adjust = "bonferroni")
BEH6_PWC

# Assupmtions 
# Sphericity
# Checked in rmANOVA. Correction applied if needed.

# Normal distibution
df_block_f0_z %>%
  ggplot(aes(x = recording_f0_z)) +
  geom_histogram(bins = 50) +
  facet_wrap(df_block_f0_z$block) +
  theme_ggstatsplot()

byf.shapiro(recording_f0_z ~ block, 
            data = df_block_f0_z)

#Balance of the design
ezDesign(df_block_f0_z, x = block, y = subjt) 

#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
# - Sphericity: Not applicable due to only 2 factor levels per factor
# - Balance of the design: OK
#------------------------------------------------------------------------------#

# BEH7
# rmANOVA (DV = Z-transformed vocal onset time, Within = Block (1,2,3,4,5,6,7,8)) 
# Analysis BEH7 concerns how the F0 of the vocal onset time is influenced by the block of the experiment.
BEH7 <- aov_ez(id = "subj",
               dv = "recording_vot_z",
               data = df_block_vot_z,
               within = c("block"),
               type = 3)
summary(BEH7)

# Plot: Vocal onset time by probe-onset and probe-type
ezPlot(
  data = df_block_vot_z 
  , dv = recording_vot_z 
  , wid = subj  
  , within= .(block)
  , x = .(block)
  , type = 1
)

# Descriptive statistics
psych::describeBy(
  df_block_vot_z$recording_vot_z,
  list(df_block_vot_z$block)
)

# Follow-Up T-Tests
BEH7_EM <- emmeans(BEH7, ~ block)
BEH7_PWC <- pairs(BEH7_EM, adjust = "bonferroni")
BEH7_PWC

# Assupmtions 
# Sphericity
# Checked in rmANOVA. Correction applied if needed.

# Normal distibution
df_block_vot_z %>%
  ggplot(aes(x = recording_vot_z)) +
  geom_histogram(bins = 50) +
  facet_wrap(df_block_vot_z$block) +
  theme_ggstatsplot()

byf.shapiro(recording_vot_z ~ block, 
            data = df_block_vot_z)

#Balance of the design
ezDesign(df_block_vot_z, x = block, y = subj) 

#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
# - Sphericity: Not applicable due to only 2 factor levels per factor
# - Balance of the design: OK
#------------------------------------------------------------------------------#

# BEH8

# BEH9
