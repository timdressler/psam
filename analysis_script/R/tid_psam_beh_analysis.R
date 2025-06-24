# tid_psam_beh_analysis.R
#
# Performs analysis of behavioural data.
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
INPATH <- file.path(MAINPATH, "data", "processed_data", "beh_preprocessed_clean")
INPATH_FAL <- file.path(MAINPATH, "data", "questionnaire_data_clean")
OUTPATH <- file.path(MAINPATH, "data", "analysis_data", "stats_beh_analysis")

FUNPATH <- file.path(MAINPATH, "functions")
source(file.path(FUNPATH, "tid_psam_check_folder_TD.R"))
source(file.path(FUNPATH, "tid_psam_clean_up_folder_TD.R"))

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH, INPATH_FAL)
tid_psam_clean_up_folder_TD(OUTPATH)

setwd(OUTPATH)

#-----------------------------------Load data-----------------------------------

# Load data
df_beh <- read_excel(file.path(INPATH, "all_subj_beh_preprocessed_clean.xlsx"))

# Load questionnaire data
df_fal <- read_excel(file.path(INPATH_FAL, "fal_data_clean.xlsx"))


#-------------------------------Create needed dfs-------------------------------

# Z-standardized probe F0 values for unaltered and altered probes for each subject
df_probe_properties_z <- df_beh %>%
  distinct(subj, probe_f0_unaltered_z, probe_f0_altered_z)
df_probe_properties_z <- df_probe_properties_z %>%
  pivot_longer(
    cols = c("probe_f0_unaltered_z", "probe_f0_altered_z"),
    names_to = "probe_type",
    values_to = "f0_z"
  )
df_probe_properties_z$probe_type <- as.factor(df_probe_properties_z$probe_type)

df_probe_properties_z_wide <- df_probe_properties_z %>% # Same df in wide format (for t.test)
  pivot_wider(names_from = probe_type, values_from = f0_z)

# Z-standardized vocal responses F0 values for probe and no-probe trials for each subject 
df_probe_f0_z <- df_beh %>% 
  group_by(probe, subj) %>%
  summarise(mean(recording_f0_z, na.rm = T))  
colnames(df_probe_f0_z) <- c("probe", "subj", "recording_f0_z")
df_probe_f0_z$probe <- as.factor(df_probe_f0_z$probe)

df_probe_f0_z_wide <- df_probe_f0_z %>% # Same df in wide format (for t.test)
  pivot_wider(names_from = probe, values_from = recording_f0_z)

# Vocal responses F0 values for probe and no-probe trials for each subject
df_probe_f0 <- df_beh %>% 
  group_by(probe, subj) %>%
  summarise(mean(recording_f0, na.rm = T))  
colnames(df_probe_f0) <- c("probe", "subj", "recording_f0")
df_probe_f0$probe <- as.factor(df_probe_f0$probe)

df_probe_f0_wide <- df_probe_f0 %>% # Same df in wide format (for t.test)
  pivot_wider(names_from = probe, values_from = recording_f0)

# Z-standardized vocal responses F0 values for probe and no-probe trials for each subject
df_probe_type_onset_f0_z <- filter(df_beh, probe == "Yes") %>% 
  group_by(probe_type, probe_onset_cat, subj) %>%
  summarise(mean(recording_f0_z, na.rm = T))  
colnames(df_probe_type_onset_f0_z) <- c("probe_type", "probe_onset_cat", "subj", "recording_f0")

# Vocal responses F0 values for probe and no-probe trials for each subject
df_probe_type_onset_f0 <- filter(df_beh, probe == "Yes") %>% 
  group_by(probe_type, probe_onset_cat, subj) %>%
  summarise(mean(recording_f0, na.rm = T))  
colnames(df_probe_type_onset_f0) <- c("probe_type", "probe_onset_cat", "subj", "recording_f0")

# Z-standardized VOTs for probe and no-probe trials for each subject
df_probe_vot_z <- df_beh %>% 
  group_by(probe, subj) %>%
  summarise(mean(recording_vot_z, na.rm = T))  
colnames(df_probe_vot_z) <- c("probe", "subj", "recording_vot_z")
df_probe_vot_z$probe <- as.factor(df_probe_vot_z$probe)

df_probe_vot_z_wide <- df_probe_vot_z %>% # Same df in wide format (for t.test)
  pivot_wider(names_from = probe, values_from = recording_vot_z)

# VOTs for probe and no-probe trials for each subject
df_probe_vot <- df_beh %>% 
  group_by(probe, subj) %>%
  summarise(mean(recording_vot, na.rm = T))  
colnames(df_probe_vot) <- c("probe", "subj", "recording_vot")
df_probe_vot$probe <- as.factor(df_probe_vot$probe)

df_probe_vot_wide <- df_probe_vot %>% # Same df in wide format (for t.test)
  pivot_wider(names_from = probe, values_from = recording_vot)

# Z-standardized VOTs for probe and no-probe trials for each subject
df_probe_type_onset_vot_z <- filter(df_beh, probe == "Yes") %>% 
  group_by(probe_type, probe_onset_cat, subj) %>%
  summarise(mean(recording_vot_z, na.rm = T))  
colnames(df_probe_type_onset_vot_z) <- c("probe_type", "probe_onset_cat", "subj", "recording_vot_z")

# VOTs for probe and no-probe trials for each subject
df_probe_type_onset_vot <- filter(df_beh, probe == "Yes") %>% 
  group_by(probe_type, probe_onset_cat, subj) %>%
  summarise(mean(recording_vot, na.rm = T))  
colnames(df_probe_type_onset_vot) <- c("probe_type", "probe_onset_cat", "subj", "recording_vot")

# Z-standardized vocal responses F0 values for each block and subject
df_block_f0_z <- df_beh %>% 
  group_by(block, subj) %>%
  summarise(mean(recording_f0_z, na.rm = T))  
colnames(df_block_f0_z) <- c("block", "subj", "recording_f0_z")
df_block_f0_z$block <- as.factor(df_block_f0_z$block)

# Vocal responses F0 values for each block and subject
df_block_f0 <- df_beh %>% 
  group_by(block, subj) %>%
  summarise(mean(recording_f0, na.rm = T))  
colnames(df_block_f0) <- c("block", "subj", "recording_f0")
df_block_f0$block <- as.factor(df_block_f0$block)

# Z-standardized VOTs for each block and subject
df_block_vot_z <- df_beh %>% 
  group_by(block, subj) %>%
  summarise(mean(recording_vot_z, na.rm = T))  
colnames(df_block_vot_z) <- c("block", "subj", "recording_vot_z")
df_block_vot_z$block <- as.factor(df_block_vot_z$block)

# VOTs for each block and subject
df_block_vot <- df_beh %>% 
  group_by(block, subj) %>%
  summarise(mean(recording_vot, na.rm = T))  
colnames(df_block_vot) <- c("block", "subj", "recording_vot")
df_block_vot$block <- as.factor(df_block_vot$block)

# Vocal responses F0 values for each sex
df_subj_f0 <- filter(df_beh, probe == "Yes") %>% 
  group_by(subj) %>%
  summarise(mean(recording_f0, na.rm = T))  
colnames(df_subj_f0) <- c("subj", "recording_f0")
df_subj_f0 <- merge(df_fal, df_subj_f0, by = "subj", all.x = TRUE)
df_subj_f0 <- df_subj_f0[c("subj", "var3_sex", "recording_f0")]
df_subj_f0$var3_sex <- as.factor(df_subj_f0$var3_sex)
df_subj_f0$var3_sex <- car::recode(df_subj_f0$var3_sex, "1='male';2='female'; 3 = 'diverse'")



#------------------------------------Analysis-----------------------------------

# BEH1
# Paired T-Test (DV = Z-transformed probe F0 value, within = Probe-type)
# Analysis BEH1 concerns the F0 of the auditory probes (unaltered and altered) relative to the distribution of the F0 of the vocal responses during the experiment.
BEH1 <- t.test(
  df_probe_properties_z_wide$probe_f0_altered_z,
  df_probe_properties_z_wide$probe_f0_unaltered_z,
  paired = TRUE
)
BEH1

BEH1_ES <- as.data.frame(df_probe_properties_z) %>% 
  cohens_d(f0_z ~ probe_type, paired = TRUE)
BEH1_ES

# Plot: Probe F0 by probe-type
ggplot(df_probe_properties_z, aes(x = probe_type, y = f0_z, fill = probe_type)) +
  geom_boxplot(show.legend = T) +  
  theme_ggstatsplot() 

# Descriptive statistics
psych::describeBy(df_probe_properties_z$f0_z,
                  group = df_probe_properties_z$probe_type)

# Assumptions 
# Normal distibution
df_probe_properties_z %>%
  ggplot(aes(x = f0_z)) +
  geom_histogram(bins = 50) +
  facet_wrap(df_probe_properties_z$probe_type) +
  theme_ggstatsplot()

byf.shapiro(f0_z ~ probe_type, 
            data = df_probe_properties_z)

#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
#------------------------------------------------------------------------------#

# BEH2
# Paired T-Test (DV = F0 value, within = Probe (Yes, No))
# Analysis BEH2 concerns how the F0 of the vocal responses during the experiment is influenced by a probe being presented.
BEH2 <- t.test(
  df_probe_f0_wide$No,
  df_probe_f0_wide$Yes,
  paired = TRUE
) 
BEH2

BEH2_ES <- as.data.frame(df_probe_f0) %>% 
  cohens_d(recording_f0 ~ probe, paired = TRUE) 
BEH2_ES

# Plot: F0 by probe
ggplot(df_probe_f0, aes(x = probe, y = recording_f0, fill = probe)) +
  geom_boxplot(show.legend = T) +  
  theme_ggstatsplot() 

# Descriptive statistics
psych::describeBy(df_probe_f0$recording_f0,
                  group = df_probe_f0$probe)

# Assumptions 
# Normal distibution
df_probe_f0 %>%
  ggplot(aes(x = recording_f0)) +
  geom_histogram(bins = 50) +
  facet_wrap(df_probe_f0$probe) +
  theme_ggstatsplot()

byf.shapiro(recording_f0 ~ probe, 
            data = df_probe_f0)

#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
#------------------------------------------------------------------------------#


# BEH2_ALT
# Wilcoxon Signed-Rank Test (DV = F0 value, within = Probe (Yes, No))
# Analysis BEH2_ALT concerns how the F0 of the vocal responses during the experiment is influenced by a probe being presented.
# Since the assumptions were not met for BEH2, a non-parametric alternative is used.
BEH2_ALT <- wilcox.test(df_probe_f0_wide$No,
                             df_probe_f0_wide$Yes,
                             paired = TRUE)
BEH2_ALT

BEH2_ES_ALT <- effectsize::rank_biserial(df_probe_f0_wide$No,
                                         df_probe_f0_wide$Yes,
                                         paired = TRUE)
BEH2_ES_ALT

# Plot: F0 by probe
ggplot(df_probe_f0, aes(x = probe, y = recording_f0, fill = probe)) +
  geom_boxplot(show.legend = T) +  
  theme_ggstatsplot() 

# Descriptive statistics
psych::describeBy(df_probe_f0$recording_f0,
                  group = df_probe_f0$probe)

#------------------------------------------------------------------------------#
#
#
#------------------------------------------------------------------------------#

# BEH3
# rmANOVA (DV = F0 value, Within = Probe-type (unaltered, altered), Probe-onset (early, late)) 
# Analysis BEH3 concerns how the F0 of the vocal responses during the experiment is influenced by probe-type and probe-onset.
BEH3 <- aov_ez(id = "subj",
                dv = "recording_f0",
                data = df_probe_type_onset_f0,
                within = c("probe_type", "probe_onset_cat"),
                type = 3)
summary(BEH3)


BEH3_ES <- BEH3$anova_table
BEH3_ES

# Plot: F0 by probe-onset and probe-type
ezPlot(
  data = df_probe_type_onset_f0 
  , dv = recording_f0 
  , wid = subj  
  , within= .(probe_type, probe_onset_cat)
  , x = .(probe_onset_cat)
  , split   = .(probe_type)
)

# Descriptive statistics
psych::describeBy(
  df_probe_type_onset_f0$recording_f0,
  list(df_probe_type_onset_f0$probe_type, df_probe_type_onset_f0$probe_onset_cat)
)

# Follow-Up T-Tests
BEH3_EM <- emmeans(BEH3, ~ probe_type | probe_onset_cat)
BEH3_PWC <- pairs(BEH3_EM, adjust = "bonferroni")
BEH3_PWC

# Assumptions 
# Sphericity
# Not applicable due to only 2 factor levels per factor

# Normal distibution
df_probe_type_onset_f0 %>%
  ggplot(aes(x = recording_f0)) +
  geom_histogram(bins = 50) +
  facet_grid(probe_type ~ probe_onset_cat) +
  theme_ggstatsplot()

byf.shapiro(recording_f0 ~ probe_type * probe_onset_cat, 
            data = df_probe_type_onset_f0)

#Balance of the design
ezDesign(df_probe_type_onset_f0, x = probe_type, y = subj, row = probe_onset_cat) 

#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
# - Sphericity: Not applicable due to only 2 factor levels per factor
# - Balance of the design: OK
#------------------------------------------------------------------------------#

# BEH4
# Paired T-Test (DV = Vocal onset time, within = Probe (Yes, No))
# Analysis BEH4 concerns how the vocal onset time is influenced by a probe being presented.
BEH4 <- t.test(
  df_probe_vot_wide$No,
  df_probe_vot_wide$Yes,
  paired = TRUE
)
BEH4

BEH4_ES <- as.data.frame(df_probe_vot) %>% 
  cohens_d(recording_vot ~ probe, paired = TRUE) 
BEH4_ES

# Plot: vot by probe
ggplot(df_probe_vot, aes(x = probe, y = recording_vot, fill = probe)) +
  geom_boxplot(show.legend = T) +  
  theme_ggstatsplot() 

# Descriptive statistics
psych::describeBy(df_probe_vot$recording_vot,
                  group = df_probe_vot$probe)

# Assumptions 
# Normal distibution
df_probe_vot %>%
  ggplot(aes(x = recording_vot)) +
  geom_histogram(bins = 50) +
  facet_wrap(df_probe_vot$probe) +
  theme_ggstatsplot()

byf.shapiro(recording_vot ~ probe, 
            data = df_probe_vot)

#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
#------------------------------------------------------------------------------#

# BEH4_ALT
# Wilcoxon Signed-Rank Test (DV = Vocal onset time, within = Probe (Yes, No))
# Analysis BEH4_ALT concerns how the vocal onset time is influenced by a probe being presented.
# Since the assumptions were not met for BEH4, a non-parametric alternative is used.
BEH4_ALT <- wilcox.test(df_probe_vot_wide$No,
                        df_probe_vot_wide$Yes,
                        paired = TRUE)
BEH4_ALT

BEH4_ES_ALT <- effectsize::rank_biserial(df_probe_vot_wide$No,
                                         df_probe_vot_wide$Yes,
                                         paired = TRUE)
BEH4_ES_ALT

# Plot: vot by probe
ggplot(df_probe_vot, aes(x = probe, y = recording_vot, fill = probe)) +
  geom_boxplot(show.legend = T) +  
  theme_ggstatsplot() 

# Descriptive statistics
psych::describeBy(df_probe_vot$recording_vot,
                  group = df_probe_vot$probe)

#------------------------------------------------------------------------------#
#
#
#------------------------------------------------------------------------------#

# BEH5
# rmANOVA (DV = Vocal onset time, Within = Probe-type (unaltered, altered), Probe-onset (early, late)) 
# Analysis BEH5 concerns how the vocal onset time is influenced by probe-type and probe-onset.
BEH5 <- aov_ez(id = "subj",
               dv = "recording_vot",
               data = df_probe_type_onset_vot,
               within = c("probe_type", "probe_onset_cat"),
               type = 3)
summary(BEH5)

BEH5_ES <- BEH5$anova_table
BEH5_ES


# Plot: Vocal onset time by probe-onset and probe-type
ezPlot(
  data = df_probe_type_onset_vot 
  , dv = recording_vot 
  , wid = subj  
  , within= .(probe_type, probe_onset_cat)
  , x = .(probe_onset_cat)
  , split   = .(probe_type)
)

# Descriptive statistics
psych::describeBy(
  df_probe_type_onset_vot$recording_vot,
  list(df_probe_type_onset_vot$probe_type, df_probe_type_onset_vot$probe_onset_cat)
)

# Follow-Up T-Tests
BEH5_EM <- emmeans(BEH5, ~ probe_onset_cat| probe_type)
BEH5_PWC <- pairs(BEH5_EM, adjust = "bonferroni")
BEH5_PWC

# Assumptions 
# Sphericity
# Not applicable due to only 2 factor levels per factor

# Normal distibution
df_probe_type_onset_vot %>%
  ggplot(aes(x = recording_vot)) +
  geom_histogram(bins = 50) +
  facet_grid(probe_type ~ probe_onset_cat) +
  theme_ggstatsplot()

byf.shapiro(recording_vot ~ probe_type * probe_onset_cat, 
            data = df_probe_type_onset_vot)

#Balance of the design
ezDesign(df_probe_type_onset_vot, x = probe_type, y = subj, row = probe_onset_cat) 

#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
# - Sphericity: Not applicable due to only 2 factor levels per factor
# - Balance of the design: OK
#------------------------------------------------------------------------------#

# BEH6
# rmANOVA (DV = F0 value, Within = Block (1,2,3,4,5,6,7,8)) 
# Analysis BEH6 concerns how the F0 of the vocal responses during the experiment is influenced by the block of the experiment.
BEH6 <- aov_ez(id = "subj",
               dv = "recording_f0",
               data = df_block_f0,
               within = c("block"),
               type = 3)
summary(BEH6)

BEH6_ES <- BEH6$anova_table
BEH6_ES

# Plot: Vocal onset time by probe-onset and probe-type
ezPlot(
  data = df_block_f0 
  , dv = recording_f0 
  , wid = subj  
  , within= .(block)
  , x = .(block)
)

# Descriptive statistics
psych::describeBy(
  df_block_f0$recording_f0,
  list(df_block_f0$block)
)

# Follow-Up T-Tests
BEH6_EM <- emmeans(BEH6, ~ block)
BEH6_PWC <- pairs(BEH6_EM, adjust = "bonferroni")
BEH6_PWC

# Assumptions 
# Sphericity
# Checked in rmANOVA. Correction applied if needed.

# Normal distibution
df_block_f0 %>%
  ggplot(aes(x = recording_f0)) +
  geom_histogram(bins = 50) +
  facet_wrap(df_block_f0$block) +
  theme_ggstatsplot()

byf.shapiro(recording_f0 ~ block, 
            data = df_block_f0)

#Balance of the design
ezDesign(df_block_f0, x = block, y = subj) 

#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
# - Sphericity: Not applicable due to only 2 factor levels per factor
# - Balance of the design: OK
#------------------------------------------------------------------------------#

# BEH7
# rmANOVA (DV = Vocal onset time, Within = Block (1,2,3,4,5,6,7,8)) 
# Analysis BEH7 concerns how the F0 of the vocal onset time is influenced by the block of the experiment.
BEH7 <- aov_ez(id = "subj",
               dv = "recording_vot",
               data = df_block_vot,
               within = c("block"),
               type = 3)
summary(BEH7)

BEH7_ES <- BEH7$anova_table
BEH7_ES

# Plot: Vocal onset time by probe-onset and probe-type
ezPlot(
  data = df_block_vot 
  , dv = recording_vot 
  , wid = subj  
  , within= .(block)
  , x = .(block)
)

# Descriptive statistics
psych::describeBy(
  df_block_vot$recording_vot,
  list(df_block_vot$block)
)

# Follow-Up T-Tests
BEH7_EM <- emmeans(BEH7, ~ block)
BEH7_PWC <- pairs(BEH7_EM, adjust = "bonferroni")
BEH7_PWC

# Assumptions 
# Sphericity
# Checked in rmANOVA. Correction applied if needed.

# Normal distibution
df_block_vot %>%
  ggplot(aes(x = recording_vot)) +
  geom_histogram(bins = 50) +
  facet_wrap(df_block_vot$block) +
  theme_ggstatsplot()

byf.shapiro(recording_vot ~ block, 
            data = df_block_vot)

#Balance of the design
ezDesign(df_block_vot, x = block, y = subj) 

#------------------------------------------------------------------------------#
#
#

# Assumptions
# - Normal distribution:  
# - Sphericity: Not applicable due to only 2 factor levels per factor
# - Balance of the design: OK
#------------------------------------------------------------------------------#

# BEH8
# Descriptive only.
# Analysis BEH8 concerns how the F0 is influenced by sex.
psych::describeBy(
  df_subj_f0$recording_f0,
  list(df_subj_f0$var3_sex)
)


#-------------------------------------Plots-------------------------------------

# Plot 1: ERP Amplitudes by Probe Type, Task Instruction, and Probe Onset
P1 <- df_probe_properties_z %>%
  ggplot(aes(x = probe_type, y = f0_z, fill = probe_type)) +
  scale_fill_manual(values = c(colors$main_red, colors$main_blue), name = "Probe Type") +
  geom_boxplot(width=0.5, alpha = 1) +
  scale_x_discrete(limits=c("probe_f0_unaltered_z","probe_f0_altered_z"), labels = c('probe_f0_unaltered_z' = 'Unaltered', 'probe_f0_altered_z' = 'Altered')) +
  scale_y_continuous(n.breaks = 8) +
  labs(x = "Probe Type", y = "F0 [z]") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.position = "none") +
  geom_signif(comparisons=list(c("probe_f0_altered_z", "probe_f0_unaltered_z")), annotations="***",
              y_position = 5.5, tip_length = 0.02,  vjust=0.4) 

P1

# Save plot
ggsave(
  filename = "tid_psam_probe_z_f0.png", 
  plot = P1,
  width = 8,      
  height = 6,     
  dpi = 300,
  bg = "white"
)






