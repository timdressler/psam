# tid_psam_questionnaire_analysis.R
#
# Performs analysis of questionnaire data.
#
# Tim Dressler, 21.04.25

#-------------------------------------Set up------------------------------------
# Load packages
library(devtools)
library(renv)
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
INPATH <- file.path(MAINPATH, "data", "questionnaire_data_clean")
OUTPATH <- file.path(MAINPATH, "data", "analysis_data", "stats_questionnaire_analysis")
FUNPATH <- file.path(MAINPATH, "functions")

source(file.path(FUNPATH, "tid_psam_check_folder_TD.R"))
source(file.path(FUNPATH, "tid_psam_clean_up_folder_TD.R"))

tid_psam_check_folder_TD(MAINPATH, INPATH, OUTPATH)
tid_psam_clean_up_folder_TD(OUTPATH)

setwd(OUTPATH)

#------------------------------Load & Modify data-------------------------------

# Load FAL data
df_fal <- read_excel(file.path(INPATH, "fal_data_clean.xlsx"))
# Load NASA-TLX data
df_nasatlx <- read_excel(file.path(INPATH, "nasatlx_data_clean.xlsx"))
# Load SAM data
df_sam <- read_excel(file.path(INPATH, "sam_data_clean.xlsx"))

# Change FAL variable types
df_fal <- df_fal %>%
  mutate(across(c(2, 10, 18, 6), as.numeric))

df_fal <- df_fal %>%
  mutate(across(c(3,4,5,6,7,9,11,12,14,16,17,19,21,23,25,27), as.factor))

# Change NASA-TLX variable types
df_nasatlx <- df_nasatlx %>%
  mutate(across(c(2:6), as.numeric))

# Change SAM variable types
df_sam <- df_sam %>%
  mutate(across(c(2:17), as.numeric))

# Recode FAL variables
df_fal$var2_handedness <- recode(df_fal$var2_handedness, '1' = "right-handed", '2' = "left-handed", '3' = "two-handed")
df_fal$var3_sex <- recode(df_fal$var3_sex, '1' = "male", '2' = "female", '3' = "diverse")
df_fal$var4_education <- recode(df_fal$var4_education, '0' = "no degree", '1' = "Hauptschule", '2' = "Mittlere-Reife", '3' = "Abitur")
df_fal$var5_occupation <- recode(df_fal$var5_occupation, '1' = "student", '2' = "employed", '3' = "unemployed")
df_fal$var5_hearing_problems <- recode(df_fal$var5_hearing_problems, '1' = "yes", '2' = "no")
df_fal$var7_ringing_ears <- recode(df_fal$var7_ringing_ears, '1' = "yes", '2' = "no")
df_fal$var9_sleep_assessment <- recode(df_fal$var9_sleep_assessment, '1' = "normal", '2' = "rather long", '3' = "way too short")
df_fal$var10_alcohol_yesterday <- recode(df_fal$var10_alcohol_yesterday, '1' = "yes", '2' = "no")
df_fal$var11_alcohol_today <- recode(df_fal$var11_alcohol_today, '1' = "yes", '2' = "no")
df_fal$var12_smoking <- recode(df_fal$var12_smoking, '0' = "non-smoker", '1' = "little", '2' = "normal", '3' = "heavy")
df_fal$var13_coffee_and_other <- recode(df_fal$var13_coffee_and_other, '0' = "no coffee",'1' = "little", '2' = "normal", '3' = "much")
df_fal$var15_currently_neuro_treatment <- recode(df_fal$var15_currently_neuro_treatment, '1' = "yes", '2' = "no")
df_fal$var16_earlier_neuro_treatment <- recode(df_fal$var16_earlier_neuro_treatment, '1' = "yes", '2' = "no")
df_fal$var17_other_treatment <- recode(df_fal$var17_other_treatment, '1' = "yes", '2' = "no")
df_fal$var18_medication <- recode(df_fal$var18_medication, '1' = "yes", '2' = "no")
df_fal$var19_drugs <- recode(df_fal$var19_drugs, '1' = "yes", '2' = "no")

#-------------------------------Create needed dfs-------------------------------

df_sam_long <- df_sam %>%
  pivot_longer(
    cols = matches("_(mood|tiredness)_break\\d+"),
    names_to = c("dimension", "break_num"),
    names_pattern = ".*_(mood|tiredness)_break(\\d+)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = dimension,
    values_from = value
  ) %>%
  mutate(break_num = as.integer(break_num))
df_sam_long$break_num <- as.factor(df_sam_long$break_num)

df_nasatlx_long <- df_nasatlx %>%
  pivot_longer(
    cols = starts_with("var"),
    names_to = "dimension",
    values_to = "rating"
  ) %>%
  mutate(
    dimension = factor(dimension,
                       levels = c("var1_mental_demand", "var2_physical_demand", "var3_performance",
                                  "var4_effort", "var5_frustration"),
                       labels = c("Mental Demand", "Physical Demand", "Performance", "Effort", "Frustration")
    )
  )


#------------------------------------Analysis-----------------------------------

# FAL
# FAL1
# Descriptive only.
# Analysis FAL1 concerns the average age of the included participants. 
describe(df_fal$var1_age)

# FAL2
# Descriptive only.
# Analysis Analysis FAL2 concerns the sex of the included participants.  
table(df_fal$var3_sex)
prop.table(table(df_fal$var3_sex))

# Plot: Barplot of number of subjects of each sex
barplot(table(df_fal$var3_sex))

# FAL3
# Descriptive only.
# Analysis FAL3 concerns the education of the included participants. 
table(df_fal$var4_education)
prop.table(table(df_fal$var4_education))

# Plot: Barplot of number of subjects of each education level
barplot(table(df_fal$var4_education))

# FAL4
# Descriptive only.
# Analysis FAL4 concerns the occupation of the included participants. 
table(df_fal$var5_occupation)
prop.table(table(df_fal$var5_occupation))

# Plot: Barplot of number of subjects of each education level
barplot(table(df_fal$var5_occupation))

# NASA-TLX
# NASATLX1
# Descriptive only.
# Analysis NASATLX1 concerns the mental demand dimension of the NASA-TLX. 
describe(df_nasatlx$var1_mental_demand)

# NASATLX2
# Descriptive only.
# Analysis NASATLX2 concerns the physical demand dimension of the NASA-TLX. 
describe(df_nasatlx$var2_physical_demand)

# NASATLX3
# Descriptive only.
# Analysis NASATLX3 concerns the performance dimension of the NASA-TLX. 
describe(df_nasatlx$var3_performance)

# NASATLX4
# Descriptive only.
# Analysis NASATLX4 concerns the effort dimension of the NASA-TLX. 
describe(df_nasatlx$var4_effort)

# NASATLX5
# Descriptive only.
# Analysis Analysis NASATLX5 concerns the frustration dimension of the NASA-TLX. 
describe(df_nasatlx$var5_frustration)

# SAM
# SAM1
# Descriptive only.
# Analysis SAM1 concerns the mood dimension of the SAM during break 1. 
describe(df_sam$var1_mood_break1)

# SAM2
# Descriptive only.
# Analysis SAM2 concerns the mood dimension of the SAM during break 2. 
describe(df_sam$var3_mood_break2)

# SAM3
# Descriptive only.
# Analysis SAM3 concerns the mood dimension of the SAM during break 3. 
describe(df_sam$var5_mood_break3)

# SAM4
# Descriptive only.
# Analysis SAM4 concerns the mood dimension of the SAM during break 4. 
describe(df_sam$var7_mood_break4)

# SAM5
# Descriptive only.
# Analysis SAM5 concerns the mood dimension of the SAM during break 5. 
describe(df_sam$var9_mood_break5)

# SAM6
# Descriptive only.
# Analysis SAM6 concerns the mood dimension of the SAM during break 6. 
describe(df_sam$var11_mood_break6)

# SAM7
# Descriptive only.
# Analysis SAM7 concerns the mood dimension of the SAM during break 7. 
describe(df_sam$var13_mood_break7)

# SAM8
# Descriptive only.
# Analysis SAM8 concerns the mood dimension of the SAM during break 8 / after the experiment. 
describe(df_sam$var15_mood_break8)

# SAM9
# Descriptive only.
# Analysis SAM9 concerns the tiredness dimension of the SAM during break 1.  
describe(df_sam$var2_tiredness_break1)

# SAM10
# Descriptive only.
# Analysis SAM10 concerns the tiredness dimension of the SAM during break 2.  
describe(df_sam$var4_tiredness_break2)

# SAM11
# Descriptive only.
# Analysis SAM11 concerns the tiredness dimension of the SAM during break 3.  
describe(df_sam$var6_tiredness_break3)

# SAM12
# Descriptive only.
# Analysis SAM12 concerns the tiredness dimension of the SAM during break 4.  
describe(df_sam$var8_tiredness_break4)

# SAM13
# Descriptive only.
# Analysis SAM13 concerns the tiredness dimension of the SAM during break 5.  
describe(df_sam$var10_tiredness_break5)

# SAM14
# Descriptive only.
# Analysis SAM14 concerns the tiredness dimension of the SAM during break 6.  
describe(df_sam$var12_tiredness_break6)

# SAM15
# Descriptive only.
# Analysis SAM15 concerns the tiredness dimension of the SAM during break 7.  
describe(df_sam$var14_tiredness_break7)

# SAM16
# Descriptive only.
# Analysis SAM16 concerns the tiredness dimension of the SAM during break 8 / after the experiment.  
describe(df_sam$var16_tiredness_break8)

# SAM17
# rmANOVA (DV = Mood, Within = Break (1,2,3,4,5,6,7,8)) 
# Analysis SAM17 concerns how the mood dimension of the SAM is influenced by the block/break.
SAM17 <- aov_ez(id = "subj",
               dv = "mood",
               data = df_sam_long,
               within = c("break_num"),
               type = 3)
summary(SAM17)

SAM17_ES <- SAM17$anova_table
SAM17_ES

# Plot: Vocal onset time by probe-onset and probe-type
ezPlot(
  data = df_sam_long 
  , dv = mood 
  , wid = subj  
  , within= .(break_num)
  , x = .(break_num)
  , type = 3
)

# Descriptive statistics
psych::describeBy(
  df_sam_long$mood,
  list(df_sam_long$break_num)
)

# Follow-Up T-Tests
SAM17_EM <- emmeans(SAM17, ~ break_num)
SAM17_PWC <- pairs(SAM17_EM, adjust = "bonferroni")
SAM17_PWC

# Assupmtions 
# Sphericity
# Checked in rmANOVA. Correction applied if needed.

# Normal distibution
df_sam_long %>%
  ggplot(aes(x = mood)) +
  geom_histogram(bins = 50) +
  facet_wrap(df_sam_long$break_num) +
  theme_ggstatsplot()

byf.shapiro(mood ~ break_num, 
            data = df_sam_long)

#Balance of the design
ezDesign(df_sam_long, x = break_num, y = subj) 

#------------------------------------------------------------------------------#
# break_num n.s.

# Assumptions
# - Normal distribution: Shapiro-Wilk sig. but due to visual inspection OK  
# - Sphericity: not OK, GG correction used
# - Balance of the design: OK
#------------------------------------------------------------------------------#

# SAM18
# rmANOVA (DV = Tiredness, Within = Break (1,2,3,4,5,6,7,8)) 
# Analysis SAM18 concerns how the tiredness dimension of the SAM is influenced by the block/break.
SAM18 <- aov_ez(id = "subj",
                dv = "tiredness",
                data = df_sam_long,
                within = c("break_num"),
                type = 3)
summary(SAM18)

SAM18_ES <- SAM18$anova_table
SAM18_ES

# Plot: Vocal onset time by probe-onset and probe-type
ezPlot(
  data = df_sam_long 
  , dv = tiredness 
  , wid = subj  
  , within= .(break_num)
  , x = .(break_num)
)

# Descriptive statistics
psych::describeBy(
  df_sam_long$tiredness,
  list(df_sam_long$break_num)
)

# Follow-Up T-Tests
SAM18_EM <- emmeans(SAM18, ~ break_num)
SAM18_PWC <- pairs(SAM18_EM, adjust = "bonferroni")
SAM18_PWC

# Assupmtions 
# Sphericity
# Checked in rmANOVA. Correction applied if needed.

# Normal distibution
df_sam_long %>%
  ggplot(aes(x = tiredness)) +
  geom_histogram(bins = 50) +
  facet_wrap(df_sam_long$break_num) +
  theme_ggstatsplot()

byf.shapiro(tiredness ~ break_num, 
             data = df_sam_long)

#Balance of the design
ezDesign(df_sam_long, x = break_num, y = subj) 

#------------------------------------------------------------------------------#
# break_num n.s.

# Assumptions
# - Normal distribution: OK
# - Sphericity: not OK, GG correction used
# - Balance of the design: OK
#------------------------------------------------------------------------------#

#-------------------------------------Plots-------------------------------------

# Plot 1: Mood by Block

P1 <- ggplot(df_sam_long, aes(x = break_num, y = mood)) +
  geom_boxplot(width = 0.3, fill = colors$main_blue, color = "black", alpha = 1) +
  scale_y_continuous(
    limits = c(1, 9),
    breaks = 1:9,
    labels = c(
      "1\nVery Negative", "2", "3", "4", "5", "6", "7", "8", "9\nVery Positive"
    )
  ) +
  labs(
    x = 'Break Number',
    y = 'Mood'
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
P1

ggsave(
  filename = "tid_psam_SAM_mood.pdf", 
  plot = P1,
  width = 8,      
  height = 6,     
  dpi = 900,
  bg = "white"
)

# Plot 2: Tiredness by Block

P2 <- ggplot(df_sam_long, aes(x = break_num, y = tiredness)) +
  geom_boxplot(width = 0.3, fill = colors$main_green, color = "black", alpha = 1) +
  scale_y_continuous(
    limits = c(1, 9),
    breaks = 1:9,
    labels = c(
      "1\nVery Awake", "2", "3", "4", "5", "6", "7", "8", "9\nVery Tired"
    )
  ) +
  labs(
    x = 'Break Number',
    y = 'Tiredness'
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
P2

ggsave(
  filename = "tid_psam_SAM_tiredness.pdf", 
  plot = P2,
  width = 8,      
  height = 6,     
  dpi = 900,
  bg = "white"
)

# Combine Plot 1 and Plot 2
P1_P2 <- ggdraw() +
  draw_plot(P1, x = 0, y = 0, width = .5, height = 1) +
  draw_plot(P2, x = .5, y = 0, width = .5, height = 1) +
  draw_plot_label(label = c("A", "B"), size = 15,
                  x = c(0, 0.5), y = c(1, 1))
P1_P2

ggsave(
  filename = "tid_psam_SAM_combined.pdf", 
  plot = P1_P2,
  width = 8,      
  height = 6,     
  dpi = 900,
  bg = "white"
)

# Plot 3: NASA-TLX

P3 <- ggplot(df_nasatlx_long, aes(x = dimension, y = rating)) +
  geom_boxplot(fill = colors$main_blue, color = "black", width = 0.5) +
  scale_x_discrete(labels = c(
    "Mental Demand" = "Mental Demand",
    "Physical Demand" = "Physical Demand",
    "Effort" = "Effort",
    "Frustration" = "Frustration",
    "Performance" = "Performance*"
  )) +
  scale_y_continuous(
    limits = c(0, 10),
    breaks = seq(0, 10, 2),
    labels = c("0\nLow/\nGood*", "2", "4", "6", "8", "10\nHigh/\nBad*")
  ) +
  labs(
    x = "Dimension",
    y = "Rating"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
P3

ggsave(
  filename = "tid_psam_NASATLX.pdf", 
  plot = P3,
  width = 8,      
  height = 6,     
  dpi = 900,
  bg = "white"
)







