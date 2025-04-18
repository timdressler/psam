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

MAINPATH <- gsub("analysis_script/R", "", SCRIPTPATH)
OUTPATH <- file.path(MAINPATH, "data", "analysis_data", "stats_beh_analysis")

setwd(OUTPATH)

#-----------------------------------Load Data-----------------------------------

#load raw main data
df_beh <- read.csv(file.path(MAINPATH, "data/analysis_data/erp_analysis/erp_analysis.csv"))








