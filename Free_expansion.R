################################################################################
# Script to conduct GO analysis of CLIP-seq & Ribo-footprinting data [052824]
# ------------------------------------------------------------------------------
## 1. Load libraries & Set workspace
# ------------------------------------------------------------------------------
library(clusterProfiler)
library(ggplot2)
library(colorRamp2)

setwd("/home/local/hoeyoungkim_000504/coursework/SNU-BI1")
# ------------------------------------------------------------------------------
## Load raw counts data generated with FeatureCounts
