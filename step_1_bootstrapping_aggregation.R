# title: 'Venn diagrams and random forest bootstrapping aggregation analysis'
# author: 'Marie-Madlen Pust'
# date: '05 December 2021'

############################################################################################################
# clean global environment
rm(list=ls())

# set working directory
setwd('C:/R')

############################################################################################################

# define global variables
quantile_range = c(0.15, 0.25, 0.35)
set.seed(1)
random_seeds <- sample(1:200, 100, replace = FALSE) 

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, 'Package'])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)}

# bootstrap the standard error of the median
bootstrap_median <- function(data, num) {
  resamples <- lapply(1:num, function(i) sample(data, replace=T))
  r.median <- sapply(resamples, median)
  std.err <- sqrt(var(r.median))
  list(std.err=std.err, resamples=resamples, medians=r.median)}

# required packages
packages <- c('readr', 'dplyr','stringr','tidyr','ggplot2', 'rstatix', 'factoextra','gmodels','ggpubr','plyr','purrr',
              'Hmisc','reshape','vegan','magrittr','scales', 'grid','reshape2','phyloseq','rcompanion','viridis',
              'randomForest','knitr','ggrepel','forcats','ggVennDiagram', 'pheatmap', 'Boruta', 'rstatix')

############################################################################################################

# load required R packages
ipak(packages)

############################################################################################################

# import datasets
# import meta data 
md <- read_delim('input_files/meta_data/spatial_metadata_2020_12_2.csv', ';', escape_double = FALSE, trim_ws = TRUE)
md <- data.frame(md)
rownames(md) <- md$Sample 
md$Sample <- NULL
# add column with BMI from Weight and Height column
md$BMI <- round(md$Weight / ((md$Height/100) * (md$Height/100)),1)

# subset metadata table into smaller tables
md_0 <- subset(md, AgeGroup == '0 years')
md_0_h <- subset(md_0, State == 'Healthy')
md_0_cf <- subset(md_0, State == 'CF')
md_1_3 <- subset(md, AgeGroup == '1-3 years')
md_1_3_h <- subset(md_1_3, State == 'Healthy')
md_1_3_cf <- subset(md_1_3, State == 'CF')
md_4_6 <- subset(md, AgeGroup == '4-6 years')
md_4_6_h <- subset(md_4_6, State == 'Healthy')
md_4_6_cf <- subset(md_4_6, State == 'CF')


# import pangenome data
# BCPHC-normalised count table
ds_pangenome <- read_delim('input_files/taxonomic_data/pangenome.bcphc.merged.csv', ';', escape_double = FALSE, trim_ws = TRUE, skip_empty_rows = TRUE)
ds_pangenome <- data.frame(ds_pangenome)
# concatenate duplicate Species rows and sum count value
ds_pangenome <- ddply(ds_pangenome, 'Species', numcolwise(sum))
# set species as row names
rownames(ds_pangenome) <- ds_pangenome$Species
# remove species column
ds_pangenome$Species <- NULL


# import pangenome data
# RLE-normalised count table
ds_pangenome_rle <- read_delim('input_files/taxonomic_data/pangenome.rle.merged.csv', ';', escape_double = FALSE, trim_ws = TRUE, skip_empty_rows = TRUE)
ds_pangenome_rle <- data.frame(ds_pangenome_rle)
# concatenate duplicate Species rows and sum count value
ds_pangenome_rle <- ddply(ds_pangenome_rle, 'Species', numcolwise(sum))
# set species as row names
rownames(ds_pangenome_rle) <- ds_pangenome_rle$Species
# remove species column
ds_pangenome_rle$Species <- NULL


# import pangenome data
# vst-normalised count table
ds_pangenome_vst <- read_delim('input_files/taxonomic_data/pangenome.vst.merged.csv', ';', escape_double = FALSE, trim_ws = TRUE, skip_empty_rows = TRUE)
ds_pangenome_vst <- data.frame(ds_pangenome_vst)
# concatenate duplicate Species rows and sum count value
ds_pangenome_vst <- ddply(ds_pangenome_vst, 'Species', numcolwise(sum))
# set species as row names
rownames(ds_pangenome_vst) <- ds_pangenome_vst$Species
# remove species column
ds_pangenome_vst$Species <- NULL


# import one-strain per species data
# bcphc-normalised count table
ds_osps <- read_delim('input_files/taxonomic_data/osps.bcphc.merged.csv', ';', escape_double = FALSE, trim_ws = TRUE, skip_empty_rows = TRUE)
ds_osps <- data.frame(ds_osps)
# concatenate duplicate Species rows and sum count value
ds_osps <- ddply(ds_osps, 'Species', numcolwise(sum))
# set species as row names
rownames(ds_osps) <- ds_osps$Species
# remove species column
ds_osps$Species <- NULL


# import one-strain per species data
# RLE normalisation
ds_osps_rle <- read_delim('input_files/taxonomic_data/osps.rle.merged.csv', ';', escape_double = FALSE, trim_ws = TRUE, skip_empty_rows = TRUE)
ds_osps_rle <- data.frame(ds_osps_rle)
# concatenate duplicate Species rows and sum count value
ds_osps_rle <- ddply(ds_osps_rle, 'Species', numcolwise(sum))
rownames(ds_osps_rle) <- ds_osps_rle$Species
ds_osps_rle$Species <- NULL


# import one-strain per species data
# vst normalisation
ds_osps_vst <- read_delim('input_files/taxonomic_data/osps.vst.merged.csv', ';', escape_double = FALSE, trim_ws = TRUE, skip_empty_rows = TRUE)
ds_osps_vst <- data.frame(ds_osps_vst)
# concatenate duplicate Species rows and sum count value
ds_osps_vst <- ddply(ds_osps_vst, 'Species', numcolwise(sum))
rownames(ds_osps_vst) <- ds_osps_vst$Species
# remove species column
ds_osps_vst$Species <- NULL

############################################################################################################
# PANGENOME
# define core and rare species based on three rarity thresholds (defined above in variable abund_quantile)
# this is repeated for all three normalisations steps (bcphc, rle, vst)

# age group (0 years), healthy, pangenome database, bcphc-normalised
c_0_h_pangenome <- select(ds_pangenome, rownames(md_0_h))
# remove species rows that sum to zero
c_0_h_pangenome <- c_0_h_pangenome[rowSums(c_0_h_pangenome[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_0_pangenome = colSums(c_0_h_pangenome)
# total sum
sum_all_healthy_0_pangenome = sum(sum_healthy_0_pangenome)
# add abundance column
c_0_h_pangenome$abundance <- (rowSums(c_0_h_pangenome[,1:ncol(c_0_h_pangenome)]) / sum_all_healthy_0_pangenome) * 100
# sort abundance decreasing
c_0_h_pangenome <- c_0_h_pangenome[with(c_0_h_pangenome, order(-abundance)), ]
# obtain cumulative sum
c_0_h_pangenome$cumsum <- cumsum(c_0_h_pangenome$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile <- quantile(c_0_h_pangenome$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_0_h_pangenome_core_a1 <- subset(c_0_h_pangenome, cumsum <= abund_quantile[1])
# remove non-numeric columns
c_0_h_pangenome_core_a1$cumsum <- NULL
c_0_h_pangenome_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_0_h_pangenome_core_a2 <- subset(c_0_h_pangenome, cumsum <= abund_quantile[2])
# remove non-numeric columns
c_0_h_pangenome_core_a2$cumsum <- NULL
c_0_h_pangenome_core_a2$abundance <- NULL
# # subset core species based on third defined abundance quantile
c_0_h_pangenome_core_a3 <- subset(c_0_h_pangenome, cumsum <= abund_quantile[3])
# remove non-numeric columns
c_0_h_pangenome_core_a3$cumsum <- NULL
c_0_h_pangenome_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_0_h_pangenome_rare_a1 <- subset(c_0_h_pangenome, cumsum > abund_quantile[1])
# remove non-numeric columns
c_0_h_pangenome_rare_a1$cumsum <- NULL
c_0_h_pangenome_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_0_h_pangenome_rare_a2 <- subset(c_0_h_pangenome, cumsum > abund_quantile[2])
# remove non-numeric columns
c_0_h_pangenome_rare_a2$cumsum <- NULL
c_0_h_pangenome_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_0_h_pangenome_rare_a3 <- subset(c_0_h_pangenome, cumsum > abund_quantile[3])
# remove non-numeric columns
c_0_h_pangenome_rare_a3$cumsum <- NULL
c_0_h_pangenome_rare_a3$abundance <- NULL


# Age group (1-3 years), healthy, pangenome database, bcphc-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_1_3_h_pangenome <- select(ds_pangenome, rownames(md_1_3_h))
# remove species rows that sum to zero
c_1_3_h_pangenome <- c_1_3_h_pangenome[rowSums(c_1_3_h_pangenome[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_1_3_pangenome = colSums(c_1_3_h_pangenome)
# obtain total sum
sum_all_healthy_1_3_pangenome = sum(sum_healthy_1_3_pangenome)
# add abundance column
c_1_3_h_pangenome$abundance <- (rowSums(c_1_3_h_pangenome[,1:ncol(c_1_3_h_pangenome)]) / sum_all_healthy_1_3_pangenome) * 100
# sort abundance decreasing
c_1_3_h_pangenome <- c_1_3_h_pangenome[with(c_1_3_h_pangenome, order(-abundance)), ]
# obtain cumulative sum
c_1_3_h_pangenome$cumsum <- cumsum(c_1_3_h_pangenome$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_1_3_h <- quantile(c_1_3_h_pangenome$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_1_3_h_pangenome_core_a1 <- subset(c_1_3_h_pangenome, cumsum <= abund_quantile_1_3_h[1])
# remove non-numeric columns
c_1_3_h_pangenome_core_a1$cumsum <- NULL
c_1_3_h_pangenome_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_1_3_h_pangenome_core_a2 <- subset(c_1_3_h_pangenome, cumsum <= abund_quantile_1_3_h[2])
# remove non-numeric columns
c_1_3_h_pangenome_core_a2$cumsum <- NULL
c_1_3_h_pangenome_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_1_3_h_pangenome_core_a3 <- subset(c_1_3_h_pangenome, cumsum <= abund_quantile_1_3_h[3])
# remove non-numeric columns
c_1_3_h_pangenome_core_a3$cumsum <- NULL
c_1_3_h_pangenome_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_1_3_h_pangenome_rare_a1 <- subset(c_1_3_h_pangenome, cumsum > abund_quantile_1_3_h[1])
# remove non-numeric columns
c_1_3_h_pangenome_rare_a1$cumsum <- NULL
c_1_3_h_pangenome_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_1_3_h_pangenome_rare_a2 <- subset(c_1_3_h_pangenome, cumsum > abund_quantile_1_3_h[2])
# remove non-numeric columns
c_1_3_h_pangenome_rare_a2$cumsum <- NULL
c_1_3_h_pangenome_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_1_3_h_pangenome_rare_a3 <- subset(c_1_3_h_pangenome, cumsum > abund_quantile_1_3_h[3])
# remove non-numeric columns
c_1_3_h_pangenome_rare_a3$cumsum <- NULL
c_1_3_h_pangenome_rare_a3$abundance <- NULL


# Age group (4-6 years), healthy, pangenome database, bcphc-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_4_6_h_pangenome <- select(ds_pangenome, rownames(md_4_6_h))
# remove species rows that sum to zero
c_4_6_h_pangenome <- c_4_6_h_pangenome[rowSums(c_4_6_h_pangenome[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_4_6_pangenome = colSums(c_4_6_h_pangenome)
# get total sum
sum_all_healthy_4_6_pangenome = sum(sum_healthy_4_6_pangenome)
# add abundance column
c_4_6_h_pangenome$abundance <- (rowSums(c_4_6_h_pangenome[,1:ncol(c_4_6_h_pangenome)]) / sum_all_healthy_4_6_pangenome) * 100
# sort abundance decreasing
c_4_6_h_pangenome <- c_4_6_h_pangenome[with(c_4_6_h_pangenome, order(-abundance)), ]
# get cumulative sum
c_4_6_h_pangenome$cumsum <- cumsum(c_4_6_h_pangenome$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_4_6_h <- quantile(c_4_6_h_pangenome$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_4_6_h_pangenome_core_a1 <- subset(c_4_6_h_pangenome, cumsum <= abund_quantile_4_6_h[1])
# remove non-numeric columns
c_4_6_h_pangenome_core_a1$cumsum <- NULL
c_4_6_h_pangenome_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_4_6_h_pangenome_core_a2 <- subset(c_4_6_h_pangenome, cumsum <= abund_quantile_4_6_h[2])
# remove non-numeric columns
c_4_6_h_pangenome_core_a2$cumsum <- NULL
c_4_6_h_pangenome_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_4_6_h_pangenome_core_a3 <- subset(c_4_6_h_pangenome, cumsum <= abund_quantile_4_6_h[3])
# remove non-numeric columns
c_4_6_h_pangenome_core_a3$cumsum <- NULL
c_4_6_h_pangenome_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_4_6_h_pangenome_rare_a1 <- subset(c_4_6_h_pangenome, cumsum > abund_quantile_4_6_h[1])
# remove non-numeric columns
c_4_6_h_pangenome_rare_a1$cumsum <- NULL
c_4_6_h_pangenome_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_4_6_h_pangenome_rare_a2 <- subset(c_4_6_h_pangenome, cumsum > abund_quantile_4_6_h[2])
# remove non-numeric columns
c_4_6_h_pangenome_rare_a2$cumsum <- NULL
c_4_6_h_pangenome_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_4_6_h_pangenome_rare_a3 <- subset(c_4_6_h_pangenome, cumsum > abund_quantile_4_6_h[3])
# remove non-numeric columns
c_4_6_h_pangenome_rare_a3$cumsum <- NULL
c_4_6_h_pangenome_rare_a3$abundance <- NULL


# Age group (0 years), CF, pangenome database, bcphc-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_0_cf_pangenome <- select(ds_pangenome, rownames(md_0_cf))
# remove species rows that sum to zero
c_0_cf_pangenome <- c_0_cf_pangenome[rowSums(c_0_cf_pangenome[, -1])>0, ]
# sum species counts per CF child
sum_cf_0_pangenome = colSums(c_0_cf_pangenome)
# get total sum
sum_all_cf_0_pangenome = sum(sum_cf_0_pangenome)
# add abundance column
c_0_cf_pangenome$abundance <- (rowSums(c_0_cf_pangenome[,1:ncol(c_0_cf_pangenome)]) / sum_all_cf_0_pangenome) * 100
# sort abundance decreasing
c_0_cf_pangenome <- c_0_cf_pangenome[with(c_0_cf_pangenome, order(-abundance)), ]
# add column with cumulative sum
c_0_cf_pangenome$cumsum <- cumsum(c_0_cf_pangenome$abundance)

# create sample quantiles in dataset
abund_quantile_0_cf <- quantile(c_0_cf_pangenome$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_0_cf_pangenome_core_a1 <- subset(c_0_cf_pangenome, cumsum <= abund_quantile_0_cf[1])
# remove non-numeric columns
c_0_cf_pangenome_core_a1$cumsum <- NULL
c_0_cf_pangenome_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_0_cf_pangenome_core_a2 <- subset(c_0_cf_pangenome, cumsum <= abund_quantile_0_cf[2])
# remove non-numeric columns
c_0_cf_pangenome_core_a2$cumsum <- NULL
c_0_cf_pangenome_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_0_cf_pangenome_core_a3 <- subset(c_0_cf_pangenome, cumsum <= abund_quantile_0_cf[3])
# remove non-numeric columns
c_0_cf_pangenome_core_a3$cumsum <- NULL
c_0_cf_pangenome_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_0_cf_pangenome_rare_a1 <- subset(c_0_cf_pangenome, cumsum > abund_quantile_0_cf[1])
# remove non-numeric columns
c_0_cf_pangenome_rare_a1$cumsum <- NULL
c_0_cf_pangenome_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_0_cf_pangenome_rare_a2 <- subset(c_0_cf_pangenome, cumsum > abund_quantile_0_cf[2])
# remove non-numeric columns
c_0_cf_pangenome_rare_a2$cumsum <- NULL
c_0_cf_pangenome_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_0_cf_pangenome_rare_a3 <- subset(c_0_cf_pangenome, cumsum > abund_quantile_0_cf[3])
# remove non-numeric columns
c_0_cf_pangenome_rare_a3$cumsum <- NULL
c_0_cf_pangenome_rare_a3$abundance <- NULL


# Age group (1-3 years), CF, pangenome database, bcphc-normalised
c_1_3_cf_pangenome <- select(ds_pangenome, rownames(md_1_3_cf))
# remove species rows that sum to zero
c_1_3_cf_pangenome <- c_1_3_cf_pangenome[rowSums(c_1_3_cf_pangenome[, -1])>0, ]
# sum species counts per CF child
sum_cf_1_3_pangenome = colSums(c_1_3_cf_pangenome)
# get total sum
sum_all_cf_1_3_pangenome = sum(sum_cf_1_3_pangenome)
# add abundance column
c_1_3_cf_pangenome$abundance <- (rowSums(c_1_3_cf_pangenome[,1:ncol(c_1_3_cf_pangenome)]) / sum_all_cf_1_3_pangenome) * 100
# sort abundance decreasing
c_1_3_cf_pangenome <- c_1_3_cf_pangenome[with(c_1_3_cf_pangenome, order(-abundance)), ]
# add column with cumulative sum
c_1_3_cf_pangenome$cumsum <- cumsum(c_1_3_cf_pangenome$abundance)

# table with abundant species
# create sample quantiles in dataset
abund_quantile_1_3_cf <- quantile(c_1_3_cf_pangenome$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_1_3_cf_pangenome_core_a1 <- subset(c_1_3_cf_pangenome, cumsum <= abund_quantile_1_3_cf[1])
# remove non-numeric columns
c_1_3_cf_pangenome_core_a1$cumsum <- NULL
c_1_3_cf_pangenome_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_1_3_cf_pangenome_core_a2 <- subset(c_1_3_cf_pangenome, cumsum <= abund_quantile_1_3_cf[2])
# remove non-numeric columns
c_1_3_cf_pangenome_core_a2$cumsum <- NULL
c_1_3_cf_pangenome_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_1_3_cf_pangenome_core_a3 <- subset(c_1_3_cf_pangenome, cumsum <= abund_quantile_1_3_cf[3])
# remove non-numeric columns
c_1_3_cf_pangenome_core_a3$cumsum <- NULL
c_1_3_cf_pangenome_core_a3$abundance <- NULL


# table with rare species
# subset rare species based on first defined abundance quantile
c_1_3_cf_pangenome_rare_a1 <- subset(c_1_3_cf_pangenome, cumsum > abund_quantile_1_3_cf[1])
# remove non-numeric columns
c_1_3_cf_pangenome_rare_a1$cumsum <- NULL
c_1_3_cf_pangenome_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_1_3_cf_pangenome_rare_a2 <- subset(c_1_3_cf_pangenome, cumsum > abund_quantile_1_3_cf[2])
# remove non-numeric columns
c_1_3_cf_pangenome_rare_a2$cumsum <- NULL
c_1_3_cf_pangenome_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_1_3_cf_pangenome_rare_a3 <- subset(c_1_3_cf_pangenome, cumsum > abund_quantile_1_3_cf[3])
# remove non-numeric columns
c_1_3_cf_pangenome_rare_a3$cumsum <- NULL
c_1_3_cf_pangenome_rare_a3$abundance <- NULL


# Age group (4-6 years), CF, pangenome database, bcphc-normalised
c_4_6_cf_pangenome <- select(ds_pangenome, rownames(md_4_6_cf))
# remove species rows that sum to zero
c_4_6_cf_pangenome <- c_4_6_cf_pangenome[rowSums(c_4_6_cf_pangenome[, -1])>0, ]
# sum species counts per CF child
sum_cf_4_6_pangenome = colSums(c_4_6_cf_pangenome)
# get total sum
sum_all_cf_4_6_pangenome = sum(sum_cf_4_6_pangenome)
# add abundance column
c_4_6_cf_pangenome$abundance <- (rowSums(c_4_6_cf_pangenome[,1:ncol(c_4_6_cf_pangenome)]) / sum_all_cf_4_6_pangenome) * 100
# sort abundance decreasing
c_4_6_cf_pangenome <- c_4_6_cf_pangenome[with(c_4_6_cf_pangenome, order(-abundance)), ]
# add column with cumulative sum
c_4_6_cf_pangenome$cumsum <- cumsum(c_4_6_cf_pangenome$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_4_6_cf <- quantile(c_4_6_cf_pangenome$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_4_6_cf_pangenome_core_a1 <- subset(c_4_6_cf_pangenome, cumsum <= abund_quantile_4_6_cf[1])
# remove non-numeric columns
c_4_6_cf_pangenome_core_a1$cumsum <- NULL
c_4_6_cf_pangenome_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_4_6_cf_pangenome_core_a2 <- subset(c_4_6_cf_pangenome, cumsum <= abund_quantile_4_6_cf[2])
# remove non-numeric columns
c_4_6_cf_pangenome_core_a2$cumsum <- NULL
c_4_6_cf_pangenome_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_4_6_cf_pangenome_core_a3 <- subset(c_4_6_cf_pangenome, cumsum <= abund_quantile_4_6_cf[3])
# remove non-numeric columns
c_4_6_cf_pangenome_core_a3$cumsum <- NULL
c_4_6_cf_pangenome_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_4_6_cf_pangenome_rare_a1 <- subset(c_4_6_cf_pangenome, cumsum > abund_quantile_4_6_cf[1])
# remove non-numeric columns
c_4_6_cf_pangenome_rare_a1$cumsum <- NULL
c_4_6_cf_pangenome_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_4_6_cf_pangenome_rare_a2 <- subset(c_4_6_cf_pangenome, cumsum > abund_quantile_4_6_cf[2])
# remove non-numeric columns
c_4_6_cf_pangenome_rare_a2$cumsum <- NULL
c_4_6_cf_pangenome_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_4_6_cf_pangenome_rare_a3 <- subset(c_4_6_cf_pangenome, cumsum > abund_quantile_4_6_cf[3])
# remove non-numeric columns
c_4_6_cf_pangenome_rare_a3$cumsum <- NULL
c_4_6_cf_pangenome_rare_a3$abundance <- NULL


# create vector with all rare species (15th abundance percentile, pangenome, healthy)
healthy_rare_a1_pangenome <- c(rownames(c_0_h_pangenome_rare_a1), 
                               rownames(c_1_3_h_pangenome_rare_a1), 
                               rownames(c_4_6_h_pangenome_rare_a1))
# remove duplicate entries
healthy_rare_a1_pangenome <- healthy_rare_a1_pangenome[!duplicated(healthy_rare_a1_pangenome)]

# create vector with all core species (15th abundance percentile, pangenome, healthy)
healthy_core_a1_pangenome <- c(rownames(c_0_h_pangenome_core_a1), 
                               rownames(c_1_3_h_pangenome_core_a1), 
                               rownames(c_4_6_h_pangenome_core_a1))
# remove duplicate entries
healthy_core_a1_pangenome <- healthy_core_a1_pangenome[!duplicated(healthy_core_a1_pangenome)]

# create vector with all rare species (25th abundance percentile, pangenome, healthy)
healthy_rare_a2_pangenome <- c(rownames(c_0_h_pangenome_rare_a2), 
                               rownames(c_1_3_h_pangenome_rare_a2), 
                               rownames(c_4_6_h_pangenome_rare_a2))
# remove duplicate entries
healthy_rare_a2_pangenome <- healthy_rare_a2_pangenome[!duplicated(healthy_rare_a2_pangenome)]

# create vector with all core species (25th abundance percentile, pangenome, healthy)
healthy_core_a2_pangenome <- c(rownames(c_0_h_pangenome_core_a2), 
                               rownames(c_1_3_h_pangenome_core_a2), 
                               rownames(c_4_6_h_pangenome_core_a2))
# remove duplicate entries
healthy_core_a2_pangenome <- healthy_core_a2_pangenome[!duplicated(healthy_core_a2_pangenome)]

# create vector with all rare species (35th abundance percentile, pangenome, healthy)
healthy_rare_a3_pangenome <- c(rownames(c_0_h_pangenome_rare_a3), 
                               rownames(c_1_3_h_pangenome_rare_a3), 
                               rownames(c_4_6_h_pangenome_rare_a3))
# remove duplicate entries
healthy_rare_a3_pangenome <- healthy_rare_a3_pangenome[!duplicated(healthy_rare_a3_pangenome)]

# create vector with all core species (35th abundance percentile, pangenome, healthy)
healthy_core_a3_pangenome <- c(rownames(c_0_h_pangenome_core_a3), 
                               rownames(c_1_3_h_pangenome_core_a3), 
                               rownames(c_4_6_h_pangenome_core_a3))
# remove duplicate entries
healthy_core_a3_pangenome <- healthy_core_a3_pangenome[!duplicated(healthy_core_a3_pangenome)]


# create vector with all rare species (15th abundance percentile, pangenome, CF)
cf_rare_a1_pangenome <- c(rownames(c_0_cf_pangenome_rare_a1), 
                          rownames(c_1_3_cf_pangenome_rare_a1), 
                          rownames(c_4_6_cf_pangenome_rare_a1))
# remove duplicate entries
cf_rare_a1_pangenome <- cf_rare_a1_pangenome[!duplicated(cf_rare_a1_pangenome)]

# create vector with all core species (15th abundance percentile, pangenome, CF)
cf_core_a1_pangenome <- c(rownames(c_0_cf_pangenome_core_a1), 
                          rownames(c_1_3_cf_pangenome_core_a1), 
                          rownames(c_4_6_cf_pangenome_core_a1))
# remove duplicate entries
cf_core_a1_pangenome <- cf_core_a1_pangenome[!duplicated(cf_core_a1_pangenome)]

# create vector with all rare species (25th abundance percentile, pangenome, CF)
cf_rare_a2_pangenome <- c(rownames(c_0_cf_pangenome_rare_a2), 
                          rownames(c_1_3_cf_pangenome_rare_a2), 
                          rownames(c_4_6_cf_pangenome_rare_a2))
# remove duplicate entries
cf_rare_a2_pangenome <- cf_rare_a2_pangenome[!duplicated(cf_rare_a2_pangenome)]

# create vector with all core species (25th abundance percentile, pangenome, CF)
cf_core_a2_pangenome <- c(rownames(c_0_cf_pangenome_core_a2), 
                          rownames(c_1_3_cf_pangenome_core_a2), 
                          rownames(c_4_6_cf_pangenome_core_a2))
# remove duplicate entries
cf_core_a2_pangenome <- cf_core_a2_pangenome[!duplicated(cf_core_a2_pangenome)]

# create vector with all rare species (35th abundance percentile, pangenome, CF)
cf_rare_a3_pangenome <- c(rownames(c_0_cf_pangenome_rare_a3), 
                          rownames(c_1_3_cf_pangenome_rare_a3), 
                          rownames(c_4_6_cf_pangenome_rare_a3))
# remove duplicate entries
cf_rare_a3_pangenome <- cf_rare_a3_pangenome[!duplicated(cf_rare_a3_pangenome)]

# create vector with all core species (35th abundance percentile, pangenome, CF)
cf_core_a3_pangenome <- c(rownames(c_0_cf_pangenome_core_a3), 
                          rownames(c_1_3_cf_pangenome_core_a3), 
                          rownames(c_4_6_cf_pangenome_core_a3))
# remove duplicate entries
cf_core_a3_pangenome <- cf_core_a3_pangenome[!duplicated(cf_core_a3_pangenome)]



# age group (0 years), healthy, pangenome database, rle-normalised
c_0_h_pangenome_rle <- select(ds_pangenome_rle, rownames(md_0_h))
# remove species rows that sum to zero
c_0_h_pangenome_rle <- c_0_h_pangenome_rle[rowSums(c_0_h_pangenome_rle[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_0_pangenome_rle = colSums(c_0_h_pangenome_rle)
# total sum
sum_all_healthy_0_pangenome_rle = sum(sum_healthy_0_pangenome_rle)
# add abundance column
c_0_h_pangenome_rle$abundance <- (rowSums(c_0_h_pangenome_rle[,1:ncol(c_0_h_pangenome_rle)]) / sum_all_healthy_0_pangenome_rle) * 100
# sort abundance decreasing
c_0_h_pangenome_rle <- c_0_h_pangenome_rle[with(c_0_h_pangenome_rle, order(-abundance)), ]
# obtain cumulative sum
c_0_h_pangenome_rle$cumsum <- cumsum(c_0_h_pangenome_rle$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile <- quantile(c_0_h_pangenome_rle$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_0_h_pangenome_rle_core_a1 <- subset(c_0_h_pangenome_rle, cumsum <= abund_quantile[1])
# remove non-numeric columns
c_0_h_pangenome_rle_core_a1$cumsum <- NULL
c_0_h_pangenome_rle_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_0_h_pangenome_rle_core_a2 <- subset(c_0_h_pangenome_rle, cumsum <= abund_quantile[2])
# remove non-numeric columns
c_0_h_pangenome_rle_core_a2$cumsum <- NULL
c_0_h_pangenome_rle_core_a2$abundance <- NULL
# # subset core species based on third defined abundance quantile
c_0_h_pangenome_rle_core_a3 <- subset(c_0_h_pangenome_rle, cumsum <= abund_quantile[3])
# remove non-numeric columns
c_0_h_pangenome_rle_core_a3$cumsum <- NULL
c_0_h_pangenome_rle_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_0_h_pangenome_rle_rare_a1 <- subset(c_0_h_pangenome_rle, cumsum > abund_quantile[1])
# remove non-numeric columns
c_0_h_pangenome_rle_rare_a1$cumsum <- NULL
c_0_h_pangenome_rle_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_0_h_pangenome_rle_rare_a2 <- subset(c_0_h_pangenome_rle, cumsum > abund_quantile[2])
# remove non-numeric columns
c_0_h_pangenome_rle_rare_a2$cumsum <- NULL
c_0_h_pangenome_rle_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_0_h_pangenome_rle_rare_a3 <- subset(c_0_h_pangenome_rle, cumsum > abund_quantile[3])
# remove non-numeric columns
c_0_h_pangenome_rle_rare_a3$cumsum <- NULL
c_0_h_pangenome_rle_rare_a3$abundance <- NULL


# Age group (1-3 years), healthy, pangenome database, rle-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_1_3_h_pangenome_rle <- select(ds_pangenome_rle, rownames(md_1_3_h))
# remove species rows that sum to zero
c_1_3_h_pangenome_rle <- c_1_3_h_pangenome_rle[rowSums(c_1_3_h_pangenome_rle[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_1_3_pangenome_rle = colSums(c_1_3_h_pangenome_rle)
# obtain total sum
sum_all_healthy_1_3_pangenome_rle = sum(sum_healthy_1_3_pangenome_rle)
# add abundance column
c_1_3_h_pangenome_rle$abundance <- (rowSums(c_1_3_h_pangenome_rle[,1:ncol(c_1_3_h_pangenome_rle)]) / sum_all_healthy_1_3_pangenome_rle) * 100
# sort abundance decreasing
c_1_3_h_pangenome_rle <- c_1_3_h_pangenome_rle[with(c_1_3_h_pangenome_rle, order(-abundance)), ]
# obtain cumulative sum
c_1_3_h_pangenome_rle$cumsum <- cumsum(c_1_3_h_pangenome_rle$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_1_3_h <- quantile(c_1_3_h_pangenome_rle$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_1_3_h_pangenome_rle_core_a1 <- subset(c_1_3_h_pangenome_rle, cumsum <= abund_quantile_1_3_h[1])
# remove non-numeric columns
c_1_3_h_pangenome_rle_core_a1$cumsum <- NULL
c_1_3_h_pangenome_rle_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_1_3_h_pangenome_rle_core_a2 <- subset(c_1_3_h_pangenome_rle, cumsum <= abund_quantile_1_3_h[2])
# remove non-numeric columns
c_1_3_h_pangenome_rle_core_a2$cumsum <- NULL
c_1_3_h_pangenome_rle_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_1_3_h_pangenome_rle_core_a3 <- subset(c_1_3_h_pangenome_rle, cumsum <= abund_quantile_1_3_h[3])
# remove non-numeric columns
c_1_3_h_pangenome_rle_core_a3$cumsum <- NULL
c_1_3_h_pangenome_rle_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_1_3_h_pangenome_rle_rare_a1 <- subset(c_1_3_h_pangenome_rle, cumsum > abund_quantile_1_3_h[1])
# remove non-numeric columns
c_1_3_h_pangenome_rle_rare_a1$cumsum <- NULL
c_1_3_h_pangenome_rle_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_1_3_h_pangenome_rle_rare_a2 <- subset(c_1_3_h_pangenome_rle, cumsum > abund_quantile_1_3_h[2])
# remove non-numeric columns
c_1_3_h_pangenome_rle_rare_a2$cumsum <- NULL
c_1_3_h_pangenome_rle_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_1_3_h_pangenome_rle_rare_a3 <- subset(c_1_3_h_pangenome_rle, cumsum > abund_quantile_1_3_h[3])
# remove non-numeric columns
c_1_3_h_pangenome_rle_rare_a3$cumsum <- NULL
c_1_3_h_pangenome_rle_rare_a3$abundance <- NULL


# Age group (4-6 years), healthy, pangenome database, rle-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_4_6_h_pangenome_rle <- select(ds_pangenome_rle, rownames(md_4_6_h))
# remove species rows that sum to zero
c_4_6_h_pangenome_rle <- c_4_6_h_pangenome_rle[rowSums(c_4_6_h_pangenome_rle[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_4_6_pangenome_rle = colSums(c_4_6_h_pangenome_rle)
# get total sum
sum_all_healthy_4_6_pangenome_rle = sum(sum_healthy_4_6_pangenome_rle)
# add abundance column
c_4_6_h_pangenome_rle$abundance <- (rowSums(c_4_6_h_pangenome_rle[,1:ncol(c_4_6_h_pangenome_rle)]) / sum_all_healthy_4_6_pangenome_rle) * 100
# sort abundance decreasing
c_4_6_h_pangenome_rle <- c_4_6_h_pangenome_rle[with(c_4_6_h_pangenome_rle, order(-abundance)), ]
# get cumulative sum
c_4_6_h_pangenome_rle$cumsum <- cumsum(c_4_6_h_pangenome_rle$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_4_6_h <- quantile(c_4_6_h_pangenome_rle$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_4_6_h_pangenome_rle_core_a1 <- subset(c_4_6_h_pangenome_rle, cumsum <= abund_quantile_4_6_h[1])
# remove non-numeric columns
c_4_6_h_pangenome_rle_core_a1$cumsum <- NULL
c_4_6_h_pangenome_rle_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_4_6_h_pangenome_rle_core_a2 <- subset(c_4_6_h_pangenome_rle, cumsum <= abund_quantile_4_6_h[2])
# remove non-numeric columns
c_4_6_h_pangenome_rle_core_a2$cumsum <- NULL
c_4_6_h_pangenome_rle_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_4_6_h_pangenome_rle_core_a3 <- subset(c_4_6_h_pangenome_rle, cumsum <= abund_quantile_4_6_h[3])
# remove non-numeric columns
c_4_6_h_pangenome_rle_core_a3$cumsum <- NULL
c_4_6_h_pangenome_rle_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_4_6_h_pangenome_rle_rare_a1 <- subset(c_4_6_h_pangenome_rle, cumsum > abund_quantile_4_6_h[1])
# remove non-numeric columns
c_4_6_h_pangenome_rle_rare_a1$cumsum <- NULL
c_4_6_h_pangenome_rle_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_4_6_h_pangenome_rle_rare_a2 <- subset(c_4_6_h_pangenome_rle, cumsum > abund_quantile_4_6_h[2])
# remove non-numeric columns
c_4_6_h_pangenome_rle_rare_a2$cumsum <- NULL
c_4_6_h_pangenome_rle_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_4_6_h_pangenome_rle_rare_a3 <- subset(c_4_6_h_pangenome_rle, cumsum > abund_quantile_4_6_h[3])
# remove non-numeric columns
c_4_6_h_pangenome_rle_rare_a3$cumsum <- NULL
c_4_6_h_pangenome_rle_rare_a3$abundance <- NULL


# Age group (0 years), CF, pangenome database, rle-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_0_cf_pangenome_rle <- select(ds_pangenome_rle, rownames(md_0_cf))
# remove species rows that sum to zero
c_0_cf_pangenome_rle <- c_0_cf_pangenome_rle[rowSums(c_0_cf_pangenome_rle[, -1])>0, ]
# sum species counts per CF child
sum_cf_0_pangenome_rle = colSums(c_0_cf_pangenome_rle)
# get total sum
sum_all_cf_0_pangenome_rle = sum(sum_cf_0_pangenome_rle)
# add abundance column
c_0_cf_pangenome_rle$abundance <- (rowSums(c_0_cf_pangenome_rle[,1:ncol(c_0_cf_pangenome_rle)]) / sum_all_cf_0_pangenome_rle) * 100
# sort abundance decreasing
c_0_cf_pangenome_rle <- c_0_cf_pangenome_rle[with(c_0_cf_pangenome_rle, order(-abundance)), ]
# add column with cumulative sum
c_0_cf_pangenome_rle$cumsum <- cumsum(c_0_cf_pangenome_rle$abundance)

# create sample quantiles in dataset
abund_quantile_0_cf <- quantile(c_0_cf_pangenome_rle$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_0_cf_pangenome_rle_core_a1 <- subset(c_0_cf_pangenome_rle, cumsum <= abund_quantile_0_cf[1])
# remove non-numeric columns
c_0_cf_pangenome_rle_core_a1$cumsum <- NULL
c_0_cf_pangenome_rle_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_0_cf_pangenome_rle_core_a2 <- subset(c_0_cf_pangenome_rle, cumsum <= abund_quantile_0_cf[2])
# remove non-numeric columns
c_0_cf_pangenome_rle_core_a2$cumsum <- NULL
c_0_cf_pangenome_rle_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_0_cf_pangenome_rle_core_a3 <- subset(c_0_cf_pangenome_rle, cumsum <= abund_quantile_0_cf[3])
# remove non-numeric columns
c_0_cf_pangenome_rle_core_a3$cumsum <- NULL
c_0_cf_pangenome_rle_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_0_cf_pangenome_rle_rare_a1 <- subset(c_0_cf_pangenome_rle, cumsum > abund_quantile_0_cf[1])
# remove non-numeric columns
c_0_cf_pangenome_rle_rare_a1$cumsum <- NULL
c_0_cf_pangenome_rle_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_0_cf_pangenome_rle_rare_a2 <- subset(c_0_cf_pangenome_rle, cumsum > abund_quantile_0_cf[2])
# remove non-numeric columns
c_0_cf_pangenome_rle_rare_a2$cumsum <- NULL
c_0_cf_pangenome_rle_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_0_cf_pangenome_rle_rare_a3 <- subset(c_0_cf_pangenome_rle, cumsum > abund_quantile_0_cf[3])
# remove non-numeric columns
c_0_cf_pangenome_rle_rare_a3$cumsum <- NULL
c_0_cf_pangenome_rle_rare_a3$abundance <- NULL


# Age group (1-3 years), CF, pangenome database, rle-normalised
c_1_3_cf_pangenome_rle <- select(ds_pangenome_rle, rownames(md_1_3_cf))
# remove species rows that sum to zero
c_1_3_cf_pangenome_rle <- c_1_3_cf_pangenome_rle[rowSums(c_1_3_cf_pangenome_rle[, -1])>0, ]
# sum species counts per CF child
sum_cf_1_3_pangenome_rle = colSums(c_1_3_cf_pangenome_rle)
# get total sum
sum_all_cf_1_3_pangenome_rle = sum(sum_cf_1_3_pangenome_rle)
# add abundance column
c_1_3_cf_pangenome_rle$abundance <- (rowSums(c_1_3_cf_pangenome_rle[,1:ncol(c_1_3_cf_pangenome_rle)]) / sum_all_cf_1_3_pangenome_rle) * 100
# sort abundance decreasing
c_1_3_cf_pangenome_rle <- c_1_3_cf_pangenome_rle[with(c_1_3_cf_pangenome_rle, order(-abundance)), ]
# add column with cumulative sum
c_1_3_cf_pangenome_rle$cumsum <- cumsum(c_1_3_cf_pangenome_rle$abundance)

# table with abundant species
# create sample quantiles in dataset
abund_quantile_1_3_cf <- quantile(c_1_3_cf_pangenome_rle$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_1_3_cf_pangenome_rle_core_a1 <- subset(c_1_3_cf_pangenome_rle, cumsum <= abund_quantile_1_3_cf[1])
# remove non-numeric columns
c_1_3_cf_pangenome_rle_core_a1$cumsum <- NULL
c_1_3_cf_pangenome_rle_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_1_3_cf_pangenome_rle_core_a2 <- subset(c_1_3_cf_pangenome_rle, cumsum <= abund_quantile_1_3_cf[2])
# remove non-numeric columns
c_1_3_cf_pangenome_rle_core_a2$cumsum <- NULL
c_1_3_cf_pangenome_rle_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_1_3_cf_pangenome_rle_core_a3 <- subset(c_1_3_cf_pangenome_rle, cumsum <= abund_quantile_1_3_cf[3])
# remove non-numeric columns
c_1_3_cf_pangenome_rle_core_a3$cumsum <- NULL
c_1_3_cf_pangenome_rle_core_a3$abundance <- NULL


# table with rare species
# subset rare species based on first defined abundance quantile
c_1_3_cf_pangenome_rle_rare_a1 <- subset(c_1_3_cf_pangenome_rle, cumsum > abund_quantile_1_3_cf[1])
# remove non-numeric columns
c_1_3_cf_pangenome_rle_rare_a1$cumsum <- NULL
c_1_3_cf_pangenome_rle_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_1_3_cf_pangenome_rle_rare_a2 <- subset(c_1_3_cf_pangenome_rle, cumsum > abund_quantile_1_3_cf[2])
# remove non-numeric columns
c_1_3_cf_pangenome_rle_rare_a2$cumsum <- NULL
c_1_3_cf_pangenome_rle_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_1_3_cf_pangenome_rle_rare_a3 <- subset(c_1_3_cf_pangenome_rle, cumsum > abund_quantile_1_3_cf[3])
# remove non-numeric columns
c_1_3_cf_pangenome_rle_rare_a3$cumsum <- NULL
c_1_3_cf_pangenome_rle_rare_a3$abundance <- NULL


# Age group (4-6 years), CF, pangenome database, rle-normalised
c_4_6_cf_pangenome_rle <- select(ds_pangenome_rle, rownames(md_4_6_cf))
# remove species rows that sum to zero
c_4_6_cf_pangenome_rle <- c_4_6_cf_pangenome_rle[rowSums(c_4_6_cf_pangenome_rle[, -1])>0, ]
# sum species counts per CF child
sum_cf_4_6_pangenome_rle = colSums(c_4_6_cf_pangenome_rle)
# get total sum
sum_all_cf_4_6_pangenome_rle = sum(sum_cf_4_6_pangenome_rle)
# add abundance column
c_4_6_cf_pangenome_rle$abundance <- (rowSums(c_4_6_cf_pangenome_rle[,1:ncol(c_4_6_cf_pangenome_rle)]) / sum_all_cf_4_6_pangenome_rle) * 100
# sort abundance decreasing
c_4_6_cf_pangenome_rle <- c_4_6_cf_pangenome_rle[with(c_4_6_cf_pangenome_rle, order(-abundance)), ]
# add column with cumulative sum
c_4_6_cf_pangenome_rle$cumsum <- cumsum(c_4_6_cf_pangenome_rle$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_4_6_cf <- quantile(c_4_6_cf_pangenome_rle$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_4_6_cf_pangenome_rle_core_a1 <- subset(c_4_6_cf_pangenome_rle, cumsum <= abund_quantile_4_6_cf[1])
# remove non-numeric columns
c_4_6_cf_pangenome_rle_core_a1$cumsum <- NULL
c_4_6_cf_pangenome_rle_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_4_6_cf_pangenome_rle_core_a2 <- subset(c_4_6_cf_pangenome_rle, cumsum <= abund_quantile_4_6_cf[2])
# remove non-numeric columns
c_4_6_cf_pangenome_rle_core_a2$cumsum <- NULL
c_4_6_cf_pangenome_rle_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_4_6_cf_pangenome_rle_core_a3 <- subset(c_4_6_cf_pangenome_rle, cumsum <= abund_quantile_4_6_cf[3])
# remove non-numeric columns
c_4_6_cf_pangenome_rle_core_a3$cumsum <- NULL
c_4_6_cf_pangenome_rle_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_4_6_cf_pangenome_rle_rare_a1 <- subset(c_4_6_cf_pangenome_rle, cumsum > abund_quantile_4_6_cf[1])
# remove non-numeric columns
c_4_6_cf_pangenome_rle_rare_a1$cumsum <- NULL
c_4_6_cf_pangenome_rle_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_4_6_cf_pangenome_rle_rare_a2 <- subset(c_4_6_cf_pangenome_rle, cumsum > abund_quantile_4_6_cf[2])
# remove non-numeric columns
c_4_6_cf_pangenome_rle_rare_a2$cumsum <- NULL
c_4_6_cf_pangenome_rle_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_4_6_cf_pangenome_rle_rare_a3 <- subset(c_4_6_cf_pangenome_rle, cumsum > abund_quantile_4_6_cf[3])
# remove non-numeric columns
c_4_6_cf_pangenome_rle_rare_a3$cumsum <- NULL
c_4_6_cf_pangenome_rle_rare_a3$abundance <- NULL


# create vector with all rare species (15th abundance percentile, pangenome, healthy)
healthy_rare_a1_pangenome_rle <- c(rownames(c_0_h_pangenome_rle_rare_a1), 
                               rownames(c_1_3_h_pangenome_rle_rare_a1), 
                               rownames(c_4_6_h_pangenome_rle_rare_a1))
# remove duplicate entries
healthy_rare_a1_pangenome_rle <- healthy_rare_a1_pangenome_rle[!duplicated(healthy_rare_a1_pangenome_rle)]

# create vector with all core species (15th abundance percentile, pangenome, healthy)
healthy_core_a1_pangenome_rle <- c(rownames(c_0_h_pangenome_rle_core_a1), 
                               rownames(c_1_3_h_pangenome_rle_core_a1), 
                               rownames(c_4_6_h_pangenome_rle_core_a1))
# remove duplicate entries
healthy_core_a1_pangenome_rle <- healthy_core_a1_pangenome_rle[!duplicated(healthy_core_a1_pangenome_rle)]

# create vector with all rare species (25th abundance percentile, pangenome, healthy)
healthy_rare_a2_pangenome_rle <- c(rownames(c_0_h_pangenome_rle_rare_a2), 
                               rownames(c_1_3_h_pangenome_rle_rare_a2), 
                               rownames(c_4_6_h_pangenome_rle_rare_a2))
# remove duplicate entries
healthy_rare_a2_pangenome_rle <- healthy_rare_a2_pangenome_rle[!duplicated(healthy_rare_a2_pangenome_rle)]

# create vector with all core species (25th abundance percentile, pangenome, healthy)
healthy_core_a2_pangenome_rle <- c(rownames(c_0_h_pangenome_rle_core_a2), 
                               rownames(c_1_3_h_pangenome_rle_core_a2), 
                               rownames(c_4_6_h_pangenome_rle_core_a2))
# remove duplicate entries
healthy_core_a2_pangenome_rle <- healthy_core_a2_pangenome_rle[!duplicated(healthy_core_a2_pangenome_rle)]

# create vector with all rare species (35th abundance percentile, pangenome, healthy)
healthy_rare_a3_pangenome_rle <- c(rownames(c_0_h_pangenome_rle_rare_a3), 
                               rownames(c_1_3_h_pangenome_rle_rare_a3), 
                               rownames(c_4_6_h_pangenome_rle_rare_a3))
# remove duplicate entries
healthy_rare_a3_pangenome_rle <- healthy_rare_a3_pangenome_rle[!duplicated(healthy_rare_a3_pangenome_rle)]

# create vector with all core species (35th abundance percentile, pangenome, healthy)
healthy_core_a3_pangenome_rle <- c(rownames(c_0_h_pangenome_rle_core_a3), 
                               rownames(c_1_3_h_pangenome_rle_core_a3), 
                               rownames(c_4_6_h_pangenome_rle_core_a3))
# remove duplicate entries
healthy_core_a3_pangenome_rle <- healthy_core_a3_pangenome_rle[!duplicated(healthy_core_a3_pangenome_rle)]


# create vector with all rare species (15th abundance percentile, pangenome, CF)
cf_rare_a1_pangenome_rle <- c(rownames(c_0_cf_pangenome_rle_rare_a1), 
                          rownames(c_1_3_cf_pangenome_rle_rare_a1), 
                          rownames(c_4_6_cf_pangenome_rle_rare_a1))
# remove duplicate entries
cf_rare_a1_pangenome_rle <- cf_rare_a1_pangenome_rle[!duplicated(cf_rare_a1_pangenome_rle)]

# create vector with all core species (15th abundance percentile, pangenome, CF)
cf_core_a1_pangenome_rle <- c(rownames(c_0_cf_pangenome_rle_core_a1), 
                          rownames(c_1_3_cf_pangenome_rle_core_a1), 
                          rownames(c_4_6_cf_pangenome_rle_core_a1))
# remove duplicate entries
cf_core_a1_pangenome_rle <- cf_core_a1_pangenome_rle[!duplicated(cf_core_a1_pangenome_rle)]

# create vector with all rare species (25th abundance percentile, pangenome, CF)
cf_rare_a2_pangenome_rle <- c(rownames(c_0_cf_pangenome_rle_rare_a2), 
                          rownames(c_1_3_cf_pangenome_rle_rare_a2), 
                          rownames(c_4_6_cf_pangenome_rle_rare_a2))
# remove duplicate entries
cf_rare_a2_pangenome_rle <- cf_rare_a2_pangenome_rle[!duplicated(cf_rare_a2_pangenome_rle)]

# create vector with all core species (25th abundance percentile, pangenome, CF)
cf_core_a2_pangenome_rle <- c(rownames(c_0_cf_pangenome_rle_core_a2), 
                          rownames(c_1_3_cf_pangenome_rle_core_a2), 
                          rownames(c_4_6_cf_pangenome_rle_core_a2))
# remove duplicate entries
cf_core_a2_pangenome_rle <- cf_core_a2_pangenome_rle[!duplicated(cf_core_a2_pangenome_rle)]

# create vector with all rare species (35th abundance percentile, pangenome, CF)
cf_rare_a3_pangenome_rle <- c(rownames(c_0_cf_pangenome_rle_rare_a3), 
                          rownames(c_1_3_cf_pangenome_rle_rare_a3), 
                          rownames(c_4_6_cf_pangenome_rle_rare_a3))
# remove duplicate entries
cf_rare_a3_pangenome_rle <- cf_rare_a3_pangenome_rle[!duplicated(cf_rare_a3_pangenome_rle)]

# create vector with all core species (35th abundance percentile, pangenome, CF)
cf_core_a3_pangenome_rle <- c(rownames(c_0_cf_pangenome_rle_core_a3), 
                          rownames(c_1_3_cf_pangenome_rle_core_a3), 
                          rownames(c_4_6_cf_pangenome_rle_core_a3))
# remove duplicate entries
cf_core_a3_pangenome_rle <- cf_core_a3_pangenome_rle[!duplicated(cf_core_a3_pangenome_rle)]


# age group (0 years), healthy, pangenome database, vst-normalised
c_0_h_pangenome_vst <- select(ds_pangenome_vst, rownames(md_0_h))
# remove species rows that sum to zero
c_0_h_pangenome_vst <- c_0_h_pangenome_vst[rowSums(c_0_h_pangenome_vst[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_0_pangenome_vst = colSums(c_0_h_pangenome_vst)
# total sum
sum_all_healthy_0_pangenome_vst = sum(sum_healthy_0_pangenome_vst)
# add abundance column
c_0_h_pangenome_vst$abundance <- (rowSums(c_0_h_pangenome_vst[,1:ncol(c_0_h_pangenome_vst)]) / sum_all_healthy_0_pangenome_vst) * 100
# sort abundance decreasing
c_0_h_pangenome_vst <- c_0_h_pangenome_vst[with(c_0_h_pangenome_vst, order(-abundance)), ]
# obtain cumulative sum
c_0_h_pangenome_vst$cumsum <- cumsum(c_0_h_pangenome_vst$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile <- quantile(c_0_h_pangenome_vst$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_0_h_pangenome_vst_core_a1 <- subset(c_0_h_pangenome_vst, cumsum <= abund_quantile[1])
# remove non-numeric columns
c_0_h_pangenome_vst_core_a1$cumsum <- NULL
c_0_h_pangenome_vst_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_0_h_pangenome_vst_core_a2 <- subset(c_0_h_pangenome_vst, cumsum <= abund_quantile[2])
# remove non-numeric columns
c_0_h_pangenome_vst_core_a2$cumsum <- NULL
c_0_h_pangenome_vst_core_a2$abundance <- NULL
# # subset core species based on third defined abundance quantile
c_0_h_pangenome_vst_core_a3 <- subset(c_0_h_pangenome_vst, cumsum <= abund_quantile[3])
# remove non-numeric columns
c_0_h_pangenome_vst_core_a3$cumsum <- NULL
c_0_h_pangenome_vst_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_0_h_pangenome_vst_rare_a1 <- subset(c_0_h_pangenome_vst, cumsum > abund_quantile[1])
# remove non-numeric columns
c_0_h_pangenome_vst_rare_a1$cumsum <- NULL
c_0_h_pangenome_vst_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_0_h_pangenome_vst_rare_a2 <- subset(c_0_h_pangenome_vst, cumsum > abund_quantile[2])
# remove non-numeric columns
c_0_h_pangenome_vst_rare_a2$cumsum <- NULL
c_0_h_pangenome_vst_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_0_h_pangenome_vst_rare_a3 <- subset(c_0_h_pangenome_vst, cumsum > abund_quantile[3])
# remove non-numeric columns
c_0_h_pangenome_vst_rare_a3$cumsum <- NULL
c_0_h_pangenome_vst_rare_a3$abundance <- NULL


# Age group (1-3 years), healthy, pangenome database, vst-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_1_3_h_pangenome_vst <- select(ds_pangenome_vst, rownames(md_1_3_h))
# remove species rows that sum to zero
c_1_3_h_pangenome_vst <- c_1_3_h_pangenome_vst[rowSums(c_1_3_h_pangenome_vst[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_1_3_pangenome_vst = colSums(c_1_3_h_pangenome_vst)
# obtain total sum
sum_all_healthy_1_3_pangenome_vst = sum(sum_healthy_1_3_pangenome_vst)
# add abundance column
c_1_3_h_pangenome_vst$abundance <- (rowSums(c_1_3_h_pangenome_vst[,1:ncol(c_1_3_h_pangenome_vst)]) / sum_all_healthy_1_3_pangenome_vst) * 100
# sort abundance decreasing
c_1_3_h_pangenome_vst <- c_1_3_h_pangenome_vst[with(c_1_3_h_pangenome_vst, order(-abundance)), ]
# obtain cumulative sum
c_1_3_h_pangenome_vst$cumsum <- cumsum(c_1_3_h_pangenome_vst$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_1_3_h <- quantile(c_1_3_h_pangenome_vst$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_1_3_h_pangenome_vst_core_a1 <- subset(c_1_3_h_pangenome_vst, cumsum <= abund_quantile_1_3_h[1])
# remove non-numeric columns
c_1_3_h_pangenome_vst_core_a1$cumsum <- NULL
c_1_3_h_pangenome_vst_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_1_3_h_pangenome_vst_core_a2 <- subset(c_1_3_h_pangenome_vst, cumsum <= abund_quantile_1_3_h[2])
# remove non-numeric columns
c_1_3_h_pangenome_vst_core_a2$cumsum <- NULL
c_1_3_h_pangenome_vst_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_1_3_h_pangenome_vst_core_a3 <- subset(c_1_3_h_pangenome_vst, cumsum <= abund_quantile_1_3_h[3])
# remove non-numeric columns
c_1_3_h_pangenome_vst_core_a3$cumsum <- NULL
c_1_3_h_pangenome_vst_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_1_3_h_pangenome_vst_rare_a1 <- subset(c_1_3_h_pangenome_vst, cumsum > abund_quantile_1_3_h[1])
# remove non-numeric columns
c_1_3_h_pangenome_vst_rare_a1$cumsum <- NULL
c_1_3_h_pangenome_vst_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_1_3_h_pangenome_vst_rare_a2 <- subset(c_1_3_h_pangenome_vst, cumsum > abund_quantile_1_3_h[2])
# remove non-numeric columns
c_1_3_h_pangenome_vst_rare_a2$cumsum <- NULL
c_1_3_h_pangenome_vst_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_1_3_h_pangenome_vst_rare_a3 <- subset(c_1_3_h_pangenome_vst, cumsum > abund_quantile_1_3_h[3])
# remove non-numeric columns
c_1_3_h_pangenome_vst_rare_a3$cumsum <- NULL
c_1_3_h_pangenome_vst_rare_a3$abundance <- NULL


# Age group (4-6 years), healthy, pangenome database, vst-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_4_6_h_pangenome_vst <- select(ds_pangenome_vst, rownames(md_4_6_h))
# remove species rows that sum to zero
c_4_6_h_pangenome_vst <- c_4_6_h_pangenome_vst[rowSums(c_4_6_h_pangenome_vst[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_4_6_pangenome_vst = colSums(c_4_6_h_pangenome_vst)
# get total sum
sum_all_healthy_4_6_pangenome_vst = sum(sum_healthy_4_6_pangenome_vst)
# add abundance column
c_4_6_h_pangenome_vst$abundance <- (rowSums(c_4_6_h_pangenome_vst[,1:ncol(c_4_6_h_pangenome_vst)]) / sum_all_healthy_4_6_pangenome_vst) * 100
# sort abundance decreasing
c_4_6_h_pangenome_vst <- c_4_6_h_pangenome_vst[with(c_4_6_h_pangenome_vst, order(-abundance)), ]
# get cumulative sum
c_4_6_h_pangenome_vst$cumsum <- cumsum(c_4_6_h_pangenome_vst$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_4_6_h <- quantile(c_4_6_h_pangenome_vst$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_4_6_h_pangenome_vst_core_a1 <- subset(c_4_6_h_pangenome_vst, cumsum <= abund_quantile_4_6_h[1])
# remove non-numeric columns
c_4_6_h_pangenome_vst_core_a1$cumsum <- NULL
c_4_6_h_pangenome_vst_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_4_6_h_pangenome_vst_core_a2 <- subset(c_4_6_h_pangenome_vst, cumsum <= abund_quantile_4_6_h[2])
# remove non-numeric columns
c_4_6_h_pangenome_vst_core_a2$cumsum <- NULL
c_4_6_h_pangenome_vst_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_4_6_h_pangenome_vst_core_a3 <- subset(c_4_6_h_pangenome_vst, cumsum <= abund_quantile_4_6_h[3])
# remove non-numeric columns
c_4_6_h_pangenome_vst_core_a3$cumsum <- NULL
c_4_6_h_pangenome_vst_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_4_6_h_pangenome_vst_rare_a1 <- subset(c_4_6_h_pangenome_vst, cumsum > abund_quantile_4_6_h[1])
# remove non-numeric columns
c_4_6_h_pangenome_vst_rare_a1$cumsum <- NULL
c_4_6_h_pangenome_vst_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_4_6_h_pangenome_vst_rare_a2 <- subset(c_4_6_h_pangenome_vst, cumsum > abund_quantile_4_6_h[2])
# remove non-numeric columns
c_4_6_h_pangenome_vst_rare_a2$cumsum <- NULL
c_4_6_h_pangenome_vst_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_4_6_h_pangenome_vst_rare_a3 <- subset(c_4_6_h_pangenome_vst, cumsum > abund_quantile_4_6_h[3])
# remove non-numeric columns
c_4_6_h_pangenome_vst_rare_a3$cumsum <- NULL
c_4_6_h_pangenome_vst_rare_a3$abundance <- NULL


# Age group (0 years), CF, pangenome database, vst-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_0_cf_pangenome_vst <- select(ds_pangenome_vst, rownames(md_0_cf))
# remove species rows that sum to zero
c_0_cf_pangenome_vst <- c_0_cf_pangenome_vst[rowSums(c_0_cf_pangenome_vst[, -1])>0, ]
# sum species counts per CF child
sum_cf_0_pangenome_vst = colSums(c_0_cf_pangenome_vst)
# get total sum
sum_all_cf_0_pangenome_vst = sum(sum_cf_0_pangenome_vst)
# add abundance column
c_0_cf_pangenome_vst$abundance <- (rowSums(c_0_cf_pangenome_vst[,1:ncol(c_0_cf_pangenome_vst)]) / sum_all_cf_0_pangenome_vst) * 100
# sort abundance decreasing
c_0_cf_pangenome_vst <- c_0_cf_pangenome_vst[with(c_0_cf_pangenome_vst, order(-abundance)), ]
# add column with cumulative sum
c_0_cf_pangenome_vst$cumsum <- cumsum(c_0_cf_pangenome_vst$abundance)

# create sample quantiles in dataset
abund_quantile_0_cf <- quantile(c_0_cf_pangenome_vst$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_0_cf_pangenome_vst_core_a1 <- subset(c_0_cf_pangenome_vst, cumsum <= abund_quantile_0_cf[1])
# remove non-numeric columns
c_0_cf_pangenome_vst_core_a1$cumsum <- NULL
c_0_cf_pangenome_vst_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_0_cf_pangenome_vst_core_a2 <- subset(c_0_cf_pangenome_vst, cumsum <= abund_quantile_0_cf[2])
# remove non-numeric columns
c_0_cf_pangenome_vst_core_a2$cumsum <- NULL
c_0_cf_pangenome_vst_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_0_cf_pangenome_vst_core_a3 <- subset(c_0_cf_pangenome_vst, cumsum <= abund_quantile_0_cf[3])
# remove non-numeric columns
c_0_cf_pangenome_vst_core_a3$cumsum <- NULL
c_0_cf_pangenome_vst_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_0_cf_pangenome_vst_rare_a1 <- subset(c_0_cf_pangenome_vst, cumsum > abund_quantile_0_cf[1])
# remove non-numeric columns
c_0_cf_pangenome_vst_rare_a1$cumsum <- NULL
c_0_cf_pangenome_vst_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_0_cf_pangenome_vst_rare_a2 <- subset(c_0_cf_pangenome_vst, cumsum > abund_quantile_0_cf[2])
# remove non-numeric columns
c_0_cf_pangenome_vst_rare_a2$cumsum <- NULL
c_0_cf_pangenome_vst_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_0_cf_pangenome_vst_rare_a3 <- subset(c_0_cf_pangenome_vst, cumsum > abund_quantile_0_cf[3])
# remove non-numeric columns
c_0_cf_pangenome_vst_rare_a3$cumsum <- NULL
c_0_cf_pangenome_vst_rare_a3$abundance <- NULL


# Age group (1-3 years), CF, pangenome database, vst-normalised
c_1_3_cf_pangenome_vst <- select(ds_pangenome_vst, rownames(md_1_3_cf))
# remove species rows that sum to zero
c_1_3_cf_pangenome_vst <- c_1_3_cf_pangenome_vst[rowSums(c_1_3_cf_pangenome_vst[, -1])>0, ]
# sum species counts per CF child
sum_cf_1_3_pangenome_vst = colSums(c_1_3_cf_pangenome_vst)
# get total sum
sum_all_cf_1_3_pangenome_vst = sum(sum_cf_1_3_pangenome_vst)
# add abundance column
c_1_3_cf_pangenome_vst$abundance <- (rowSums(c_1_3_cf_pangenome_vst[,1:ncol(c_1_3_cf_pangenome_vst)]) / sum_all_cf_1_3_pangenome_vst) * 100
# sort abundance decreasing
c_1_3_cf_pangenome_vst <- c_1_3_cf_pangenome_vst[with(c_1_3_cf_pangenome_vst, order(-abundance)), ]
# add column with cumulative sum
c_1_3_cf_pangenome_vst$cumsum <- cumsum(c_1_3_cf_pangenome_vst$abundance)

# table with abundant species
# create sample quantiles in dataset
abund_quantile_1_3_cf <- quantile(c_1_3_cf_pangenome_vst$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_1_3_cf_pangenome_vst_core_a1 <- subset(c_1_3_cf_pangenome_vst, cumsum <= abund_quantile_1_3_cf[1])
# remove non-numeric columns
c_1_3_cf_pangenome_vst_core_a1$cumsum <- NULL
c_1_3_cf_pangenome_vst_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_1_3_cf_pangenome_vst_core_a2 <- subset(c_1_3_cf_pangenome_vst, cumsum <= abund_quantile_1_3_cf[2])
# remove non-numeric columns
c_1_3_cf_pangenome_vst_core_a2$cumsum <- NULL
c_1_3_cf_pangenome_vst_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_1_3_cf_pangenome_vst_core_a3 <- subset(c_1_3_cf_pangenome_vst, cumsum <= abund_quantile_1_3_cf[3])
# remove non-numeric columns
c_1_3_cf_pangenome_vst_core_a3$cumsum <- NULL
c_1_3_cf_pangenome_vst_core_a3$abundance <- NULL


# table with rare species
# subset rare species based on first defined abundance quantile
c_1_3_cf_pangenome_vst_rare_a1 <- subset(c_1_3_cf_pangenome_vst, cumsum > abund_quantile_1_3_cf[1])
# remove non-numeric columns
c_1_3_cf_pangenome_vst_rare_a1$cumsum <- NULL
c_1_3_cf_pangenome_vst_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_1_3_cf_pangenome_vst_rare_a2 <- subset(c_1_3_cf_pangenome_vst, cumsum > abund_quantile_1_3_cf[2])
# remove non-numeric columns
c_1_3_cf_pangenome_vst_rare_a2$cumsum <- NULL
c_1_3_cf_pangenome_vst_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_1_3_cf_pangenome_vst_rare_a3 <- subset(c_1_3_cf_pangenome_vst, cumsum > abund_quantile_1_3_cf[3])
# remove non-numeric columns
c_1_3_cf_pangenome_vst_rare_a3$cumsum <- NULL
c_1_3_cf_pangenome_vst_rare_a3$abundance <- NULL


# Age group (4-6 years), CF, pangenome database, vst-normalised
c_4_6_cf_pangenome_vst <- select(ds_pangenome_vst, rownames(md_4_6_cf))
# remove species rows that sum to zero
c_4_6_cf_pangenome_vst <- c_4_6_cf_pangenome_vst[rowSums(c_4_6_cf_pangenome_vst[, -1])>0, ]
# sum species counts per CF child
sum_cf_4_6_pangenome_vst = colSums(c_4_6_cf_pangenome_vst)
# get total sum
sum_all_cf_4_6_pangenome_vst = sum(sum_cf_4_6_pangenome_vst)
# add abundance column
c_4_6_cf_pangenome_vst$abundance <- (rowSums(c_4_6_cf_pangenome_vst[,1:ncol(c_4_6_cf_pangenome_vst)]) / sum_all_cf_4_6_pangenome_vst) * 100
# sort abundance decreasing
c_4_6_cf_pangenome_vst <- c_4_6_cf_pangenome_vst[with(c_4_6_cf_pangenome_vst, order(-abundance)), ]
# add column with cumulative sum
c_4_6_cf_pangenome_vst$cumsum <- cumsum(c_4_6_cf_pangenome_vst$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_4_6_cf <- quantile(c_4_6_cf_pangenome_vst$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_4_6_cf_pangenome_vst_core_a1 <- subset(c_4_6_cf_pangenome_vst, cumsum <= abund_quantile_4_6_cf[1])
# remove non-numeric columns
c_4_6_cf_pangenome_vst_core_a1$cumsum <- NULL
c_4_6_cf_pangenome_vst_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_4_6_cf_pangenome_vst_core_a2 <- subset(c_4_6_cf_pangenome_vst, cumsum <= abund_quantile_4_6_cf[2])
# remove non-numeric columns
c_4_6_cf_pangenome_vst_core_a2$cumsum <- NULL
c_4_6_cf_pangenome_vst_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_4_6_cf_pangenome_vst_core_a3 <- subset(c_4_6_cf_pangenome_vst, cumsum <= abund_quantile_4_6_cf[3])
# remove non-numeric columns
c_4_6_cf_pangenome_vst_core_a3$cumsum <- NULL
c_4_6_cf_pangenome_vst_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_4_6_cf_pangenome_vst_rare_a1 <- subset(c_4_6_cf_pangenome_vst, cumsum > abund_quantile_4_6_cf[1])
# remove non-numeric columns
c_4_6_cf_pangenome_vst_rare_a1$cumsum <- NULL
c_4_6_cf_pangenome_vst_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_4_6_cf_pangenome_vst_rare_a2 <- subset(c_4_6_cf_pangenome_vst, cumsum > abund_quantile_4_6_cf[2])
# remove non-numeric columns
c_4_6_cf_pangenome_vst_rare_a2$cumsum <- NULL
c_4_6_cf_pangenome_vst_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_4_6_cf_pangenome_vst_rare_a3 <- subset(c_4_6_cf_pangenome_vst, cumsum > abund_quantile_4_6_cf[3])
# remove non-numeric columns
c_4_6_cf_pangenome_vst_rare_a3$cumsum <- NULL
c_4_6_cf_pangenome_vst_rare_a3$abundance <- NULL


# create vector with all rare species (15th abundance percentile, pangenome, healthy)
healthy_rare_a1_pangenome_vst <- c(rownames(c_0_h_pangenome_vst_rare_a1), 
                                   rownames(c_1_3_h_pangenome_vst_rare_a1), 
                                   rownames(c_4_6_h_pangenome_vst_rare_a1))
# remove duplicate entries
healthy_rare_a1_pangenome_vst <- healthy_rare_a1_pangenome_vst[!duplicated(healthy_rare_a1_pangenome_vst)]

# create vector with all core species (15th abundance percentile, pangenome, healthy)
healthy_core_a1_pangenome_vst <- c(rownames(c_0_h_pangenome_vst_core_a1), 
                                   rownames(c_1_3_h_pangenome_vst_core_a1), 
                                   rownames(c_4_6_h_pangenome_vst_core_a1))
# remove duplicate entries
healthy_core_a1_pangenome_vst <- healthy_core_a1_pangenome_vst[!duplicated(healthy_core_a1_pangenome_vst)]

# create vector with all rare species (25th abundance percentile, pangenome, healthy)
healthy_rare_a2_pangenome_vst <- c(rownames(c_0_h_pangenome_vst_rare_a2), 
                                   rownames(c_1_3_h_pangenome_vst_rare_a2), 
                                   rownames(c_4_6_h_pangenome_vst_rare_a2))
# remove duplicate entries
healthy_rare_a2_pangenome_vst <- healthy_rare_a2_pangenome_vst[!duplicated(healthy_rare_a2_pangenome_vst)]

# create vector with all core species (25th abundance percentile, pangenome, healthy)
healthy_core_a2_pangenome_vst <- c(rownames(c_0_h_pangenome_vst_core_a2), 
                                   rownames(c_1_3_h_pangenome_vst_core_a2), 
                                   rownames(c_4_6_h_pangenome_vst_core_a2))
# remove duplicate entries
healthy_core_a2_pangenome_vst <- healthy_core_a2_pangenome_vst[!duplicated(healthy_core_a2_pangenome_vst)]

# create vector with all rare species (35th abundance percentile, pangenome, healthy)
healthy_rare_a3_pangenome_vst <- c(rownames(c_0_h_pangenome_vst_rare_a3), 
                                   rownames(c_1_3_h_pangenome_vst_rare_a3), 
                                   rownames(c_4_6_h_pangenome_vst_rare_a3))
# remove duplicate entries
healthy_rare_a3_pangenome_vst <- healthy_rare_a3_pangenome_vst[!duplicated(healthy_rare_a3_pangenome_vst)]

# create vector with all core species (35th abundance percentile, pangenome, healthy)
healthy_core_a3_pangenome_vst <- c(rownames(c_0_h_pangenome_vst_core_a3), 
                                   rownames(c_1_3_h_pangenome_vst_core_a3), 
                                   rownames(c_4_6_h_pangenome_vst_core_a3))
# remove duplicate entries
healthy_core_a3_pangenome_vst <- healthy_core_a3_pangenome_vst[!duplicated(healthy_core_a3_pangenome_vst)]


# create vector with all rare species (15th abundance percentile, pangenome, CF)
cf_rare_a1_pangenome_vst <- c(rownames(c_0_cf_pangenome_vst_rare_a1), 
                              rownames(c_1_3_cf_pangenome_vst_rare_a1), 
                              rownames(c_4_6_cf_pangenome_vst_rare_a1))
# remove duplicate entries
cf_rare_a1_pangenome_vst <- cf_rare_a1_pangenome_vst[!duplicated(cf_rare_a1_pangenome_vst)]

# create vector with all core species (15th abundance percentile, pangenome, CF)
cf_core_a1_pangenome_vst <- c(rownames(c_0_cf_pangenome_vst_core_a1), 
                              rownames(c_1_3_cf_pangenome_vst_core_a1), 
                              rownames(c_4_6_cf_pangenome_vst_core_a1))
# remove duplicate entries
cf_core_a1_pangenome_vst <- cf_core_a1_pangenome_vst[!duplicated(cf_core_a1_pangenome_vst)]

# create vector with all rare species (25th abundance percentile, pangenome, CF)
cf_rare_a2_pangenome_vst <- c(rownames(c_0_cf_pangenome_vst_rare_a2), 
                              rownames(c_1_3_cf_pangenome_vst_rare_a2), 
                              rownames(c_4_6_cf_pangenome_vst_rare_a2))
# remove duplicate entries
cf_rare_a2_pangenome_vst <- cf_rare_a2_pangenome_vst[!duplicated(cf_rare_a2_pangenome_vst)]

# create vector with all core species (25th abundance percentile, pangenome, CF)
cf_core_a2_pangenome_vst <- c(rownames(c_0_cf_pangenome_vst_core_a2), 
                              rownames(c_1_3_cf_pangenome_vst_core_a2), 
                              rownames(c_4_6_cf_pangenome_vst_core_a2))
# remove duplicate entries
cf_core_a2_pangenome_vst <- cf_core_a2_pangenome_vst[!duplicated(cf_core_a2_pangenome_vst)]

# create vector with all rare species (35th abundance percentile, pangenome, CF)
cf_rare_a3_pangenome_vst <- c(rownames(c_0_cf_pangenome_vst_rare_a3), 
                              rownames(c_1_3_cf_pangenome_vst_rare_a3), 
                              rownames(c_4_6_cf_pangenome_vst_rare_a3))
# remove duplicate entries
cf_rare_a3_pangenome_vst <- cf_rare_a3_pangenome_vst[!duplicated(cf_rare_a3_pangenome_vst)]

# create vector with all core species (35th abundance percentile, pangenome, CF)
cf_core_a3_pangenome_vst <- c(rownames(c_0_cf_pangenome_vst_core_a3), 
                              rownames(c_1_3_cf_pangenome_vst_core_a3), 
                              rownames(c_4_6_cf_pangenome_vst_core_a3))
# remove duplicate entries
cf_core_a3_pangenome_vst <- cf_core_a3_pangenome_vst[!duplicated(cf_core_a3_pangenome_vst)]

############################################################################################################
# One-strain per species
# define core and rare species based on three rarity thresholds (defined above in variable abund_quantile)
# this is repeated for all three normalisations steps (bcphc, rle, vst)

# age group (0 years), healthy, pangenome database, bcphc-normalised
c_0_h_osps <- select(ds_osps, rownames(md_0_h))
# remove species rows that sum to zero
c_0_h_osps <- c_0_h_osps[rowSums(c_0_h_osps[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_0_osps = colSums(c_0_h_osps)
# total sum
sum_all_healthy_0_osps = sum(sum_healthy_0_osps)
# add abundance column
c_0_h_osps$abundance <- (rowSums(c_0_h_osps[,1:ncol(c_0_h_osps)]) / sum_all_healthy_0_osps) * 100
# sort abundance decreasing
c_0_h_osps <- c_0_h_osps[with(c_0_h_osps, order(-abundance)), ]
# obtain cumulative sum
c_0_h_osps$cumsum <- cumsum(c_0_h_osps$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile <- quantile(c_0_h_osps$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_0_h_osps_core_a1 <- subset(c_0_h_osps, cumsum <= abund_quantile[1])
# remove non-numeric columns
c_0_h_osps_core_a1$cumsum <- NULL
c_0_h_osps_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_0_h_osps_core_a2 <- subset(c_0_h_osps, cumsum <= abund_quantile[2])
# remove non-numeric columns
c_0_h_osps_core_a2$cumsum <- NULL
c_0_h_osps_core_a2$abundance <- NULL
# # subset core species based on third defined abundance quantile
c_0_h_osps_core_a3 <- subset(c_0_h_osps, cumsum <= abund_quantile[3])
# remove non-numeric columns
c_0_h_osps_core_a3$cumsum <- NULL
c_0_h_osps_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_0_h_osps_rare_a1 <- subset(c_0_h_osps, cumsum > abund_quantile[1])
# remove non-numeric columns
c_0_h_osps_rare_a1$cumsum <- NULL
c_0_h_osps_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_0_h_osps_rare_a2 <- subset(c_0_h_osps, cumsum > abund_quantile[2])
# remove non-numeric columns
c_0_h_osps_rare_a2$cumsum <- NULL
c_0_h_osps_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_0_h_osps_rare_a3 <- subset(c_0_h_osps, cumsum > abund_quantile[3])
# remove non-numeric columns
c_0_h_osps_rare_a3$cumsum <- NULL
c_0_h_osps_rare_a3$abundance <- NULL


# Age group (1-3 years), healthy, osps database, bcphc-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_1_3_h_osps <- select(ds_osps, rownames(md_1_3_h))
# remove species rows that sum to zero
c_1_3_h_osps <- c_1_3_h_osps[rowSums(c_1_3_h_osps[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_1_3_osps = colSums(c_1_3_h_osps)
# obtain total sum
sum_all_healthy_1_3_osps = sum(sum_healthy_1_3_osps)
# add abundance column
c_1_3_h_osps$abundance <- (rowSums(c_1_3_h_osps[,1:ncol(c_1_3_h_osps)]) / sum_all_healthy_1_3_osps) * 100
# sort abundance decreasing
c_1_3_h_osps <- c_1_3_h_osps[with(c_1_3_h_osps, order(-abundance)), ]
# obtain cumulative sum
c_1_3_h_osps$cumsum <- cumsum(c_1_3_h_osps$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_1_3_h <- quantile(c_1_3_h_osps$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_1_3_h_osps_core_a1 <- subset(c_1_3_h_osps, cumsum <= abund_quantile_1_3_h[1])
# remove non-numeric columns
c_1_3_h_osps_core_a1$cumsum <- NULL
c_1_3_h_osps_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_1_3_h_osps_core_a2 <- subset(c_1_3_h_osps, cumsum <= abund_quantile_1_3_h[2])
# remove non-numeric columns
c_1_3_h_osps_core_a2$cumsum <- NULL
c_1_3_h_osps_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_1_3_h_osps_core_a3 <- subset(c_1_3_h_osps, cumsum <= abund_quantile_1_3_h[3])
# remove non-numeric columns
c_1_3_h_osps_core_a3$cumsum <- NULL
c_1_3_h_osps_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_1_3_h_osps_rare_a1 <- subset(c_1_3_h_osps, cumsum > abund_quantile_1_3_h[1])
# remove non-numeric columns
c_1_3_h_osps_rare_a1$cumsum <- NULL
c_1_3_h_osps_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_1_3_h_osps_rare_a2 <- subset(c_1_3_h_osps, cumsum > abund_quantile_1_3_h[2])
# remove non-numeric columns
c_1_3_h_osps_rare_a2$cumsum <- NULL
c_1_3_h_osps_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_1_3_h_osps_rare_a3 <- subset(c_1_3_h_osps, cumsum > abund_quantile_1_3_h[3])
# remove non-numeric columns
c_1_3_h_osps_rare_a3$cumsum <- NULL
c_1_3_h_osps_rare_a3$abundance <- NULL


# Age group (4-6 years), healthy, osps database, bcphc-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_4_6_h_osps <- select(ds_osps, rownames(md_4_6_h))
# remove species rows that sum to zero
c_4_6_h_osps <- c_4_6_h_osps[rowSums(c_4_6_h_osps[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_4_6_osps = colSums(c_4_6_h_osps)
# get total sum
sum_all_healthy_4_6_osps = sum(sum_healthy_4_6_osps)
# add abundance column
c_4_6_h_osps$abundance <- (rowSums(c_4_6_h_osps[,1:ncol(c_4_6_h_osps)]) / sum_all_healthy_4_6_osps) * 100
# sort abundance decreasing
c_4_6_h_osps <- c_4_6_h_osps[with(c_4_6_h_osps, order(-abundance)), ]
# get cumulative sum
c_4_6_h_osps$cumsum <- cumsum(c_4_6_h_osps$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_4_6_h <- quantile(c_4_6_h_osps$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_4_6_h_osps_core_a1 <- subset(c_4_6_h_osps, cumsum <= abund_quantile_4_6_h[1])
# remove non-numeric columns
c_4_6_h_osps_core_a1$cumsum <- NULL
c_4_6_h_osps_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_4_6_h_osps_core_a2 <- subset(c_4_6_h_osps, cumsum <= abund_quantile_4_6_h[2])
# remove non-numeric columns
c_4_6_h_osps_core_a2$cumsum <- NULL
c_4_6_h_osps_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_4_6_h_osps_core_a3 <- subset(c_4_6_h_osps, cumsum <= abund_quantile_4_6_h[3])
# remove non-numeric columns
c_4_6_h_osps_core_a3$cumsum <- NULL
c_4_6_h_osps_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_4_6_h_osps_rare_a1 <- subset(c_4_6_h_osps, cumsum > abund_quantile_4_6_h[1])
# remove non-numeric columns
c_4_6_h_osps_rare_a1$cumsum <- NULL
c_4_6_h_osps_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_4_6_h_osps_rare_a2 <- subset(c_4_6_h_osps, cumsum > abund_quantile_4_6_h[2])
# remove non-numeric columns
c_4_6_h_osps_rare_a2$cumsum <- NULL
c_4_6_h_osps_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_4_6_h_osps_rare_a3 <- subset(c_4_6_h_osps, cumsum > abund_quantile_4_6_h[3])
# remove non-numeric columns
c_4_6_h_osps_rare_a3$cumsum <- NULL
c_4_6_h_osps_rare_a3$abundance <- NULL


# Age group (0 years), CF, osps database, bcphc-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_0_cf_osps <- select(ds_osps, rownames(md_0_cf))
# remove species rows that sum to zero
c_0_cf_osps <- c_0_cf_osps[rowSums(c_0_cf_osps[, -1])>0, ]
# sum species counts per CF child
sum_cf_0_osps = colSums(c_0_cf_osps)
# get total sum
sum_all_cf_0_osps = sum(sum_cf_0_osps)
# add abundance column
c_0_cf_osps$abundance <- (rowSums(c_0_cf_osps[,1:ncol(c_0_cf_osps)]) / sum_all_cf_0_osps) * 100
# sort abundance decreasing
c_0_cf_osps <- c_0_cf_osps[with(c_0_cf_osps, order(-abundance)), ]
# add column with cumulative sum
c_0_cf_osps$cumsum <- cumsum(c_0_cf_osps$abundance)

# create sample quantiles in dataset
abund_quantile_0_cf <- quantile(c_0_cf_osps$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_0_cf_osps_core_a1 <- subset(c_0_cf_osps, cumsum <= abund_quantile_0_cf[1])
# remove non-numeric columns
c_0_cf_osps_core_a1$cumsum <- NULL
c_0_cf_osps_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_0_cf_osps_core_a2 <- subset(c_0_cf_osps, cumsum <= abund_quantile_0_cf[2])
# remove non-numeric columns
c_0_cf_osps_core_a2$cumsum <- NULL
c_0_cf_osps_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_0_cf_osps_core_a3 <- subset(c_0_cf_osps, cumsum <= abund_quantile_0_cf[3])
# remove non-numeric columns
c_0_cf_osps_core_a3$cumsum <- NULL
c_0_cf_osps_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_0_cf_osps_rare_a1 <- subset(c_0_cf_osps, cumsum > abund_quantile_0_cf[1])
# remove non-numeric columns
c_0_cf_osps_rare_a1$cumsum <- NULL
c_0_cf_osps_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_0_cf_osps_rare_a2 <- subset(c_0_cf_osps, cumsum > abund_quantile_0_cf[2])
# remove non-numeric columns
c_0_cf_osps_rare_a2$cumsum <- NULL
c_0_cf_osps_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_0_cf_osps_rare_a3 <- subset(c_0_cf_osps, cumsum > abund_quantile_0_cf[3])
# remove non-numeric columns
c_0_cf_osps_rare_a3$cumsum <- NULL
c_0_cf_osps_rare_a3$abundance <- NULL


# Age group (1-3 years), CF, osps database, bcphc-normalised
c_1_3_cf_osps <- select(ds_osps, rownames(md_1_3_cf))
# remove species rows that sum to zero
c_1_3_cf_osps <- c_1_3_cf_osps[rowSums(c_1_3_cf_osps[, -1])>0, ]
# sum species counts per CF child
sum_cf_1_3_osps = colSums(c_1_3_cf_osps)
# get total sum
sum_all_cf_1_3_osps = sum(sum_cf_1_3_osps)
# add abundance column
c_1_3_cf_osps$abundance <- (rowSums(c_1_3_cf_osps[,1:ncol(c_1_3_cf_osps)]) / sum_all_cf_1_3_osps) * 100
# sort abundance decreasing
c_1_3_cf_osps <- c_1_3_cf_osps[with(c_1_3_cf_osps, order(-abundance)), ]
# add column with cumulative sum
c_1_3_cf_osps$cumsum <- cumsum(c_1_3_cf_osps$abundance)

# table with abundant species
# create sample quantiles in dataset
abund_quantile_1_3_cf <- quantile(c_1_3_cf_osps$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_1_3_cf_osps_core_a1 <- subset(c_1_3_cf_osps, cumsum <= abund_quantile_1_3_cf[1])
# remove non-numeric columns
c_1_3_cf_osps_core_a1$cumsum <- NULL
c_1_3_cf_osps_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_1_3_cf_osps_core_a2 <- subset(c_1_3_cf_osps, cumsum <= abund_quantile_1_3_cf[2])
# remove non-numeric columns
c_1_3_cf_osps_core_a2$cumsum <- NULL
c_1_3_cf_osps_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_1_3_cf_osps_core_a3 <- subset(c_1_3_cf_osps, cumsum <= abund_quantile_1_3_cf[3])
# remove non-numeric columns
c_1_3_cf_osps_core_a3$cumsum <- NULL
c_1_3_cf_osps_core_a3$abundance <- NULL


# table with rare species
# subset rare species based on first defined abundance quantile
c_1_3_cf_osps_rare_a1 <- subset(c_1_3_cf_osps, cumsum > abund_quantile_1_3_cf[1])
# remove non-numeric columns
c_1_3_cf_osps_rare_a1$cumsum <- NULL
c_1_3_cf_osps_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_1_3_cf_osps_rare_a2 <- subset(c_1_3_cf_osps, cumsum > abund_quantile_1_3_cf[2])
# remove non-numeric columns
c_1_3_cf_osps_rare_a2$cumsum <- NULL
c_1_3_cf_osps_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_1_3_cf_osps_rare_a3 <- subset(c_1_3_cf_osps, cumsum > abund_quantile_1_3_cf[3])
# remove non-numeric columns
c_1_3_cf_osps_rare_a3$cumsum <- NULL
c_1_3_cf_osps_rare_a3$abundance <- NULL


# Age group (4-6 years), CF, osps database, bcphc-normalised
c_4_6_cf_osps <- select(ds_osps, rownames(md_4_6_cf))
# remove species rows that sum to zero
c_4_6_cf_osps <- c_4_6_cf_osps[rowSums(c_4_6_cf_osps[, -1])>0, ]
# sum species counts per CF child
sum_cf_4_6_osps = colSums(c_4_6_cf_osps)
# get total sum
sum_all_cf_4_6_osps = sum(sum_cf_4_6_osps)
# add abundance column
c_4_6_cf_osps$abundance <- (rowSums(c_4_6_cf_osps[,1:ncol(c_4_6_cf_osps)]) / sum_all_cf_4_6_osps) * 100
# sort abundance decreasing
c_4_6_cf_osps <- c_4_6_cf_osps[with(c_4_6_cf_osps, order(-abundance)), ]
# add column with cumulative sum
c_4_6_cf_osps$cumsum <- cumsum(c_4_6_cf_osps$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_4_6_cf <- quantile(c_4_6_cf_osps$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_4_6_cf_osps_core_a1 <- subset(c_4_6_cf_osps, cumsum <= abund_quantile_4_6_cf[1])
# remove non-numeric columns
c_4_6_cf_osps_core_a1$cumsum <- NULL
c_4_6_cf_osps_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_4_6_cf_osps_core_a2 <- subset(c_4_6_cf_osps, cumsum <= abund_quantile_4_6_cf[2])
# remove non-numeric columns
c_4_6_cf_osps_core_a2$cumsum <- NULL
c_4_6_cf_osps_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_4_6_cf_osps_core_a3 <- subset(c_4_6_cf_osps, cumsum <= abund_quantile_4_6_cf[3])
# remove non-numeric columns
c_4_6_cf_osps_core_a3$cumsum <- NULL
c_4_6_cf_osps_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_4_6_cf_osps_rare_a1 <- subset(c_4_6_cf_osps, cumsum > abund_quantile_4_6_cf[1])
# remove non-numeric columns
c_4_6_cf_osps_rare_a1$cumsum <- NULL
c_4_6_cf_osps_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_4_6_cf_osps_rare_a2 <- subset(c_4_6_cf_osps, cumsum > abund_quantile_4_6_cf[2])
# remove non-numeric columns
c_4_6_cf_osps_rare_a2$cumsum <- NULL
c_4_6_cf_osps_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_4_6_cf_osps_rare_a3 <- subset(c_4_6_cf_osps, cumsum > abund_quantile_4_6_cf[3])
# remove non-numeric columns
c_4_6_cf_osps_rare_a3$cumsum <- NULL
c_4_6_cf_osps_rare_a3$abundance <- NULL


# create vector with all rare species (15th abundance percentile, osps, healthy)
healthy_rare_a1_osps <- c(rownames(c_0_h_osps_rare_a1), 
                               rownames(c_1_3_h_osps_rare_a1), 
                               rownames(c_4_6_h_osps_rare_a1))
# remove duplicate entries
healthy_rare_a1_osps <- healthy_rare_a1_osps[!duplicated(healthy_rare_a1_osps)]

# create vector with all core species (15th abundance percentile, osps, healthy)
healthy_core_a1_osps <- c(rownames(c_0_h_osps_core_a1), 
                               rownames(c_1_3_h_osps_core_a1), 
                               rownames(c_4_6_h_osps_core_a1))
# remove duplicate entries
healthy_core_a1_osps <- healthy_core_a1_osps[!duplicated(healthy_core_a1_osps)]

# create vector with all rare species (25th abundance percentile, osps, healthy)
healthy_rare_a2_osps <- c(rownames(c_0_h_osps_rare_a2), 
                               rownames(c_1_3_h_osps_rare_a2), 
                               rownames(c_4_6_h_osps_rare_a2))
# remove duplicate entries
healthy_rare_a2_osps <- healthy_rare_a2_osps[!duplicated(healthy_rare_a2_osps)]

# create vector with all core species (25th abundance percentile, osps, healthy)
healthy_core_a2_osps <- c(rownames(c_0_h_osps_core_a2), 
                               rownames(c_1_3_h_osps_core_a2), 
                               rownames(c_4_6_h_osps_core_a2))
# remove duplicate entries
healthy_core_a2_osps <- healthy_core_a2_osps[!duplicated(healthy_core_a2_osps)]

# create vector with all rare species (35th abundance percentile, osps, healthy)
healthy_rare_a3_osps <- c(rownames(c_0_h_osps_rare_a3), 
                               rownames(c_1_3_h_osps_rare_a3), 
                               rownames(c_4_6_h_osps_rare_a3))
# remove duplicate entries
healthy_rare_a3_osps <- healthy_rare_a3_osps[!duplicated(healthy_rare_a3_osps)]

# create vector with all core species (35th abundance percentile, osps, healthy)
healthy_core_a3_osps <- c(rownames(c_0_h_osps_core_a3), 
                               rownames(c_1_3_h_osps_core_a3), 
                               rownames(c_4_6_h_osps_core_a3))
# remove duplicate entries
healthy_core_a3_osps <- healthy_core_a3_osps[!duplicated(healthy_core_a3_osps)]


# create vector with all rare species (15th abundance percentile, osps, CF)
cf_rare_a1_osps <- c(rownames(c_0_cf_osps_rare_a1), 
                          rownames(c_1_3_cf_osps_rare_a1), 
                          rownames(c_4_6_cf_osps_rare_a1))
# remove duplicate entries
cf_rare_a1_osps <- cf_rare_a1_osps[!duplicated(cf_rare_a1_osps)]

# create vector with all core species (15th abundance percentile, osps, CF)
cf_core_a1_osps <- c(rownames(c_0_cf_osps_core_a1), 
                          rownames(c_1_3_cf_osps_core_a1), 
                          rownames(c_4_6_cf_osps_core_a1))
# remove duplicate entries
cf_core_a1_osps <- cf_core_a1_osps[!duplicated(cf_core_a1_osps)]

# create vector with all rare species (25th abundance percentile, osps, CF)
cf_rare_a2_osps <- c(rownames(c_0_cf_osps_rare_a2), 
                          rownames(c_1_3_cf_osps_rare_a2), 
                          rownames(c_4_6_cf_osps_rare_a2))
# remove duplicate entries
cf_rare_a2_osps <- cf_rare_a2_osps[!duplicated(cf_rare_a2_osps)]

# create vector with all core species (25th abundance percentile, osps, CF)
cf_core_a2_osps <- c(rownames(c_0_cf_osps_core_a2), 
                          rownames(c_1_3_cf_osps_core_a2), 
                          rownames(c_4_6_cf_osps_core_a2))
# remove duplicate entries
cf_core_a2_osps <- cf_core_a2_osps[!duplicated(cf_core_a2_osps)]

# create vector with all rare species (35th abundance percentile, osps, CF)
cf_rare_a3_osps <- c(rownames(c_0_cf_osps_rare_a3), 
                          rownames(c_1_3_cf_osps_rare_a3), 
                          rownames(c_4_6_cf_osps_rare_a3))
# remove duplicate entries
cf_rare_a3_osps <- cf_rare_a3_osps[!duplicated(cf_rare_a3_osps)]

# create vector with all core species (35th abundance percentile, osps, CF)
cf_core_a3_osps <- c(rownames(c_0_cf_osps_core_a3), 
                          rownames(c_1_3_cf_osps_core_a3), 
                          rownames(c_4_6_cf_osps_core_a3))
# remove duplicate entries
cf_core_a3_osps <- cf_core_a3_osps[!duplicated(cf_core_a3_osps)]



# age group (0 years), healthy, osps database, rle-normalised
c_0_h_osps_rle <- select(ds_osps_rle, rownames(md_0_h))
# remove species rows that sum to zero
c_0_h_osps_rle <- c_0_h_osps_rle[rowSums(c_0_h_osps_rle[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_0_osps_rle = colSums(c_0_h_osps_rle)
# total sum
sum_all_healthy_0_osps_rle = sum(sum_healthy_0_osps_rle)
# add abundance column
c_0_h_osps_rle$abundance <- (rowSums(c_0_h_osps_rle[,1:ncol(c_0_h_osps_rle)]) / sum_all_healthy_0_osps_rle) * 100
# sort abundance decreasing
c_0_h_osps_rle <- c_0_h_osps_rle[with(c_0_h_osps_rle, order(-abundance)), ]
# obtain cumulative sum
c_0_h_osps_rle$cumsum <- cumsum(c_0_h_osps_rle$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile <- quantile(c_0_h_osps_rle$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_0_h_osps_rle_core_a1 <- subset(c_0_h_osps_rle, cumsum <= abund_quantile[1])
# remove non-numeric columns
c_0_h_osps_rle_core_a1$cumsum <- NULL
c_0_h_osps_rle_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_0_h_osps_rle_core_a2 <- subset(c_0_h_osps_rle, cumsum <= abund_quantile[2])
# remove non-numeric columns
c_0_h_osps_rle_core_a2$cumsum <- NULL
c_0_h_osps_rle_core_a2$abundance <- NULL
# # subset core species based on third defined abundance quantile
c_0_h_osps_rle_core_a3 <- subset(c_0_h_osps_rle, cumsum <= abund_quantile[3])
# remove non-numeric columns
c_0_h_osps_rle_core_a3$cumsum <- NULL
c_0_h_osps_rle_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_0_h_osps_rle_rare_a1 <- subset(c_0_h_osps_rle, cumsum > abund_quantile[1])
# remove non-numeric columns
c_0_h_osps_rle_rare_a1$cumsum <- NULL
c_0_h_osps_rle_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_0_h_osps_rle_rare_a2 <- subset(c_0_h_osps_rle, cumsum > abund_quantile[2])
# remove non-numeric columns
c_0_h_osps_rle_rare_a2$cumsum <- NULL
c_0_h_osps_rle_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_0_h_osps_rle_rare_a3 <- subset(c_0_h_osps_rle, cumsum > abund_quantile[3])
# remove non-numeric columns
c_0_h_osps_rle_rare_a3$cumsum <- NULL
c_0_h_osps_rle_rare_a3$abundance <- NULL


# Age group (1-3 years), healthy, osps database, rle-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_1_3_h_osps_rle <- select(ds_osps_rle, rownames(md_1_3_h))
# remove species rows that sum to zero
c_1_3_h_osps_rle <- c_1_3_h_osps_rle[rowSums(c_1_3_h_osps_rle[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_1_3_osps_rle = colSums(c_1_3_h_osps_rle)
# obtain total sum
sum_all_healthy_1_3_osps_rle = sum(sum_healthy_1_3_osps_rle)
# add abundance column
c_1_3_h_osps_rle$abundance <- (rowSums(c_1_3_h_osps_rle[,1:ncol(c_1_3_h_osps_rle)]) / sum_all_healthy_1_3_osps_rle) * 100
# sort abundance decreasing
c_1_3_h_osps_rle <- c_1_3_h_osps_rle[with(c_1_3_h_osps_rle, order(-abundance)), ]
# obtain cumulative sum
c_1_3_h_osps_rle$cumsum <- cumsum(c_1_3_h_osps_rle$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_1_3_h <- quantile(c_1_3_h_osps_rle$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_1_3_h_osps_rle_core_a1 <- subset(c_1_3_h_osps_rle, cumsum <= abund_quantile_1_3_h[1])
# remove non-numeric columns
c_1_3_h_osps_rle_core_a1$cumsum <- NULL
c_1_3_h_osps_rle_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_1_3_h_osps_rle_core_a2 <- subset(c_1_3_h_osps_rle, cumsum <= abund_quantile_1_3_h[2])
# remove non-numeric columns
c_1_3_h_osps_rle_core_a2$cumsum <- NULL
c_1_3_h_osps_rle_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_1_3_h_osps_rle_core_a3 <- subset(c_1_3_h_osps_rle, cumsum <= abund_quantile_1_3_h[3])
# remove non-numeric columns
c_1_3_h_osps_rle_core_a3$cumsum <- NULL
c_1_3_h_osps_rle_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_1_3_h_osps_rle_rare_a1 <- subset(c_1_3_h_osps_rle, cumsum > abund_quantile_1_3_h[1])
# remove non-numeric columns
c_1_3_h_osps_rle_rare_a1$cumsum <- NULL
c_1_3_h_osps_rle_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_1_3_h_osps_rle_rare_a2 <- subset(c_1_3_h_osps_rle, cumsum > abund_quantile_1_3_h[2])
# remove non-numeric columns
c_1_3_h_osps_rle_rare_a2$cumsum <- NULL
c_1_3_h_osps_rle_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_1_3_h_osps_rle_rare_a3 <- subset(c_1_3_h_osps_rle, cumsum > abund_quantile_1_3_h[3])
# remove non-numeric columns
c_1_3_h_osps_rle_rare_a3$cumsum <- NULL
c_1_3_h_osps_rle_rare_a3$abundance <- NULL


# Age group (4-6 years), healthy, osps database, rle-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_4_6_h_osps_rle <- select(ds_osps_rle, rownames(md_4_6_h))
# remove species rows that sum to zero
c_4_6_h_osps_rle <- c_4_6_h_osps_rle[rowSums(c_4_6_h_osps_rle[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_4_6_osps_rle = colSums(c_4_6_h_osps_rle)
# get total sum
sum_all_healthy_4_6_osps_rle = sum(sum_healthy_4_6_osps_rle)
# add abundance column
c_4_6_h_osps_rle$abundance <- (rowSums(c_4_6_h_osps_rle[,1:ncol(c_4_6_h_osps_rle)]) / sum_all_healthy_4_6_osps_rle) * 100
# sort abundance decreasing
c_4_6_h_osps_rle <- c_4_6_h_osps_rle[with(c_4_6_h_osps_rle, order(-abundance)), ]
# get cumulative sum
c_4_6_h_osps_rle$cumsum <- cumsum(c_4_6_h_osps_rle$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_4_6_h <- quantile(c_4_6_h_osps_rle$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_4_6_h_osps_rle_core_a1 <- subset(c_4_6_h_osps_rle, cumsum <= abund_quantile_4_6_h[1])
# remove non-numeric columns
c_4_6_h_osps_rle_core_a1$cumsum <- NULL
c_4_6_h_osps_rle_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_4_6_h_osps_rle_core_a2 <- subset(c_4_6_h_osps_rle, cumsum <= abund_quantile_4_6_h[2])
# remove non-numeric columns
c_4_6_h_osps_rle_core_a2$cumsum <- NULL
c_4_6_h_osps_rle_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_4_6_h_osps_rle_core_a3 <- subset(c_4_6_h_osps_rle, cumsum <= abund_quantile_4_6_h[3])
# remove non-numeric columns
c_4_6_h_osps_rle_core_a3$cumsum <- NULL
c_4_6_h_osps_rle_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_4_6_h_osps_rle_rare_a1 <- subset(c_4_6_h_osps_rle, cumsum > abund_quantile_4_6_h[1])
# remove non-numeric columns
c_4_6_h_osps_rle_rare_a1$cumsum <- NULL
c_4_6_h_osps_rle_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_4_6_h_osps_rle_rare_a2 <- subset(c_4_6_h_osps_rle, cumsum > abund_quantile_4_6_h[2])
# remove non-numeric columns
c_4_6_h_osps_rle_rare_a2$cumsum <- NULL
c_4_6_h_osps_rle_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_4_6_h_osps_rle_rare_a3 <- subset(c_4_6_h_osps_rle, cumsum > abund_quantile_4_6_h[3])
# remove non-numeric columns
c_4_6_h_osps_rle_rare_a3$cumsum <- NULL
c_4_6_h_osps_rle_rare_a3$abundance <- NULL


# Age group (0 years), CF, osps database, rle-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_0_cf_osps_rle <- select(ds_osps_rle, rownames(md_0_cf))
# remove species rows that sum to zero
c_0_cf_osps_rle <- c_0_cf_osps_rle[rowSums(c_0_cf_osps_rle[, -1])>0, ]
# sum species counts per CF child
sum_cf_0_osps_rle = colSums(c_0_cf_osps_rle)
# get total sum
sum_all_cf_0_osps_rle = sum(sum_cf_0_osps_rle)
# add abundance column
c_0_cf_osps_rle$abundance <- (rowSums(c_0_cf_osps_rle[,1:ncol(c_0_cf_osps_rle)]) / sum_all_cf_0_osps_rle) * 100
# sort abundance decreasing
c_0_cf_osps_rle <- c_0_cf_osps_rle[with(c_0_cf_osps_rle, order(-abundance)), ]
# add column with cumulative sum
c_0_cf_osps_rle$cumsum <- cumsum(c_0_cf_osps_rle$abundance)

# create sample quantiles in dataset
abund_quantile_0_cf <- quantile(c_0_cf_osps_rle$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_0_cf_osps_rle_core_a1 <- subset(c_0_cf_osps_rle, cumsum <= abund_quantile_0_cf[1])
# remove non-numeric columns
c_0_cf_osps_rle_core_a1$cumsum <- NULL
c_0_cf_osps_rle_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_0_cf_osps_rle_core_a2 <- subset(c_0_cf_osps_rle, cumsum <= abund_quantile_0_cf[2])
# remove non-numeric columns
c_0_cf_osps_rle_core_a2$cumsum <- NULL
c_0_cf_osps_rle_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_0_cf_osps_rle_core_a3 <- subset(c_0_cf_osps_rle, cumsum <= abund_quantile_0_cf[3])
# remove non-numeric columns
c_0_cf_osps_rle_core_a3$cumsum <- NULL
c_0_cf_osps_rle_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_0_cf_osps_rle_rare_a1 <- subset(c_0_cf_osps_rle, cumsum > abund_quantile_0_cf[1])
# remove non-numeric columns
c_0_cf_osps_rle_rare_a1$cumsum <- NULL
c_0_cf_osps_rle_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_0_cf_osps_rle_rare_a2 <- subset(c_0_cf_osps_rle, cumsum > abund_quantile_0_cf[2])
# remove non-numeric columns
c_0_cf_osps_rle_rare_a2$cumsum <- NULL
c_0_cf_osps_rle_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_0_cf_osps_rle_rare_a3 <- subset(c_0_cf_osps_rle, cumsum > abund_quantile_0_cf[3])
# remove non-numeric columns
c_0_cf_osps_rle_rare_a3$cumsum <- NULL
c_0_cf_osps_rle_rare_a3$abundance <- NULL


# Age group (1-3 years), CF, osps database, rle-normalised
c_1_3_cf_osps_rle <- select(ds_osps_rle, rownames(md_1_3_cf))
# remove species rows that sum to zero
c_1_3_cf_osps_rle <- c_1_3_cf_osps_rle[rowSums(c_1_3_cf_osps_rle[, -1])>0, ]
# sum species counts per CF child
sum_cf_1_3_osps_rle = colSums(c_1_3_cf_osps_rle)
# get total sum
sum_all_cf_1_3_osps_rle = sum(sum_cf_1_3_osps_rle)
# add abundance column
c_1_3_cf_osps_rle$abundance <- (rowSums(c_1_3_cf_osps_rle[,1:ncol(c_1_3_cf_osps_rle)]) / sum_all_cf_1_3_osps_rle) * 100
# sort abundance decreasing
c_1_3_cf_osps_rle <- c_1_3_cf_osps_rle[with(c_1_3_cf_osps_rle, order(-abundance)), ]
# add column with cumulative sum
c_1_3_cf_osps_rle$cumsum <- cumsum(c_1_3_cf_osps_rle$abundance)

# table with abundant species
# create sample quantiles in dataset
abund_quantile_1_3_cf <- quantile(c_1_3_cf_osps_rle$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_1_3_cf_osps_rle_core_a1 <- subset(c_1_3_cf_osps_rle, cumsum <= abund_quantile_1_3_cf[1])
# remove non-numeric columns
c_1_3_cf_osps_rle_core_a1$cumsum <- NULL
c_1_3_cf_osps_rle_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_1_3_cf_osps_rle_core_a2 <- subset(c_1_3_cf_osps_rle, cumsum <= abund_quantile_1_3_cf[2])
# remove non-numeric columns
c_1_3_cf_osps_rle_core_a2$cumsum <- NULL
c_1_3_cf_osps_rle_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_1_3_cf_osps_rle_core_a3 <- subset(c_1_3_cf_osps_rle, cumsum <= abund_quantile_1_3_cf[3])
# remove non-numeric columns
c_1_3_cf_osps_rle_core_a3$cumsum <- NULL
c_1_3_cf_osps_rle_core_a3$abundance <- NULL


# table with rare species
# subset rare species based on first defined abundance quantile
c_1_3_cf_osps_rle_rare_a1 <- subset(c_1_3_cf_osps_rle, cumsum > abund_quantile_1_3_cf[1])
# remove non-numeric columns
c_1_3_cf_osps_rle_rare_a1$cumsum <- NULL
c_1_3_cf_osps_rle_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_1_3_cf_osps_rle_rare_a2 <- subset(c_1_3_cf_osps_rle, cumsum > abund_quantile_1_3_cf[2])
# remove non-numeric columns
c_1_3_cf_osps_rle_rare_a2$cumsum <- NULL
c_1_3_cf_osps_rle_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_1_3_cf_osps_rle_rare_a3 <- subset(c_1_3_cf_osps_rle, cumsum > abund_quantile_1_3_cf[3])
# remove non-numeric columns
c_1_3_cf_osps_rle_rare_a3$cumsum <- NULL
c_1_3_cf_osps_rle_rare_a3$abundance <- NULL


# Age group (4-6 years), CF, osps database, rle-normalised
c_4_6_cf_osps_rle <- select(ds_osps_rle, rownames(md_4_6_cf))
# remove species rows that sum to zero
c_4_6_cf_osps_rle <- c_4_6_cf_osps_rle[rowSums(c_4_6_cf_osps_rle[, -1])>0, ]
# sum species counts per CF child
sum_cf_4_6_osps_rle = colSums(c_4_6_cf_osps_rle)
# get total sum
sum_all_cf_4_6_osps_rle = sum(sum_cf_4_6_osps_rle)
# add abundance column
c_4_6_cf_osps_rle$abundance <- (rowSums(c_4_6_cf_osps_rle[,1:ncol(c_4_6_cf_osps_rle)]) / sum_all_cf_4_6_osps_rle) * 100
# sort abundance decreasing
c_4_6_cf_osps_rle <- c_4_6_cf_osps_rle[with(c_4_6_cf_osps_rle, order(-abundance)), ]
# add column with cumulative sum
c_4_6_cf_osps_rle$cumsum <- cumsum(c_4_6_cf_osps_rle$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_4_6_cf <- quantile(c_4_6_cf_osps_rle$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_4_6_cf_osps_rle_core_a1 <- subset(c_4_6_cf_osps_rle, cumsum <= abund_quantile_4_6_cf[1])
# remove non-numeric columns
c_4_6_cf_osps_rle_core_a1$cumsum <- NULL
c_4_6_cf_osps_rle_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_4_6_cf_osps_rle_core_a2 <- subset(c_4_6_cf_osps_rle, cumsum <= abund_quantile_4_6_cf[2])
# remove non-numeric columns
c_4_6_cf_osps_rle_core_a2$cumsum <- NULL
c_4_6_cf_osps_rle_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_4_6_cf_osps_rle_core_a3 <- subset(c_4_6_cf_osps_rle, cumsum <= abund_quantile_4_6_cf[3])
# remove non-numeric columns
c_4_6_cf_osps_rle_core_a3$cumsum <- NULL
c_4_6_cf_osps_rle_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_4_6_cf_osps_rle_rare_a1 <- subset(c_4_6_cf_osps_rle, cumsum > abund_quantile_4_6_cf[1])
# remove non-numeric columns
c_4_6_cf_osps_rle_rare_a1$cumsum <- NULL
c_4_6_cf_osps_rle_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_4_6_cf_osps_rle_rare_a2 <- subset(c_4_6_cf_osps_rle, cumsum > abund_quantile_4_6_cf[2])
# remove non-numeric columns
c_4_6_cf_osps_rle_rare_a2$cumsum <- NULL
c_4_6_cf_osps_rle_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_4_6_cf_osps_rle_rare_a3 <- subset(c_4_6_cf_osps_rle, cumsum > abund_quantile_4_6_cf[3])
# remove non-numeric columns
c_4_6_cf_osps_rle_rare_a3$cumsum <- NULL
c_4_6_cf_osps_rle_rare_a3$abundance <- NULL


# create vector with all rare species (15th abundance percentile, osps, healthy)
healthy_rare_a1_osps_rle <- c(rownames(c_0_h_osps_rle_rare_a1), 
                                   rownames(c_1_3_h_osps_rle_rare_a1), 
                                   rownames(c_4_6_h_osps_rle_rare_a1))
# remove duplicate entries
healthy_rare_a1_osps_rle <- healthy_rare_a1_osps_rle[!duplicated(healthy_rare_a1_osps_rle)]

# create vector with all core species (15th abundance percentile, osps, healthy)
healthy_core_a1_osps_rle <- c(rownames(c_0_h_osps_rle_core_a1), 
                                   rownames(c_1_3_h_osps_rle_core_a1), 
                                   rownames(c_4_6_h_osps_rle_core_a1))
# remove duplicate entries
healthy_core_a1_osps_rle <- healthy_core_a1_osps_rle[!duplicated(healthy_core_a1_osps_rle)]

# create vector with all rare species (25th abundance percentile, osps, healthy)
healthy_rare_a2_osps_rle <- c(rownames(c_0_h_osps_rle_rare_a2), 
                                   rownames(c_1_3_h_osps_rle_rare_a2), 
                                   rownames(c_4_6_h_osps_rle_rare_a2))
# remove duplicate entries
healthy_rare_a2_osps_rle <- healthy_rare_a2_osps_rle[!duplicated(healthy_rare_a2_osps_rle)]

# create vector with all core species (25th abundance percentile, osps, healthy)
healthy_core_a2_osps_rle <- c(rownames(c_0_h_osps_rle_core_a2), 
                                   rownames(c_1_3_h_osps_rle_core_a2), 
                                   rownames(c_4_6_h_osps_rle_core_a2))
# remove duplicate entries
healthy_core_a2_osps_rle <- healthy_core_a2_osps_rle[!duplicated(healthy_core_a2_osps_rle)]

# create vector with all rare species (35th abundance percentile, osps, healthy)
healthy_rare_a3_osps_rle <- c(rownames(c_0_h_osps_rle_rare_a3), 
                                   rownames(c_1_3_h_osps_rle_rare_a3), 
                                   rownames(c_4_6_h_osps_rle_rare_a3))
# remove duplicate entries
healthy_rare_a3_osps_rle <- healthy_rare_a3_osps_rle[!duplicated(healthy_rare_a3_osps_rle)]

# create vector with all core species (35th abundance percentile, osps, healthy)
healthy_core_a3_osps_rle <- c(rownames(c_0_h_osps_rle_core_a3), 
                                   rownames(c_1_3_h_osps_rle_core_a3), 
                                   rownames(c_4_6_h_osps_rle_core_a3))
# remove duplicate entries
healthy_core_a3_osps_rle <- healthy_core_a3_osps_rle[!duplicated(healthy_core_a3_osps_rle)]


# create vector with all rare species (15th abundance percentile, osps, CF)
cf_rare_a1_osps_rle <- c(rownames(c_0_cf_osps_rle_rare_a1), 
                              rownames(c_1_3_cf_osps_rle_rare_a1), 
                              rownames(c_4_6_cf_osps_rle_rare_a1))
# remove duplicate entries
cf_rare_a1_osps_rle <- cf_rare_a1_osps_rle[!duplicated(cf_rare_a1_osps_rle)]

# create vector with all core species (15th abundance percentile, osps, CF)
cf_core_a1_osps_rle <- c(rownames(c_0_cf_osps_rle_core_a1), 
                              rownames(c_1_3_cf_osps_rle_core_a1), 
                              rownames(c_4_6_cf_osps_rle_core_a1))
# remove duplicate entries
cf_core_a1_osps_rle <- cf_core_a1_osps_rle[!duplicated(cf_core_a1_osps_rle)]

# create vector with all rare species (25th abundance percentile, osps, CF)
cf_rare_a2_osps_rle <- c(rownames(c_0_cf_osps_rle_rare_a2), 
                              rownames(c_1_3_cf_osps_rle_rare_a2), 
                              rownames(c_4_6_cf_osps_rle_rare_a2))
# remove duplicate entries
cf_rare_a2_osps_rle <- cf_rare_a2_osps_rle[!duplicated(cf_rare_a2_osps_rle)]

# create vector with all core species (25th abundance percentile, osps, CF)
cf_core_a2_osps_rle <- c(rownames(c_0_cf_osps_rle_core_a2), 
                              rownames(c_1_3_cf_osps_rle_core_a2), 
                              rownames(c_4_6_cf_osps_rle_core_a2))
# remove duplicate entries
cf_core_a2_osps_rle <- cf_core_a2_osps_rle[!duplicated(cf_core_a2_osps_rle)]

# create vector with all rare species (35th abundance percentile, osps, CF)
cf_rare_a3_osps_rle <- c(rownames(c_0_cf_osps_rle_rare_a3), 
                              rownames(c_1_3_cf_osps_rle_rare_a3), 
                              rownames(c_4_6_cf_osps_rle_rare_a3))
# remove duplicate entries
cf_rare_a3_osps_rle <- cf_rare_a3_osps_rle[!duplicated(cf_rare_a3_osps_rle)]

# create vector with all core species (35th abundance percentile, osps, CF)
cf_core_a3_osps_rle <- c(rownames(c_0_cf_osps_rle_core_a3), 
                              rownames(c_1_3_cf_osps_rle_core_a3), 
                              rownames(c_4_6_cf_osps_rle_core_a3))
# remove duplicate entries
cf_core_a3_osps_rle <- cf_core_a3_osps_rle[!duplicated(cf_core_a3_osps_rle)]

# age group (0 years), healthy, osps database, vst-normalised
c_0_h_osps_vst <- select(ds_osps_vst, rownames(md_0_h))
# remove species rows that sum to zero
c_0_h_osps_vst <- c_0_h_osps_vst[rowSums(c_0_h_osps_vst[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_0_osps_vst = colSums(c_0_h_osps_vst)
# total sum
sum_all_healthy_0_osps_vst = sum(sum_healthy_0_osps_vst)
# add abundance column
c_0_h_osps_vst$abundance <- (rowSums(c_0_h_osps_vst[,1:ncol(c_0_h_osps_vst)]) / sum_all_healthy_0_osps_vst) * 100
# sort abundance decreasing
c_0_h_osps_vst <- c_0_h_osps_vst[with(c_0_h_osps_vst, order(-abundance)), ]
# obtain cumulative sum
c_0_h_osps_vst$cumsum <- cumsum(c_0_h_osps_vst$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile <- quantile(c_0_h_osps_vst$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_0_h_osps_vst_core_a1 <- subset(c_0_h_osps_vst, cumsum <= abund_quantile[1])
# remove non-numeric columns
c_0_h_osps_vst_core_a1$cumsum <- NULL
c_0_h_osps_vst_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_0_h_osps_vst_core_a2 <- subset(c_0_h_osps_vst, cumsum <= abund_quantile[2])
# remove non-numeric columns
c_0_h_osps_vst_core_a2$cumsum <- NULL
c_0_h_osps_vst_core_a2$abundance <- NULL
# # subset core species based on third defined abundance quantile
c_0_h_osps_vst_core_a3 <- subset(c_0_h_osps_vst, cumsum <= abund_quantile[3])
# remove non-numeric columns
c_0_h_osps_vst_core_a3$cumsum <- NULL
c_0_h_osps_vst_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_0_h_osps_vst_rare_a1 <- subset(c_0_h_osps_vst, cumsum > abund_quantile[1])
# remove non-numeric columns
c_0_h_osps_vst_rare_a1$cumsum <- NULL
c_0_h_osps_vst_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_0_h_osps_vst_rare_a2 <- subset(c_0_h_osps_vst, cumsum > abund_quantile[2])
# remove non-numeric columns
c_0_h_osps_vst_rare_a2$cumsum <- NULL
c_0_h_osps_vst_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_0_h_osps_vst_rare_a3 <- subset(c_0_h_osps_vst, cumsum > abund_quantile[3])
# remove non-numeric columns
c_0_h_osps_vst_rare_a3$cumsum <- NULL
c_0_h_osps_vst_rare_a3$abundance <- NULL


# Age group (1-3 years), healthy, osps database, vst-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_1_3_h_osps_vst <- select(ds_osps_vst, rownames(md_1_3_h))
# remove species rows that sum to zero
c_1_3_h_osps_vst <- c_1_3_h_osps_vst[rowSums(c_1_3_h_osps_vst[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_1_3_osps_vst = colSums(c_1_3_h_osps_vst)
# obtain total sum
sum_all_healthy_1_3_osps_vst = sum(sum_healthy_1_3_osps_vst)
# add abundance column
c_1_3_h_osps_vst$abundance <- (rowSums(c_1_3_h_osps_vst[,1:ncol(c_1_3_h_osps_vst)]) / sum_all_healthy_1_3_osps_vst) * 100
# sort abundance decreasing
c_1_3_h_osps_vst <- c_1_3_h_osps_vst[with(c_1_3_h_osps_vst, order(-abundance)), ]
# obtain cumulative sum
c_1_3_h_osps_vst$cumsum <- cumsum(c_1_3_h_osps_vst$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_1_3_h <- quantile(c_1_3_h_osps_vst$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_1_3_h_osps_vst_core_a1 <- subset(c_1_3_h_osps_vst, cumsum <= abund_quantile_1_3_h[1])
# remove non-numeric columns
c_1_3_h_osps_vst_core_a1$cumsum <- NULL
c_1_3_h_osps_vst_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_1_3_h_osps_vst_core_a2 <- subset(c_1_3_h_osps_vst, cumsum <= abund_quantile_1_3_h[2])
# remove non-numeric columns
c_1_3_h_osps_vst_core_a2$cumsum <- NULL
c_1_3_h_osps_vst_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_1_3_h_osps_vst_core_a3 <- subset(c_1_3_h_osps_vst, cumsum <= abund_quantile_1_3_h[3])
# remove non-numeric columns
c_1_3_h_osps_vst_core_a3$cumsum <- NULL
c_1_3_h_osps_vst_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_1_3_h_osps_vst_rare_a1 <- subset(c_1_3_h_osps_vst, cumsum > abund_quantile_1_3_h[1])
# remove non-numeric columns
c_1_3_h_osps_vst_rare_a1$cumsum <- NULL
c_1_3_h_osps_vst_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_1_3_h_osps_vst_rare_a2 <- subset(c_1_3_h_osps_vst, cumsum > abund_quantile_1_3_h[2])
# remove non-numeric columns
c_1_3_h_osps_vst_rare_a2$cumsum <- NULL
c_1_3_h_osps_vst_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_1_3_h_osps_vst_rare_a3 <- subset(c_1_3_h_osps_vst, cumsum > abund_quantile_1_3_h[3])
# remove non-numeric columns
c_1_3_h_osps_vst_rare_a3$cumsum <- NULL
c_1_3_h_osps_vst_rare_a3$abundance <- NULL


# Age group (4-6 years), healthy, osps database, vst-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_4_6_h_osps_vst <- select(ds_osps_vst, rownames(md_4_6_h))
# remove species rows that sum to zero
c_4_6_h_osps_vst <- c_4_6_h_osps_vst[rowSums(c_4_6_h_osps_vst[, -1])>0, ]
# sum species counts per healthy child
sum_healthy_4_6_osps_vst = colSums(c_4_6_h_osps_vst)
# get total sum
sum_all_healthy_4_6_osps_vst = sum(sum_healthy_4_6_osps_vst)
# add abundance column
c_4_6_h_osps_vst$abundance <- (rowSums(c_4_6_h_osps_vst[,1:ncol(c_4_6_h_osps_vst)]) / sum_all_healthy_4_6_osps_vst) * 100
# sort abundance decreasing
c_4_6_h_osps_vst <- c_4_6_h_osps_vst[with(c_4_6_h_osps_vst, order(-abundance)), ]
# get cumulative sum
c_4_6_h_osps_vst$cumsum <- cumsum(c_4_6_h_osps_vst$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_4_6_h <- quantile(c_4_6_h_osps_vst$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_4_6_h_osps_vst_core_a1 <- subset(c_4_6_h_osps_vst, cumsum <= abund_quantile_4_6_h[1])
# remove non-numeric columns
c_4_6_h_osps_vst_core_a1$cumsum <- NULL
c_4_6_h_osps_vst_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_4_6_h_osps_vst_core_a2 <- subset(c_4_6_h_osps_vst, cumsum <= abund_quantile_4_6_h[2])
# remove non-numeric columns
c_4_6_h_osps_vst_core_a2$cumsum <- NULL
c_4_6_h_osps_vst_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_4_6_h_osps_vst_core_a3 <- subset(c_4_6_h_osps_vst, cumsum <= abund_quantile_4_6_h[3])
# remove non-numeric columns
c_4_6_h_osps_vst_core_a3$cumsum <- NULL
c_4_6_h_osps_vst_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_4_6_h_osps_vst_rare_a1 <- subset(c_4_6_h_osps_vst, cumsum > abund_quantile_4_6_h[1])
# remove non-numeric columns
c_4_6_h_osps_vst_rare_a1$cumsum <- NULL
c_4_6_h_osps_vst_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_4_6_h_osps_vst_rare_a2 <- subset(c_4_6_h_osps_vst, cumsum > abund_quantile_4_6_h[2])
# remove non-numeric columns
c_4_6_h_osps_vst_rare_a2$cumsum <- NULL
c_4_6_h_osps_vst_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_4_6_h_osps_vst_rare_a3 <- subset(c_4_6_h_osps_vst, cumsum > abund_quantile_4_6_h[3])
# remove non-numeric columns
c_4_6_h_osps_vst_rare_a3$cumsum <- NULL
c_4_6_h_osps_vst_rare_a3$abundance <- NULL


# Age group (0 years), CF, osps database, vst-normalised
# extract core and rare species based on the three rarity thresholds, defined above (abund_quantile)
c_0_cf_osps_vst <- select(ds_osps_vst, rownames(md_0_cf))
# remove species rows that sum to zero
c_0_cf_osps_vst <- c_0_cf_osps_vst[rowSums(c_0_cf_osps_vst[, -1])>0, ]
# sum species counts per CF child
sum_cf_0_osps_vst = colSums(c_0_cf_osps_vst)
# get total sum
sum_all_cf_0_osps_vst = sum(sum_cf_0_osps_vst)
# add abundance column
c_0_cf_osps_vst$abundance <- (rowSums(c_0_cf_osps_vst[,1:ncol(c_0_cf_osps_vst)]) / sum_all_cf_0_osps_vst) * 100
# sort abundance decreasing
c_0_cf_osps_vst <- c_0_cf_osps_vst[with(c_0_cf_osps_vst, order(-abundance)), ]
# add column with cumulative sum
c_0_cf_osps_vst$cumsum <- cumsum(c_0_cf_osps_vst$abundance)

# create sample quantiles in dataset
abund_quantile_0_cf <- quantile(c_0_cf_osps_vst$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_0_cf_osps_vst_core_a1 <- subset(c_0_cf_osps_vst, cumsum <= abund_quantile_0_cf[1])
# remove non-numeric columns
c_0_cf_osps_vst_core_a1$cumsum <- NULL
c_0_cf_osps_vst_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_0_cf_osps_vst_core_a2 <- subset(c_0_cf_osps_vst, cumsum <= abund_quantile_0_cf[2])
# remove non-numeric columns
c_0_cf_osps_vst_core_a2$cumsum <- NULL
c_0_cf_osps_vst_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_0_cf_osps_vst_core_a3 <- subset(c_0_cf_osps_vst, cumsum <= abund_quantile_0_cf[3])
# remove non-numeric columns
c_0_cf_osps_vst_core_a3$cumsum <- NULL
c_0_cf_osps_vst_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_0_cf_osps_vst_rare_a1 <- subset(c_0_cf_osps_vst, cumsum > abund_quantile_0_cf[1])
# remove non-numeric columns
c_0_cf_osps_vst_rare_a1$cumsum <- NULL
c_0_cf_osps_vst_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_0_cf_osps_vst_rare_a2 <- subset(c_0_cf_osps_vst, cumsum > abund_quantile_0_cf[2])
# remove non-numeric columns
c_0_cf_osps_vst_rare_a2$cumsum <- NULL
c_0_cf_osps_vst_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_0_cf_osps_vst_rare_a3 <- subset(c_0_cf_osps_vst, cumsum > abund_quantile_0_cf[3])
# remove non-numeric columns
c_0_cf_osps_vst_rare_a3$cumsum <- NULL
c_0_cf_osps_vst_rare_a3$abundance <- NULL


# Age group (1-3 years), CF, osps database, vst-normalised
c_1_3_cf_osps_vst <- select(ds_osps_vst, rownames(md_1_3_cf))
# remove species rows that sum to zero
c_1_3_cf_osps_vst <- c_1_3_cf_osps_vst[rowSums(c_1_3_cf_osps_vst[, -1])>0, ]
# sum species counts per CF child
sum_cf_1_3_osps_vst = colSums(c_1_3_cf_osps_vst)
# get total sum
sum_all_cf_1_3_osps_vst = sum(sum_cf_1_3_osps_vst)
# add abundance column
c_1_3_cf_osps_vst$abundance <- (rowSums(c_1_3_cf_osps_vst[,1:ncol(c_1_3_cf_osps_vst)]) / sum_all_cf_1_3_osps_vst) * 100
# sort abundance decreasing
c_1_3_cf_osps_vst <- c_1_3_cf_osps_vst[with(c_1_3_cf_osps_vst, order(-abundance)), ]
# add column with cumulative sum
c_1_3_cf_osps_vst$cumsum <- cumsum(c_1_3_cf_osps_vst$abundance)

# table with abundant species
# create sample quantiles in dataset
abund_quantile_1_3_cf <- quantile(c_1_3_cf_osps_vst$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_1_3_cf_osps_vst_core_a1 <- subset(c_1_3_cf_osps_vst, cumsum <= abund_quantile_1_3_cf[1])
# remove non-numeric columns
c_1_3_cf_osps_vst_core_a1$cumsum <- NULL
c_1_3_cf_osps_vst_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_1_3_cf_osps_vst_core_a2 <- subset(c_1_3_cf_osps_vst, cumsum <= abund_quantile_1_3_cf[2])
# remove non-numeric columns
c_1_3_cf_osps_vst_core_a2$cumsum <- NULL
c_1_3_cf_osps_vst_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_1_3_cf_osps_vst_core_a3 <- subset(c_1_3_cf_osps_vst, cumsum <= abund_quantile_1_3_cf[3])
# remove non-numeric columns
c_1_3_cf_osps_vst_core_a3$cumsum <- NULL
c_1_3_cf_osps_vst_core_a3$abundance <- NULL


# table with rare species
# subset rare species based on first defined abundance quantile
c_1_3_cf_osps_vst_rare_a1 <- subset(c_1_3_cf_osps_vst, cumsum > abund_quantile_1_3_cf[1])
# remove non-numeric columns
c_1_3_cf_osps_vst_rare_a1$cumsum <- NULL
c_1_3_cf_osps_vst_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_1_3_cf_osps_vst_rare_a2 <- subset(c_1_3_cf_osps_vst, cumsum > abund_quantile_1_3_cf[2])
# remove non-numeric columns
c_1_3_cf_osps_vst_rare_a2$cumsum <- NULL
c_1_3_cf_osps_vst_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_1_3_cf_osps_vst_rare_a3 <- subset(c_1_3_cf_osps_vst, cumsum > abund_quantile_1_3_cf[3])
# remove non-numeric columns
c_1_3_cf_osps_vst_rare_a3$cumsum <- NULL
c_1_3_cf_osps_vst_rare_a3$abundance <- NULL


# Age group (4-6 years), CF, osps database, vst-normalised
c_4_6_cf_osps_vst <- select(ds_osps_vst, rownames(md_4_6_cf))
# remove species rows that sum to zero
c_4_6_cf_osps_vst <- c_4_6_cf_osps_vst[rowSums(c_4_6_cf_osps_vst[, -1])>0, ]
# sum species counts per CF child
sum_cf_4_6_osps_vst = colSums(c_4_6_cf_osps_vst)
# get total sum
sum_all_cf_4_6_osps_vst = sum(sum_cf_4_6_osps_vst)
# add abundance column
c_4_6_cf_osps_vst$abundance <- (rowSums(c_4_6_cf_osps_vst[,1:ncol(c_4_6_cf_osps_vst)]) / sum_all_cf_4_6_osps_vst) * 100
# sort abundance decreasing
c_4_6_cf_osps_vst <- c_4_6_cf_osps_vst[with(c_4_6_cf_osps_vst, order(-abundance)), ]
# add column with cumulative sum
c_4_6_cf_osps_vst$cumsum <- cumsum(c_4_6_cf_osps_vst$abundance)

# table with core species
# create sample quantiles in dataset
abund_quantile_4_6_cf <- quantile(c_4_6_cf_osps_vst$cumsum, probs = quantile_range)
# subset core species based on first defined abundance quantile
c_4_6_cf_osps_vst_core_a1 <- subset(c_4_6_cf_osps_vst, cumsum <= abund_quantile_4_6_cf[1])
# remove non-numeric columns
c_4_6_cf_osps_vst_core_a1$cumsum <- NULL
c_4_6_cf_osps_vst_core_a1$abundance <- NULL
# subset core species based on second defined abundance quantile
c_4_6_cf_osps_vst_core_a2 <- subset(c_4_6_cf_osps_vst, cumsum <= abund_quantile_4_6_cf[2])
# remove non-numeric columns
c_4_6_cf_osps_vst_core_a2$cumsum <- NULL
c_4_6_cf_osps_vst_core_a2$abundance <- NULL
# subset core species based on third defined abundance quantile
c_4_6_cf_osps_vst_core_a3 <- subset(c_4_6_cf_osps_vst, cumsum <= abund_quantile_4_6_cf[3])
# remove non-numeric columns
c_4_6_cf_osps_vst_core_a3$cumsum <- NULL
c_4_6_cf_osps_vst_core_a3$abundance <- NULL

# table with rare species
# subset rare species based on first defined abundance quantile
c_4_6_cf_osps_vst_rare_a1 <- subset(c_4_6_cf_osps_vst, cumsum > abund_quantile_4_6_cf[1])
# remove non-numeric columns
c_4_6_cf_osps_vst_rare_a1$cumsum <- NULL
c_4_6_cf_osps_vst_rare_a1$abundance <- NULL
# subset rare species based on second defined abundance quantile
c_4_6_cf_osps_vst_rare_a2 <- subset(c_4_6_cf_osps_vst, cumsum > abund_quantile_4_6_cf[2])
# remove non-numeric columns
c_4_6_cf_osps_vst_rare_a2$cumsum <- NULL
c_4_6_cf_osps_vst_rare_a2$abundance <- NULL
# subset rare species based on third defined abundance quantile
c_4_6_cf_osps_vst_rare_a3 <- subset(c_4_6_cf_osps_vst, cumsum > abund_quantile_4_6_cf[3])
# remove non-numeric columns
c_4_6_cf_osps_vst_rare_a3$cumsum <- NULL
c_4_6_cf_osps_vst_rare_a3$abundance <- NULL


# create vector with all rare species (15th abundance percentile, osps, healthy)
healthy_rare_a1_osps_vst <- c(rownames(c_0_h_osps_vst_rare_a1), 
                                   rownames(c_1_3_h_osps_vst_rare_a1), 
                                   rownames(c_4_6_h_osps_vst_rare_a1))
# remove duplicate entries
healthy_rare_a1_osps_vst <- healthy_rare_a1_osps_vst[!duplicated(healthy_rare_a1_osps_vst)]

# create vector with all core species (15th abundance percentile, osps, healthy)
healthy_core_a1_osps_vst <- c(rownames(c_0_h_osps_vst_core_a1), 
                                   rownames(c_1_3_h_osps_vst_core_a1), 
                                   rownames(c_4_6_h_osps_vst_core_a1))
# remove duplicate entries
healthy_core_a1_osps_vst <- healthy_core_a1_osps_vst[!duplicated(healthy_core_a1_osps_vst)]

# create vector with all rare species (25th abundance percentile, osps, healthy)
healthy_rare_a2_osps_vst <- c(rownames(c_0_h_osps_vst_rare_a2), 
                                   rownames(c_1_3_h_osps_vst_rare_a2), 
                                   rownames(c_4_6_h_osps_vst_rare_a2))
# remove duplicate entries
healthy_rare_a2_osps_vst <- healthy_rare_a2_osps_vst[!duplicated(healthy_rare_a2_osps_vst)]

# create vector with all core species (25th abundance percentile, osps, healthy)
healthy_core_a2_osps_vst <- c(rownames(c_0_h_osps_vst_core_a2), 
                                   rownames(c_1_3_h_osps_vst_core_a2), 
                                   rownames(c_4_6_h_osps_vst_core_a2))
# remove duplicate entries
healthy_core_a2_osps_vst <- healthy_core_a2_osps_vst[!duplicated(healthy_core_a2_osps_vst)]

# create vector with all rare species (35th abundance percentile, osps, healthy)
healthy_rare_a3_osps_vst <- c(rownames(c_0_h_osps_vst_rare_a3), 
                                   rownames(c_1_3_h_osps_vst_rare_a3), 
                                   rownames(c_4_6_h_osps_vst_rare_a3))
# remove duplicate entries
healthy_rare_a3_osps_vst <- healthy_rare_a3_osps_vst[!duplicated(healthy_rare_a3_osps_vst)]

# create vector with all core species (35th abundance percentile, osps, healthy)
healthy_core_a3_osps_vst <- c(rownames(c_0_h_osps_vst_core_a3), 
                                   rownames(c_1_3_h_osps_vst_core_a3), 
                                   rownames(c_4_6_h_osps_vst_core_a3))
# remove duplicate entries
healthy_core_a3_osps_vst <- healthy_core_a3_osps_vst[!duplicated(healthy_core_a3_osps_vst)]


# create vector with all rare species (15th abundance percentile, osps, CF)
cf_rare_a1_osps_vst <- c(rownames(c_0_cf_osps_vst_rare_a1), 
                              rownames(c_1_3_cf_osps_vst_rare_a1), 
                              rownames(c_4_6_cf_osps_vst_rare_a1))
# remove duplicate entries
cf_rare_a1_osps_vst <- cf_rare_a1_osps_vst[!duplicated(cf_rare_a1_osps_vst)]

# create vector with all core species (15th abundance percentile, osps, CF)
cf_core_a1_osps_vst <- c(rownames(c_0_cf_osps_vst_core_a1), 
                              rownames(c_1_3_cf_osps_vst_core_a1), 
                              rownames(c_4_6_cf_osps_vst_core_a1))
# remove duplicate entries
cf_core_a1_osps_vst <- cf_core_a1_osps_vst[!duplicated(cf_core_a1_osps_vst)]

# create vector with all rare species (25th abundance percentile, osps, CF)
cf_rare_a2_osps_vst <- c(rownames(c_0_cf_osps_vst_rare_a2), 
                              rownames(c_1_3_cf_osps_vst_rare_a2), 
                              rownames(c_4_6_cf_osps_vst_rare_a2))
# remove duplicate entries
cf_rare_a2_osps_vst <- cf_rare_a2_osps_vst[!duplicated(cf_rare_a2_osps_vst)]

# create vector with all core species (25th abundance percentile, osps, CF)
cf_core_a2_osps_vst <- c(rownames(c_0_cf_osps_vst_core_a2), 
                              rownames(c_1_3_cf_osps_vst_core_a2), 
                              rownames(c_4_6_cf_osps_vst_core_a2))
# remove duplicate entries
cf_core_a2_osps_vst <- cf_core_a2_osps_vst[!duplicated(cf_core_a2_osps_vst)]

# create vector with all rare species (35th abundance percentile, osps, CF)
cf_rare_a3_osps_vst <- c(rownames(c_0_cf_osps_vst_rare_a3), 
                              rownames(c_1_3_cf_osps_vst_rare_a3), 
                              rownames(c_4_6_cf_osps_vst_rare_a3))
# remove duplicate entries
cf_rare_a3_osps_vst <- cf_rare_a3_osps_vst[!duplicated(cf_rare_a3_osps_vst)]

# create vector with all core species (35th abundance percentile, osps, CF)
cf_core_a3_osps_vst <- c(rownames(c_0_cf_osps_vst_core_a3), 
                              rownames(c_1_3_cf_osps_vst_core_a3), 
                              rownames(c_4_6_cf_osps_vst_core_a3))
# remove duplicate entries
cf_core_a3_osps_vst <- cf_core_a3_osps_vst[!duplicated(cf_core_a3_osps_vst)]

############################################################################################################
# Venn diagram analysis, pangenome, BCPHC
# make vector with all rare species in healthy and CF children, 15th abundance percentile
all_rare_a1_pangenome <- c(healthy_rare_a1_pangenome, cf_rare_a1_pangenome)
# make vector with all core species in healthy and CF children, 15th abundance percentile
all_core_a1_pangenome <- c(healthy_core_a1_pangenome, cf_core_a1_pangenome)
# make vector with all rare species in healthy and CF children, 25th abundance percentile
all_rare_a2_pangenome <- c(healthy_rare_a2_pangenome, cf_rare_a2_pangenome)
# make vector with all core species in healthy and CF children, 25th abundance percentile
all_core_a2_pangenome <- c(healthy_core_a2_pangenome, cf_core_a2_pangenome)
# make vector with all rare species in healthy and CF children, 35th abundance percentile
all_rare_a3_pangenome <- c(healthy_rare_a3_pangenome, cf_rare_a3_pangenome)
# make vector with all core species in healthy and CF children, 35th abundance percentile
all_core_a3_pangenome <- c(healthy_core_a3_pangenome, cf_core_a3_pangenome)

# rename original dataframe (BCPHC)
venn_rare_h_pangenome <- data.frame(ds_pangenome)
# rare species, healthy, BCPHC
# 15th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_pangenome$h_0_summary_a1 <- ifelse(rownames(venn_rare_h_pangenome) %in% rownames(c_0_h_pangenome_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_pangenome$h_1_3_summary_a1 <- ifelse(rownames(venn_rare_h_pangenome) %in% rownames(c_1_3_h_pangenome_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_pangenome$h_4_6_summary_a1 <- ifelse(rownames(venn_rare_h_pangenome) %in% rownames(c_4_6_h_pangenome_rare_a1), 1, 0)
# 25th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_pangenome$h_0_summary_a2 <- ifelse(rownames(venn_rare_h_pangenome) %in% rownames(c_0_h_pangenome_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_pangenome$h_1_3_summary_a2 <- ifelse(rownames(venn_rare_h_pangenome) %in% rownames(c_1_3_h_pangenome_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_pangenome$h_4_6_summary_a2 <- ifelse(rownames(venn_rare_h_pangenome) %in% rownames(c_4_6_h_pangenome_rare_a2), 1, 0)
# 35th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_pangenome$h_0_summary_a3 <- ifelse(rownames(venn_rare_h_pangenome) %in% rownames(c_0_h_pangenome_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_pangenome$h_1_3_summary_a3 <- ifelse(rownames(venn_rare_h_pangenome) %in% rownames(c_1_3_h_pangenome_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_pangenome$h_4_6_summary_a3 <- ifelse(rownames(venn_rare_h_pangenome) %in% rownames(c_4_6_h_pangenome_rare_a3), 1, 0)

# select the first three columns (rare species presence/absence per age, 15th abundance percentile)
venn_rare_h_pangenome_a1 <- venn_rare_h_pangenome %>% select(h_0_summary_a1, h_1_3_summary_a1, h_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_pangenome_a1 <- venn_rare_h_pangenome_a1[rowSums(venn_rare_h_pangenome_a1) > 0,]
# define the background species, present in all age groups
background_rare_h_pangenome_a1 <- venn_rare_h_pangenome_a1[rowSums(venn_rare_h_pangenome_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_h_pangenome_a1 <- venn_rare_h_pangenome_a1[rowSums(venn_rare_h_pangenome_a1)<3,]

# select the next three columns (rare species presence/absence per age, 25th abundance percentile)
venn_rare_h_pangenome_a2 <- venn_rare_h_pangenome %>% select(h_0_summary_a2, h_1_3_summary_a2, h_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_pangenome_a2 <- venn_rare_h_pangenome_a2[rowSums(venn_rare_h_pangenome_a2) > 0,]
background_rare_h_pangenome_a2 <- venn_rare_h_pangenome_a2[rowSums(venn_rare_h_pangenome_a2)>2,]
non_persistent_rare_h_pangenome_a2 <- venn_rare_h_pangenome_a2[rowSums(venn_rare_h_pangenome_a2)<3,]

# select the next three columns (rare species presence/absence per age, 35th abundance percentile)
venn_rare_h_pangenome_a3 <- venn_rare_h_pangenome %>% select(h_0_summary_a3, h_1_3_summary_a3, h_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_pangenome_a3 <- venn_rare_h_pangenome_a3[rowSums(venn_rare_h_pangenome_a3) > 0,]
# define the background species, present in all age groups
background_rare_h_pangenome_a3 <- venn_rare_h_pangenome_a3[rowSums(venn_rare_h_pangenome_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_h_pangenome_a3 <- venn_rare_h_pangenome_a3[rowSums(venn_rare_h_pangenome_a3)<3,]

# core species, healthy, BCPHC
venn_core_h_pangenome <- data.frame(ds_pangenome)
# 15th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_pangenome$h_0_summary_a1 <- ifelse(rownames(venn_core_h_pangenome) %in% rownames(c_0_h_pangenome_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_pangenome$h_1_3_summary_a1 <- ifelse(rownames(venn_core_h_pangenome) %in% rownames(c_1_3_h_pangenome_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_pangenome$h_4_6_summary_a1 <- ifelse(rownames(venn_core_h_pangenome) %in% rownames(c_4_6_h_pangenome_core_a1), 1, 0)
# 25th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_pangenome$h_0_summary_a2 <- ifelse(rownames(venn_core_h_pangenome) %in% rownames(c_0_h_pangenome_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_pangenome$h_1_3_summary_a2 <- ifelse(rownames(venn_core_h_pangenome) %in% rownames(c_1_3_h_pangenome_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_pangenome$h_4_6_summary_a2 <- ifelse(rownames(venn_core_h_pangenome) %in% rownames(c_4_6_h_pangenome_core_a2), 1, 0)
# 35th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_pangenome$h_0_summary_a3 <- ifelse(rownames(venn_core_h_pangenome) %in% rownames(c_0_h_pangenome_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_pangenome$h_1_3_summary_a3 <- ifelse(rownames(venn_core_h_pangenome) %in% rownames(c_1_3_h_pangenome_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_pangenome$h_4_6_summary_a3 <- ifelse(rownames(venn_core_h_pangenome) %in% rownames(c_4_6_h_pangenome_core_a3), 1, 0)

# select the first three columns (core species presence/absence per age, 15th abundance percentile)
venn_core_h_pangenome_a1 <- venn_core_h_pangenome %>% select(h_0_summary_a1, h_1_3_summary_a1, h_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_core_h_pangenome_a1 <- venn_core_h_pangenome_a1[rowSums(venn_core_h_pangenome_a1) > 0,]
# define the background species, present in all age groups
background_core_h_pangenome_a1 <- venn_core_h_pangenome_a1[rowSums(venn_core_h_pangenome_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_h_pangenome_a1 <- venn_core_h_pangenome_a1[rowSums(venn_core_h_pangenome_a1)<3,]

# select the next three columns (core species presence/absence per age, 25th abundance percentile)
venn_core_h_pangenome_a2 <- venn_core_h_pangenome %>% select(h_0_summary_a2, h_1_3_summary_a2, h_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_core_h_pangenome_a2 <- venn_core_h_pangenome_a2[rowSums(venn_core_h_pangenome_a2) > 0,]
background_core_h_pangenome_a2 <- venn_core_h_pangenome_a2[rowSums(venn_core_h_pangenome_a2)>2,]
non_persistent_core_h_pangenome_a2 <- venn_core_h_pangenome_a2[rowSums(venn_core_h_pangenome_a2)<3,]

# select the next three columns (core species presence/absence per age, 35th abundance percentile)
venn_core_h_pangenome_a3 <- venn_core_h_pangenome %>% select(h_0_summary_a3, h_1_3_summary_a3, h_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_core_h_pangenome_a3 <- venn_core_h_pangenome_a3[rowSums(venn_core_h_pangenome_a3) > 0,]
# define the background species, present in all age groups
background_core_h_pangenome_a3 <- venn_core_h_pangenome_a3[rowSums(venn_core_h_pangenome_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_h_pangenome_a3 <- venn_core_h_pangenome_a3[rowSums(venn_core_h_pangenome_a3)<3,]



# repeat with CF children, rare, BCPHC-normalised, pangenome
# rename original dataframe (BCPHC)
venn_rare_cf_pangenome <- data.frame(ds_pangenome)
# rare species, cf, BCPHC
# 15th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_pangenome$cf_0_summary_a1 <- ifelse(rownames(venn_rare_cf_pangenome) %in% rownames(c_0_cf_pangenome_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_pangenome$cf_1_3_summary_a1 <- ifelse(rownames(venn_rare_cf_pangenome) %in% rownames(c_1_3_cf_pangenome_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_pangenome$cf_4_6_summary_a1 <- ifelse(rownames(venn_rare_cf_pangenome) %in% rownames(c_4_6_cf_pangenome_rare_a1), 1, 0)
# 25th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_pangenome$cf_0_summary_a2 <- ifelse(rownames(venn_rare_cf_pangenome) %in% rownames(c_0_cf_pangenome_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_pangenome$cf_1_3_summary_a2 <- ifelse(rownames(venn_rare_cf_pangenome) %in% rownames(c_1_3_cf_pangenome_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_pangenome$cf_4_6_summary_a2 <- ifelse(rownames(venn_rare_cf_pangenome) %in% rownames(c_4_6_cf_pangenome_rare_a2), 1, 0)
# 35th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_pangenome$cf_0_summary_a3 <- ifelse(rownames(venn_rare_cf_pangenome) %in% rownames(c_0_cf_pangenome_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_pangenome$cf_1_3_summary_a3 <- ifelse(rownames(venn_rare_cf_pangenome) %in% rownames(c_1_3_cf_pangenome_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_pangenome$cf_4_6_summary_a3 <- ifelse(rownames(venn_rare_cf_pangenome) %in% rownames(c_4_6_cf_pangenome_rare_a3), 1, 0)

# select the first three columns (rare species presence/absence per age, 15th abundance percentile)
venn_rare_cf_pangenome_a1 <- venn_rare_cf_pangenome %>% select(cf_0_summary_a1, cf_1_3_summary_a1, cf_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_pangenome_a1 <- venn_rare_cf_pangenome_a1[rowSums(venn_rare_cf_pangenome_a1) > 0,]
# define the background species, present in all age groups
background_rare_cf_pangenome_a1 <- venn_rare_cf_pangenome_a1[rowSums(venn_rare_cf_pangenome_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_cf_pangenome_a1 <- venn_rare_cf_pangenome_a1[rowSums(venn_rare_cf_pangenome_a1)<3,]

# select the next three columns (rare species presence/absence per age, 25th abundance percentile)
venn_rare_cf_pangenome_a2 <- venn_rare_cf_pangenome %>% select(cf_0_summary_a2, cf_1_3_summary_a2, cf_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_pangenome_a2 <- venn_rare_cf_pangenome_a2[rowSums(venn_rare_cf_pangenome_a2) > 0,]
background_rare_cf_pangenome_a2 <- venn_rare_cf_pangenome_a2[rowSums(venn_rare_cf_pangenome_a2)>2,]
non_persistent_rare_cf_pangenome_a2 <- venn_rare_cf_pangenome_a2[rowSums(venn_rare_cf_pangenome_a2)<3,]

# select the next three columns (rare species presence/absence per age, 35th abundance percentile)
venn_rare_cf_pangenome_a3 <- venn_rare_cf_pangenome %>% select(cf_0_summary_a3, cf_1_3_summary_a3, cf_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_pangenome_a3 <- venn_rare_cf_pangenome_a3[rowSums(venn_rare_cf_pangenome_a3) > 0,]
# define the background species, present in all age groups
background_rare_cf_pangenome_a3 <- venn_rare_cf_pangenome_a3[rowSums(venn_rare_cf_pangenome_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_cf_pangenome_a3 <- venn_rare_cf_pangenome_a3[rowSums(venn_rare_cf_pangenome_a3)<3,]

# core species, cf, BCPHC
venn_core_cf_pangenome <- data.frame(ds_pangenome)
# 15th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_pangenome$cf_0_summary_a1 <- ifelse(rownames(venn_core_cf_pangenome) %in% rownames(c_0_cf_pangenome_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_pangenome$cf_1_3_summary_a1 <- ifelse(rownames(venn_core_cf_pangenome) %in% rownames(c_1_3_cf_pangenome_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_pangenome$cf_4_6_summary_a1 <- ifelse(rownames(venn_core_cf_pangenome) %in% rownames(c_4_6_cf_pangenome_core_a1), 1, 0)
# 25th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_pangenome$cf_0_summary_a2 <- ifelse(rownames(venn_core_cf_pangenome) %in% rownames(c_0_cf_pangenome_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_pangenome$cf_1_3_summary_a2 <- ifelse(rownames(venn_core_cf_pangenome) %in% rownames(c_1_3_cf_pangenome_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_pangenome$cf_4_6_summary_a2 <- ifelse(rownames(venn_core_cf_pangenome) %in% rownames(c_4_6_cf_pangenome_core_a2), 1, 0)
# 35th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_pangenome$cf_0_summary_a3 <- ifelse(rownames(venn_core_cf_pangenome) %in% rownames(c_0_cf_pangenome_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_pangenome$cf_1_3_summary_a3 <- ifelse(rownames(venn_core_cf_pangenome) %in% rownames(c_1_3_cf_pangenome_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_pangenome$cf_4_6_summary_a3 <- ifelse(rownames(venn_core_cf_pangenome) %in% rownames(c_4_6_cf_pangenome_core_a3), 1, 0)

# select the first three columns (core species presence/absence per age, 15th abundance percentile)
venn_core_cf_pangenome_a1 <- venn_core_cf_pangenome %>% select(cf_0_summary_a1, cf_1_3_summary_a1, cf_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_pangenome_a1 <- venn_core_cf_pangenome_a1[rowSums(venn_core_cf_pangenome_a1) > 0,]
# define the background species, present in all age groups
background_core_cf_pangenome_a1 <- venn_core_cf_pangenome_a1[rowSums(venn_core_cf_pangenome_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_cf_pangenome_a1 <- venn_core_cf_pangenome_a1[rowSums(venn_core_cf_pangenome_a1)<3,]

# select the next three columns (core species presence/absence per age, 25th abundance percentile)
venn_core_cf_pangenome_a2 <- venn_core_cf_pangenome %>% select(cf_0_summary_a2, cf_1_3_summary_a2, cf_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_pangenome_a2 <- venn_core_cf_pangenome_a2[rowSums(venn_core_cf_pangenome_a2) > 0,]
background_core_cf_pangenome_a2 <- venn_core_cf_pangenome_a2[rowSums(venn_core_cf_pangenome_a2)>2,]
non_persistent_core_cf_pangenome_a2 <- venn_core_cf_pangenome_a2[rowSums(venn_core_cf_pangenome_a2)<3,]

# select the next three columns (core species presence/absence per age, 35th abundance percentile)
venn_core_cf_pangenome_a3 <- venn_core_cf_pangenome %>% select(cf_0_summary_a3, cf_1_3_summary_a3, cf_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_pangenome_a3 <- venn_core_cf_pangenome_a3[rowSums(venn_core_cf_pangenome_a3) > 0,]
# define the background species, present in all age groups
background_core_cf_pangenome_a3 <- venn_core_cf_pangenome_a3[rowSums(venn_core_cf_pangenome_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_cf_pangenome_a3 <- venn_core_cf_pangenome_a3[rowSums(venn_core_cf_pangenome_a3)<3,]


# make named list (similar as before)
# rare species, pangenome, 15th abundance percentile, healthy
h_x_rare_pangenome_a1 <- list(A = rownames(c_0_h_pangenome_rare_a1), 
                              B = rownames(c_1_3_h_pangenome_rare_a1), 
                              C = rownames(c_4_6_h_pangenome_rare_a1))
# rare species, pangenome, 15th abundance percentile, CF
cf_x_rare_pangenome_a1 <- list(A = rownames(c_0_cf_pangenome_rare_a1), 
                               B = rownames(c_1_3_cf_pangenome_rare_a1), 
                               C = rownames(c_4_6_cf_pangenome_rare_a1))
# rare species, pangenome, 25th abundance percentile, healthy
h_x_rare_pangenome_a2 <- list(A = rownames(c_0_h_pangenome_rare_a2), 
                              B = rownames(c_1_3_h_pangenome_rare_a2), 
                              C = rownames(c_4_6_h_pangenome_rare_a2))
# rare species, pangenome, 25th abundance percentile, CF
cf_x_rare_pangenome_a2 <- list(A = rownames(c_0_cf_pangenome_rare_a2), 
                               B = rownames(c_1_3_cf_pangenome_rare_a2), 
                               C = rownames(c_4_6_cf_pangenome_rare_a2))
# rare species, pangenome, 35th abundance percentile, healthy
h_x_rare_pangenome_a3 <- list(A = rownames(c_0_h_pangenome_rare_a3), 
                              B = rownames(c_1_3_h_pangenome_rare_a3), 
                              C = rownames(c_4_6_h_pangenome_rare_a3))
# rare species, pangenome, 35th abundance percentile, CF
cf_x_rare_pangenome_a3 <- list(A = rownames(c_0_cf_pangenome_rare_a3), 
                               B = rownames(c_1_3_cf_pangenome_rare_a3), 
                               C = rownames(c_4_6_cf_pangenome_rare_a3))

# core species, pangenome, 15th abundance percentile, healthy
h_x_core_pangenome_a1 <- list(A = rownames(c_0_h_pangenome_core_a1), 
                              B = rownames(c_1_3_h_pangenome_core_a1), 
                              C = rownames(c_4_6_h_pangenome_core_a1))
# core species, pangenome, 15th abundance percentile, CF
cf_x_core_pangenome_a1 <- list(A = rownames(c_0_cf_pangenome_core_a1), 
                               B = rownames(c_1_3_cf_pangenome_core_a1), 
                               C = rownames(c_4_6_cf_pangenome_core_a1))
# core species, pangenome, 25th abundance percentile, healthy
h_x_core_pangenome_a2 <- list(A = rownames(c_0_h_pangenome_core_a2), 
                              B = rownames(c_1_3_h_pangenome_core_a2), 
                              C = rownames(c_4_6_h_pangenome_core_a2))
# core species, pangenome, 25th abundance percentile, CF
cf_x_core_pangenome_a2 <- list(A = rownames(c_0_cf_pangenome_core_a2), 
                               B = rownames(c_1_3_cf_pangenome_core_a2), 
                               C = rownames(c_4_6_cf_pangenome_core_a2))
# core species, pangenome, 35th abundance percentile, healthy
h_x_core_pangenome_a3 <- list(A = rownames(c_0_h_pangenome_core_a3), 
                              B = rownames(c_1_3_h_pangenome_core_a3), 
                              C = rownames(c_4_6_h_pangenome_core_a3))
# core species, pangenome, 35th abundance percentile, CF
cf_x_core_pangenome_a3 <- list(A = rownames(c_0_cf_pangenome_core_a3), 
                               B = rownames(c_1_3_cf_pangenome_core_a3), 
                               C = rownames(c_4_6_cf_pangenome_core_a3))


# Venn diagram analysis, pangenome, RLE
# make vector with all rare species in healthy and CF children, 15th abundance percentile
all_rare_a1_pangenome_rle <- c(healthy_rare_a1_pangenome_rle, cf_rare_a1_pangenome_rle)
# make vector with all core species in healthy and CF children, 15th abundance percentile
all_core_a1_pangenome_rle <- c(healthy_core_a1_pangenome_rle, cf_core_a1_pangenome_rle)
# make vector with all rare species in healthy and CF children, 25th abundance percentile
all_rare_a2_pangenome_rle <- c(healthy_rare_a2_pangenome_rle, cf_rare_a2_pangenome_rle)
# make vector with all core species in healthy and CF children, 25th abundance percentile
all_core_a2_pangenome_rle <- c(healthy_core_a2_pangenome_rle, cf_core_a2_pangenome_rle)
# make vector with all rare species in healthy and CF children, 35th abundance percentile
all_rare_a3_pangenome_rle <- c(healthy_rare_a3_pangenome_rle, cf_rare_a3_pangenome_rle)
# make vector with all core species in healthy and CF children, 35th abundance percentile
all_core_a3_pangenome_rle <- c(healthy_core_a3_pangenome_rle, cf_core_a3_pangenome_rle)

# rename original dataframe (RLE)
venn_rare_h_pangenome_rle <- data.frame(ds_pangenome_rle)
# rare species, healthy, RLE
# 15th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_pangenome_rle$h_0_summary_a1 <- ifelse(rownames(venn_rare_h_pangenome_rle) %in% rownames(c_0_h_pangenome_rle_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_pangenome_rle$h_1_3_summary_a1 <- ifelse(rownames(venn_rare_h_pangenome_rle) %in% rownames(c_1_3_h_pangenome_rle_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_pangenome_rle$h_4_6_summary_a1 <- ifelse(rownames(venn_rare_h_pangenome_rle) %in% rownames(c_4_6_h_pangenome_rle_rare_a1), 1, 0)
# 25th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_pangenome_rle$h_0_summary_a2 <- ifelse(rownames(venn_rare_h_pangenome_rle) %in% rownames(c_0_h_pangenome_rle_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_pangenome_rle$h_1_3_summary_a2 <- ifelse(rownames(venn_rare_h_pangenome_rle) %in% rownames(c_1_3_h_pangenome_rle_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_pangenome_rle$h_4_6_summary_a2 <- ifelse(rownames(venn_rare_h_pangenome_rle) %in% rownames(c_4_6_h_pangenome_rle_rare_a2), 1, 0)
# 35th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_pangenome_rle$h_0_summary_a3 <- ifelse(rownames(venn_rare_h_pangenome_rle) %in% rownames(c_0_h_pangenome_rle_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_pangenome_rle$h_1_3_summary_a3 <- ifelse(rownames(venn_rare_h_pangenome_rle) %in% rownames(c_1_3_h_pangenome_rle_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_pangenome_rle$h_4_6_summary_a3 <- ifelse(rownames(venn_rare_h_pangenome_rle) %in% rownames(c_4_6_h_pangenome_rle_rare_a3), 1, 0)

# select the first three columns (rare species presence/absence per age, 15th abundance percentile)
venn_rare_h_pangenome_rle_a1 <- venn_rare_h_pangenome_rle %>% select(h_0_summary_a1, h_1_3_summary_a1, h_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_pangenome_rle_a1 <- venn_rare_h_pangenome_rle_a1[rowSums(venn_rare_h_pangenome_rle_a1) > 0,]
# define the background species, present in all age groups
background_rare_h_pangenome_rle_a1 <- venn_rare_h_pangenome_rle_a1[rowSums(venn_rare_h_pangenome_rle_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_h_pangenome_rle_a1 <- venn_rare_h_pangenome_rle_a1[rowSums(venn_rare_h_pangenome_rle_a1)<3,]

# select the next three columns (rare species presence/absence per age, 25th abundance percentile)
venn_rare_h_pangenome_rle_a2 <- venn_rare_h_pangenome_rle %>% select(h_0_summary_a2, h_1_3_summary_a2, h_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_pangenome_rle_a2 <- venn_rare_h_pangenome_rle_a2[rowSums(venn_rare_h_pangenome_rle_a2) > 0,]
background_rare_h_pangenome_rle_a2 <- venn_rare_h_pangenome_rle_a2[rowSums(venn_rare_h_pangenome_rle_a2)>2,]
non_persistent_rare_h_pangenome_rle_a2 <- venn_rare_h_pangenome_rle_a2[rowSums(venn_rare_h_pangenome_rle_a2)<3,]

# select the next three columns (rare species presence/absence per age, 35th abundance percentile)
venn_rare_h_pangenome_rle_a3 <- venn_rare_h_pangenome_rle %>% select(h_0_summary_a3, h_1_3_summary_a3, h_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_pangenome_rle_a3 <- venn_rare_h_pangenome_rle_a3[rowSums(venn_rare_h_pangenome_rle_a3) > 0,]
# define the background species, present in all age groups
background_rare_h_pangenome_rle_a3 <- venn_rare_h_pangenome_rle_a3[rowSums(venn_rare_h_pangenome_rle_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_h_pangenome_rle_a3 <- venn_rare_h_pangenome_rle_a3[rowSums(venn_rare_h_pangenome_rle_a3)<3,]

# core species, healthy, RLE
venn_core_h_pangenome_rle <- data.frame(ds_pangenome_rle)
# 15th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_pangenome_rle$h_0_summary_a1 <- ifelse(rownames(venn_core_h_pangenome_rle) %in% rownames(c_0_h_pangenome_rle_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_pangenome_rle$h_1_3_summary_a1 <- ifelse(rownames(venn_core_h_pangenome_rle) %in% rownames(c_1_3_h_pangenome_rle_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_pangenome_rle$h_4_6_summary_a1 <- ifelse(rownames(venn_core_h_pangenome_rle) %in% rownames(c_4_6_h_pangenome_rle_core_a1), 1, 0)
# 25th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_pangenome_rle$h_0_summary_a2 <- ifelse(rownames(venn_core_h_pangenome_rle) %in% rownames(c_0_h_pangenome_rle_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_pangenome_rle$h_1_3_summary_a2 <- ifelse(rownames(venn_core_h_pangenome_rle) %in% rownames(c_1_3_h_pangenome_rle_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_pangenome_rle$h_4_6_summary_a2 <- ifelse(rownames(venn_core_h_pangenome_rle) %in% rownames(c_4_6_h_pangenome_rle_core_a2), 1, 0)
# 35th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_pangenome_rle$h_0_summary_a3 <- ifelse(rownames(venn_core_h_pangenome_rle) %in% rownames(c_0_h_pangenome_rle_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_pangenome_rle$h_1_3_summary_a3 <- ifelse(rownames(venn_core_h_pangenome_rle) %in% rownames(c_1_3_h_pangenome_rle_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_pangenome_rle$h_4_6_summary_a3 <- ifelse(rownames(venn_core_h_pangenome_rle) %in% rownames(c_4_6_h_pangenome_rle_core_a3), 1, 0)

# select the first three columns (core species presence/absence per age, 15th abundance percentile)
venn_core_h_pangenome_rle_a1 <- venn_core_h_pangenome_rle %>% select(h_0_summary_a1, h_1_3_summary_a1, h_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_core_h_pangenome_rle_a1 <- venn_core_h_pangenome_rle_a1[rowSums(venn_core_h_pangenome_rle_a1) > 0,]
# define the background species, present in all age groups
background_core_h_pangenome_rle_a1 <- venn_core_h_pangenome_rle_a1[rowSums(venn_core_h_pangenome_rle_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_h_pangenome_rle_a1 <- venn_core_h_pangenome_rle_a1[rowSums(venn_core_h_pangenome_rle_a1)<3,]

# select the next three columns (core species presence/absence per age, 25th abundance percentile)
venn_core_h_pangenome_rle_a2 <- venn_core_h_pangenome_rle %>% select(h_0_summary_a2, h_1_3_summary_a2, h_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_core_h_pangenome_rle_a2 <- venn_core_h_pangenome_rle_a2[rowSums(venn_core_h_pangenome_rle_a2) > 0,]
background_core_h_pangenome_rle_a2 <- venn_core_h_pangenome_rle_a2[rowSums(venn_core_h_pangenome_rle_a2)>2,]
non_persistent_core_h_pangenome_rle_a2 <- venn_core_h_pangenome_rle_a2[rowSums(venn_core_h_pangenome_rle_a2)<3,]

# select the next three columns (core species presence/absence per age, 35th abundance percentile)
venn_core_h_pangenome_rle_a3 <- venn_core_h_pangenome_rle %>% select(h_0_summary_a3, h_1_3_summary_a3, h_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_core_h_pangenome_rle_a3 <- venn_core_h_pangenome_rle_a3[rowSums(venn_core_h_pangenome_rle_a3) > 0,]
# define the background species, present in all age groups
background_core_h_pangenome_rle_a3 <- venn_core_h_pangenome_rle_a3[rowSums(venn_core_h_pangenome_rle_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_h_pangenome_rle_a3 <- venn_core_h_pangenome_rle_a3[rowSums(venn_core_h_pangenome_rle_a3)<3,]



# repeat with CF children, rare, RLE-normalised, pangenome
# rename original dataframe (RLE)
venn_rare_cf_pangenome_rle <- data.frame(ds_pangenome_rle)
# rare species, cf, RLE
# 15th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_pangenome_rle$cf_0_summary_a1 <- ifelse(rownames(venn_rare_cf_pangenome_rle) %in% rownames(c_0_cf_pangenome_rle_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_pangenome_rle$cf_1_3_summary_a1 <- ifelse(rownames(venn_rare_cf_pangenome_rle) %in% rownames(c_1_3_cf_pangenome_rle_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_pangenome_rle$cf_4_6_summary_a1 <- ifelse(rownames(venn_rare_cf_pangenome_rle) %in% rownames(c_4_6_cf_pangenome_rle_rare_a1), 1, 0)
# 25th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_pangenome_rle$cf_0_summary_a2 <- ifelse(rownames(venn_rare_cf_pangenome_rle) %in% rownames(c_0_cf_pangenome_rle_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_pangenome_rle$cf_1_3_summary_a2 <- ifelse(rownames(venn_rare_cf_pangenome_rle) %in% rownames(c_1_3_cf_pangenome_rle_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_pangenome_rle$cf_4_6_summary_a2 <- ifelse(rownames(venn_rare_cf_pangenome_rle) %in% rownames(c_4_6_cf_pangenome_rle_rare_a2), 1, 0)
# 35th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_pangenome_rle$cf_0_summary_a3 <- ifelse(rownames(venn_rare_cf_pangenome_rle) %in% rownames(c_0_cf_pangenome_rle_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_pangenome_rle$cf_1_3_summary_a3 <- ifelse(rownames(venn_rare_cf_pangenome_rle) %in% rownames(c_1_3_cf_pangenome_rle_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_pangenome_rle$cf_4_6_summary_a3 <- ifelse(rownames(venn_rare_cf_pangenome_rle) %in% rownames(c_4_6_cf_pangenome_rle_rare_a3), 1, 0)

# select the first three columns (rare species presence/absence per age, 15th abundance percentile)
venn_rare_cf_pangenome_rle_a1 <- venn_rare_cf_pangenome_rle %>% select(cf_0_summary_a1, cf_1_3_summary_a1, cf_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_pangenome_rle_a1 <- venn_rare_cf_pangenome_rle_a1[rowSums(venn_rare_cf_pangenome_rle_a1) > 0,]
# define the background species, present in all age groups
background_rare_cf_pangenome_rle_a1 <- venn_rare_cf_pangenome_rle_a1[rowSums(venn_rare_cf_pangenome_rle_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_cf_pangenome_rle_a1 <- venn_rare_cf_pangenome_rle_a1[rowSums(venn_rare_cf_pangenome_rle_a1)<3,]

# select the next three columns (rare species presence/absence per age, 25th abundance percentile)
venn_rare_cf_pangenome_rle_a2 <- venn_rare_cf_pangenome_rle %>% select(cf_0_summary_a2, cf_1_3_summary_a2, cf_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_pangenome_rle_a2 <- venn_rare_cf_pangenome_rle_a2[rowSums(venn_rare_cf_pangenome_rle_a2) > 0,]
background_rare_cf_pangenome_rle_a2 <- venn_rare_cf_pangenome_rle_a2[rowSums(venn_rare_cf_pangenome_rle_a2)>2,]
non_persistent_rare_cf_pangenome_rle_a2 <- venn_rare_cf_pangenome_rle_a2[rowSums(venn_rare_cf_pangenome_rle_a2)<3,]

# select the next three columns (rare species presence/absence per age, 35th abundance percentile)
venn_rare_cf_pangenome_rle_a3 <- venn_rare_cf_pangenome_rle %>% select(cf_0_summary_a3, cf_1_3_summary_a3, cf_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_pangenome_rle_a3 <- venn_rare_cf_pangenome_rle_a3[rowSums(venn_rare_cf_pangenome_rle_a3) > 0,]
# define the background species, present in all age groups
background_rare_cf_pangenome_rle_a3 <- venn_rare_cf_pangenome_rle_a3[rowSums(venn_rare_cf_pangenome_rle_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_cf_pangenome_rle_a3 <- venn_rare_cf_pangenome_rle_a3[rowSums(venn_rare_cf_pangenome_rle_a3)<3,]

# core species, cf, RLE
venn_core_cf_pangenome_rle <- data.frame(ds_pangenome_rle)
# 15th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_pangenome_rle$cf_0_summary_a1 <- ifelse(rownames(venn_core_cf_pangenome_rle) %in% rownames(c_0_cf_pangenome_rle_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_pangenome_rle$cf_1_3_summary_a1 <- ifelse(rownames(venn_core_cf_pangenome_rle) %in% rownames(c_1_3_cf_pangenome_rle_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_pangenome_rle$cf_4_6_summary_a1 <- ifelse(rownames(venn_core_cf_pangenome_rle) %in% rownames(c_4_6_cf_pangenome_rle_core_a1), 1, 0)
# 25th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_pangenome_rle$cf_0_summary_a2 <- ifelse(rownames(venn_core_cf_pangenome_rle) %in% rownames(c_0_cf_pangenome_rle_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_pangenome_rle$cf_1_3_summary_a2 <- ifelse(rownames(venn_core_cf_pangenome_rle) %in% rownames(c_1_3_cf_pangenome_rle_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_pangenome_rle$cf_4_6_summary_a2 <- ifelse(rownames(venn_core_cf_pangenome_rle) %in% rownames(c_4_6_cf_pangenome_rle_core_a2), 1, 0)
# 35th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_pangenome_rle$cf_0_summary_a3 <- ifelse(rownames(venn_core_cf_pangenome_rle) %in% rownames(c_0_cf_pangenome_rle_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_pangenome_rle$cf_1_3_summary_a3 <- ifelse(rownames(venn_core_cf_pangenome_rle) %in% rownames(c_1_3_cf_pangenome_rle_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_pangenome_rle$cf_4_6_summary_a3 <- ifelse(rownames(venn_core_cf_pangenome_rle) %in% rownames(c_4_6_cf_pangenome_rle_core_a3), 1, 0)

# select the first three columns (core species presence/absence per age, 15th abundance percentile)
venn_core_cf_pangenome_rle_a1 <- venn_core_cf_pangenome_rle %>% select(cf_0_summary_a1, cf_1_3_summary_a1, cf_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_pangenome_rle_a1 <- venn_core_cf_pangenome_rle_a1[rowSums(venn_core_cf_pangenome_rle_a1) > 0,]
# define the background species, present in all age groups
background_core_cf_pangenome_rle_a1 <- venn_core_cf_pangenome_rle_a1[rowSums(venn_core_cf_pangenome_rle_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_cf_pangenome_rle_a1 <- venn_core_cf_pangenome_rle_a1[rowSums(venn_core_cf_pangenome_rle_a1)<3,]

# select the next three columns (core species presence/absence per age, 25th abundance percentile)
venn_core_cf_pangenome_rle_a2 <- venn_core_cf_pangenome_rle %>% select(cf_0_summary_a2, cf_1_3_summary_a2, cf_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_pangenome_rle_a2 <- venn_core_cf_pangenome_rle_a2[rowSums(venn_core_cf_pangenome_rle_a2) > 0,]
background_core_cf_pangenome_rle_a2 <- venn_core_cf_pangenome_rle_a2[rowSums(venn_core_cf_pangenome_rle_a2)>2,]
non_persistent_core_cf_pangenome_rle_a2 <- venn_core_cf_pangenome_rle_a2[rowSums(venn_core_cf_pangenome_rle_a2)<3,]

# select the next three columns (core species presence/absence per age, 35th abundance percentile)
venn_core_cf_pangenome_rle_a3 <- venn_core_cf_pangenome_rle %>% select(cf_0_summary_a3, cf_1_3_summary_a3, cf_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_pangenome_rle_a3 <- venn_core_cf_pangenome_rle_a3[rowSums(venn_core_cf_pangenome_rle_a3) > 0,]
# define the background species, present in all age groups
background_core_cf_pangenome_rle_a3 <- venn_core_cf_pangenome_rle_a3[rowSums(venn_core_cf_pangenome_rle_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_cf_pangenome_rle_a3 <- venn_core_cf_pangenome_rle_a3[rowSums(venn_core_cf_pangenome_rle_a3)<3,]


# make named list (similar as before)
# rare species, pangenome, 15th abundance percentile, healthy
h_x_rare_pangenome_rle_a1 <- list(A = rownames(c_0_h_pangenome_rle_rare_a1), 
                              B = rownames(c_1_3_h_pangenome_rle_rare_a1), 
                              C = rownames(c_4_6_h_pangenome_rle_rare_a1))
# rare species, pangenome, 15th abundance percentile, CF
cf_x_rare_pangenome_rle_a1 <- list(A = rownames(c_0_cf_pangenome_rle_rare_a1), 
                               B = rownames(c_1_3_cf_pangenome_rle_rare_a1), 
                               C = rownames(c_4_6_cf_pangenome_rle_rare_a1))
# rare species, pangenome, 25th abundance percentile, healthy
h_x_rare_pangenome_rle_a2 <- list(A = rownames(c_0_h_pangenome_rle_rare_a2), 
                              B = rownames(c_1_3_h_pangenome_rle_rare_a2), 
                              C = rownames(c_4_6_h_pangenome_rle_rare_a2))
# rare species, pangenome, 25th abundance percentile, CF
cf_x_rare_pangenome_rle_a2 <- list(A = rownames(c_0_cf_pangenome_rle_rare_a2), 
                               B = rownames(c_1_3_cf_pangenome_rle_rare_a2), 
                               C = rownames(c_4_6_cf_pangenome_rle_rare_a2))
# rare species, pangenome, 35th abundance percentile, healthy
h_x_rare_pangenome_rle_a3 <- list(A = rownames(c_0_h_pangenome_rle_rare_a3), 
                              B = rownames(c_1_3_h_pangenome_rle_rare_a3), 
                              C = rownames(c_4_6_h_pangenome_rle_rare_a3))
# rare species, pangenome, 35th abundance percentile, CF
cf_x_rare_pangenome_rle_a3 <- list(A = rownames(c_0_cf_pangenome_rle_rare_a3), 
                               B = rownames(c_1_3_cf_pangenome_rle_rare_a3), 
                               C = rownames(c_4_6_cf_pangenome_rle_rare_a3))

# core species, pangenome, 15th abundance percentile, healthy
h_x_core_pangenome_rle_a1 <- list(A = rownames(c_0_h_pangenome_rle_core_a1), 
                              B = rownames(c_1_3_h_pangenome_rle_core_a1), 
                              C = rownames(c_4_6_h_pangenome_rle_core_a1))
# core species, pangenome, 15th abundance percentile, CF
cf_x_core_pangenome_rle_a1 <- list(A = rownames(c_0_cf_pangenome_rle_core_a1), 
                               B = rownames(c_1_3_cf_pangenome_rle_core_a1), 
                               C = rownames(c_4_6_cf_pangenome_rle_core_a1))
# core species, pangenome, 25th abundance percentile, healthy
h_x_core_pangenome_rle_a2 <- list(A = rownames(c_0_h_pangenome_rle_core_a2), 
                              B = rownames(c_1_3_h_pangenome_rle_core_a2), 
                              C = rownames(c_4_6_h_pangenome_rle_core_a2))
# core species, pangenome, 25th abundance percentile, CF
cf_x_core_pangenome_rle_a2 <- list(A = rownames(c_0_cf_pangenome_rle_core_a2), 
                               B = rownames(c_1_3_cf_pangenome_rle_core_a2), 
                               C = rownames(c_4_6_cf_pangenome_rle_core_a2))
# core species, pangenome, 35th abundance percentile, healthy
h_x_core_pangenome_rle_a3 <- list(A = rownames(c_0_h_pangenome_rle_core_a3), 
                              B = rownames(c_1_3_h_pangenome_rle_core_a3), 
                              C = rownames(c_4_6_h_pangenome_rle_core_a3))
# core species, pangenome, 35th abundance percentile, CF
cf_x_core_pangenome_rle_a3 <- list(A = rownames(c_0_cf_pangenome_rle_core_a3), 
                               B = rownames(c_1_3_cf_pangenome_rle_core_a3), 
                               C = rownames(c_4_6_cf_pangenome_rle_core_a3))


# Venn diagram analysis, pangenome, VST
# make vector with all rare species in healthy and CF children, 15th abundance percentile
all_rare_a1_pangenome_vst <- c(healthy_rare_a1_pangenome_vst, cf_rare_a1_pangenome_vst)
# make vector with all core species in healthy and CF children, 15th abundance percentile
all_core_a1_pangenome_vst <- c(healthy_core_a1_pangenome_vst, cf_core_a1_pangenome_vst)
# make vector with all rare species in healthy and CF children, 25th abundance percentile
all_rare_a2_pangenome_vst <- c(healthy_rare_a2_pangenome_vst, cf_rare_a2_pangenome_vst)
# make vector with all core species in healthy and CF children, 25th abundance percentile
all_core_a2_pangenome_vst <- c(healthy_core_a2_pangenome_vst, cf_core_a2_pangenome_vst)
# make vector with all rare species in healthy and CF children, 35th abundance percentile
all_rare_a3_pangenome_vst <- c(healthy_rare_a3_pangenome_vst, cf_rare_a3_pangenome_vst)
# make vector with all core species in healthy and CF children, 35th abundance percentile
all_core_a3_pangenome_vst <- c(healthy_core_a3_pangenome_vst, cf_core_a3_pangenome_vst)

# rename original dataframe (vst)
venn_rare_h_pangenome_vst <- data.frame(ds_pangenome_vst)
# rare species, healthy, vst
# 15th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_pangenome_vst$h_0_summary_a1 <- ifelse(rownames(venn_rare_h_pangenome_vst) %in% rownames(c_0_h_pangenome_vst_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_pangenome_vst$h_1_3_summary_a1 <- ifelse(rownames(venn_rare_h_pangenome_vst) %in% rownames(c_1_3_h_pangenome_vst_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_pangenome_vst$h_4_6_summary_a1 <- ifelse(rownames(venn_rare_h_pangenome_vst) %in% rownames(c_4_6_h_pangenome_vst_rare_a1), 1, 0)
# 25th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_pangenome_vst$h_0_summary_a2 <- ifelse(rownames(venn_rare_h_pangenome_vst) %in% rownames(c_0_h_pangenome_vst_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_pangenome_vst$h_1_3_summary_a2 <- ifelse(rownames(venn_rare_h_pangenome_vst) %in% rownames(c_1_3_h_pangenome_vst_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_pangenome_vst$h_4_6_summary_a2 <- ifelse(rownames(venn_rare_h_pangenome_vst) %in% rownames(c_4_6_h_pangenome_vst_rare_a2), 1, 0)
# 35th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_pangenome_vst$h_0_summary_a3 <- ifelse(rownames(venn_rare_h_pangenome_vst) %in% rownames(c_0_h_pangenome_vst_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_pangenome_vst$h_1_3_summary_a3 <- ifelse(rownames(venn_rare_h_pangenome_vst) %in% rownames(c_1_3_h_pangenome_vst_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_pangenome_vst$h_4_6_summary_a3 <- ifelse(rownames(venn_rare_h_pangenome_vst) %in% rownames(c_4_6_h_pangenome_vst_rare_a3), 1, 0)

# select the first three columns (rare species presence/absence per age, 15th abundance percentile)
venn_rare_h_pangenome_vst_a1 <- venn_rare_h_pangenome_vst %>% select(h_0_summary_a1, h_1_3_summary_a1, h_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_pangenome_vst_a1 <- venn_rare_h_pangenome_vst_a1[rowSums(venn_rare_h_pangenome_vst_a1) > 0,]
# define the background species, present in all age groups
background_rare_h_pangenome_vst_a1 <- venn_rare_h_pangenome_vst_a1[rowSums(venn_rare_h_pangenome_vst_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_h_pangenome_vst_a1 <- venn_rare_h_pangenome_vst_a1[rowSums(venn_rare_h_pangenome_vst_a1)<3,]

# select the next three columns (rare species presence/absence per age, 25th abundance percentile)
venn_rare_h_pangenome_vst_a2 <- venn_rare_h_pangenome_vst %>% select(h_0_summary_a2, h_1_3_summary_a2, h_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_pangenome_vst_a2 <- venn_rare_h_pangenome_vst_a2[rowSums(venn_rare_h_pangenome_vst_a2) > 0,]
background_rare_h_pangenome_vst_a2 <- venn_rare_h_pangenome_vst_a2[rowSums(venn_rare_h_pangenome_vst_a2)>2,]
non_persistent_rare_h_pangenome_vst_a2 <- venn_rare_h_pangenome_vst_a2[rowSums(venn_rare_h_pangenome_vst_a2)<3,]

# select the next three columns (rare species presence/absence per age, 35th abundance percentile)
venn_rare_h_pangenome_vst_a3 <- venn_rare_h_pangenome_vst %>% select(h_0_summary_a3, h_1_3_summary_a3, h_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_pangenome_vst_a3 <- venn_rare_h_pangenome_vst_a3[rowSums(venn_rare_h_pangenome_vst_a3) > 0,]
# define the background species, present in all age groups
background_rare_h_pangenome_vst_a3 <- venn_rare_h_pangenome_vst_a3[rowSums(venn_rare_h_pangenome_vst_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_h_pangenome_vst_a3 <- venn_rare_h_pangenome_vst_a3[rowSums(venn_rare_h_pangenome_vst_a3)<3,]

# core species, healthy, vst
venn_core_h_pangenome_vst <- data.frame(ds_pangenome_vst)
# 15th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_pangenome_vst$h_0_summary_a1 <- ifelse(rownames(venn_core_h_pangenome_vst) %in% rownames(c_0_h_pangenome_vst_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_pangenome_vst$h_1_3_summary_a1 <- ifelse(rownames(venn_core_h_pangenome_vst) %in% rownames(c_1_3_h_pangenome_vst_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_pangenome_vst$h_4_6_summary_a1 <- ifelse(rownames(venn_core_h_pangenome_vst) %in% rownames(c_4_6_h_pangenome_vst_core_a1), 1, 0)
# 25th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_pangenome_vst$h_0_summary_a2 <- ifelse(rownames(venn_core_h_pangenome_vst) %in% rownames(c_0_h_pangenome_vst_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_pangenome_vst$h_1_3_summary_a2 <- ifelse(rownames(venn_core_h_pangenome_vst) %in% rownames(c_1_3_h_pangenome_vst_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_pangenome_vst$h_4_6_summary_a2 <- ifelse(rownames(venn_core_h_pangenome_vst) %in% rownames(c_4_6_h_pangenome_vst_core_a2), 1, 0)
# 35th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_pangenome_vst$h_0_summary_a3 <- ifelse(rownames(venn_core_h_pangenome_vst) %in% rownames(c_0_h_pangenome_vst_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_pangenome_vst$h_1_3_summary_a3 <- ifelse(rownames(venn_core_h_pangenome_vst) %in% rownames(c_1_3_h_pangenome_vst_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_pangenome_vst$h_4_6_summary_a3 <- ifelse(rownames(venn_core_h_pangenome_vst) %in% rownames(c_4_6_h_pangenome_vst_core_a3), 1, 0)

# select the first three columns (core species presence/absence per age, 15th abundance percentile)
venn_core_h_pangenome_vst_a1 <- venn_core_h_pangenome_vst %>% select(h_0_summary_a1, h_1_3_summary_a1, h_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_core_h_pangenome_vst_a1 <- venn_core_h_pangenome_vst_a1[rowSums(venn_core_h_pangenome_vst_a1) > 0,]
# define the background species, present in all age groups
background_core_h_pangenome_vst_a1 <- venn_core_h_pangenome_vst_a1[rowSums(venn_core_h_pangenome_vst_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_h_pangenome_vst_a1 <- venn_core_h_pangenome_vst_a1[rowSums(venn_core_h_pangenome_vst_a1)<3,]

# select the next three columns (core species presence/absence per age, 25th abundance percentile)
venn_core_h_pangenome_vst_a2 <- venn_core_h_pangenome_vst %>% select(h_0_summary_a2, h_1_3_summary_a2, h_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_core_h_pangenome_vst_a2 <- venn_core_h_pangenome_vst_a2[rowSums(venn_core_h_pangenome_vst_a2) > 0,]
background_core_h_pangenome_vst_a2 <- venn_core_h_pangenome_vst_a2[rowSums(venn_core_h_pangenome_vst_a2)>2,]
non_persistent_core_h_pangenome_vst_a2 <- venn_core_h_pangenome_vst_a2[rowSums(venn_core_h_pangenome_vst_a2)<3,]

# select the next three columns (core species presence/absence per age, 35th abundance percentile)
venn_core_h_pangenome_vst_a3 <- venn_core_h_pangenome_vst %>% select(h_0_summary_a3, h_1_3_summary_a3, h_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_core_h_pangenome_vst_a3 <- venn_core_h_pangenome_vst_a3[rowSums(venn_core_h_pangenome_vst_a3) > 0,]
# define the background species, present in all age groups
background_core_h_pangenome_vst_a3 <- venn_core_h_pangenome_vst_a3[rowSums(venn_core_h_pangenome_vst_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_h_pangenome_vst_a3 <- venn_core_h_pangenome_vst_a3[rowSums(venn_core_h_pangenome_vst_a3)<3,]



# repeat with CF children, rare, vst-normalised, pangenome
# rename original dataframe (vst)
venn_rare_cf_pangenome_vst <- data.frame(ds_pangenome_vst)
# rare species, cf, vst
# 15th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_pangenome_vst$cf_0_summary_a1 <- ifelse(rownames(venn_rare_cf_pangenome_vst) %in% rownames(c_0_cf_pangenome_vst_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_pangenome_vst$cf_1_3_summary_a1 <- ifelse(rownames(venn_rare_cf_pangenome_vst) %in% rownames(c_1_3_cf_pangenome_vst_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_pangenome_vst$cf_4_6_summary_a1 <- ifelse(rownames(venn_rare_cf_pangenome_vst) %in% rownames(c_4_6_cf_pangenome_vst_rare_a1), 1, 0)
# 25th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_pangenome_vst$cf_0_summary_a2 <- ifelse(rownames(venn_rare_cf_pangenome_vst) %in% rownames(c_0_cf_pangenome_vst_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_pangenome_vst$cf_1_3_summary_a2 <- ifelse(rownames(venn_rare_cf_pangenome_vst) %in% rownames(c_1_3_cf_pangenome_vst_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_pangenome_vst$cf_4_6_summary_a2 <- ifelse(rownames(venn_rare_cf_pangenome_vst) %in% rownames(c_4_6_cf_pangenome_vst_rare_a2), 1, 0)
# 35th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_pangenome_vst$cf_0_summary_a3 <- ifelse(rownames(venn_rare_cf_pangenome_vst) %in% rownames(c_0_cf_pangenome_vst_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_pangenome_vst$cf_1_3_summary_a3 <- ifelse(rownames(venn_rare_cf_pangenome_vst) %in% rownames(c_1_3_cf_pangenome_vst_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_pangenome_vst$cf_4_6_summary_a3 <- ifelse(rownames(venn_rare_cf_pangenome_vst) %in% rownames(c_4_6_cf_pangenome_vst_rare_a3), 1, 0)

# select the first three columns (rare species presence/absence per age, 15th abundance percentile)
venn_rare_cf_pangenome_vst_a1 <- venn_rare_cf_pangenome_vst %>% select(cf_0_summary_a1, cf_1_3_summary_a1, cf_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_pangenome_vst_a1 <- venn_rare_cf_pangenome_vst_a1[rowSums(venn_rare_cf_pangenome_vst_a1) > 0,]
# define the background species, present in all age groups
background_rare_cf_pangenome_vst_a1 <- venn_rare_cf_pangenome_vst_a1[rowSums(venn_rare_cf_pangenome_vst_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_cf_pangenome_vst_a1 <- venn_rare_cf_pangenome_vst_a1[rowSums(venn_rare_cf_pangenome_vst_a1)<3,]

# select the next three columns (rare species presence/absence per age, 25th abundance percentile)
venn_rare_cf_pangenome_vst_a2 <- venn_rare_cf_pangenome_vst %>% select(cf_0_summary_a2, cf_1_3_summary_a2, cf_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_pangenome_vst_a2 <- venn_rare_cf_pangenome_vst_a2[rowSums(venn_rare_cf_pangenome_vst_a2) > 0,]
background_rare_cf_pangenome_vst_a2 <- venn_rare_cf_pangenome_vst_a2[rowSums(venn_rare_cf_pangenome_vst_a2)>2,]
non_persistent_rare_cf_pangenome_vst_a2 <- venn_rare_cf_pangenome_vst_a2[rowSums(venn_rare_cf_pangenome_vst_a2)<3,]

# select the next three columns (rare species presence/absence per age, 35th abundance percentile)
venn_rare_cf_pangenome_vst_a3 <- venn_rare_cf_pangenome_vst %>% select(cf_0_summary_a3, cf_1_3_summary_a3, cf_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_pangenome_vst_a3 <- venn_rare_cf_pangenome_vst_a3[rowSums(venn_rare_cf_pangenome_vst_a3) > 0,]
# define the background species, present in all age groups
background_rare_cf_pangenome_vst_a3 <- venn_rare_cf_pangenome_vst_a3[rowSums(venn_rare_cf_pangenome_vst_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_cf_pangenome_vst_a3 <- venn_rare_cf_pangenome_vst_a3[rowSums(venn_rare_cf_pangenome_vst_a3)<3,]

# core species, cf, vst
venn_core_cf_pangenome_vst <- data.frame(ds_pangenome_vst)
# 15th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_pangenome_vst$cf_0_summary_a1 <- ifelse(rownames(venn_core_cf_pangenome_vst) %in% rownames(c_0_cf_pangenome_vst_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_pangenome_vst$cf_1_3_summary_a1 <- ifelse(rownames(venn_core_cf_pangenome_vst) %in% rownames(c_1_3_cf_pangenome_vst_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_pangenome_vst$cf_4_6_summary_a1 <- ifelse(rownames(venn_core_cf_pangenome_vst) %in% rownames(c_4_6_cf_pangenome_vst_core_a1), 1, 0)
# 25th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_pangenome_vst$cf_0_summary_a2 <- ifelse(rownames(venn_core_cf_pangenome_vst) %in% rownames(c_0_cf_pangenome_vst_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_pangenome_vst$cf_1_3_summary_a2 <- ifelse(rownames(venn_core_cf_pangenome_vst) %in% rownames(c_1_3_cf_pangenome_vst_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_pangenome_vst$cf_4_6_summary_a2 <- ifelse(rownames(venn_core_cf_pangenome_vst) %in% rownames(c_4_6_cf_pangenome_vst_core_a2), 1, 0)
# 35th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_pangenome_vst$cf_0_summary_a3 <- ifelse(rownames(venn_core_cf_pangenome_vst) %in% rownames(c_0_cf_pangenome_vst_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_pangenome_vst$cf_1_3_summary_a3 <- ifelse(rownames(venn_core_cf_pangenome_vst) %in% rownames(c_1_3_cf_pangenome_vst_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_pangenome_vst$cf_4_6_summary_a3 <- ifelse(rownames(venn_core_cf_pangenome_vst) %in% rownames(c_4_6_cf_pangenome_vst_core_a3), 1, 0)

# select the first three columns (core species presence/absence per age, 15th abundance percentile)
venn_core_cf_pangenome_vst_a1 <- venn_core_cf_pangenome_vst %>% select(cf_0_summary_a1, cf_1_3_summary_a1, cf_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_pangenome_vst_a1 <- venn_core_cf_pangenome_vst_a1[rowSums(venn_core_cf_pangenome_vst_a1) > 0,]
# define the background species, present in all age groups
background_core_cf_pangenome_vst_a1 <- venn_core_cf_pangenome_vst_a1[rowSums(venn_core_cf_pangenome_vst_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_cf_pangenome_vst_a1 <- venn_core_cf_pangenome_vst_a1[rowSums(venn_core_cf_pangenome_vst_a1)<3,]

# select the next three columns (core species presence/absence per age, 25th abundance percentile)
venn_core_cf_pangenome_vst_a2 <- venn_core_cf_pangenome_vst %>% select(cf_0_summary_a2, cf_1_3_summary_a2, cf_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_pangenome_vst_a2 <- venn_core_cf_pangenome_vst_a2[rowSums(venn_core_cf_pangenome_vst_a2) > 0,]
background_core_cf_pangenome_vst_a2 <- venn_core_cf_pangenome_vst_a2[rowSums(venn_core_cf_pangenome_vst_a2)>2,]
non_persistent_core_cf_pangenome_vst_a2 <- venn_core_cf_pangenome_vst_a2[rowSums(venn_core_cf_pangenome_vst_a2)<3,]

# select the next three columns (core species presence/absence per age, 35th abundance percentile)
venn_core_cf_pangenome_vst_a3 <- venn_core_cf_pangenome_vst %>% select(cf_0_summary_a3, cf_1_3_summary_a3, cf_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_pangenome_vst_a3 <- venn_core_cf_pangenome_vst_a3[rowSums(venn_core_cf_pangenome_vst_a3) > 0,]
# define the background species, present in all age groups
background_core_cf_pangenome_vst_a3 <- venn_core_cf_pangenome_vst_a3[rowSums(venn_core_cf_pangenome_vst_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_cf_pangenome_vst_a3 <- venn_core_cf_pangenome_vst_a3[rowSums(venn_core_cf_pangenome_vst_a3)<3,]


# make named list (similar as before)
# rare species, pangenome, 15th abundance percentile, healthy
h_x_rare_pangenome_vst_a1 <- list(A = rownames(c_0_h_pangenome_vst_rare_a1), 
                                  B = rownames(c_1_3_h_pangenome_vst_rare_a1), 
                                  C = rownames(c_4_6_h_pangenome_vst_rare_a1))
# rare species, pangenome, 15th abundance percentile, CF
cf_x_rare_pangenome_vst_a1 <- list(A = rownames(c_0_cf_pangenome_vst_rare_a1), 
                                   B = rownames(c_1_3_cf_pangenome_vst_rare_a1), 
                                   C = rownames(c_4_6_cf_pangenome_vst_rare_a1))
# rare species, pangenome, 25th abundance percentile, healthy
h_x_rare_pangenome_vst_a2 <- list(A = rownames(c_0_h_pangenome_vst_rare_a2), 
                                  B = rownames(c_1_3_h_pangenome_vst_rare_a2), 
                                  C = rownames(c_4_6_h_pangenome_vst_rare_a2))
# rare species, pangenome, 25th abundance percentile, CF
cf_x_rare_pangenome_vst_a2 <- list(A = rownames(c_0_cf_pangenome_vst_rare_a2), 
                                   B = rownames(c_1_3_cf_pangenome_vst_rare_a2), 
                                   C = rownames(c_4_6_cf_pangenome_vst_rare_a2))
# rare species, pangenome, 35th abundance percentile, healthy
h_x_rare_pangenome_vst_a3 <- list(A = rownames(c_0_h_pangenome_vst_rare_a3), 
                                  B = rownames(c_1_3_h_pangenome_vst_rare_a3), 
                                  C = rownames(c_4_6_h_pangenome_vst_rare_a3))
# rare species, pangenome, 35th abundance percentile, CF
cf_x_rare_pangenome_vst_a3 <- list(A = rownames(c_0_cf_pangenome_vst_rare_a3), 
                                   B = rownames(c_1_3_cf_pangenome_vst_rare_a3), 
                                   C = rownames(c_4_6_cf_pangenome_vst_rare_a3))

# core species, pangenome, 15th abundance percentile, healthy
h_x_core_pangenome_vst_a1 <- list(A = rownames(c_0_h_pangenome_vst_core_a1), 
                                  B = rownames(c_1_3_h_pangenome_vst_core_a1), 
                                  C = rownames(c_4_6_h_pangenome_vst_core_a1))
# core species, pangenome, 15th abundance percentile, CF
cf_x_core_pangenome_vst_a1 <- list(A = rownames(c_0_cf_pangenome_vst_core_a1), 
                                   B = rownames(c_1_3_cf_pangenome_vst_core_a1), 
                                   C = rownames(c_4_6_cf_pangenome_vst_core_a1))
# core species, pangenome, 25th abundance percentile, healthy
h_x_core_pangenome_vst_a2 <- list(A = rownames(c_0_h_pangenome_vst_core_a2), 
                                  B = rownames(c_1_3_h_pangenome_vst_core_a2), 
                                  C = rownames(c_4_6_h_pangenome_vst_core_a2))
# core species, pangenome, 25th abundance percentile, CF
cf_x_core_pangenome_vst_a2 <- list(A = rownames(c_0_cf_pangenome_vst_core_a2), 
                                   B = rownames(c_1_3_cf_pangenome_vst_core_a2), 
                                   C = rownames(c_4_6_cf_pangenome_vst_core_a2))
# core species, pangenome, 35th abundance percentile, healthy
h_x_core_pangenome_vst_a3 <- list(A = rownames(c_0_h_pangenome_vst_core_a3), 
                                  B = rownames(c_1_3_h_pangenome_vst_core_a3), 
                                  C = rownames(c_4_6_h_pangenome_vst_core_a3))
# core species, pangenome, 35th abundance percentile, CF
cf_x_core_pangenome_vst_a3 <- list(A = rownames(c_0_cf_pangenome_vst_core_a3), 
                                   B = rownames(c_1_3_cf_pangenome_vst_core_a3), 
                                   C = rownames(c_4_6_cf_pangenome_vst_core_a3))


############################################################################################################
# Venn diagram analysis, repeat with one strain per species database, BCPHC
# make vector with all rare species in healthy and CF children, 15th abundance percentile
all_rare_a1_osps <- c(healthy_rare_a1_osps, cf_rare_a1_osps)
# make vector with all core species in healthy and CF children, 15th abundance percentile
all_core_a1_osps <- c(healthy_core_a1_osps, cf_core_a1_osps)
# make vector with all rare species in healthy and CF children, 25th abundance percentile
all_rare_a2_osps <- c(healthy_rare_a2_osps, cf_rare_a2_osps)
# make vector with all core species in healthy and CF children, 25th abundance percentile
all_core_a2_osps <- c(healthy_core_a2_osps, cf_core_a2_osps)
# make vector with all rare species in healthy and CF children, 35th abundance percentile
all_rare_a3_osps <- c(healthy_rare_a3_osps, cf_rare_a3_osps)
# make vector with all core species in healthy and CF children, 35th abundance percentile
all_core_a3_osps <- c(healthy_core_a3_osps, cf_core_a3_osps)

# rename original dataframe (BCPHC)
venn_rare_h_osps <- data.frame(ds_osps)
# rare species, healthy, BCPHC
# 15th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_osps$h_0_summary_a1 <- ifelse(rownames(venn_rare_h_osps) %in% rownames(c_0_h_osps_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_osps$h_1_3_summary_a1 <- ifelse(rownames(venn_rare_h_osps) %in% rownames(c_1_3_h_osps_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_osps$h_4_6_summary_a1 <- ifelse(rownames(venn_rare_h_osps) %in% rownames(c_4_6_h_osps_rare_a1), 1, 0)
# 25th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_osps$h_0_summary_a2 <- ifelse(rownames(venn_rare_h_osps) %in% rownames(c_0_h_osps_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_osps$h_1_3_summary_a2 <- ifelse(rownames(venn_rare_h_osps) %in% rownames(c_1_3_h_osps_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_osps$h_4_6_summary_a2 <- ifelse(rownames(venn_rare_h_osps) %in% rownames(c_4_6_h_osps_rare_a2), 1, 0)
# 35th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_osps$h_0_summary_a3 <- ifelse(rownames(venn_rare_h_osps) %in% rownames(c_0_h_osps_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_osps$h_1_3_summary_a3 <- ifelse(rownames(venn_rare_h_osps) %in% rownames(c_1_3_h_osps_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_osps$h_4_6_summary_a3 <- ifelse(rownames(venn_rare_h_osps) %in% rownames(c_4_6_h_osps_rare_a3), 1, 0)

# select the first three columns (rare species presence/absence per age, 15th abundance percentile)
venn_rare_h_osps_a1 <- venn_rare_h_osps %>% select(h_0_summary_a1, h_1_3_summary_a1, h_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_osps_a1 <- venn_rare_h_osps_a1[rowSums(venn_rare_h_osps_a1) > 0,]
# define the background species, present in all age groups
background_rare_h_osps_a1 <- venn_rare_h_osps_a1[rowSums(venn_rare_h_osps_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_h_osps_a1 <- venn_rare_h_osps_a1[rowSums(venn_rare_h_osps_a1)<3,]

# select the next three columns (rare species presence/absence per age, 25th abundance percentile)
venn_rare_h_osps_a2 <- venn_rare_h_osps %>% select(h_0_summary_a2, h_1_3_summary_a2, h_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_osps_a2 <- venn_rare_h_osps_a2[rowSums(venn_rare_h_osps_a2) > 0,]
background_rare_h_osps_a2 <- venn_rare_h_osps_a2[rowSums(venn_rare_h_osps_a2)>2,]
non_persistent_rare_h_osps_a2 <- venn_rare_h_osps_a2[rowSums(venn_rare_h_osps_a2)<3,]

# select the next three columns (rare species presence/absence per age, 35th abundance percentile)
venn_rare_h_osps_a3 <- venn_rare_h_osps %>% select(h_0_summary_a3, h_1_3_summary_a3, h_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_osps_a3 <- venn_rare_h_osps_a3[rowSums(venn_rare_h_osps_a3) > 0,]
# define the background species, present in all age groups
background_rare_h_osps_a3 <- venn_rare_h_osps_a3[rowSums(venn_rare_h_osps_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_h_osps_a3 <- venn_rare_h_osps_a3[rowSums(venn_rare_h_osps_a3)<3,]

# core species, healthy, BCPHC
venn_core_h_osps <- data.frame(ds_osps)
# 15th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_osps$h_0_summary_a1 <- ifelse(rownames(venn_core_h_osps) %in% rownames(c_0_h_osps_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_osps$h_1_3_summary_a1 <- ifelse(rownames(venn_core_h_osps) %in% rownames(c_1_3_h_osps_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_osps$h_4_6_summary_a1 <- ifelse(rownames(venn_core_h_osps) %in% rownames(c_4_6_h_osps_core_a1), 1, 0)
# 25th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_osps$h_0_summary_a2 <- ifelse(rownames(venn_core_h_osps) %in% rownames(c_0_h_osps_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_osps$h_1_3_summary_a2 <- ifelse(rownames(venn_core_h_osps) %in% rownames(c_1_3_h_osps_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_osps$h_4_6_summary_a2 <- ifelse(rownames(venn_core_h_osps) %in% rownames(c_4_6_h_osps_core_a2), 1, 0)
# 35th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_osps$h_0_summary_a3 <- ifelse(rownames(venn_core_h_osps) %in% rownames(c_0_h_osps_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_osps$h_1_3_summary_a3 <- ifelse(rownames(venn_core_h_osps) %in% rownames(c_1_3_h_osps_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_osps$h_4_6_summary_a3 <- ifelse(rownames(venn_core_h_osps) %in% rownames(c_4_6_h_osps_core_a3), 1, 0)

# select the first three columns (core species presence/absence per age, 15th abundance percentile)
venn_core_h_osps_a1 <- venn_core_h_osps %>% select(h_0_summary_a1, h_1_3_summary_a1, h_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_core_h_osps_a1 <- venn_core_h_osps_a1[rowSums(venn_core_h_osps_a1) > 0,]
# define the background species, present in all age groups
background_core_h_osps_a1 <- venn_core_h_osps_a1[rowSums(venn_core_h_osps_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_h_osps_a1 <- venn_core_h_osps_a1[rowSums(venn_core_h_osps_a1)<3,]

# select the next three columns (core species presence/absence per age, 25th abundance percentile)
venn_core_h_osps_a2 <- venn_core_h_osps %>% select(h_0_summary_a2, h_1_3_summary_a2, h_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_core_h_osps_a2 <- venn_core_h_osps_a2[rowSums(venn_core_h_osps_a2) > 0,]
background_core_h_osps_a2 <- venn_core_h_osps_a2[rowSums(venn_core_h_osps_a2)>2,]
non_persistent_core_h_osps_a2 <- venn_core_h_osps_a2[rowSums(venn_core_h_osps_a2)<3,]

# select the next three columns (core species presence/absence per age, 35th abundance percentile)
venn_core_h_osps_a3 <- venn_core_h_osps %>% select(h_0_summary_a3, h_1_3_summary_a3, h_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_core_h_osps_a3 <- venn_core_h_osps_a3[rowSums(venn_core_h_osps_a3) > 0,]
# define the background species, present in all age groups
background_core_h_osps_a3 <- venn_core_h_osps_a3[rowSums(venn_core_h_osps_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_h_osps_a3 <- venn_core_h_osps_a3[rowSums(venn_core_h_osps_a3)<3,]


# repeat with CF children, rare, BCPHC-normalised, osps
# rename original dataframe (BCPHC)
venn_rare_cf_osps <- data.frame(ds_osps)
# rare species, cf, BCPHC
# 15th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_osps$cf_0_summary_a1 <- ifelse(rownames(venn_rare_cf_osps) %in% rownames(c_0_cf_osps_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_osps$cf_1_3_summary_a1 <- ifelse(rownames(venn_rare_cf_osps) %in% rownames(c_1_3_cf_osps_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_osps$cf_4_6_summary_a1 <- ifelse(rownames(venn_rare_cf_osps) %in% rownames(c_4_6_cf_osps_rare_a1), 1, 0)
# 25th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_osps$cf_0_summary_a2 <- ifelse(rownames(venn_rare_cf_osps) %in% rownames(c_0_cf_osps_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_osps$cf_1_3_summary_a2 <- ifelse(rownames(venn_rare_cf_osps) %in% rownames(c_1_3_cf_osps_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_osps$cf_4_6_summary_a2 <- ifelse(rownames(venn_rare_cf_osps) %in% rownames(c_4_6_cf_osps_rare_a2), 1, 0)
# 35th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_osps$cf_0_summary_a3 <- ifelse(rownames(venn_rare_cf_osps) %in% rownames(c_0_cf_osps_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_osps$cf_1_3_summary_a3 <- ifelse(rownames(venn_rare_cf_osps) %in% rownames(c_1_3_cf_osps_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_osps$cf_4_6_summary_a3 <- ifelse(rownames(venn_rare_cf_osps) %in% rownames(c_4_6_cf_osps_rare_a3), 1, 0)

# select the first three columns (rare species presence/absence per age, 15th abundance percentile)
venn_rare_cf_osps_a1 <- venn_rare_cf_osps %>% select(cf_0_summary_a1, cf_1_3_summary_a1, cf_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_osps_a1 <- venn_rare_cf_osps_a1[rowSums(venn_rare_cf_osps_a1) > 0,]
# define the background species, present in all age groups
background_rare_cf_osps_a1 <- venn_rare_cf_osps_a1[rowSums(venn_rare_cf_osps_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_cf_osps_a1 <- venn_rare_cf_osps_a1[rowSums(venn_rare_cf_osps_a1)<3,]

# select the next three columns (rare species presence/absence per age, 25th abundance percentile)
venn_rare_cf_osps_a2 <- venn_rare_cf_osps %>% select(cf_0_summary_a2, cf_1_3_summary_a2, cf_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_osps_a2 <- venn_rare_cf_osps_a2[rowSums(venn_rare_cf_osps_a2) > 0,]
background_rare_cf_osps_a2 <- venn_rare_cf_osps_a2[rowSums(venn_rare_cf_osps_a2)>2,]
non_persistent_rare_cf_osps_a2 <- venn_rare_cf_osps_a2[rowSums(venn_rare_cf_osps_a2)<3,]

# select the next three columns (rare species presence/absence per age, 35th abundance percentile)
venn_rare_cf_osps_a3 <- venn_rare_cf_osps %>% select(cf_0_summary_a3, cf_1_3_summary_a3, cf_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_osps_a3 <- venn_rare_cf_osps_a3[rowSums(venn_rare_cf_osps_a3) > 0,]
# define the background species, present in all age groups
background_rare_cf_osps_a3 <- venn_rare_cf_osps_a3[rowSums(venn_rare_cf_osps_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_cf_osps_a3 <- venn_rare_cf_osps_a3[rowSums(venn_rare_cf_osps_a3)<3,]

# core species, cf, BCPHC
venn_core_cf_osps <- data.frame(ds_osps)
# 15th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_osps$cf_0_summary_a1 <- ifelse(rownames(venn_core_cf_osps) %in% rownames(c_0_cf_osps_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_osps$cf_1_3_summary_a1 <- ifelse(rownames(venn_core_cf_osps) %in% rownames(c_1_3_cf_osps_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_osps$cf_4_6_summary_a1 <- ifelse(rownames(venn_core_cf_osps) %in% rownames(c_4_6_cf_osps_core_a1), 1, 0)
# 25th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_osps$cf_0_summary_a2 <- ifelse(rownames(venn_core_cf_osps) %in% rownames(c_0_cf_osps_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_osps$cf_1_3_summary_a2 <- ifelse(rownames(venn_core_cf_osps) %in% rownames(c_1_3_cf_osps_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_osps$cf_4_6_summary_a2 <- ifelse(rownames(venn_core_cf_osps) %in% rownames(c_4_6_cf_osps_core_a2), 1, 0)
# 35th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_osps$cf_0_summary_a3 <- ifelse(rownames(venn_core_cf_osps) %in% rownames(c_0_cf_osps_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_osps$cf_1_3_summary_a3 <- ifelse(rownames(venn_core_cf_osps) %in% rownames(c_1_3_cf_osps_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_osps$cf_4_6_summary_a3 <- ifelse(rownames(venn_core_cf_osps) %in% rownames(c_4_6_cf_osps_core_a3), 1, 0)

# select the first three columns (core species presence/absence per age, 15th abundance percentile)
venn_core_cf_osps_a1 <- venn_core_cf_osps %>% select(cf_0_summary_a1, cf_1_3_summary_a1, cf_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_osps_a1 <- venn_core_cf_osps_a1[rowSums(venn_core_cf_osps_a1) > 0,]
# define the background species, present in all age groups
background_core_cf_osps_a1 <- venn_core_cf_osps_a1[rowSums(venn_core_cf_osps_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_cf_osps_a1 <- venn_core_cf_osps_a1[rowSums(venn_core_cf_osps_a1)<3,]

# select the next three columns (core species presence/absence per age, 25th abundance percentile)
venn_core_cf_osps_a2 <- venn_core_cf_osps %>% select(cf_0_summary_a2, cf_1_3_summary_a2, cf_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_osps_a2 <- venn_core_cf_osps_a2[rowSums(venn_core_cf_osps_a2) > 0,]
background_core_cf_osps_a2 <- venn_core_cf_osps_a2[rowSums(venn_core_cf_osps_a2)>2,]
non_persistent_core_cf_osps_a2 <- venn_core_cf_osps_a2[rowSums(venn_core_cf_osps_a2)<3,]

# select the next three columns (core species presence/absence per age, 35th abundance percentile)
venn_core_cf_osps_a3 <- venn_core_cf_osps %>% select(cf_0_summary_a3, cf_1_3_summary_a3, cf_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_osps_a3 <- venn_core_cf_osps_a3[rowSums(venn_core_cf_osps_a3) > 0,]
# define the background species, present in all age groups
background_core_cf_osps_a3 <- venn_core_cf_osps_a3[rowSums(venn_core_cf_osps_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_cf_osps_a3 <- venn_core_cf_osps_a3[rowSums(venn_core_cf_osps_a3)<3,]


# make named list (similar as before)
# rare species, osps, 15th abundance percentile, healthy
h_x_rare_osps_a1 <- list(A = rownames(c_0_h_osps_rare_a1), 
                              B = rownames(c_1_3_h_osps_rare_a1), 
                              C = rownames(c_4_6_h_osps_rare_a1))
# rare species, osps, 15th abundance percentile, CF
cf_x_rare_osps_a1 <- list(A = rownames(c_0_cf_osps_rare_a1), 
                               B = rownames(c_1_3_cf_osps_rare_a1), 
                               C = rownames(c_4_6_cf_osps_rare_a1))
# rare species, osps, 25th abundance percentile, healthy
h_x_rare_osps_a2 <- list(A = rownames(c_0_h_osps_rare_a2), 
                              B = rownames(c_1_3_h_osps_rare_a2), 
                              C = rownames(c_4_6_h_osps_rare_a2))
# rare species, osps, 25th abundance percentile, CF
cf_x_rare_osps_a2 <- list(A = rownames(c_0_cf_osps_rare_a2), 
                               B = rownames(c_1_3_cf_osps_rare_a2), 
                               C = rownames(c_4_6_cf_osps_rare_a2))
# rare species, osps, 35th abundance percentile, healthy
h_x_rare_osps_a3 <- list(A = rownames(c_0_h_osps_rare_a3), 
                              B = rownames(c_1_3_h_osps_rare_a3), 
                              C = rownames(c_4_6_h_osps_rare_a3))
# rare species, osps, 35th abundance percentile, CF
cf_x_rare_osps_a3 <- list(A = rownames(c_0_cf_osps_rare_a3), 
                               B = rownames(c_1_3_cf_osps_rare_a3), 
                               C = rownames(c_4_6_cf_osps_rare_a3))

# core species, osps, 15th abundance percentile, healthy
h_x_core_osps_a1 <- list(A = rownames(c_0_h_osps_core_a1), 
                              B = rownames(c_1_3_h_osps_core_a1), 
                              C = rownames(c_4_6_h_osps_core_a1))
# core species, osps, 15th abundance percentile, CF
cf_x_core_osps_a1 <- list(A = rownames(c_0_cf_osps_core_a1), 
                               B = rownames(c_1_3_cf_osps_core_a1), 
                               C = rownames(c_4_6_cf_osps_core_a1))
# core species, osps, 25th abundance percentile, healthy
h_x_core_osps_a2 <- list(A = rownames(c_0_h_osps_core_a2), 
                              B = rownames(c_1_3_h_osps_core_a2), 
                              C = rownames(c_4_6_h_osps_core_a2))
# core species, osps, 25th abundance percentile, CF
cf_x_core_osps_a2 <- list(A = rownames(c_0_cf_osps_core_a2), 
                               B = rownames(c_1_3_cf_osps_core_a2), 
                               C = rownames(c_4_6_cf_osps_core_a2))
# core species, osps, 35th abundance percentile, healthy
h_x_core_osps_a3 <- list(A = rownames(c_0_h_osps_core_a3), 
                              B = rownames(c_1_3_h_osps_core_a3), 
                              C = rownames(c_4_6_h_osps_core_a3))
# core species, osps, 35th abundance percentile, CF
cf_x_core_osps_a3 <- list(A = rownames(c_0_cf_osps_core_a3), 
                               B = rownames(c_1_3_cf_osps_core_a3), 
                               C = rownames(c_4_6_cf_osps_core_a3))


# Venn diagram analysis, osps, RLE
# make vector with all rare species in healthy and CF children, 15th abundance percentile
all_rare_a1_osps_rle <- c(healthy_rare_a1_osps_rle, cf_rare_a1_osps_rle)
# make vector with all core species in healthy and CF children, 15th abundance percentile
all_core_a1_osps_rle <- c(healthy_core_a1_osps_rle, cf_core_a1_osps_rle)
# make vector with all rare species in healthy and CF children, 25th abundance percentile
all_rare_a2_osps_rle <- c(healthy_rare_a2_osps_rle, cf_rare_a2_osps_rle)
# make vector with all core species in healthy and CF children, 25th abundance percentile
all_core_a2_osps_rle <- c(healthy_core_a2_osps_rle, cf_core_a2_osps_rle)
# make vector with all rare species in healthy and CF children, 35th abundance percentile
all_rare_a3_osps_rle <- c(healthy_rare_a3_osps_rle, cf_rare_a3_osps_rle)
# make vector with all core species in healthy and CF children, 35th abundance percentile
all_core_a3_osps_rle <- c(healthy_core_a3_osps_rle, cf_core_a3_osps_rle)

# rename original dataframe (RLE)
venn_rare_h_osps_rle <- data.frame(ds_osps_rle)
# rare species, healthy, RLE
# 15th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_osps_rle$h_0_summary_a1 <- ifelse(rownames(venn_rare_h_osps_rle) %in% rownames(c_0_h_osps_rle_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_osps_rle$h_1_3_summary_a1 <- ifelse(rownames(venn_rare_h_osps_rle) %in% rownames(c_1_3_h_osps_rle_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_osps_rle$h_4_6_summary_a1 <- ifelse(rownames(venn_rare_h_osps_rle) %in% rownames(c_4_6_h_osps_rle_rare_a1), 1, 0)
# 25th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_osps_rle$h_0_summary_a2 <- ifelse(rownames(venn_rare_h_osps_rle) %in% rownames(c_0_h_osps_rle_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_osps_rle$h_1_3_summary_a2 <- ifelse(rownames(venn_rare_h_osps_rle) %in% rownames(c_1_3_h_osps_rle_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_osps_rle$h_4_6_summary_a2 <- ifelse(rownames(venn_rare_h_osps_rle) %in% rownames(c_4_6_h_osps_rle_rare_a2), 1, 0)
# 35th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_osps_rle$h_0_summary_a3 <- ifelse(rownames(venn_rare_h_osps_rle) %in% rownames(c_0_h_osps_rle_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_osps_rle$h_1_3_summary_a3 <- ifelse(rownames(venn_rare_h_osps_rle) %in% rownames(c_1_3_h_osps_rle_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_osps_rle$h_4_6_summary_a3 <- ifelse(rownames(venn_rare_h_osps_rle) %in% rownames(c_4_6_h_osps_rle_rare_a3), 1, 0)

# select the first three columns (rare species presence/absence per age, 15th abundance percentile)
venn_rare_h_osps_rle_a1 <- venn_rare_h_osps_rle %>% select(h_0_summary_a1, h_1_3_summary_a1, h_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_osps_rle_a1 <- venn_rare_h_osps_rle_a1[rowSums(venn_rare_h_osps_rle_a1) > 0,]
# define the background species, present in all age groups
background_rare_h_osps_rle_a1 <- venn_rare_h_osps_rle_a1[rowSums(venn_rare_h_osps_rle_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_h_osps_rle_a1 <- venn_rare_h_osps_rle_a1[rowSums(venn_rare_h_osps_rle_a1)<3,]

# select the next three columns (rare species presence/absence per age, 25th abundance percentile)
venn_rare_h_osps_rle_a2 <- venn_rare_h_osps_rle %>% select(h_0_summary_a2, h_1_3_summary_a2, h_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_osps_rle_a2 <- venn_rare_h_osps_rle_a2[rowSums(venn_rare_h_osps_rle_a2) > 0,]
background_rare_h_osps_rle_a2 <- venn_rare_h_osps_rle_a2[rowSums(venn_rare_h_osps_rle_a2)>2,]
non_persistent_rare_h_osps_rle_a2 <- venn_rare_h_osps_rle_a2[rowSums(venn_rare_h_osps_rle_a2)<3,]

# select the next three columns (rare species presence/absence per age, 35th abundance percentile)
venn_rare_h_osps_rle_a3 <- venn_rare_h_osps_rle %>% select(h_0_summary_a3, h_1_3_summary_a3, h_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_osps_rle_a3 <- venn_rare_h_osps_rle_a3[rowSums(venn_rare_h_osps_rle_a3) > 0,]
# define the background species, present in all age groups
background_rare_h_osps_rle_a3 <- venn_rare_h_osps_rle_a3[rowSums(venn_rare_h_osps_rle_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_h_osps_rle_a3 <- venn_rare_h_osps_rle_a3[rowSums(venn_rare_h_osps_rle_a3)<3,]

# core species, healthy, RLE
venn_core_h_osps_rle <- data.frame(ds_osps_rle)
# 15th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_osps_rle$h_0_summary_a1 <- ifelse(rownames(venn_core_h_osps_rle) %in% rownames(c_0_h_osps_rle_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_osps_rle$h_1_3_summary_a1 <- ifelse(rownames(venn_core_h_osps_rle) %in% rownames(c_1_3_h_osps_rle_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_osps_rle$h_4_6_summary_a1 <- ifelse(rownames(venn_core_h_osps_rle) %in% rownames(c_4_6_h_osps_rle_core_a1), 1, 0)
# 25th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_osps_rle$h_0_summary_a2 <- ifelse(rownames(venn_core_h_osps_rle) %in% rownames(c_0_h_osps_rle_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_osps_rle$h_1_3_summary_a2 <- ifelse(rownames(venn_core_h_osps_rle) %in% rownames(c_1_3_h_osps_rle_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_osps_rle$h_4_6_summary_a2 <- ifelse(rownames(venn_core_h_osps_rle) %in% rownames(c_4_6_h_osps_rle_core_a2), 1, 0)
# 35th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_osps_rle$h_0_summary_a3 <- ifelse(rownames(venn_core_h_osps_rle) %in% rownames(c_0_h_osps_rle_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_osps_rle$h_1_3_summary_a3 <- ifelse(rownames(venn_core_h_osps_rle) %in% rownames(c_1_3_h_osps_rle_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_osps_rle$h_4_6_summary_a3 <- ifelse(rownames(venn_core_h_osps_rle) %in% rownames(c_4_6_h_osps_rle_core_a3), 1, 0)

# select the first three columns (core species presence/absence per age, 15th abundance percentile)
venn_core_h_osps_rle_a1 <- venn_core_h_osps_rle %>% select(h_0_summary_a1, h_1_3_summary_a1, h_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_core_h_osps_rle_a1 <- venn_core_h_osps_rle_a1[rowSums(venn_core_h_osps_rle_a1) > 0,]
# define the background species, present in all age groups
background_core_h_osps_rle_a1 <- venn_core_h_osps_rle_a1[rowSums(venn_core_h_osps_rle_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_h_osps_rle_a1 <- venn_core_h_osps_rle_a1[rowSums(venn_core_h_osps_rle_a1)<3,]

# select the next three columns (core species presence/absence per age, 25th abundance percentile)
venn_core_h_osps_rle_a2 <- venn_core_h_osps_rle %>% select(h_0_summary_a2, h_1_3_summary_a2, h_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_core_h_osps_rle_a2 <- venn_core_h_osps_rle_a2[rowSums(venn_core_h_osps_rle_a2) > 0,]
background_core_h_osps_rle_a2 <- venn_core_h_osps_rle_a2[rowSums(venn_core_h_osps_rle_a2)>2,]
non_persistent_core_h_osps_rle_a2 <- venn_core_h_osps_rle_a2[rowSums(venn_core_h_osps_rle_a2)<3,]

# select the next three columns (core species presence/absence per age, 35th abundance percentile)
venn_core_h_osps_rle_a3 <- venn_core_h_osps_rle %>% select(h_0_summary_a3, h_1_3_summary_a3, h_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_core_h_osps_rle_a3 <- venn_core_h_osps_rle_a3[rowSums(venn_core_h_osps_rle_a3) > 0,]
# define the background species, present in all age groups
background_core_h_osps_rle_a3 <- venn_core_h_osps_rle_a3[rowSums(venn_core_h_osps_rle_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_h_osps_rle_a3 <- venn_core_h_osps_rle_a3[rowSums(venn_core_h_osps_rle_a3)<3,]



# repeat with CF children, rare, RLE-normalised, osps
# rename original dataframe (RLE)
venn_rare_cf_osps_rle <- data.frame(ds_osps_rle)
# rare species, cf, RLE
# 15th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_osps_rle$cf_0_summary_a1 <- ifelse(rownames(venn_rare_cf_osps_rle) %in% rownames(c_0_cf_osps_rle_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_osps_rle$cf_1_3_summary_a1 <- ifelse(rownames(venn_rare_cf_osps_rle) %in% rownames(c_1_3_cf_osps_rle_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_osps_rle$cf_4_6_summary_a1 <- ifelse(rownames(venn_rare_cf_osps_rle) %in% rownames(c_4_6_cf_osps_rle_rare_a1), 1, 0)
# 25th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_osps_rle$cf_0_summary_a2 <- ifelse(rownames(venn_rare_cf_osps_rle) %in% rownames(c_0_cf_osps_rle_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_osps_rle$cf_1_3_summary_a2 <- ifelse(rownames(venn_rare_cf_osps_rle) %in% rownames(c_1_3_cf_osps_rle_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_osps_rle$cf_4_6_summary_a2 <- ifelse(rownames(venn_rare_cf_osps_rle) %in% rownames(c_4_6_cf_osps_rle_rare_a2), 1, 0)
# 35th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_osps_rle$cf_0_summary_a3 <- ifelse(rownames(venn_rare_cf_osps_rle) %in% rownames(c_0_cf_osps_rle_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_osps_rle$cf_1_3_summary_a3 <- ifelse(rownames(venn_rare_cf_osps_rle) %in% rownames(c_1_3_cf_osps_rle_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_osps_rle$cf_4_6_summary_a3 <- ifelse(rownames(venn_rare_cf_osps_rle) %in% rownames(c_4_6_cf_osps_rle_rare_a3), 1, 0)

# select the first three columns (rare species presence/absence per age, 15th abundance percentile)
venn_rare_cf_osps_rle_a1 <- venn_rare_cf_osps_rle %>% select(cf_0_summary_a1, cf_1_3_summary_a1, cf_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_osps_rle_a1 <- venn_rare_cf_osps_rle_a1[rowSums(venn_rare_cf_osps_rle_a1) > 0,]
# define the background species, present in all age groups
background_rare_cf_osps_rle_a1 <- venn_rare_cf_osps_rle_a1[rowSums(venn_rare_cf_osps_rle_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_cf_osps_rle_a1 <- venn_rare_cf_osps_rle_a1[rowSums(venn_rare_cf_osps_rle_a1)<3,]

# select the next three columns (rare species presence/absence per age, 25th abundance percentile)
venn_rare_cf_osps_rle_a2 <- venn_rare_cf_osps_rle %>% select(cf_0_summary_a2, cf_1_3_summary_a2, cf_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_osps_rle_a2 <- venn_rare_cf_osps_rle_a2[rowSums(venn_rare_cf_osps_rle_a2) > 0,]
background_rare_cf_osps_rle_a2 <- venn_rare_cf_osps_rle_a2[rowSums(venn_rare_cf_osps_rle_a2)>2,]
non_persistent_rare_cf_osps_rle_a2 <- venn_rare_cf_osps_rle_a2[rowSums(venn_rare_cf_osps_rle_a2)<3,]

# select the next three columns (rare species presence/absence per age, 35th abundance percentile)
venn_rare_cf_osps_rle_a3 <- venn_rare_cf_osps_rle %>% select(cf_0_summary_a3, cf_1_3_summary_a3, cf_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_osps_rle_a3 <- venn_rare_cf_osps_rle_a3[rowSums(venn_rare_cf_osps_rle_a3) > 0,]
# define the background species, present in all age groups
background_rare_cf_osps_rle_a3 <- venn_rare_cf_osps_rle_a3[rowSums(venn_rare_cf_osps_rle_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_cf_osps_rle_a3 <- venn_rare_cf_osps_rle_a3[rowSums(venn_rare_cf_osps_rle_a3)<3,]

# core species, cf, RLE
venn_core_cf_osps_rle <- data.frame(ds_osps_rle)
# 15th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_osps_rle$cf_0_summary_a1 <- ifelse(rownames(venn_core_cf_osps_rle) %in% rownames(c_0_cf_osps_rle_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_osps_rle$cf_1_3_summary_a1 <- ifelse(rownames(venn_core_cf_osps_rle) %in% rownames(c_1_3_cf_osps_rle_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_osps_rle$cf_4_6_summary_a1 <- ifelse(rownames(venn_core_cf_osps_rle) %in% rownames(c_4_6_cf_osps_rle_core_a1), 1, 0)
# 25th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_osps_rle$cf_0_summary_a2 <- ifelse(rownames(venn_core_cf_osps_rle) %in% rownames(c_0_cf_osps_rle_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_osps_rle$cf_1_3_summary_a2 <- ifelse(rownames(venn_core_cf_osps_rle) %in% rownames(c_1_3_cf_osps_rle_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_osps_rle$cf_4_6_summary_a2 <- ifelse(rownames(venn_core_cf_osps_rle) %in% rownames(c_4_6_cf_osps_rle_core_a2), 1, 0)
# 35th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_osps_rle$cf_0_summary_a3 <- ifelse(rownames(venn_core_cf_osps_rle) %in% rownames(c_0_cf_osps_rle_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_osps_rle$cf_1_3_summary_a3 <- ifelse(rownames(venn_core_cf_osps_rle) %in% rownames(c_1_3_cf_osps_rle_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_osps_rle$cf_4_6_summary_a3 <- ifelse(rownames(venn_core_cf_osps_rle) %in% rownames(c_4_6_cf_osps_rle_core_a3), 1, 0)

# select the first three columns (core species presence/absence per age, 15th abundance percentile)
venn_core_cf_osps_rle_a1 <- venn_core_cf_osps_rle %>% select(cf_0_summary_a1, cf_1_3_summary_a1, cf_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_osps_rle_a1 <- venn_core_cf_osps_rle_a1[rowSums(venn_core_cf_osps_rle_a1) > 0,]
# define the background species, present in all age groups
background_core_cf_osps_rle_a1 <- venn_core_cf_osps_rle_a1[rowSums(venn_core_cf_osps_rle_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_cf_osps_rle_a1 <- venn_core_cf_osps_rle_a1[rowSums(venn_core_cf_osps_rle_a1)<3,]

# select the next three columns (core species presence/absence per age, 25th abundance percentile)
venn_core_cf_osps_rle_a2 <- venn_core_cf_osps_rle %>% select(cf_0_summary_a2, cf_1_3_summary_a2, cf_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_osps_rle_a2 <- venn_core_cf_osps_rle_a2[rowSums(venn_core_cf_osps_rle_a2) > 0,]
background_core_cf_osps_rle_a2 <- venn_core_cf_osps_rle_a2[rowSums(venn_core_cf_osps_rle_a2)>2,]
non_persistent_core_cf_osps_rle_a2 <- venn_core_cf_osps_rle_a2[rowSums(venn_core_cf_osps_rle_a2)<3,]

# select the next three columns (core species presence/absence per age, 35th abundance percentile)
venn_core_cf_osps_rle_a3 <- venn_core_cf_osps_rle %>% select(cf_0_summary_a3, cf_1_3_summary_a3, cf_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_osps_rle_a3 <- venn_core_cf_osps_rle_a3[rowSums(venn_core_cf_osps_rle_a3) > 0,]
# define the background species, present in all age groups
background_core_cf_osps_rle_a3 <- venn_core_cf_osps_rle_a3[rowSums(venn_core_cf_osps_rle_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_cf_osps_rle_a3 <- venn_core_cf_osps_rle_a3[rowSums(venn_core_cf_osps_rle_a3)<3,]


# make named list (similar as before)
# rare species, osps, 15th abundance percentile, healthy
h_x_rare_osps_rle_a1 <- list(A = rownames(c_0_h_osps_rle_rare_a1), 
                                  B = rownames(c_1_3_h_osps_rle_rare_a1), 
                                  C = rownames(c_4_6_h_osps_rle_rare_a1))
# rare species, osps, 15th abundance percentile, CF
cf_x_rare_osps_rle_a1 <- list(A = rownames(c_0_cf_osps_rle_rare_a1), 
                                   B = rownames(c_1_3_cf_osps_rle_rare_a1), 
                                   C = rownames(c_4_6_cf_osps_rle_rare_a1))
# rare species, osps, 25th abundance percentile, healthy
h_x_rare_osps_rle_a2 <- list(A = rownames(c_0_h_osps_rle_rare_a2), 
                                  B = rownames(c_1_3_h_osps_rle_rare_a2), 
                                  C = rownames(c_4_6_h_osps_rle_rare_a2))
# rare species, osps, 25th abundance percentile, CF
cf_x_rare_osps_rle_a2 <- list(A = rownames(c_0_cf_osps_rle_rare_a2), 
                                   B = rownames(c_1_3_cf_osps_rle_rare_a2), 
                                   C = rownames(c_4_6_cf_osps_rle_rare_a2))
# rare species, osps, 35th abundance percentile, healthy
h_x_rare_osps_rle_a3 <- list(A = rownames(c_0_h_osps_rle_rare_a3), 
                                  B = rownames(c_1_3_h_osps_rle_rare_a3), 
                                  C = rownames(c_4_6_h_osps_rle_rare_a3))
# rare species, osps, 35th abundance percentile, CF
cf_x_rare_osps_rle_a3 <- list(A = rownames(c_0_cf_osps_rle_rare_a3), 
                                   B = rownames(c_1_3_cf_osps_rle_rare_a3), 
                                   C = rownames(c_4_6_cf_osps_rle_rare_a3))

# core species, osps, 15th abundance percentile, healthy
h_x_core_osps_rle_a1 <- list(A = rownames(c_0_h_osps_rle_core_a1), 
                                  B = rownames(c_1_3_h_osps_rle_core_a1), 
                                  C = rownames(c_4_6_h_osps_rle_core_a1))
# core species, osps, 15th abundance percentile, CF
cf_x_core_osps_rle_a1 <- list(A = rownames(c_0_cf_osps_rle_core_a1), 
                                   B = rownames(c_1_3_cf_osps_rle_core_a1), 
                                   C = rownames(c_4_6_cf_osps_rle_core_a1))
# core species, osps, 25th abundance percentile, healthy
h_x_core_osps_rle_a2 <- list(A = rownames(c_0_h_osps_rle_core_a2), 
                                  B = rownames(c_1_3_h_osps_rle_core_a2), 
                                  C = rownames(c_4_6_h_osps_rle_core_a2))
# core species, osps, 25th abundance percentile, CF
cf_x_core_osps_rle_a2 <- list(A = rownames(c_0_cf_osps_rle_core_a2), 
                                   B = rownames(c_1_3_cf_osps_rle_core_a2), 
                                   C = rownames(c_4_6_cf_osps_rle_core_a2))
# core species, osps, 35th abundance percentile, healthy
h_x_core_osps_rle_a3 <- list(A = rownames(c_0_h_osps_rle_core_a3), 
                                  B = rownames(c_1_3_h_osps_rle_core_a3), 
                                  C = rownames(c_4_6_h_osps_rle_core_a3))
# core species, osps, 35th abundance percentile, CF
cf_x_core_osps_rle_a3 <- list(A = rownames(c_0_cf_osps_rle_core_a3), 
                                   B = rownames(c_1_3_cf_osps_rle_core_a3), 
                                   C = rownames(c_4_6_cf_osps_rle_core_a3))


# Venn diagram analysis, osps, VST
# make vector with all rare species in healthy and CF children, 15th abundance percentile
all_rare_a1_osps_vst <- c(healthy_rare_a1_osps_vst, cf_rare_a1_osps_vst)
# make vector with all core species in healthy and CF children, 15th abundance percentile
all_core_a1_osps_vst <- c(healthy_core_a1_osps_vst, cf_core_a1_osps_vst)
# make vector with all rare species in healthy and CF children, 25th abundance percentile
all_rare_a2_osps_vst <- c(healthy_rare_a2_osps_vst, cf_rare_a2_osps_vst)
# make vector with all core species in healthy and CF children, 25th abundance percentile
all_core_a2_osps_vst <- c(healthy_core_a2_osps_vst, cf_core_a2_osps_vst)
# make vector with all rare species in healthy and CF children, 35th abundance percentile
all_rare_a3_osps_vst <- c(healthy_rare_a3_osps_vst, cf_rare_a3_osps_vst)
# make vector with all core species in healthy and CF children, 35th abundance percentile
all_core_a3_osps_vst <- c(healthy_core_a3_osps_vst, cf_core_a3_osps_vst)

# rename original dataframe (vst)
venn_rare_h_osps_vst <- data.frame(ds_osps_vst)
# rare species, healthy, vst
# 15th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_osps_vst$h_0_summary_a1 <- ifelse(rownames(venn_rare_h_osps_vst) %in% rownames(c_0_h_osps_vst_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_osps_vst$h_1_3_summary_a1 <- ifelse(rownames(venn_rare_h_osps_vst) %in% rownames(c_1_3_h_osps_vst_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_osps_vst$h_4_6_summary_a1 <- ifelse(rownames(venn_rare_h_osps_vst) %in% rownames(c_4_6_h_osps_vst_rare_a1), 1, 0)
# 25th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_osps_vst$h_0_summary_a2 <- ifelse(rownames(venn_rare_h_osps_vst) %in% rownames(c_0_h_osps_vst_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_osps_vst$h_1_3_summary_a2 <- ifelse(rownames(venn_rare_h_osps_vst) %in% rownames(c_1_3_h_osps_vst_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_osps_vst$h_4_6_summary_a2 <- ifelse(rownames(venn_rare_h_osps_vst) %in% rownames(c_4_6_h_osps_vst_rare_a2), 1, 0)
# 35th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_h_osps_vst$h_0_summary_a3 <- ifelse(rownames(venn_rare_h_osps_vst) %in% rownames(c_0_h_osps_vst_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_h_osps_vst$h_1_3_summary_a3 <- ifelse(rownames(venn_rare_h_osps_vst) %in% rownames(c_1_3_h_osps_vst_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_h_osps_vst$h_4_6_summary_a3 <- ifelse(rownames(venn_rare_h_osps_vst) %in% rownames(c_4_6_h_osps_vst_rare_a3), 1, 0)

# select the first three columns (rare species presence/absence per age, 15th abundance percentile)
venn_rare_h_osps_vst_a1 <- venn_rare_h_osps_vst %>% select(h_0_summary_a1, h_1_3_summary_a1, h_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_osps_vst_a1 <- venn_rare_h_osps_vst_a1[rowSums(venn_rare_h_osps_vst_a1) > 0,]
# define the background species, present in all age groups
background_rare_h_osps_vst_a1 <- venn_rare_h_osps_vst_a1[rowSums(venn_rare_h_osps_vst_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_h_osps_vst_a1 <- venn_rare_h_osps_vst_a1[rowSums(venn_rare_h_osps_vst_a1)<3,]

# select the next three columns (rare species presence/absence per age, 25th abundance percentile)
venn_rare_h_osps_vst_a2 <- venn_rare_h_osps_vst %>% select(h_0_summary_a2, h_1_3_summary_a2, h_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_osps_vst_a2 <- venn_rare_h_osps_vst_a2[rowSums(venn_rare_h_osps_vst_a2) > 0,]
background_rare_h_osps_vst_a2 <- venn_rare_h_osps_vst_a2[rowSums(venn_rare_h_osps_vst_a2)>2,]
non_persistent_rare_h_osps_vst_a2 <- venn_rare_h_osps_vst_a2[rowSums(venn_rare_h_osps_vst_a2)<3,]

# select the next three columns (rare species presence/absence per age, 35th abundance percentile)
venn_rare_h_osps_vst_a3 <- venn_rare_h_osps_vst %>% select(h_0_summary_a3, h_1_3_summary_a3, h_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_rare_h_osps_vst_a3 <- venn_rare_h_osps_vst_a3[rowSums(venn_rare_h_osps_vst_a3) > 0,]
# define the background species, present in all age groups
background_rare_h_osps_vst_a3 <- venn_rare_h_osps_vst_a3[rowSums(venn_rare_h_osps_vst_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_h_osps_vst_a3 <- venn_rare_h_osps_vst_a3[rowSums(venn_rare_h_osps_vst_a3)<3,]

# core species, healthy, vst
venn_core_h_osps_vst <- data.frame(ds_osps_vst)
# 15th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_osps_vst$h_0_summary_a1 <- ifelse(rownames(venn_core_h_osps_vst) %in% rownames(c_0_h_osps_vst_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_osps_vst$h_1_3_summary_a1 <- ifelse(rownames(venn_core_h_osps_vst) %in% rownames(c_1_3_h_osps_vst_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_osps_vst$h_4_6_summary_a1 <- ifelse(rownames(venn_core_h_osps_vst) %in% rownames(c_4_6_h_osps_vst_core_a1), 1, 0)
# 25th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_osps_vst$h_0_summary_a2 <- ifelse(rownames(venn_core_h_osps_vst) %in% rownames(c_0_h_osps_vst_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_osps_vst$h_1_3_summary_a2 <- ifelse(rownames(venn_core_h_osps_vst) %in% rownames(c_1_3_h_osps_vst_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_osps_vst$h_4_6_summary_a2 <- ifelse(rownames(venn_core_h_osps_vst) %in% rownames(c_4_6_h_osps_vst_core_a2), 1, 0)
# 35th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_h_osps_vst$h_0_summary_a3 <- ifelse(rownames(venn_core_h_osps_vst) %in% rownames(c_0_h_osps_vst_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_h_osps_vst$h_1_3_summary_a3 <- ifelse(rownames(venn_core_h_osps_vst) %in% rownames(c_1_3_h_osps_vst_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_h_osps_vst$h_4_6_summary_a3 <- ifelse(rownames(venn_core_h_osps_vst) %in% rownames(c_4_6_h_osps_vst_core_a3), 1, 0)

# select the first three columns (core species presence/absence per age, 15th abundance percentile)
venn_core_h_osps_vst_a1 <- venn_core_h_osps_vst %>% select(h_0_summary_a1, h_1_3_summary_a1, h_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_core_h_osps_vst_a1 <- venn_core_h_osps_vst_a1[rowSums(venn_core_h_osps_vst_a1) > 0,]
# define the background species, present in all age groups
background_core_h_osps_vst_a1 <- venn_core_h_osps_vst_a1[rowSums(venn_core_h_osps_vst_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_h_osps_vst_a1 <- venn_core_h_osps_vst_a1[rowSums(venn_core_h_osps_vst_a1)<3,]

# select the next three columns (core species presence/absence per age, 25th abundance percentile)
venn_core_h_osps_vst_a2 <- venn_core_h_osps_vst %>% select(h_0_summary_a2, h_1_3_summary_a2, h_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_core_h_osps_vst_a2 <- venn_core_h_osps_vst_a2[rowSums(venn_core_h_osps_vst_a2) > 0,]
background_core_h_osps_vst_a2 <- venn_core_h_osps_vst_a2[rowSums(venn_core_h_osps_vst_a2)>2,]
non_persistent_core_h_osps_vst_a2 <- venn_core_h_osps_vst_a2[rowSums(venn_core_h_osps_vst_a2)<3,]

# select the next three columns (core species presence/absence per age, 35th abundance percentile)
venn_core_h_osps_vst_a3 <- venn_core_h_osps_vst %>% select(h_0_summary_a3, h_1_3_summary_a3, h_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_core_h_osps_vst_a3 <- venn_core_h_osps_vst_a3[rowSums(venn_core_h_osps_vst_a3) > 0,]
# define the background species, present in all age groups
background_core_h_osps_vst_a3 <- venn_core_h_osps_vst_a3[rowSums(venn_core_h_osps_vst_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_h_osps_vst_a3 <- venn_core_h_osps_vst_a3[rowSums(venn_core_h_osps_vst_a3)<3,]


# repeat with CF children, rare, vst-normalised, osps
# rename original dataframe (vst)
venn_rare_cf_osps_vst <- data.frame(ds_osps_vst)
# rare species, cf, vst
# 15th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_osps_vst$cf_0_summary_a1 <- ifelse(rownames(venn_rare_cf_osps_vst) %in% rownames(c_0_cf_osps_vst_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_osps_vst$cf_1_3_summary_a1 <- ifelse(rownames(venn_rare_cf_osps_vst) %in% rownames(c_1_3_cf_osps_vst_rare_a1), 1, 0)
# 15th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_osps_vst$cf_4_6_summary_a1 <- ifelse(rownames(venn_rare_cf_osps_vst) %in% rownames(c_4_6_cf_osps_vst_rare_a1), 1, 0)
# 25th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_osps_vst$cf_0_summary_a2 <- ifelse(rownames(venn_rare_cf_osps_vst) %in% rownames(c_0_cf_osps_vst_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_osps_vst$cf_1_3_summary_a2 <- ifelse(rownames(venn_rare_cf_osps_vst) %in% rownames(c_1_3_cf_osps_vst_rare_a2), 1, 0)
# 25th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_osps_vst$cf_4_6_summary_a2 <- ifelse(rownames(venn_rare_cf_osps_vst) %in% rownames(c_4_6_cf_osps_vst_rare_a2), 1, 0)
# 35th abundance percentile, if rare species is present in first year of life, add 1, otherwise add 0
venn_rare_cf_osps_vst$cf_0_summary_a3 <- ifelse(rownames(venn_rare_cf_osps_vst) %in% rownames(c_0_cf_osps_vst_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in toddler age group, add 1, otherwise add 0
venn_rare_cf_osps_vst$cf_1_3_summary_a3 <- ifelse(rownames(venn_rare_cf_osps_vst) %in% rownames(c_1_3_cf_osps_vst_rare_a3), 1, 0)
# 35th abundance percentile, if rare species is present in preschool age group, add 1, otherwise add 0
venn_rare_cf_osps_vst$cf_4_6_summary_a3 <- ifelse(rownames(venn_rare_cf_osps_vst) %in% rownames(c_4_6_cf_osps_vst_rare_a3), 1, 0)

# select the first three columns (rare species presence/absence per age, 15th abundance percentile)
venn_rare_cf_osps_vst_a1 <- venn_rare_cf_osps_vst %>% select(cf_0_summary_a1, cf_1_3_summary_a1, cf_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_osps_vst_a1 <- venn_rare_cf_osps_vst_a1[rowSums(venn_rare_cf_osps_vst_a1) > 0,]
# define the background species, present in all age groups
background_rare_cf_osps_vst_a1 <- venn_rare_cf_osps_vst_a1[rowSums(venn_rare_cf_osps_vst_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_cf_osps_vst_a1 <- venn_rare_cf_osps_vst_a1[rowSums(venn_rare_cf_osps_vst_a1)<3,]

# select the next three columns (rare species presence/absence per age, 25th abundance percentile)
venn_rare_cf_osps_vst_a2 <- venn_rare_cf_osps_vst %>% select(cf_0_summary_a2, cf_1_3_summary_a2, cf_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_osps_vst_a2 <- venn_rare_cf_osps_vst_a2[rowSums(venn_rare_cf_osps_vst_a2) > 0,]
background_rare_cf_osps_vst_a2 <- venn_rare_cf_osps_vst_a2[rowSums(venn_rare_cf_osps_vst_a2)>2,]
non_persistent_rare_cf_osps_vst_a2 <- venn_rare_cf_osps_vst_a2[rowSums(venn_rare_cf_osps_vst_a2)<3,]

# select the next three columns (rare species presence/absence per age, 35th abundance percentile)
venn_rare_cf_osps_vst_a3 <- venn_rare_cf_osps_vst %>% select(cf_0_summary_a3, cf_1_3_summary_a3, cf_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_rare_cf_osps_vst_a3 <- venn_rare_cf_osps_vst_a3[rowSums(venn_rare_cf_osps_vst_a3) > 0,]
# define the background species, present in all age groups
background_rare_cf_osps_vst_a3 <- venn_rare_cf_osps_vst_a3[rowSums(venn_rare_cf_osps_vst_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_rare_cf_osps_vst_a3 <- venn_rare_cf_osps_vst_a3[rowSums(venn_rare_cf_osps_vst_a3)<3,]

# core species, cf, vst
venn_core_cf_osps_vst <- data.frame(ds_osps_vst)
# 15th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_osps_vst$cf_0_summary_a1 <- ifelse(rownames(venn_core_cf_osps_vst) %in% rownames(c_0_cf_osps_vst_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_osps_vst$cf_1_3_summary_a1 <- ifelse(rownames(venn_core_cf_osps_vst) %in% rownames(c_1_3_cf_osps_vst_core_a1), 1, 0)
# 15th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_osps_vst$cf_4_6_summary_a1 <- ifelse(rownames(venn_core_cf_osps_vst) %in% rownames(c_4_6_cf_osps_vst_core_a1), 1, 0)
# 25th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_osps_vst$cf_0_summary_a2 <- ifelse(rownames(venn_core_cf_osps_vst) %in% rownames(c_0_cf_osps_vst_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_osps_vst$cf_1_3_summary_a2 <- ifelse(rownames(venn_core_cf_osps_vst) %in% rownames(c_1_3_cf_osps_vst_core_a2), 1, 0)
# 25th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_osps_vst$cf_4_6_summary_a2 <- ifelse(rownames(venn_core_cf_osps_vst) %in% rownames(c_4_6_cf_osps_vst_core_a2), 1, 0)
# 35th abundance percentile, if core species is present in first year of life, add 1, otherwise add 0
venn_core_cf_osps_vst$cf_0_summary_a3 <- ifelse(rownames(venn_core_cf_osps_vst) %in% rownames(c_0_cf_osps_vst_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in toddler age group, add 1, otherwise add 0
venn_core_cf_osps_vst$cf_1_3_summary_a3 <- ifelse(rownames(venn_core_cf_osps_vst) %in% rownames(c_1_3_cf_osps_vst_core_a3), 1, 0)
# 35th abundance percentile, if core species is present in preschool age group, add 1, otherwise add 0
venn_core_cf_osps_vst$cf_4_6_summary_a3 <- ifelse(rownames(venn_core_cf_osps_vst) %in% rownames(c_4_6_cf_osps_vst_core_a3), 1, 0)

# select the first three columns (core species presence/absence per age, 15th abundance percentile)
venn_core_cf_osps_vst_a1 <- venn_core_cf_osps_vst %>% select(cf_0_summary_a1, cf_1_3_summary_a1, cf_4_6_summary_a1)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_osps_vst_a1 <- venn_core_cf_osps_vst_a1[rowSums(venn_core_cf_osps_vst_a1) > 0,]
# define the background species, present in all age groups
background_core_cf_osps_vst_a1 <- venn_core_cf_osps_vst_a1[rowSums(venn_core_cf_osps_vst_a1)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_cf_osps_vst_a1 <- venn_core_cf_osps_vst_a1[rowSums(venn_core_cf_osps_vst_a1)<3,]

# select the next three columns (core species presence/absence per age, 25th abundance percentile)
venn_core_cf_osps_vst_a2 <- venn_core_cf_osps_vst %>% select(cf_0_summary_a2, cf_1_3_summary_a2, cf_4_6_summary_a2)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_osps_vst_a2 <- venn_core_cf_osps_vst_a2[rowSums(venn_core_cf_osps_vst_a2) > 0,]
background_core_cf_osps_vst_a2 <- venn_core_cf_osps_vst_a2[rowSums(venn_core_cf_osps_vst_a2)>2,]
non_persistent_core_cf_osps_vst_a2 <- venn_core_cf_osps_vst_a2[rowSums(venn_core_cf_osps_vst_a2)<3,]

# select the next three columns (core species presence/absence per age, 35th abundance percentile)
venn_core_cf_osps_vst_a3 <- venn_core_cf_osps_vst %>% select(cf_0_summary_a3, cf_1_3_summary_a3, cf_4_6_summary_a3)
# remove all rows which sum to zero (species is never detected)
venn_core_cf_osps_vst_a3 <- venn_core_cf_osps_vst_a3[rowSums(venn_core_cf_osps_vst_a3) > 0,]
# define the background species, present in all age groups
background_core_cf_osps_vst_a3 <- venn_core_cf_osps_vst_a3[rowSums(venn_core_cf_osps_vst_a3)>2,]
# define the non-persisting species, present in some age groups but not all
non_persistent_core_cf_osps_vst_a3 <- venn_core_cf_osps_vst_a3[rowSums(venn_core_cf_osps_vst_a3)<3,]


# make named list (similar as before)
# rare species, osps, 15th abundance percentile, healthy
h_x_rare_osps_vst_a1 <- list(A = rownames(c_0_h_osps_vst_rare_a1), 
                                  B = rownames(c_1_3_h_osps_vst_rare_a1), 
                                  C = rownames(c_4_6_h_osps_vst_rare_a1))
# rare species, osps, 15th abundance percentile, CF
cf_x_rare_osps_vst_a1 <- list(A = rownames(c_0_cf_osps_vst_rare_a1), 
                                   B = rownames(c_1_3_cf_osps_vst_rare_a1), 
                                   C = rownames(c_4_6_cf_osps_vst_rare_a1))
# rare species, osps, 25th abundance percentile, healthy
h_x_rare_osps_vst_a2 <- list(A = rownames(c_0_h_osps_vst_rare_a2), 
                                  B = rownames(c_1_3_h_osps_vst_rare_a2), 
                                  C = rownames(c_4_6_h_osps_vst_rare_a2))

# rare species, osps, 25th abundance percentile, CF
cf_x_rare_osps_vst_a2 <- list(A = rownames(c_0_cf_osps_vst_rare_a2), 
                                   B = rownames(c_1_3_cf_osps_vst_rare_a2), 
                                   C = rownames(c_4_6_cf_osps_vst_rare_a2))
# rare species, osps, 35th abundance percentile, healthy
h_x_rare_osps_vst_a3 <- list(A = rownames(c_0_h_osps_vst_rare_a3), 
                                  B = rownames(c_1_3_h_osps_vst_rare_a3), 
                                  C = rownames(c_4_6_h_osps_vst_rare_a3))
# rare species, osps, 35th abundance percentile, CF
cf_x_rare_osps_vst_a3 <- list(A = rownames(c_0_cf_osps_vst_rare_a3), 
                                   B = rownames(c_1_3_cf_osps_vst_rare_a3), 
                                   C = rownames(c_4_6_cf_osps_vst_rare_a3))

# core species, osps, 15th abundance percentile, healthy
h_x_core_osps_vst_a1 <- list(A = rownames(c_0_h_osps_vst_core_a1), 
                                  B = rownames(c_1_3_h_osps_vst_core_a1), 
                                  C = rownames(c_4_6_h_osps_vst_core_a1))
# core species, osps, 15th abundance percentile, CF
cf_x_core_osps_vst_a1 <- list(A = rownames(c_0_cf_osps_vst_core_a1), 
                                   B = rownames(c_1_3_cf_osps_vst_core_a1), 
                                   C = rownames(c_4_6_cf_osps_vst_core_a1))
# core species, osps, 25th abundance percentile, healthy
h_x_core_osps_vst_a2 <- list(A = rownames(c_0_h_osps_vst_core_a2), 
                                  B = rownames(c_1_3_h_osps_vst_core_a2), 
                                  C = rownames(c_4_6_h_osps_vst_core_a2))
# core species, osps, 25th abundance percentile, CF
cf_x_core_osps_vst_a2 <- list(A = rownames(c_0_cf_osps_vst_core_a2), 
                                   B = rownames(c_1_3_cf_osps_vst_core_a2), 
                                   C = rownames(c_4_6_cf_osps_vst_core_a2))
# core species, osps, 35th abundance percentile, healthy
h_x_core_osps_vst_a3 <- list(A = rownames(c_0_h_osps_vst_core_a3), 
                                  B = rownames(c_1_3_h_osps_vst_core_a3), 
                                  C = rownames(c_4_6_h_osps_vst_core_a3))
# core species, osps, 35th abundance percentile, CF
cf_x_core_osps_vst_a3 <- list(A = rownames(c_0_cf_osps_vst_core_a3), 
                                   B = rownames(c_1_3_cf_osps_vst_core_a3), 
                                   C = rownames(c_4_6_cf_osps_vst_core_a3))

############################################################################################################
# Make list of background species per database, normalisation method, disease state
# the list can be exported for further analysis later on
# CF, make list of background core species, for 25th percentile, BCPHC-normalisation, overlaps between both databases
background_cf_core_a2_bcphc <- intersect(rownames(background_core_cf_osps_a2), rownames(background_core_cf_pangenome_a2))
background_cf_core_a2_bcphc
# convert to data frame
background_a2_1 <- data.frame(background_cf_core_a2_bcphc)
# add meta data
background_a2_1$State <- 'CF'
background_a2_1$Species_type <- 'background_core'
background_a2_1$Threshold <- '25th percentile'
background_a2_1$Normalisation <- 'BCPHC'
# rename colnames
colnames(background_a2_1) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# CF, make list of background core species, for 25th percentile, VST-normalisation, overlaps between both databases
background_cf_core_a2_vst <- intersect(rownames(background_core_cf_osps_vst_a2), rownames(background_core_cf_pangenome_vst_a2))
# convert to data frame
background_a2_2 <- data.frame(background_cf_core_a2_vst)
# add meta data
background_a2_2$State <- 'CF'
background_a2_2$Species_type <- 'background_core'
background_a2_2$Threshold <- '25th percentile'
background_a2_2$Normalisation <- 'VST'
# rename colnames
colnames(background_a2_2) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# CF, make list of background core species, for 25th percentile, RLE-normalisation, overlaps between both databases
background_cf_core_a2_rle <- intersect(rownames(background_core_cf_osps_rle_a2), rownames(background_core_cf_pangenome_rle_a2))
# convert to data frame
background_a2_3 <- data.frame(background_cf_core_a2_rle)
# add meta data
background_a2_3$State <- 'CF'
background_a2_3$Species_type <- 'background_core'
background_a2_3$Threshold <- '25th percentile'
background_a2_3$Normalisation <- 'RLE'
# rename colnames
colnames(background_a2_3) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# CF, make list of background rare species, for 25th percentile, BCPHC-normalisation, overlaps between both databases
background_cf_rare_a2_bcphc <- intersect(rownames(background_rare_cf_osps_a2), rownames(background_rare_cf_pangenome_a2))
# convert to data frame
background_a2_4 <- data.frame(background_cf_rare_a2_bcphc)
# add meta data
background_a2_4$State <- 'CF'
background_a2_4$Species_type <- 'background_rare'
background_a2_4$Threshold <- '25th percentile'
background_a2_4$Normalisation <- 'BCPHC'
# rename colnames
colnames(background_a2_4) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# CF, make list of background rare species, for 25th percentile, VST-normalisation, overlaps between both databases
background_cf_rare_a2_vst <- c(rownames(background_rare_cf_osps_vst_a2), rownames(background_rare_cf_pangenome_vst_a2))
# no rare species were detected, no further meta data generation possible

# make list of background rare species, for 25th percentile, RLE-normalisation, overlaps between both databases
background_cf_rare_a2_rle <- intersect(rownames(background_rare_cf_osps_rle_a2), rownames(background_rare_cf_pangenome_rle_a2))
# convert to data frame
background_a2_6 <- data.frame(background_cf_rare_a2_rle)
# add meta data
background_a2_6$State <- 'CF'
background_a2_6$Species_type <- 'background_rare'
background_a2_6$Threshold <- '25th percentile'
background_a2_6$Normalisation <- 'RLE'
# rename colnames
colnames(background_a2_6) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background core species, for 25th percentile, BCPHC-normalisation, overlaps between both databases
background_healthy_core_a2_bcphc <- intersect(rownames(background_core_h_osps_a2), rownames(background_core_h_pangenome_a2))
# convert to data frame
background_a2_7 <- data.frame(background_healthy_core_a2_bcphc)
# add meta data
background_a2_7$State <- 'Healthy'
background_a2_7$Species_type <- 'background_core'
background_a2_7$Threshold <- '25th percentile'
background_a2_7$Normalisation <- 'BCPHC'
# rename colnames
colnames(background_a2_7) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background core species, for 25th percentile, VST-normalisation, overlaps between both databases
background_healthy_core_a2_vst <- intersect(rownames(background_core_h_osps_vst_a2), rownames(background_core_h_pangenome_vst_a2))
# convert to data frame
background_a2_8 <- data.frame(background_healthy_core_a2_vst)
# add meta data
background_a2_8$State <- 'Healthy'
background_a2_8$Species_type <- 'background_core'
background_a2_8$Threshold <- '25th percentile'
background_a2_8$Normalisation <- 'VST'
# rename colnames
colnames(background_a2_8) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background core species, for 25th percentile, RLE-normalisation, overlaps between both databases
background_healthy_core_a2_rle <- intersect(rownames(background_core_h_osps_rle_a2), rownames(background_core_h_pangenome_rle_a2))
# convert to data frame
background_a2_9 <- data.frame(background_healthy_core_a2_rle)
# add meta data
background_a2_9$State <- 'Healthy'
background_a2_9$Species_type <- 'background_core'
background_a2_9$Threshold <- '25th percentile'
background_a2_9$Normalisation <- 'RLE'
# rename colnames
colnames(background_a2_9) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background rare species, for 25th percentile, BCPHC-normalisation, overlaps between both databases
background_healthy_rare_a2_bcphc <- intersect(rownames(background_rare_h_osps_a2), rownames(background_rare_h_pangenome_a2))
# convert to data frame
background_a2_10 <- data.frame(background_healthy_rare_a2_bcphc)
# add meta data
background_a2_10$State <- 'Healthy'
background_a2_10$Species_type <- 'background_rare'
background_a2_10$Threshold <- '25th percentile'
background_a2_10$Normalisation <- 'BCPHC'
# rename colnames
colnames(background_a2_10) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background rare species, for 25th percentile, VST-normalisation, overlaps between both databases
background_healthy_rare_a2_vst <- intersect(rownames(background_rare_h_osps_vst_a2),rownames(background_rare_h_pangenome_vst_a2))
# convert to data frame
background_a2_11 <- data.frame(background_healthy_rare_a2_vst)
# add meta data
background_a2_11$State <- 'Healthy'
background_a2_11$Species_type <- 'background_rare'
background_a2_11$Threshold <- '25th percentile'
background_a2_11$Normalisation <- 'VST'
# rename colnames
colnames(background_a2_11) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background rare species, for 25th percentile, RLE-normalisation, overlaps between both databases
background_healthy_rare_a2_rle <- intersect(rownames(background_rare_h_osps_rle_a2), rownames(background_rare_h_pangenome_rle_a2))
# convert to data frame
background_a2_12 <- data.frame(background_healthy_rare_a2_rle)
# add meta data
background_a2_12$State <- 'Healthy'
background_a2_12$Species_type <- 'background_rare'
background_a2_12$Threshold <- '25th percentile'
background_a2_12$Normalisation <- 'RLE'
# rename colnames
colnames(background_a2_12) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')


# make list of background species, for 15th percentile
# CF, make list of background core species, for 15th percentile, BCPHC-normalisation, overlaps between both databases
background_cf_core_a1_bcphc <- intersect(rownames(background_core_cf_osps_a1), rownames(background_core_cf_pangenome_a1))
# convert to data frame
background_a1_1 <- data.frame(background_cf_core_a1_bcphc)
# add meta data
background_a1_1$State <- 'CF'
background_a1_1$Species_type <- 'background_core'
background_a1_1$Threshold <- '15th percentile'
background_a1_1$Normalisation <- 'BCPHC'
# rename column
colnames(background_a1_1) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# CF, make list of background core species, for 15th percentile, VST-normalisation, overlaps between both databases
background_cf_core_a1_vst <- intersect(rownames(background_core_cf_osps_vst_a1),rownames(background_core_cf_pangenome_vst_a1))
# convert to data frame
background_a1_2 <- data.frame(background_cf_core_a1_vst)
# add meta data
background_a1_2$State <- 'CF'
background_a1_2$Species_type <- 'background_core'
background_a1_2$Threshold <- '15th percentile'
background_a1_2$Normalisation <- 'VST'
# rename column
colnames(background_a1_2) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')


# CF, make list of background core species, for 15th percentile, RLE-normalisation, overlaps between both databases
background_cf_core_a1_rle <- c(rownames(background_core_cf_osps_rle_a1), rownames(background_core_cf_pangenome_rle_a1))
# convert to data frame
background_a1_3 <- data.frame(background_cf_core_a1_rle)
# add meta data
background_a1_3$State <- 'CF'
background_a1_3$Species_type <- 'background_core'
background_a1_3$Threshold <- '15th percentile'
background_a1_3$Normalisation <- 'RLE'
# rename columns
colnames(background_a1_3) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# CF, make list of background rare species, for 15th percentile, BCPHC-normalisation, overlaps between both databases
background_cf_rare_a1_bcphc <- intersect(rownames(background_rare_cf_osps_a1), rownames(background_rare_cf_pangenome_a1))
# convert to data frame
background_a1_4 <- data.frame(background_cf_rare_a1_bcphc)
# add meta data
background_a1_4$State <- 'CF'
background_a1_4$Species_type <- 'background_rare'
background_a1_4$Threshold <- '15th percentile'
background_a1_4$Normalisation <- 'BCPHC'
# rename columns
colnames(background_a1_4) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# CF, make list of background rare species, for 15th percentile, VST-normalisation, overlaps between both databases
background_cf_rare_a1_vst <- intersect(rownames(background_rare_cf_osps_vst_a1), rownames(background_rare_cf_pangenome_vst_a1))
# convert to data frame
background_a1_5 <- data.frame(background_cf_rare_a1_vst)
# add meta data
background_a1_5$State <- 'CF'
background_a1_5$Species_type <- 'background_rare'
background_a1_5$Threshold <- '15th percentile'
background_a1_5$Normalisation <- 'VST'
# rename columns
colnames(background_a1_5) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')


# CF, make list of background rare species, for 15th percentile, RLE-normalisation, overlaps between both databases
background_cf_rare_a1_rle <- intersect(rownames(background_rare_cf_osps_rle_a1), rownames(background_rare_cf_pangenome_rle_a1))
# convert to data frame
background_a1_6 <- data.frame(background_cf_rare_a1_rle)
# add meta data
background_a1_6$State <- 'CF'
background_a1_6$Species_type <- 'background_rare'
background_a1_6$Threshold <- '15th percentile'
background_a1_6$Normalisation <- 'RLE'
# rename columns
colnames(background_a1_6) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background core species, for 15th percentile, BCPHC-normalisation, overlaps between both databases
background_healthy_core_a1_bcphc <- intersect(rownames(background_core_h_osps_a1), rownames(background_core_h_pangenome_a1))
# convert to data frame
background_a1_7 <- data.frame(background_healthy_core_a1_bcphc)
# add meta data
background_a1_7$State <- 'Healthy'
background_a1_7$Species_type <- 'background_core'
background_a1_7$Threshold <- '15th percentile'
background_a1_7$Normalisation <- 'BCPHC'
# rename columns
colnames(background_a1_7) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background core species, for 15th percentile, VST-normalisation, overlaps between both databases
background_healthy_core_a1_vst <- intersect(rownames(background_core_h_osps_vst_a1), rownames(background_core_h_pangenome_vst_a1))
# convert to data frame
background_a1_8 <- data.frame(background_healthy_core_a1_vst)
# add meta data
background_a1_8$State <- 'Healthy'
background_a1_8$Species_type <- 'background_core'
background_a1_8$Threshold <- '15th percentile'
background_a1_8$Normalisation <- 'VST'
# rename columns
colnames(background_a1_8) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background core species, for 15th percentile, RLE-normalisation, overlaps between both databases
background_healthy_core_a1_rle <- intersect(rownames(background_core_h_osps_rle_a1), rownames(background_core_h_pangenome_rle_a1))
# convert to data frame
background_a1_9 <- data.frame(background_healthy_core_a1_rle)
# add meta data
background_a1_9$State <- 'Healthy'
background_a1_9$Species_type <- 'background_core'
background_a1_9$Threshold <- '15th percentile'
background_a1_9$Normalisation <- 'RLE'
# rename columns
colnames(background_a1_9) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background rare species, for 15th percentile, BCPHC-normalisation, overlaps between both databases
background_healthy_rare_a1_bcphc <- intersect(rownames(background_rare_h_osps_a1), rownames(background_rare_h_pangenome_a1))
# convert to data frame
background_a1_10 <- data.frame(background_healthy_rare_a1_bcphc)
# add meta data
background_a1_10$State <- 'Healthy'
background_a1_10$Species_type <- 'background_rare'
background_a1_10$Threshold <- '15th percentile'
background_a1_10$Normalisation <- 'BCPHC'
# rename columns
colnames(background_a1_10) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background rare species, for 15th percentile, VST-normalisation, overlaps between both databases
background_healthy_rare_a1_vst <- intersect(rownames(background_rare_h_osps_vst_a1), rownames(background_rare_h_pangenome_vst_a1))
# convert to data frame
background_a1_11 <- data.frame(background_healthy_rare_a1_vst)
# add meta data
background_a1_11$State <- 'Healthy'
background_a1_11$Species_type <- 'background_rare'
background_a1_11$Threshold <- '15th percentile'
background_a1_11$Normalisation <- 'VST'
# rename columns
colnames(background_a1_11) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background rare species, for 15th percentile, RLE-normalisation, overlaps between both databases
background_healthy_rare_a1_rle <- intersect(rownames(background_rare_h_osps_rle_a1), rownames(background_rare_h_pangenome_rle_a1))
# convert to data frame
background_a1_12 <- data.frame(background_healthy_rare_a1_rle)
# add meta data
background_a1_12$State <- 'Healthy'
background_a1_12$Species_type <- 'background_rare'
background_a1_12$Threshold <- '15th percentile'
background_a1_12$Normalisation <- 'RLE'
# rename columns
colnames(background_a1_12) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')


# make list of background species, for 35th percentile
# CF, make list of background core species, for 35th percentile, BCPHC-normalisation, overlaps between both databases
background_cf_core_a3_bcphc <- intersect(rownames(background_core_cf_osps_a3), rownames(background_core_cf_pangenome_a3))
# convert to data frame
background_a3_1 <- data.frame(background_cf_core_a3_bcphc)
# add meta data
background_a3_1$State <- 'CF'
background_a3_1$Species_type <- 'background_core'
background_a3_1$Threshold <- '15th percentile'
background_a3_1$Normalisation <- 'BCPHC'
# rename columns
colnames(background_a3_1) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# CF, make list of background core species, for 35th percentile, VST-normalisation, overlaps between both databases
background_cf_core_a3_vst <- intersect(rownames(background_core_cf_osps_vst_a3), rownames(background_core_cf_pangenome_vst_a3))
# convert to data frame
background_a3_2 <- data.frame(background_cf_core_a3_vst)
# add meta data
background_a3_2$State <- 'CF'
background_a3_2$Species_type <- 'background_core'
background_a3_2$Threshold <- '15th percentile'
background_a3_2$Normalisation <- 'VST'
# rename columns
colnames(background_a3_2) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# CF, make list of background core species, for 35th percentile, RLE-normalisation, overlaps between both databases
background_cf_core_a3_rle <- intersect(rownames(background_core_cf_osps_rle_a3), rownames(background_core_cf_pangenome_rle_a3))
# convert to data frame
background_a3_3 <- data.frame(background_cf_core_a3_rle)
# add meta data
background_a3_3$State <- 'CF'
background_a3_3$Species_type <- 'background_core'
background_a3_3$Threshold <- '15th percentile'
background_a3_3$Normalisation <- 'RLE'
# rename columns
colnames(background_a3_3) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# CF, make list of background rare species, for 35th percentile, BCPHC-normalisation, overlaps between both databases
background_cf_rare_a3_bcphc <- intersect(rownames(background_rare_cf_osps_a3), rownames(background_rare_cf_pangenome_a3))
# convert to data frame
background_a3_4 <- data.frame(background_cf_rare_a3_bcphc)
# add meta data
background_a3_4$State <- 'CF'
background_a3_4$Species_type <- 'background_rare'
background_a3_4$Threshold <- '15th percentile'
background_a3_4$Normalisation <- 'BCPHC'
# rename columns
colnames(background_a3_4) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# CF, make list of background rare species, for 35th percentile, VST-normalisation, overlaps between both databases
background_cf_rare_a3_vst <- c(rownames(background_rare_cf_osps_vst_a3), rownames(background_rare_cf_pangenome_vst_a3))
# convert to data frame
background_a3_5 <- data.frame(background_cf_rare_a3_vst)
# add meta data
background_a3_5$State <- 'CF'
background_a3_5$Species_type <- 'background_rare'
background_a3_5$Threshold <- '15th percentile'
background_a3_5$Normalisation <- 'VST'
# rename columns
colnames(background_a3_5) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# CF, make list of background rare species, for 35th percentile, RLE-normalisation, overlaps between both databases
background_cf_rare_a3_rle <- c(rownames(background_rare_cf_osps_rle_a3), rownames(background_rare_cf_pangenome_rle_a3))
# convert to data frame
background_a3_6 <- data.frame(background_cf_rare_a3_rle)
# add meta data
background_a3_6$State <- 'CF'
background_a3_6$Species_type <- 'background_rare'
background_a3_6$Threshold <- '15th percentile'
background_a3_6$Normalisation <- 'RLE'
# rename columns
colnames(background_a3_6) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background core species, for 35th percentile, BCPHC-normalisation, overlaps between both databases
background_healthy_core_a3_bcphc <- intersect(rownames(background_core_h_osps_a3), rownames(background_core_h_pangenome_a3))
# convert to data frame
background_a3_7 <- data.frame(background_healthy_core_a3_bcphc)
# add meta data
background_a3_7$State <- 'Healthy'
background_a3_7$Species_type <- 'background_core'
background_a3_7$Threshold <- '15th percentile'
background_a3_7$Normalisation <- 'BCPHC'
# rename columns
colnames(background_a3_7) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background core species, for 35th percentile, VST-normalisation, overlaps between both databases
background_healthy_core_a3_vst <- intersect(rownames(background_core_h_osps_vst_a3), rownames(background_core_h_pangenome_vst_a3))
# convert to data frame
background_a3_8 <- data.frame(background_healthy_core_a3_vst)
# add meta data
background_a3_8$State <- 'Healthy'
background_a3_8$Species_type <- 'background_core'
background_a3_8$Threshold <- '15th percentile'
background_a3_8$Normalisation <- 'VST'
# rename columns
colnames(background_a3_8) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background core species, for 35th percentile, RLE-normalisation, overlaps between both databases
background_healthy_core_a3_rle <- intersect(rownames(background_core_h_osps_rle_a3), rownames(background_core_h_pangenome_rle_a3))
# convert to data frame
background_a3_9 <- data.frame(background_healthy_core_a3_rle)
# add meta data
background_a3_9$State <- 'Healthy'
background_a3_9$Species_type <- 'background_core'
background_a3_9$Threshold <- '15th percentile'
background_a3_9$Normalisation <- 'RLE'
# rename columns
colnames(background_a3_9) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background rare species, for 35th percentile, BCPHC-normalisation, overlaps between both databases
background_healthy_rare_a3_bcphc <- intersect(rownames(background_rare_h_osps_a3), rownames(background_rare_h_pangenome_a3))
# convert to data frame
background_a3_10 <- data.frame(background_healthy_rare_a3_bcphc)
# add meta data
background_a3_10$State <- 'Healthy'
background_a3_10$Species_type <- 'background_rare'
background_a3_10$Threshold <- '15th percentile'
background_a3_10$Normalisation <- 'BCPHC'
# rename columns
colnames(background_a3_10) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background rare species, for 35th percentile, VST-normalisation, overlaps between both databases
background_healthy_rare_a3_vst <- intersect(rownames(background_rare_h_osps_vst_a3), rownames(background_rare_h_pangenome_vst_a3))
# convert to data frame
background_a3_11 <- data.frame(background_healthy_rare_a3_vst)
# add meta data
background_a3_11$State <- 'Healthy'
background_a3_11$Species_type <- 'background_rare'
background_a3_11$Threshold <- '15th percentile'
background_a3_11$Normalisation <- 'VST'
# rename columns
colnames(background_a3_11) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# Healthy, make list of background rare species, for 35th percentile, RLE-normalisation, overlaps between both databases
background_healthy_rare_a3_rle <- intersect(rownames(background_rare_h_osps_rle_a3), rownames(background_rare_h_pangenome_rle_a3))
# convert to data frame
background_a3_12 <- data.frame(background_healthy_rare_a3_rle)
# add meta data
background_a3_12$State <- 'Healthy'
background_a3_12$Species_type <- 'background_rare'
background_a3_12$Threshold <- '15th percentile'
background_a3_12$Normalisation <- 'RLE'
# rename columns
colnames(background_a3_12) <- c('Species', 'State', 'Species_type', 'Threshold', 'Normalisation')

# bind all tables into one
background_table <- data.frame(rbind(background_a1_1,background_a1_2,background_a1_4,background_a1_5, background_a1_6, background_a1_7,background_a1_8,background_a1_9,
                                     background_a1_10,background_a1_11,background_a1_12, background_a2_1, background_a2_2,background_a2_3,background_a2_4,background_a2_6,
                                     background_a2_7,background_a2_8,background_a2_9,background_a2_10, background_a2_11,background_a2_12,background_a3_1,background_a3_2,
                                     background_a3_3,background_a3_4,background_a3_7,background_a3_8, background_a3_9,background_a3_10,background_a3_11,background_a3_12))

############################################################################################################
# export table into working directory
write.table(background_table, file = 'input_files/taxonomic_data/background_species.csv', row.names = FALSE, col.names = TRUE, sep=';')
############################################################################################################

# Generate Venn Diagram plot
# rare species, pangenome, first year of life, 15th, 25th, 35th abundance percentile, CF
merge_venn_cf_y_rare_A <- c(cf_x_rare_pangenome_a2$A, cf_x_rare_pangenome_a1$A, cf_x_rare_pangenome_a3$A)
# extract overlapping species
merge_venn_cf_y_rare_A <- merge_venn_cf_y_rare_A[duplicated(merge_venn_cf_y_rare_A)]
# remove duplicate entries
merge_venn_cf_y_rare_A <- merge_venn_cf_y_rare_A[!duplicated(merge_venn_cf_y_rare_A)]

# rare species, pangenome, second year of life, 15th, 25th, 35th abundance percentile, CF
merge_venn_cf_y_rare_B <- c(cf_x_rare_pangenome_a2$B, cf_x_rare_pangenome_a1$B, cf_x_rare_pangenome_a3$B)
# extract overlapping species
merge_venn_cf_y_rare_B <- merge_venn_cf_y_rare_B[duplicated(merge_venn_cf_y_rare_B)]
# remove duplicates
merge_venn_cf_y_rare_B <- merge_venn_cf_y_rare_B[!duplicated(merge_venn_cf_y_rare_B)]

# rare species, pangenome, preschool, 15th, 25th, 35th abundance percentile, CF
merge_venn_cf_y_rare_C <- c(cf_x_rare_pangenome_a2$C, cf_x_rare_pangenome_a1$C, cf_x_rare_pangenome_a3$C)
# extract overlapping species
merge_venn_cf_y_rare_C <- merge_venn_cf_y_rare_C[duplicated(merge_venn_cf_y_rare_C)]
# remove duplicates
merge_venn_cf_y_rare_C <- merge_venn_cf_y_rare_C[!duplicated(merge_venn_cf_y_rare_C)]

# merge information, CF, all age groups, rare
merge_venn_cf_y_rare <- list(A = merge_venn_cf_y_rare_A, B = merge_venn_cf_y_rare_B, C = merge_venn_cf_y_rare_C)

# core species, pangenome, first year of life, 15th, 25th, 35th abundance percentile, CF
merge_venn_cf_y_core_A <- c(cf_x_core_pangenome_a2$A, cf_x_core_pangenome_a1$A, cf_x_core_pangenome_a3$A)
# extract overlapping species
merge_venn_cf_y_core_A <- merge_venn_cf_y_core_A[duplicated(merge_venn_cf_y_core_A)]
# remove duplicates
merge_venn_cf_y_core_A <- merge_venn_cf_y_core_A[!duplicated(merge_venn_cf_y_core_A)]

# core species, pangenome, second year of life, 15th, 25th, 35th abundance percentile, CF
merge_venn_cf_y_core_B <- c(cf_x_core_pangenome_a2$B, cf_x_core_pangenome_a1$B, cf_x_core_pangenome_a3$B)
# extract overlapping species
merge_venn_cf_y_core_B <- merge_venn_cf_y_core_B[duplicated(merge_venn_cf_y_core_B)]
# remove duplicates
merge_venn_cf_y_core_B <- merge_venn_cf_y_core_B[!duplicated(merge_venn_cf_y_core_B)]

# core species, pangenome,preschool, 15th, 25th, 35th abundance percentile, CF
merge_venn_cf_y_core_C <- c(cf_x_core_pangenome_a2$C, cf_x_core_pangenome_a1$C, cf_x_core_pangenome_a3$C)
# extract overlapping species
merge_venn_cf_y_core_C <- merge_venn_cf_y_core_C[duplicated(merge_venn_cf_y_core_C)]
# remove duplicates
merge_venn_cf_y_core_C <- merge_venn_cf_y_core_C[!duplicated(merge_venn_cf_y_core_C)]
# merge information, CF, all age groups, core
merge_venn_cf_y_core <- list(A = merge_venn_cf_y_core_A, B = merge_venn_cf_y_core_B, C = merge_venn_cf_y_core_C)


# rare species, pangenome, first year of life, 15th, 25th, 35th abundance percentile, healthy
merge_venn_h_y_rare_A <- c(h_x_rare_pangenome_a2$A, h_x_rare_pangenome_a1$A, h_x_rare_pangenome_a3$A)
# extract overlapping species
merge_venn_h_y_rare_A <- merge_venn_h_y_rare_A[duplicated(merge_venn_h_y_rare_A)]
# remove duplicates
merge_venn_h_y_rare_A <- merge_venn_h_y_rare_A[!duplicated(merge_venn_h_y_rare_A)]

# rare species, pangenome, toddler, 15th, 25th, 35th abundance percentile, healthy
merge_venn_h_y_rare_B <- c(h_x_rare_pangenome_a2$B, h_x_rare_pangenome_a1$B, h_x_rare_pangenome_a3$B)
# extract overlapping species
merge_venn_h_y_rare_B <- merge_venn_h_y_rare_B[duplicated(merge_venn_h_y_rare_B)]
# remove duplicates
merge_venn_h_y_rare_B <- merge_venn_h_y_rare_B[!duplicated(merge_venn_h_y_rare_B)]

# rare species, pangenome, preschool, 15th, 25th, 35th abundance percentile, healthy
merge_venn_h_y_rare_C <- c(h_x_rare_pangenome_a2$C, h_x_rare_pangenome_a1$C,h_x_rare_pangenome_a3$C)
# extract overlapping species
merge_venn_h_y_rare_C <- merge_venn_h_y_rare_C[duplicated(merge_venn_h_y_rare_C)]
# remove duplicates
merge_venn_h_y_rare_C <- merge_venn_h_y_rare_C[!duplicated(merge_venn_h_y_rare_C)]

# merge information, healthy, all age groups, rare
merge_venn_h_y_rare <- list(A = merge_venn_h_y_rare_A, B = merge_venn_h_y_rare_B, C = merge_venn_h_y_rare_C)

# core species, pangenome, first year of life, 15th, 25th, 35th abundance percentile, healthy
merge_venn_h_y_core_A <- c(h_x_core_pangenome_a2$A, h_x_core_pangenome_a1$A,h_x_core_pangenome_a3$A)
# extract overlapping species
merge_venn_h_y_core_A <- merge_venn_h_y_core_A[duplicated(merge_venn_h_y_core_A)]
# remove duplicates
merge_venn_h_y_core_A <- merge_venn_h_y_core_A[!duplicated(merge_venn_h_y_core_A)]

# core species, pangenome, toddler, 15th, 25th, 35th abundance percentile, healthy
merge_venn_h_y_core_B <- c(h_x_core_pangenome_a2$B, h_x_core_pangenome_a1$B,h_x_core_pangenome_a3$B)
# extract overlapping species
merge_venn_h_y_core_B <- merge_venn_h_y_core_B[duplicated(merge_venn_h_y_core_B)]
# remove duplicates
merge_venn_h_y_core_B <- merge_venn_h_y_core_B[!duplicated(merge_venn_h_y_core_B)]

# core species, pangenome,preschool, 15th, 25th, 35th abundance percentile, healthy
merge_venn_h_y_core_C <- c(h_x_core_pangenome_a2$C, h_x_core_pangenome_a1$C,h_x_core_pangenome_a3$C)
# extract overlapping species
merge_venn_h_y_core_C <- merge_venn_h_y_core_C[duplicated(merge_venn_h_y_core_C)]
# remove duplicates
merge_venn_h_y_core_C <- merge_venn_h_y_core_C[!duplicated(merge_venn_h_y_core_C)]
# merge information, healthy, all age groups, core
merge_venn_h_y_core <- list(A = merge_venn_h_y_core_A, B = merge_venn_h_y_core_B, C = merge_venn_h_y_core_C)


# Plot 1, CF, rare
cf_rare_venn_pangenome_a2 <- ggVennDiagram(merge_venn_cf_y_rare, size=1.5, category.names = c('0', '1-3', '4-6'), label = 'both', label_alpha=0, percent_digits = 0) +
  theme_bw(base_size = 6) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(),
        legend.text = element_text(size=10), legend.title = element_text(color="black", size=10), legend.position = "bottom") +
  scale_fill_gradientn(name="Species count", colours = c("white", "beige","cadetblue1", "pink", "pink2", "red"), limits=c(0,32)) +
  scale_colour_manual(values=c('seashell2', 'seashell2', 'seashell2')) +
  geom_label(aes(x=-1.3, y=9), label='CF, rare species', size=3.5)


# Plot 2, Healthy, rare
h_rare_venn_pangenome_a2 <- ggVennDiagram(merge_venn_h_y_rare, size=1.5, category.names = c('0', '1-3', '4-6'), label = 'both', label_alpha=0, percent_digits = 0) +
  theme_bw(base_size = 6) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(),
        legend.text = element_text(size=10), legend.title = element_text(color="black", size=10), legend.position = "bottom") +
  scale_fill_gradientn(name="Species count", colours = c("white", "beige","cadetblue1", "pink", "pink2", "red"), limits=c(0,32)) +
  scale_colour_manual(values=c('seashell2', 'seashell2', 'seashell2')) +
  geom_label(aes(x=-1.3, y=9), label='H, rare species', size=3.5)


# Plot 3, CF, core
cf_core_venn_pangenome_a2 <- ggVennDiagram(merge_venn_cf_y_core, size=1.5, category.names = c('0', '1-3', '4-6'), label = 'both', label_alpha=0, percent_digits = 0) +
  theme_bw(base_size = 6) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(),
        legend.text = element_text(size=10), legend.title = element_text(color="black", size=10), legend.position = "bottom") +
  scale_fill_gradientn(name="Species count", colours = c("white", "beige","cadetblue1", "pink", "pink2", "red"), limits=c(0,32)) +
  scale_colour_manual(values=c('seashell2', 'seashell2', 'seashell2')) +
  geom_label(aes(x=-1.3, y=9), label='CF, core species', size=3.5)


# Plot 4, Healthy, core
h_core_venn_pangenome_a2 <- ggVennDiagram(merge_venn_h_y_core, size=1.5, category.names = c('0', '1-3', '4-6'), label = 'both', label_alpha=0, percent_digits = 0) +
  theme_bw(base_size = 6) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), panel.border = element_blank(),
        legend.text = element_text(size=10), legend.title = element_text(color="black", size=10), legend.position = "bottom") +
  scale_fill_gradientn(name="Species count", colours = c("white", "beige","cadetblue1", "pink", "pink2", "red"), limits=c(0,32)) +
  scale_colour_manual(values=c('seashell2', 'seashell2', 'seashell2')) +
  geom_label(aes(x=-1.3, y=9), label='H, core species', size=3.5)


# merge plots into one
vennDiagramsPlot_pangenome <- ggarrange(h_rare_venn_pangenome_a2, h_core_venn_pangenome_a2, cf_rare_venn_pangenome_a2, cf_core_venn_pangenome_a2, 
                                        labels = c('A', 'B', 'C', 'D'), nrow=2, ncol=2, 
                                        font.label = list(size = 14, color = "black"), 
                                        common.legend = TRUE, legend = "right")

############################################################################################################
# make line plots
# start with generating input data frames (pangenome)
# number of core species, CF, year 0, 15th abundance percentile, BCPHC, pangenome
c_0_cf_pangenome_core_a1_n <- length(rownames(c_0_cf_pangenome_core_a1))
# number of core species, CF, year 1-3, 15th abundance percentile, BCPHC, pangenome
c_1_3_cf_pangenome_core_a1_n <- length(rownames(c_1_3_cf_pangenome_core_a1))
# number of core species, CF, year 4-6, 15th abundance percentile, BCPHC, pangenome
c_4_6_cf_pangenome_core_a1_n <- length(rownames(c_4_6_cf_pangenome_core_a1))
# merge information
c_cf_core_pangenome_a1_list <- c(c_0_cf_pangenome_core_a1_n, c_1_3_cf_pangenome_core_a1_n, c_4_6_cf_pangenome_core_a1_n)

# number of core species, CF, year 0, 25th abundance percentile, BCPHC, pangenome
c_0_cf_pangenome_core_a2_n <- length(rownames(c_0_cf_pangenome_core_a2))
# number of core species, CF, year 1-3, 25th abundance percentile, BCPHC, pangenome
c_1_3_cf_pangenome_core_a2_n <- length(rownames(c_1_3_cf_pangenome_core_a2))
# number of core species, CF, year 4-6, 25th abundance percentile, BCPHC, pangenome
c_4_6_cf_pangenome_core_a2_n <- length(rownames(c_4_6_cf_pangenome_core_a2))
# merge information
c_cf_core_pangenome_a2_list <- c(c_0_cf_pangenome_core_a2_n, c_1_3_cf_pangenome_core_a2_n, c_4_6_cf_pangenome_core_a2_n)

# number of core species, CF, year 0, 35th abundance percentile, BCPHC, pangenome
c_0_cf_pangenome_core_a3_n <- length(rownames(c_0_cf_pangenome_core_a3))
# number of core species, CF, year 1-3, 35th abundance percentile, BCPHC, pangenome
c_1_3_cf_pangenome_core_a3_n <- length(rownames(c_1_3_cf_pangenome_core_a3))
# number of core species, CF, year 4-6, 35th abundance percentile, BCPHC, pangenome
c_4_6_cf_pangenome_core_a3_n <- length(rownames(c_4_6_cf_pangenome_core_a3))
# merge information
c_cf_core_pangenome_a3_list <- c(c_0_cf_pangenome_core_a3_n, c_1_3_cf_pangenome_core_a3_n, c_4_6_cf_pangenome_core_a3_n)


# number of rare species, CF, year 0, 15th abundance percentile, BCPHC, pangenome
c_0_cf_pangenome_rare_a1_n <- length(rownames(c_0_cf_pangenome_rare_a1))
# number of rare species, CF, year 1-3, 15th abundance percentile, BCPHC, pangenome
c_1_3_cf_pangenome_rare_a1_n <- length(rownames(c_1_3_cf_pangenome_rare_a1))
# number of rare species, CF, year 4-6, 15th abundance percentile, BCPHC, pangenome
c_4_6_cf_pangenome_rare_a1_n <- length(rownames(c_4_6_cf_pangenome_rare_a1))
# merge information
c_cf_rare_pangenome_a1_list <- c(c_0_cf_pangenome_rare_a1_n, c_1_3_cf_pangenome_rare_a1_n, c_4_6_cf_pangenome_rare_a1_n)

# number of rare species, CF, year 0, 25th abundance percentile, BCPHC, pangenome
c_0_cf_pangenome_rare_a2_n <- length(rownames(c_0_cf_pangenome_rare_a2))
# number of rare species, CF, year 1-3, 25th abundance percentile, BCPHC, pangenome
c_1_3_cf_pangenome_rare_a2_n <- length(rownames(c_1_3_cf_pangenome_rare_a2))
# number of rare species, CF, year 4-6, 25th abundance percentile, BCPHC, pangenome
c_4_6_cf_pangenome_rare_a2_n <- length(rownames(c_4_6_cf_pangenome_rare_a2))
# merge information
c_cf_rare_pangenome_a2_list <- c(c_0_cf_pangenome_rare_a2_n, c_1_3_cf_pangenome_rare_a2_n, c_4_6_cf_pangenome_rare_a2_n)

# number of rare species, CF, year 0, 35th abundance percentile, BCPHC, pangenome
c_0_cf_pangenome_rare_a3_n <- length(rownames(c_0_cf_pangenome_rare_a3))
# number of rare species, CF, year 1-3, 35th abundance percentile, BCPHC, pangenome
c_1_3_cf_pangenome_rare_a3_n <- length(rownames(c_1_3_cf_pangenome_rare_a3))
# number of rare species, CF, year 4-6, 35th abundance percentile, BCPHC, pangenome
c_4_6_cf_pangenome_rare_a3_n <- length(rownames(c_4_6_cf_pangenome_rare_a3))
# merge information
c_cf_rare_pangenome_a3_list <- c(c_0_cf_pangenome_rare_a3_n, c_1_3_cf_pangenome_rare_a3_n, c_4_6_cf_pangenome_rare_a3_n)


# number of core species, CF, year 0, 15th abundance percentile, BCPHC, pangenome
c_0_h_pangenome_core_a1_n <- length(rownames(c_0_h_pangenome_core_a1))
# number of core species, CF, year 1-3, 15th abundance percentile, BCPHC, pangenome
c_1_3_h_pangenome_core_a1_n <- length(rownames(c_1_3_h_pangenome_core_a1))
# number of core species, CF, year 4-6, 15th abundance percentile, BCPHC, pangenome
c_4_6_h_pangenome_core_a1_n <- length(rownames(c_4_6_h_pangenome_core_a1))
# merge information
c_h_core_pangenome_a1_list <- c(c_0_h_pangenome_core_a1_n, c_1_3_h_pangenome_core_a1_n, c_4_6_h_pangenome_core_a1_n)

# number of core species, CF, year 0, 25th abundance percentile, BCPHC, pangenome
c_0_h_pangenome_core_a2_n <- length(rownames(c_0_h_pangenome_core_a2))
# number of core species, CF, year 1-3, 25th abundance percentile, BCPHC, pangenome
c_1_3_h_pangenome_core_a2_n <- length(rownames(c_1_3_h_pangenome_core_a2))
# number of core species, CF, year 4-6, 25th abundance percentile, BCPHC, pangenome
c_4_6_h_pangenome_core_a2_n <- length(rownames(c_4_6_h_pangenome_core_a2))
# merge information
c_h_core_pangenome_a2_list <- c(c_0_h_pangenome_core_a2_n, c_1_3_h_pangenome_core_a2_n, c_4_6_h_pangenome_core_a2_n)

# number of core species, CF, year 0, 35th abundance percentile, BCPHC, pangenome
c_0_h_pangenome_core_a3_n <- length(rownames(c_0_h_pangenome_core_a3))
# number of core species, CF, year 1-3, 35th abundance percentile, BCPHC, pangenome
c_1_3_h_pangenome_core_a3_n <- length(rownames(c_1_3_h_pangenome_core_a3))
# number of core species, CF, year 4-6, 35th abundance percentile, BCPHC, pangenome
c_4_6_h_pangenome_core_a3_n <- length(rownames(c_4_6_h_pangenome_core_a3))
# merge information
c_h_core_pangenome_a3_list <- c(c_0_h_pangenome_core_a3_n, c_1_3_h_pangenome_core_a3_n, c_4_6_h_pangenome_core_a3_n)

# number of rare species, healthy, year 0, 15th abundance percentile, BCPHC, pangenome
c_0_h_pangenome_rare_a1_n <- length(rownames(c_0_h_pangenome_rare_a1))
# number of rare species, healthy, year 1-3, 15th abundance percentile, BCPHC, pangenome
c_1_3_h_pangenome_rare_a1_n <- length(rownames(c_1_3_h_pangenome_rare_a1))
# number of rare species, healthy, year 4-6, 15th abundance percentile, BCPHC, pangenome
c_4_6_h_pangenome_rare_a1_n <- length(rownames(c_4_6_h_pangenome_rare_a1))
# merge information
c_h_rare_pangenome_a1_list <- c(c_0_h_pangenome_rare_a1_n, c_1_3_h_pangenome_rare_a1_n, c_4_6_h_pangenome_rare_a1_n)

# number of rare species, healthy, year 0, 25th abundance percentile, BCPHC, pangenome
c_0_h_pangenome_rare_a2_n <- length(rownames(c_0_h_pangenome_rare_a2))
# number of rare species, healthy, year 1-3, 25th abundance percentile, BCPHC, pangenome
c_1_3_h_pangenome_rare_a2_n <- length(rownames(c_1_3_h_pangenome_rare_a2))
# number of rare species, healthy, year 4-6, 25th abundance percentile, BCPHC, pangenome
c_4_6_h_pangenome_rare_a2_n <- length(rownames(c_4_6_h_pangenome_rare_a2))
# merge information
c_h_rare_pangenome_a2_list <- c(c_0_h_pangenome_rare_a2_n, c_1_3_h_pangenome_rare_a2_n, c_4_6_h_pangenome_rare_a2_n)

# number of rare species, healthy, year 0, 35th abundance percentile, BCPHC, pangenome
c_0_h_pangenome_rare_a3_n <- length(rownames(c_0_h_pangenome_rare_a3))
# number of rare species, healthy, year 1-3, 35th abundance percentile, BCPHC, pangenome
c_1_3_h_pangenome_rare_a3_n <- length(rownames(c_1_3_h_pangenome_rare_a3))
# number of rare species, healthy, year 4-6, 35th abundance percentile, BCPHC, pangenome
c_4_6_h_pangenome_rare_a3_n <- length(rownames(c_4_6_h_pangenome_rare_a3))
# merge information
c_h_rare_pangenome_a3_list <- c(c_0_h_pangenome_rare_a3_n, c_1_3_h_pangenome_rare_a3_n, c_4_6_h_pangenome_rare_a3_n)

# number of core species, healthy, year 0, 15th abundance percentile, RLE, pangenome
c_0_cf_pangenome_rle_core_a1_n <- length(rownames(c_0_cf_pangenome_rle_core_a1))
# number of core species, healthy, year 1-3, 15th abundance percentile, RLE, pangenome
c_1_3_cf_pangenome_rle_core_a1_n <- length(rownames(c_1_3_cf_pangenome_rle_core_a1))
# number of core species, healthy, year 1-3, 15th abundance percentile, RLE, pangenome
c_4_6_cf_pangenome_rle_core_a1_n <- length(rownames(c_4_6_cf_pangenome_rle_core_a1))
# merge information
c_cf_core_pangenome_rle_a1_list <- c(c_0_cf_pangenome_rle_core_a1_n, c_1_3_cf_pangenome_rle_core_a1_n, c_4_6_cf_pangenome_rle_core_a1_n)

# number of core species, healthy, year 0, 25th abundance percentile, RLE, pangenome
c_0_cf_pangenome_rle_core_a2_n <- length(rownames(c_0_cf_pangenome_rle_core_a2))
# number of core species, healthy, year 1-3, 25th abundance percentile, RLE, pangenome
c_1_3_cf_pangenome_rle_core_a2_n <- length(rownames(c_1_3_cf_pangenome_rle_core_a2))
# number of core species, healthy, year 4-6, 25th abundance percentile, RLE, pangenome
c_4_6_cf_pangenome_rle_core_a2_n <- length(rownames(c_4_6_cf_pangenome_rle_core_a2))
# merge information
c_cf_core_pangenome_rle_a2_list <- c(c_0_cf_pangenome_rle_core_a2_n, c_1_3_cf_pangenome_rle_core_a2_n, c_4_6_cf_pangenome_rle_core_a2_n)

# number of core species, healthy, year 0, 35th abundance percentile, RLE, pangenome
c_0_cf_pangenome_rle_core_a3_n <- length(rownames(c_0_cf_pangenome_rle_core_a3))
# number of core species, healthy, year 1-3, 35th abundance percentile, RLE, pangenome
c_1_3_cf_pangenome_rle_core_a3_n <- length(rownames(c_1_3_cf_pangenome_rle_core_a3))
# number of core species, healthy, year 4-6, 35th abundance percentile, RLE, pangenome
c_4_6_cf_pangenome_rle_core_a3_n <- length(rownames(c_4_6_cf_pangenome_rle_core_a3))
# merge information
c_cf_core_pangenome_rle_a3_list <- c(c_0_cf_pangenome_rle_core_a3_n, c_1_3_cf_pangenome_rle_core_a3_n, c_4_6_cf_pangenome_rle_core_a3_n)

# number of rare species, healthy, year 0, 15th abundance percentile, RLE, pangenome
c_0_cf_pangenome_rle_rare_a1_n <- length(rownames(c_0_cf_pangenome_rle_rare_a1))
# number of rare species, healthy, year 1-3, 15th abundance percentile, RLE, pangenome
c_1_3_cf_pangenome_rle_rare_a1_n <- length(rownames(c_1_3_cf_pangenome_rle_rare_a1))
# number of rare species, healthy, year 4-6, 15th abundance percentile, RLE, pangenome
c_4_6_cf_pangenome_rle_rare_a1_n <- length(rownames(c_4_6_cf_pangenome_rle_rare_a1))
# merge information
c_cf_rare_pangenome_rle_a1_list <- c(c_0_cf_pangenome_rle_rare_a1_n, c_1_3_cf_pangenome_rle_rare_a1_n, c_4_6_cf_pangenome_rle_rare_a1_n)

# number of rare species, healthy, year 0, 25th abundance percentile, RLE, pangenome
c_0_cf_pangenome_rle_rare_a2_n <- length(rownames(c_0_cf_pangenome_rle_rare_a2))
# number of rare species, healthy, year 1-3, 25th abundance percentile, RLE, pangenome
c_1_3_cf_pangenome_rle_rare_a2_n <- length(rownames(c_1_3_cf_pangenome_rle_rare_a2))
# number of rare species, healthy, year 4-6, 25th abundance percentile, RLE, pangenome
c_4_6_cf_pangenome_rle_rare_a2_n <- length(rownames(c_4_6_cf_pangenome_rle_rare_a2))
# merge information
c_cf_rare_pangenome_rle_a2_list <- c(c_0_cf_pangenome_rle_rare_a2_n, c_1_3_cf_pangenome_rle_rare_a2_n, c_4_6_cf_pangenome_rle_rare_a2_n)

# number of rare species, healthy, year 0, 35th abundance percentile, RLE, pangenome
c_0_cf_pangenome_rle_rare_a3_n <- length(rownames(c_0_cf_pangenome_rle_rare_a3))
# number of rare species, healthy, year 1-3, 35th abundance percentile, RLE, pangenome
c_1_3_cf_pangenome_rle_rare_a3_n <- length(rownames(c_1_3_cf_pangenome_rle_rare_a3))
# number of rare species, healthy, year 4-6, 35th abundance percentile, RLE, pangenome
c_4_6_cf_pangenome_rle_rare_a3_n <- length(rownames(c_4_6_cf_pangenome_rle_rare_a3))
# merge information
c_cf_rare_pangenome_rle_a3_list <- c(c_0_cf_pangenome_rle_rare_a3_n, c_1_3_cf_pangenome_rle_rare_a3_n, c_4_6_cf_pangenome_rle_rare_a3_n)

# number of core species, healthy, year 0, 15th abundance percentile, RLE, pangenome
c_0_h_pangenome_rle_core_a1_n <- length(rownames(c_0_h_pangenome_rle_core_a1))
# number of core species, healthy, year 1-3, 15th abundance percentile, RLE, pangenome
c_1_3_h_pangenome_rle_core_a1_n <- length(rownames(c_1_3_h_pangenome_rle_core_a1))
# number of core species, healthy, year 4-6, 15th abundance percentile, RLE, pangenome
c_4_6_h_pangenome_rle_core_a1_n <- length(rownames(c_4_6_h_pangenome_rle_core_a1))
# merge information
c_h_core_pangenome_rle_a1_list <- c(c_0_h_pangenome_rle_core_a1_n, c_1_3_h_pangenome_rle_core_a1_n, c_4_6_h_pangenome_rle_core_a1_n)

# number of core species, healthy, year 0, 25th abundance percentile, RLE, pangenome
c_0_h_pangenome_rle_core_a2_n <- length(rownames(c_0_h_pangenome_rle_core_a2))
# number of core species, healthy, year 1-3, 25th abundance percentile, RLE, pangenome
c_1_3_h_pangenome_rle_core_a2_n <- length(rownames(c_1_3_h_pangenome_rle_core_a2))
# number of core species, healthy, year 4-6, 25th abundance percentile, RLE, pangenome
c_4_6_h_pangenome_rle_core_a2_n <- length(rownames(c_4_6_h_pangenome_rle_core_a2))
# merge information
c_h_core_pangenome_rle_a2_list <- c(c_0_h_pangenome_rle_core_a2_n, c_1_3_h_pangenome_rle_core_a2_n, c_4_6_h_pangenome_rle_core_a2_n)

# number of core species, healthy, year 0, 35th abundance percentile, RLE, pangenome
c_0_h_pangenome_rle_core_a3_n <- length(rownames(c_0_h_pangenome_rle_core_a3))
# number of core species, healthy, year 1-3, 35th abundance percentile, RLE, pangenome
c_1_3_h_pangenome_rle_core_a3_n <- length(rownames(c_1_3_h_pangenome_rle_core_a3))
# number of core species, healthy, year 4-6, 35th abundance percentile, RLE, pangenome
c_4_6_h_pangenome_rle_core_a3_n <- length(rownames(c_4_6_h_pangenome_rle_core_a3))
# merge information
c_h_core_pangenome_rle_a3_list <- c(c_0_h_pangenome_rle_core_a3_n, c_1_3_h_pangenome_rle_core_a3_n, c_4_6_h_pangenome_rle_core_a3_n)

# number of rare species, healthy, year 0, 15th abundance percentile, RLE, pangenome
c_0_h_pangenome_rle_rare_a1_n <- length(rownames(c_0_h_pangenome_rle_rare_a1))
# number of rare species, healthy, year 1-3, 15th abundance percentile, RLE, pangenome
c_1_3_h_pangenome_rle_rare_a1_n <- length(rownames(c_1_3_h_pangenome_rle_rare_a1))
# number of rare species, healthy, year 4-6, 15th abundance percentile, RLE, pangenome
c_4_6_h_pangenome_rle_rare_a1_n <- length(rownames(c_4_6_h_pangenome_rle_rare_a1))
# merge information
c_h_rare_pangenome_rle_a1_list <- c(c_0_h_pangenome_rle_rare_a1_n, c_1_3_h_pangenome_rle_rare_a1_n, c_4_6_h_pangenome_rle_rare_a1_n)

# number of rare species, healthy, year 0, 25th abundance percentile, RLE, pangenome
c_0_h_pangenome_rle_rare_a2_n <- length(rownames(c_0_h_pangenome_rle_rare_a2))
# number of rare species, healthy, year 1-3, 25th abundance percentile, RLE, pangenome
c_1_3_h_pangenome_rle_rare_a2_n <- length(rownames(c_1_3_h_pangenome_rle_rare_a2))
# number of rare species, healthy, year 4-6, 25th abundance percentile, RLE, pangenome
c_4_6_h_pangenome_rle_rare_a2_n <- length(rownames(c_4_6_h_pangenome_rle_rare_a2))
# merge information
c_h_rare_pangenome_rle_a2_list <- c(c_0_h_pangenome_rle_rare_a2_n, c_1_3_h_pangenome_rle_rare_a2_n, c_4_6_h_pangenome_rle_rare_a2_n)

# number of rare species, healthy, year 0, 35th abundance percentile, RLE, pangenome
c_0_h_pangenome_rle_rare_a3_n <- length(rownames(c_0_h_pangenome_rle_rare_a3))
# number of rare species, healthy, year 1-3, 35th abundance percentile, RLE, pangenome
c_1_3_h_pangenome_rle_rare_a3_n <- length(rownames(c_1_3_h_pangenome_rle_rare_a3))
# number of rare species, healthy, year 4-6, 35th abundance percentile, RLE, pangenome
c_4_6_h_pangenome_rle_rare_a3_n <- length(rownames(c_4_6_h_pangenome_rle_rare_a3))
# merge information
c_h_rare_pangenome_rle_a3_list <- c(c_0_h_pangenome_rle_rare_a3_n, c_1_3_h_pangenome_rle_rare_a3_n, c_4_6_h_pangenome_rle_rare_a3_n)

# number of core species, CF, year 0, 15th abundance percentile, VST, pangenome
c_0_cf_pangenome_vst_core_a1_n <- length(rownames(c_0_cf_pangenome_vst_core_a1))
# number of core species, CF, year 1-3, 15th abundance percentile, VST, pangenome
c_1_3_cf_pangenome_vst_core_a1_n <- length(rownames(c_1_3_cf_pangenome_vst_core_a1))
# number of core species, CF, year 4-6, 15th abundance percentile, VST, pangenome
c_4_6_cf_pangenome_vst_core_a1_n <- length(rownames(c_4_6_cf_pangenome_vst_core_a1))
# merge information
c_cf_core_pangenome_vst_a1_list <- c(c_0_cf_pangenome_vst_core_a1_n, c_1_3_cf_pangenome_vst_core_a1_n, c_4_6_cf_pangenome_vst_core_a1_n)

# number of core species, CF, year 0, 25th abundance percentile, VST, pangenome
c_0_cf_pangenome_vst_core_a2_n <- length(rownames(c_0_cf_pangenome_vst_core_a2))
# number of core species, CF, year 1-3, 25th abundance percentile, VST, pangenome
c_1_3_cf_pangenome_vst_core_a2_n <- length(rownames(c_1_3_cf_pangenome_vst_core_a2))
# number of core species, CF, year 4-6, 25th abundance percentile, VST, pangenome
c_4_6_cf_pangenome_vst_core_a2_n <- length(rownames(c_4_6_cf_pangenome_vst_core_a2))
# merge information
c_cf_core_pangenome_vst_a2_list <- c(c_0_cf_pangenome_vst_core_a2_n, c_1_3_cf_pangenome_vst_core_a2_n, c_4_6_cf_pangenome_vst_core_a2_n)

# number of core species, CF, year 0, 35th abundance percentile, VST, pangenome
c_0_cf_pangenome_vst_core_a3_n <- length(rownames(c_0_cf_pangenome_vst_core_a3))
# number of core species, CF, year 1-3, 35th abundance percentile, VST, pangenome
c_1_3_cf_pangenome_vst_core_a3_n <- length(rownames(c_1_3_cf_pangenome_vst_core_a3))
# number of core species, CF, year 4-6, 35th abundance percentile, VST, pangenome
c_4_6_cf_pangenome_vst_core_a3_n <- length(rownames(c_4_6_cf_pangenome_vst_core_a3))
# merge information
c_cf_core_pangenome_vst_a3_list <- c(c_0_cf_pangenome_vst_core_a3_n, c_1_3_cf_pangenome_vst_core_a3_n, c_4_6_cf_pangenome_vst_core_a3_n)

# number of rare species, CF, year 0, 15th abundance percentile, VST, pangenome
c_0_cf_pangenome_vst_rare_a1_n <- length(rownames(c_0_cf_pangenome_vst_rare_a1))
# number of rare species, CF, year 1-3, 15th abundance percentile, VST, pangenome
c_1_3_cf_pangenome_vst_rare_a1_n <- length(rownames(c_1_3_cf_pangenome_vst_rare_a1))
# number of rare species, CF, year 4-6, 15th abundance percentile, VST, pangenome
c_4_6_cf_pangenome_vst_rare_a1_n <- length(rownames(c_4_6_cf_pangenome_vst_rare_a1))
# merge information
c_cf_rare_pangenome_vst_a1_list <- c(c_0_cf_pangenome_vst_rare_a1_n, c_1_3_cf_pangenome_vst_rare_a1_n, c_4_6_cf_pangenome_vst_rare_a1_n)

# number of rare species, CF, year 0, 25th abundance percentile, VST, pangenome
c_0_cf_pangenome_vst_rare_a2_n <- length(rownames(c_0_cf_pangenome_vst_rare_a2))
# number of rare species, CF, year 1-3, 25th abundance percentile, VST, pangenome
c_1_3_cf_pangenome_vst_rare_a2_n <- length(rownames(c_1_3_cf_pangenome_vst_rare_a2))
# number of rare species, CF, year 4-6, 25th abundance percentile, VST, pangenome
c_4_6_cf_pangenome_vst_rare_a2_n <- length(rownames(c_4_6_cf_pangenome_vst_rare_a2))
# merge information
c_cf_rare_pangenome_vst_a2_list <- c(c_0_cf_pangenome_vst_rare_a2_n, c_1_3_cf_pangenome_vst_rare_a2_n, c_4_6_cf_pangenome_vst_rare_a2_n)

# number of rare species, CF, year 0, 35th abundance percentile, VST, pangenome
c_0_cf_pangenome_vst_rare_a3_n <- length(rownames(c_0_cf_pangenome_vst_rare_a3))
# number of rare species, CF, year 1-3, 35th abundance percentile, VST, pangenome
c_1_3_cf_pangenome_vst_rare_a3_n <- length(rownames(c_1_3_cf_pangenome_vst_rare_a3))
# number of rare species, CF, year 4-6, 35th abundance percentile, VST, pangenome
c_4_6_cf_pangenome_vst_rare_a3_n <- length(rownames(c_4_6_cf_pangenome_vst_rare_a3))
# merge information
c_cf_rare_pangenome_vst_a3_list <- c(c_0_cf_pangenome_vst_rare_a3_n, c_1_3_cf_pangenome_vst_rare_a3_n, c_4_6_cf_pangenome_vst_rare_a3_n)

# number of core species, Healthy, year 0, 15th abundance percentile, VST, pangenome
c_0_h_pangenome_vst_core_a1_n <- length(rownames(c_0_h_pangenome_vst_core_a1))
# number of core species, Healthy, year 1-3, 15th abundance percentile, VST, pangenome
c_1_3_h_pangenome_vst_core_a1_n <- length(rownames(c_1_3_h_pangenome_vst_core_a1))
# number of core species, Healthy, year 4-6, 15th abundance percentile, VST, pangenome
c_4_6_h_pangenome_vst_core_a1_n <- length(rownames(c_4_6_h_pangenome_vst_core_a1))
# merge information
c_h_core_pangenome_vst_a1_list <- c(c_0_h_pangenome_vst_core_a1_n, c_1_3_h_pangenome_vst_core_a1_n, c_4_6_h_pangenome_vst_core_a1_n)

# number of core species, Healthy, year 0, 25th abundance percentile, VST, pangenome
c_0_h_pangenome_vst_core_a2_n <- length(rownames(c_0_h_pangenome_vst_core_a2))
# number of core species, Healthy, year 1-3, 25th abundance percentile, VST, pangenome
c_1_3_h_pangenome_vst_core_a2_n <- length(rownames(c_1_3_h_pangenome_vst_core_a2))
# number of core species, Healthy, year 4-6, 25th abundance percentile, VST, pangenome
c_4_6_h_pangenome_vst_core_a2_n <- length(rownames(c_4_6_h_pangenome_vst_core_a2))
# merge information
c_h_core_pangenome_vst_a2_list <- c(c_0_h_pangenome_vst_core_a2_n, c_1_3_h_pangenome_vst_core_a2_n, c_4_6_h_pangenome_vst_core_a2_n)

# number of core species, Healthy, year 0, 35th abundance percentile, VST, pangenome
c_0_h_pangenome_vst_core_a3_n <- length(rownames(c_0_h_pangenome_vst_core_a3))
# number of core species, Healthy, year 1-3, 35th abundance percentile, VST, pangenome
c_1_3_h_pangenome_vst_core_a3_n <- length(rownames(c_1_3_h_pangenome_vst_core_a3))
# number of core species, Healthy, year 4-6, 35th abundance percentile, VST, pangenome
c_4_6_h_pangenome_vst_core_a3_n <- length(rownames(c_4_6_h_pangenome_vst_core_a3))
# merge information
c_h_core_pangenome_vst_a3_list <- c(c_0_h_pangenome_vst_core_a3_n, c_1_3_h_pangenome_vst_core_a3_n, c_4_6_h_pangenome_vst_core_a3_n)

# number of rare species, Healthy, year 0, 15th abundance percentile, VST, pangenome
c_0_h_pangenome_vst_rare_a1_n <- length(rownames(c_0_h_pangenome_vst_rare_a1))
# number of rare species, Healthy, year 0, 15th abundance percentile, VST, pangenome
c_1_3_h_pangenome_vst_rare_a1_n <- length(rownames(c_1_3_h_pangenome_vst_rare_a1))
# number of rare species, Healthy, year 0, 15th abundance percentile, VST, pangenome
c_4_6_h_pangenome_vst_rare_a1_n <- length(rownames(c_4_6_h_pangenome_vst_rare_a1))
# merge information
c_h_rare_pangenome_vst_a1_list <- c(c_0_h_pangenome_vst_rare_a1_n, c_1_3_h_pangenome_vst_rare_a1_n, c_4_6_h_pangenome_vst_rare_a1_n)

# number of rare species, Healthy, year 0, 25th abundance percentile, VST, pangenome
c_0_h_pangenome_vst_rare_a2_n <- length(rownames(c_0_h_pangenome_vst_rare_a2))
# number of rare species, Healthy, year 1-3, 25th abundance percentile, VST, pangenome
c_1_3_h_pangenome_vst_rare_a2_n <- length(rownames(c_1_3_h_pangenome_vst_rare_a2))
# number of rare species, Healthy, year 4-6, 25th abundance percentile, VST, pangenome
c_4_6_h_pangenome_vst_rare_a2_n <- length(rownames(c_4_6_h_pangenome_vst_rare_a2))
# merge information
c_h_rare_pangenome_vst_a2_list <- c(c_0_h_pangenome_vst_rare_a2_n, c_1_3_h_pangenome_vst_rare_a2_n, c_4_6_h_pangenome_vst_rare_a2_n)

# number of rare species, Healthy, year 0, 35th abundance percentile, VST, pangenome
c_0_h_pangenome_vst_rare_a3_n <- length(rownames(c_0_h_pangenome_vst_rare_a3))
# number of rare species, Healthy, year 1-3, 35th abundance percentile, VST, pangenome
c_1_3_h_pangenome_vst_rare_a3_n <- length(rownames(c_1_3_h_pangenome_vst_rare_a3))
# number of rare species, Healthy, year 4-6, 35th abundance percentile, VST, pangenome
c_4_6_h_pangenome_vst_rare_a3_n <- length(rownames(c_4_6_h_pangenome_vst_rare_a3))
# merge information
c_h_rare_pangenome_vst_a3_list <- c(c_0_h_pangenome_vst_rare_a3_n, c_1_3_h_pangenome_vst_rare_a3_n, c_4_6_h_pangenome_vst_rare_a3_n)


# merge all information (pangenome)
# 15th abundance percentile, BCPHC 
age_pangenome_a1 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_pangenome_a1 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_pangenome_a1 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_pangenome_a1 <- c(c_cf_core_pangenome_a1_list, c_cf_rare_pangenome_a1_list, c_h_core_pangenome_a1_list, c_h_rare_pangenome_a1_list)

merge_data_venn_pangenome_a1 <- data.frame(cbind(age_pangenome_a1, state_pangenome_a1, species_type_pangenome_a1, c_pangenome_a1))
# add metadata
colnames(merge_data_venn_pangenome_a1) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_pangenome_a1$Threshold <- '15th percentile'
merge_data_venn_pangenome_a1$Normalisation <- 'BCPHC'
# convert age group to class factor and specify order
merge_data_venn_pangenome_a1$Age <- factor(merge_data_venn_pangenome_a1$Age, levels = c('0', '1-3', '4-6'))


# 25th abundance percentile, BCPHC 
age_pangenome_a2 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_pangenome_a2 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_pangenome_a2 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_pangenome_a2 <- c(c_cf_core_pangenome_a2_list, c_cf_rare_pangenome_a2_list, c_h_core_pangenome_a2_list, c_h_rare_pangenome_a2_list)

merge_data_venn_pangenome_a2 <- data.frame(cbind(age_pangenome_a2, state_pangenome_a2, species_type_pangenome_a2, c_pangenome_a2))
# add metadata
colnames(merge_data_venn_pangenome_a2) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_pangenome_a2$Threshold <- '25th percentile'
merge_data_venn_pangenome_a2$Normalisation <- 'BCPHC'
# convert age group to class factor and specify order
merge_data_venn_pangenome_a2$Age <- factor(merge_data_venn_pangenome_a2$Age, levels = c('0', '1-3', '4-6'))

# 35th abundance percentile, BCPHC 
age_pangenome_a3 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_pangenome_a3 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_pangenome_a3 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_pangenome_a3 <- c(c_cf_core_pangenome_a3_list, c_cf_rare_pangenome_a3_list, c_h_core_pangenome_a3_list, c_h_rare_pangenome_a3_list)
# add metadata 
merge_data_venn_pangenome_a3 <- data.frame(cbind(age_pangenome_a3, state_pangenome_a3, species_type_pangenome_a3, c_pangenome_a3))
colnames(merge_data_venn_pangenome_a3) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_pangenome_a3$Threshold <- '35th percentile'
merge_data_venn_pangenome_a3$Normalisation <- 'BCPHC'
# convert age group to class factor and specify order
merge_data_venn_pangenome_a3$Age <- factor(merge_data_venn_pangenome_a3$Age, levels = c('0', '1-3', '4-6'))


# 15th abundance percentile, RLE 
age_pangenome_rle_a1 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_pangenome_rle_a1 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_pangenome_rle_a1 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_pangenome_rle_a1 <- c(c_cf_core_pangenome_rle_a1_list, c_cf_rare_pangenome_rle_a1_list, c_h_core_pangenome_rle_a1_list, c_h_rare_pangenome_rle_a1_list)
# add metadata 
merge_data_venn_pangenome_rle_a1 <- data.frame(cbind(age_pangenome_rle_a1, state_pangenome_rle_a1, species_type_pangenome_rle_a1, c_pangenome_rle_a1))
colnames(merge_data_venn_pangenome_rle_a1) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_pangenome_rle_a1$Threshold <- '15th percentile'
merge_data_venn_pangenome_rle_a1$Normalisation <- 'RLE'
# convert age group to class factor and specify order
merge_data_venn_pangenome_rle_a1$Age <- factor(merge_data_venn_pangenome_rle_a1$Age, levels = c('0', '1-3', '4-6'))

# 25th abundance percentile, RLE 
age_pangenome_rle_a2 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_pangenome_rle_a2 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_pangenome_rle_a2 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_pangenome_rle_a2 <- c(c_cf_core_pangenome_rle_a2_list, c_cf_rare_pangenome_rle_a2_list, c_h_core_pangenome_rle_a2_list, c_h_rare_pangenome_rle_a2_list)
# add metadata 
merge_data_venn_pangenome_rle_a2 <- data.frame(cbind(age_pangenome_rle_a2, state_pangenome_rle_a2, species_type_pangenome_rle_a2, c_pangenome_rle_a2))
colnames(merge_data_venn_pangenome_rle_a2) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_pangenome_rle_a2$Threshold <- '25th percentile'
merge_data_venn_pangenome_rle_a2$Normalisation <- 'RLE'
# convert age group to class factor and specify order
merge_data_venn_pangenome_rle_a2$Age <- factor(merge_data_venn_pangenome_rle_a2$Age, levels = c('0', '1-3', '4-6'))

# 35th abundance percentile, RLE 
age_pangenome_rle_a3 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_pangenome_rle_a3 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_pangenome_rle_a3 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_pangenome_rle_a3 <- c(c_cf_core_pangenome_rle_a3_list, c_cf_rare_pangenome_rle_a3_list, c_h_core_pangenome_rle_a3_list, c_h_rare_pangenome_rle_a3_list)
# add metadata 
merge_data_venn_pangenome_rle_a3 <- data.frame(cbind(age_pangenome_rle_a3, state_pangenome_rle_a3, species_type_pangenome_rle_a3, c_pangenome_rle_a3))
colnames(merge_data_venn_pangenome_rle_a3) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_pangenome_rle_a3$Threshold <- '35th percentile'
merge_data_venn_pangenome_rle_a3$Normalisation <- 'RLE'
# convert age group to class factor and specify order
merge_data_venn_pangenome_rle_a3$Age <- factor(merge_data_venn_pangenome_rle_a3$Age, levels = c('0', '1-3', '4-6'))


# 15th abundance percentile, VST 
age_pangenome_vst_a1 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_pangenome_vst_a1 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_pangenome_vst_a1 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_pangenome_vst_a1 <- c(c_cf_core_pangenome_vst_a1_list, c_cf_rare_pangenome_vst_a1_list, c_h_core_pangenome_vst_a1_list, c_h_rare_pangenome_vst_a1_list)
# add metadata 
merge_data_venn_pangenome_vst_a1 <- data.frame(cbind(age_pangenome_vst_a1, state_pangenome_vst_a1, species_type_pangenome_vst_a1, c_pangenome_vst_a1))
colnames(merge_data_venn_pangenome_vst_a1) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_pangenome_vst_a1$Threshold <- '15th percentile'
merge_data_venn_pangenome_vst_a1$Normalisation <- 'VST'
# convert age group to class factor and specify order
merge_data_venn_pangenome_vst_a1$Age <- factor(merge_data_venn_pangenome_vst_a1$Age, levels = c('0', '1-3', '4-6'))

# 25th abundance percentile, VST 
age_pangenome_vst_a2 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_pangenome_vst_a2 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_pangenome_vst_a2 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_pangenome_vst_a2 <- c(c_cf_core_pangenome_vst_a2_list, c_cf_rare_pangenome_vst_a2_list, c_h_core_pangenome_vst_a2_list, c_h_rare_pangenome_vst_a2_list)
# add metadata 
merge_data_venn_pangenome_vst_a2 <- data.frame(cbind(age_pangenome_vst_a2, state_pangenome_vst_a2, species_type_pangenome_vst_a2, c_pangenome_vst_a2))
colnames(merge_data_venn_pangenome_vst_a2) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_pangenome_vst_a2$Threshold <- '25th percentile'
merge_data_venn_pangenome_vst_a2$Normalisation <- 'VST'
# convert age group to class factor and specify order
merge_data_venn_pangenome_vst_a2$Age <- factor(merge_data_venn_pangenome_vst_a2$Age, levels = c('0', '1-3', '4-6'))

# 35th abundance percentile, VST 
age_pangenome_vst_a3 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_pangenome_vst_a3 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_pangenome_vst_a3 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_pangenome_vst_a3 <- c(c_cf_core_pangenome_vst_a3_list, c_cf_rare_pangenome_vst_a3_list, c_h_core_pangenome_vst_a3_list, c_h_rare_pangenome_vst_a3_list)
# add metadata 
merge_data_venn_pangenome_vst_a3 <- data.frame(cbind(age_pangenome_vst_a3, state_pangenome_vst_a3, species_type_pangenome_vst_a3, c_pangenome_vst_a3))
colnames(merge_data_venn_pangenome_vst_a3) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_pangenome_vst_a3$Threshold <- '35th percentile'
merge_data_venn_pangenome_vst_a3$Normalisation <- 'VST'
# convert age group to class factor and specify order
merge_data_venn_pangenome_vst_a3$Age <- factor(merge_data_venn_pangenome_vst_a3$Age, levels = c('0', '1-3', '4-6'))

# merge all information in one data frame
merge_data_venn_pangenome <- data.frame(rbind(merge_data_venn_pangenome_a1, merge_data_venn_pangenome_a2, merge_data_venn_pangenome_a3,
                                              merge_data_venn_pangenome_rle_a1, merge_data_venn_pangenome_rle_a2, merge_data_venn_pangenome_rle_a3,
                                              merge_data_venn_pangenome_vst_a1, merge_data_venn_pangenome_vst_a2, merge_data_venn_pangenome_vst_a3))
# add column with pangenome information
merge_data_venn_pangenome$Database <- 'Pangenome'


# continuing with generating input data frames (osps)
# number of core species, CF, year 0, 15th abundance percentile, BCPHC, osps
c_0_cf_osps_core_a1_n <- length(rownames(c_0_cf_osps_core_a1))
# number of core species, CF, year 1-3, 15th abundance percentile, BCPHC, osps
c_1_3_cf_osps_core_a1_n <- length(rownames(c_1_3_cf_osps_core_a1))
# number of core species, CF, year 4-6, 15th abundance percentile, BCPHC, osps
c_4_6_cf_osps_core_a1_n <- length(rownames(c_4_6_cf_osps_core_a1))
# merge information
c_cf_core_osps_a1_list <- c(c_0_cf_osps_core_a1_n, c_1_3_cf_osps_core_a1_n, c_4_6_cf_osps_core_a1_n)

# number of core species, CF, year 0, 25th abundance percentile, BCPHC, osps
c_0_cf_osps_core_a2_n <- length(rownames(c_0_cf_osps_core_a2))
# number of core species, CF, year 1-3, 25th abundance percentile, BCPHC, osps
c_1_3_cf_osps_core_a2_n <- length(rownames(c_1_3_cf_osps_core_a2))
# number of core species, CF, year 4-6, 25th abundance percentile, BCPHC, osps
c_4_6_cf_osps_core_a2_n <- length(rownames(c_4_6_cf_osps_core_a2))
# merge information
c_cf_core_osps_a2_list <- c(c_0_cf_osps_core_a2_n, c_1_3_cf_osps_core_a2_n, c_4_6_cf_osps_core_a2_n)

# number of core species, CF, year 0, 35th abundance percentile, BCPHC, osps
c_0_cf_osps_core_a3_n <- length(rownames(c_0_cf_osps_core_a3))
# number of core species, CF, year 1-3, 35th abundance percentile, BCPHC, osps
c_1_3_cf_osps_core_a3_n <- length(rownames(c_1_3_cf_osps_core_a3))
# number of core species, CF, year 4-6, 35th abundance percentile, BCPHC, osps
c_4_6_cf_osps_core_a3_n <- length(rownames(c_4_6_cf_osps_core_a3))
# merge information
c_cf_core_osps_a3_list <- c(c_0_cf_osps_core_a3_n, c_1_3_cf_osps_core_a3_n, c_4_6_cf_osps_core_a3_n)


# number of rare species, CF, year 0, 15th abundance percentile, BCPHC, osps
c_0_cf_osps_rare_a1_n <- length(rownames(c_0_cf_osps_rare_a1))
# number of rare species, CF, year 1-3, 15th abundance percentile, BCPHC, osps
c_1_3_cf_osps_rare_a1_n <- length(rownames(c_1_3_cf_osps_rare_a1))
# number of rare species, CF, year 4-6, 15th abundance percentile, BCPHC, osps
c_4_6_cf_osps_rare_a1_n <- length(rownames(c_4_6_cf_osps_rare_a1))
# merge information
c_cf_rare_osps_a1_list <- c(c_0_cf_osps_rare_a1_n, c_1_3_cf_osps_rare_a1_n, c_4_6_cf_osps_rare_a1_n)

# number of rare species, CF, year 0, 25th abundance percentile, BCPHC, osps
c_0_cf_osps_rare_a2_n <- length(rownames(c_0_cf_osps_rare_a2))
# number of rare species, CF, year 1-3, 25th abundance percentile, BCPHC, osps
c_1_3_cf_osps_rare_a2_n <- length(rownames(c_1_3_cf_osps_rare_a2))
# number of rare species, CF, year 4-6, 25th abundance percentile, BCPHC, osps
c_4_6_cf_osps_rare_a2_n <- length(rownames(c_4_6_cf_osps_rare_a2))
# merge information
c_cf_rare_osps_a2_list <- c(c_0_cf_osps_rare_a2_n, c_1_3_cf_osps_rare_a2_n, c_4_6_cf_osps_rare_a2_n)

# number of rare species, CF, year 0, 35th abundance percentile, BCPHC, osps
c_0_cf_osps_rare_a3_n <- length(rownames(c_0_cf_osps_rare_a3))
# number of rare species, CF, year 1-3, 35th abundance percentile, BCPHC, osps
c_1_3_cf_osps_rare_a3_n <- length(rownames(c_1_3_cf_osps_rare_a3))
# number of rare species, CF, year 4-6, 35th abundance percentile, BCPHC, osps
c_4_6_cf_osps_rare_a3_n <- length(rownames(c_4_6_cf_osps_rare_a3))
# merge information
c_cf_rare_osps_a3_list <- c(c_0_cf_osps_rare_a3_n, c_1_3_cf_osps_rare_a3_n, c_4_6_cf_osps_rare_a3_n)


# number of core species, CF, year 0, 15th abundance percentile, BCPHC, osps
c_0_h_osps_core_a1_n <- length(rownames(c_0_h_osps_core_a1))
# number of core species, CF, year 1-3, 15th abundance percentile, BCPHC, osps
c_1_3_h_osps_core_a1_n <- length(rownames(c_1_3_h_osps_core_a1))
# number of core species, CF, year 4-6, 15th abundance percentile, BCPHC, osps
c_4_6_h_osps_core_a1_n <- length(rownames(c_4_6_h_osps_core_a1))
# merge information
c_h_core_osps_a1_list <- c(c_0_h_osps_core_a1_n, c_1_3_h_osps_core_a1_n, c_4_6_h_osps_core_a1_n)

# number of core species, CF, year 0, 25th abundance percentile, BCPHC, osps
c_0_h_osps_core_a2_n <- length(rownames(c_0_h_osps_core_a2))
# number of core species, CF, year 1-3, 25th abundance percentile, BCPHC, osps
c_1_3_h_osps_core_a2_n <- length(rownames(c_1_3_h_osps_core_a2))
# number of core species, CF, year 4-6, 25th abundance percentile, BCPHC, osps
c_4_6_h_osps_core_a2_n <- length(rownames(c_4_6_h_osps_core_a2))
# merge information
c_h_core_osps_a2_list <- c(c_0_h_osps_core_a2_n, c_1_3_h_osps_core_a2_n, c_4_6_h_osps_core_a2_n)

# number of core species, CF, year 0, 35th abundance percentile, BCPHC, osps
c_0_h_osps_core_a3_n <- length(rownames(c_0_h_osps_core_a3))
# number of core species, CF, year 1-3, 35th abundance percentile, BCPHC, osps
c_1_3_h_osps_core_a3_n <- length(rownames(c_1_3_h_osps_core_a3))
# number of core species, CF, year 4-6, 35th abundance percentile, BCPHC, osps
c_4_6_h_osps_core_a3_n <- length(rownames(c_4_6_h_osps_core_a3))
# merge information
c_h_core_osps_a3_list <- c(c_0_h_osps_core_a3_n, c_1_3_h_osps_core_a3_n, c_4_6_h_osps_core_a3_n)

# number of rare species, healthy, year 0, 15th abundance percentile, BCPHC, osps
c_0_h_osps_rare_a1_n <- length(rownames(c_0_h_osps_rare_a1))
# number of rare species, healthy, year 1-3, 15th abundance percentile, BCPHC, osps
c_1_3_h_osps_rare_a1_n <- length(rownames(c_1_3_h_osps_rare_a1))
# number of rare species, healthy, year 4-6, 15th abundance percentile, BCPHC, osps
c_4_6_h_osps_rare_a1_n <- length(rownames(c_4_6_h_osps_rare_a1))
# merge information
c_h_rare_osps_a1_list <- c(c_0_h_osps_rare_a1_n, c_1_3_h_osps_rare_a1_n, c_4_6_h_osps_rare_a1_n)

# number of rare species, healthy, year 0, 25th abundance percentile, BCPHC, osps
c_0_h_osps_rare_a2_n <- length(rownames(c_0_h_osps_rare_a2))
# number of rare species, healthy, year 1-3, 25th abundance percentile, BCPHC, osps
c_1_3_h_osps_rare_a2_n <- length(rownames(c_1_3_h_osps_rare_a2))
# number of rare species, healthy, year 4-6, 25th abundance percentile, BCPHC, osps
c_4_6_h_osps_rare_a2_n <- length(rownames(c_4_6_h_osps_rare_a2))
# merge information
c_h_rare_osps_a2_list <- c(c_0_h_osps_rare_a2_n, c_1_3_h_osps_rare_a2_n, c_4_6_h_osps_rare_a2_n)

# number of rare species, healthy, year 0, 35th abundance percentile, BCPHC, osps
c_0_h_osps_rare_a3_n <- length(rownames(c_0_h_osps_rare_a3))
# number of rare species, healthy, year 1-3, 35th abundance percentile, BCPHC, osps
c_1_3_h_osps_rare_a3_n <- length(rownames(c_1_3_h_osps_rare_a3))
# number of rare species, healthy, year 4-6, 35th abundance percentile, BCPHC, osps
c_4_6_h_osps_rare_a3_n <- length(rownames(c_4_6_h_osps_rare_a3))
# merge information
c_h_rare_osps_a3_list <- c(c_0_h_osps_rare_a3_n, c_1_3_h_osps_rare_a3_n, c_4_6_h_osps_rare_a3_n)

# number of core species, healthy, year 0, 15th abundance percentile, RLE, osps
c_0_cf_osps_rle_core_a1_n <- length(rownames(c_0_cf_osps_rle_core_a1))
# number of core species, healthy, year 1-3, 15th abundance percentile, RLE, osps
c_1_3_cf_osps_rle_core_a1_n <- length(rownames(c_1_3_cf_osps_rle_core_a1))
# number of core species, healthy, year 1-3, 15th abundance percentile, RLE, osps
c_4_6_cf_osps_rle_core_a1_n <- length(rownames(c_4_6_cf_osps_rle_core_a1))
# merge information
c_cf_core_osps_rle_a1_list <- c(c_0_cf_osps_rle_core_a1_n, c_1_3_cf_osps_rle_core_a1_n, c_4_6_cf_osps_rle_core_a1_n)

# number of core species, healthy, year 0, 25th abundance percentile, RLE, osps
c_0_cf_osps_rle_core_a2_n <- length(rownames(c_0_cf_osps_rle_core_a2))
# number of core species, healthy, year 1-3, 25th abundance percentile, RLE, osps
c_1_3_cf_osps_rle_core_a2_n <- length(rownames(c_1_3_cf_osps_rle_core_a2))
# number of core species, healthy, year 4-6, 25th abundance percentile, RLE, osps
c_4_6_cf_osps_rle_core_a2_n <- length(rownames(c_4_6_cf_osps_rle_core_a2))
# merge information
c_cf_core_osps_rle_a2_list <- c(c_0_cf_osps_rle_core_a2_n, c_1_3_cf_osps_rle_core_a2_n, c_4_6_cf_osps_rle_core_a2_n)

# number of core species, healthy, year 0, 35th abundance percentile, RLE, osps
c_0_cf_osps_rle_core_a3_n <- length(rownames(c_0_cf_osps_rle_core_a3))
# number of core species, healthy, year 1-3, 35th abundance percentile, RLE, osps
c_1_3_cf_osps_rle_core_a3_n <- length(rownames(c_1_3_cf_osps_rle_core_a3))
# number of core species, healthy, year 4-6, 35th abundance percentile, RLE, osps
c_4_6_cf_osps_rle_core_a3_n <- length(rownames(c_4_6_cf_osps_rle_core_a3))
# merge information
c_cf_core_osps_rle_a3_list <- c(c_0_cf_osps_rle_core_a3_n, c_1_3_cf_osps_rle_core_a3_n, c_4_6_cf_osps_rle_core_a3_n)

# number of rare species, healthy, year 0, 15th abundance percentile, RLE, osps
c_0_cf_osps_rle_rare_a1_n <- length(rownames(c_0_cf_osps_rle_rare_a1))
# number of rare species, healthy, year 1-3, 15th abundance percentile, RLE, osps
c_1_3_cf_osps_rle_rare_a1_n <- length(rownames(c_1_3_cf_osps_rle_rare_a1))
# number of rare species, healthy, year 4-6, 15th abundance percentile, RLE, osps
c_4_6_cf_osps_rle_rare_a1_n <- length(rownames(c_4_6_cf_osps_rle_rare_a1))
# merge information
c_cf_rare_osps_rle_a1_list <- c(c_0_cf_osps_rle_rare_a1_n, c_1_3_cf_osps_rle_rare_a1_n, c_4_6_cf_osps_rle_rare_a1_n)

# number of rare species, healthy, year 0, 25th abundance percentile, RLE, osps
c_0_cf_osps_rle_rare_a2_n <- length(rownames(c_0_cf_osps_rle_rare_a2))
# number of rare species, healthy, year 1-3, 25th abundance percentile, RLE, osps
c_1_3_cf_osps_rle_rare_a2_n <- length(rownames(c_1_3_cf_osps_rle_rare_a2))
# number of rare species, healthy, year 4-6, 25th abundance percentile, RLE, osps
c_4_6_cf_osps_rle_rare_a2_n <- length(rownames(c_4_6_cf_osps_rle_rare_a2))
# merge information
c_cf_rare_osps_rle_a2_list <- c(c_0_cf_osps_rle_rare_a2_n, c_1_3_cf_osps_rle_rare_a2_n, c_4_6_cf_osps_rle_rare_a2_n)

# number of rare species, healthy, year 0, 35th abundance percentile, RLE, osps
c_0_cf_osps_rle_rare_a3_n <- length(rownames(c_0_cf_osps_rle_rare_a3))
# number of rare species, healthy, year 1-3, 35th abundance percentile, RLE, osps
c_1_3_cf_osps_rle_rare_a3_n <- length(rownames(c_1_3_cf_osps_rle_rare_a3))
# number of rare species, healthy, year 4-6, 35th abundance percentile, RLE, osps
c_4_6_cf_osps_rle_rare_a3_n <- length(rownames(c_4_6_cf_osps_rle_rare_a3))
# merge information
c_cf_rare_osps_rle_a3_list <- c(c_0_cf_osps_rle_rare_a3_n, c_1_3_cf_osps_rle_rare_a3_n, c_4_6_cf_osps_rle_rare_a3_n)

# number of core species, healthy, year 0, 15th abundance percentile, RLE, osps
c_0_h_osps_rle_core_a1_n <- length(rownames(c_0_h_osps_rle_core_a1))
# number of core species, healthy, year 1-3, 15th abundance percentile, RLE, osps
c_1_3_h_osps_rle_core_a1_n <- length(rownames(c_1_3_h_osps_rle_core_a1))
# number of core species, healthy, year 4-6, 15th abundance percentile, RLE, osps
c_4_6_h_osps_rle_core_a1_n <- length(rownames(c_4_6_h_osps_rle_core_a1))
# merge information
c_h_core_osps_rle_a1_list <- c(c_0_h_osps_rle_core_a1_n, c_1_3_h_osps_rle_core_a1_n, c_4_6_h_osps_rle_core_a1_n)

# number of core species, healthy, year 0, 25th abundance percentile, RLE, osps
c_0_h_osps_rle_core_a2_n <- length(rownames(c_0_h_osps_rle_core_a2))
# number of core species, healthy, year 1-3, 25th abundance percentile, RLE, osps
c_1_3_h_osps_rle_core_a2_n <- length(rownames(c_1_3_h_osps_rle_core_a2))
# number of core species, healthy, year 4-6, 25th abundance percentile, RLE, osps
c_4_6_h_osps_rle_core_a2_n <- length(rownames(c_4_6_h_osps_rle_core_a2))
# merge information
c_h_core_osps_rle_a2_list <- c(c_0_h_osps_rle_core_a2_n, c_1_3_h_osps_rle_core_a2_n, c_4_6_h_osps_rle_core_a2_n)

# number of core species, healthy, year 0, 35th abundance percentile, RLE, osps
c_0_h_osps_rle_core_a3_n <- length(rownames(c_0_h_osps_rle_core_a3))
# number of core species, healthy, year 1-3, 35th abundance percentile, RLE, osps
c_1_3_h_osps_rle_core_a3_n <- length(rownames(c_1_3_h_osps_rle_core_a3))
# number of core species, healthy, year 4-6, 35th abundance percentile, RLE, osps
c_4_6_h_osps_rle_core_a3_n <- length(rownames(c_4_6_h_osps_rle_core_a3))
# merge information
c_h_core_osps_rle_a3_list <- c(c_0_h_osps_rle_core_a3_n, c_1_3_h_osps_rle_core_a3_n, c_4_6_h_osps_rle_core_a3_n)

# number of rare species, healthy, year 0, 15th abundance percentile, RLE, osps
c_0_h_osps_rle_rare_a1_n <- length(rownames(c_0_h_osps_rle_rare_a1))
# number of rare species, healthy, year 1-3, 15th abundance percentile, RLE, osps
c_1_3_h_osps_rle_rare_a1_n <- length(rownames(c_1_3_h_osps_rle_rare_a1))
# number of rare species, healthy, year 4-6, 15th abundance percentile, RLE, osps
c_4_6_h_osps_rle_rare_a1_n <- length(rownames(c_4_6_h_osps_rle_rare_a1))
# merge information
c_h_rare_osps_rle_a1_list <- c(c_0_h_osps_rle_rare_a1_n, c_1_3_h_osps_rle_rare_a1_n, c_4_6_h_osps_rle_rare_a1_n)

# number of rare species, healthy, year 0, 25th abundance percentile, RLE, osps
c_0_h_osps_rle_rare_a2_n <- length(rownames(c_0_h_osps_rle_rare_a2))
# number of rare species, healthy, year 1-3, 25th abundance percentile, RLE, osps
c_1_3_h_osps_rle_rare_a2_n <- length(rownames(c_1_3_h_osps_rle_rare_a2))
# number of rare species, healthy, year 4-6, 25th abundance percentile, RLE, osps
c_4_6_h_osps_rle_rare_a2_n <- length(rownames(c_4_6_h_osps_rle_rare_a2))
# merge information
c_h_rare_osps_rle_a2_list <- c(c_0_h_osps_rle_rare_a2_n, c_1_3_h_osps_rle_rare_a2_n, c_4_6_h_osps_rle_rare_a2_n)

# number of rare species, healthy, year 0, 35th abundance percentile, RLE, osps
c_0_h_osps_rle_rare_a3_n <- length(rownames(c_0_h_osps_rle_rare_a3))
# number of rare species, healthy, year 1-3, 35th abundance percentile, RLE, osps
c_1_3_h_osps_rle_rare_a3_n <- length(rownames(c_1_3_h_osps_rle_rare_a3))
# number of rare species, healthy, year 4-6, 35th abundance percentile, RLE, osps
c_4_6_h_osps_rle_rare_a3_n <- length(rownames(c_4_6_h_osps_rle_rare_a3))
# merge information
c_h_rare_osps_rle_a3_list <- c(c_0_h_osps_rle_rare_a3_n, c_1_3_h_osps_rle_rare_a3_n, c_4_6_h_osps_rle_rare_a3_n)

# number of core species, CF, year 0, 15th abundance percentile, VST, osps
c_0_cf_osps_vst_core_a1_n <- length(rownames(c_0_cf_osps_vst_core_a1))
# number of core species, CF, year 1-3, 15th abundance percentile, VST, osps
c_1_3_cf_osps_vst_core_a1_n <- length(rownames(c_1_3_cf_osps_vst_core_a1))
# number of core species, CF, year 4-6, 15th abundance percentile, VST, osps
c_4_6_cf_osps_vst_core_a1_n <- length(rownames(c_4_6_cf_osps_vst_core_a1))
# merge information
c_cf_core_osps_vst_a1_list <- c(c_0_cf_osps_vst_core_a1_n, c_1_3_cf_osps_vst_core_a1_n, c_4_6_cf_osps_vst_core_a1_n)

# number of core species, CF, year 0, 25th abundance percentile, VST, osps
c_0_cf_osps_vst_core_a2_n <- length(rownames(c_0_cf_osps_vst_core_a2))
# number of core species, CF, year 1-3, 25th abundance percentile, VST, osps
c_1_3_cf_osps_vst_core_a2_n <- length(rownames(c_1_3_cf_osps_vst_core_a2))
# number of core species, CF, year 4-6, 25th abundance percentile, VST, osps
c_4_6_cf_osps_vst_core_a2_n <- length(rownames(c_4_6_cf_osps_vst_core_a2))
# merge information
c_cf_core_osps_vst_a2_list <- c(c_0_cf_osps_vst_core_a2_n, c_1_3_cf_osps_vst_core_a2_n, c_4_6_cf_osps_vst_core_a2_n)

# number of core species, CF, year 0, 35th abundance percentile, VST, osps
c_0_cf_osps_vst_core_a3_n <- length(rownames(c_0_cf_osps_vst_core_a3))
# number of core species, CF, year 1-3, 35th abundance percentile, VST, osps
c_1_3_cf_osps_vst_core_a3_n <- length(rownames(c_1_3_cf_osps_vst_core_a3))
# number of core species, CF, year 4-6, 35th abundance percentile, VST, osps
c_4_6_cf_osps_vst_core_a3_n <- length(rownames(c_4_6_cf_osps_vst_core_a3))
# merge information
c_cf_core_osps_vst_a3_list <- c(c_0_cf_osps_vst_core_a3_n, c_1_3_cf_osps_vst_core_a3_n, c_4_6_cf_osps_vst_core_a3_n)

# number of rare species, CF, year 0, 15th abundance percentile, VST, osps
c_0_cf_osps_vst_rare_a1_n <- length(rownames(c_0_cf_osps_vst_rare_a1))
# number of rare species, CF, year 1-3, 15th abundance percentile, VST, osps
c_1_3_cf_osps_vst_rare_a1_n <- length(rownames(c_1_3_cf_osps_vst_rare_a1))
# number of rare species, CF, year 4-6, 15th abundance percentile, VST, osps
c_4_6_cf_osps_vst_rare_a1_n <- length(rownames(c_4_6_cf_osps_vst_rare_a1))
# merge information
c_cf_rare_osps_vst_a1_list <- c(c_0_cf_osps_vst_rare_a1_n, c_1_3_cf_osps_vst_rare_a1_n, c_4_6_cf_osps_vst_rare_a1_n)

# number of rare species, CF, year 0, 25th abundance percentile, VST, osps
c_0_cf_osps_vst_rare_a2_n <- length(rownames(c_0_cf_osps_vst_rare_a2))
# number of rare species, CF, year 1-3, 25th abundance percentile, VST, osps
c_1_3_cf_osps_vst_rare_a2_n <- length(rownames(c_1_3_cf_osps_vst_rare_a2))
# number of rare species, CF, year 4-6, 25th abundance percentile, VST, osps
c_4_6_cf_osps_vst_rare_a2_n <- length(rownames(c_4_6_cf_osps_vst_rare_a2))
# merge information
c_cf_rare_osps_vst_a2_list <- c(c_0_cf_osps_vst_rare_a2_n, c_1_3_cf_osps_vst_rare_a2_n, c_4_6_cf_osps_vst_rare_a2_n)

# number of rare species, CF, year 0, 35th abundance percentile, VST, osps
c_0_cf_osps_vst_rare_a3_n <- length(rownames(c_0_cf_osps_vst_rare_a3))
# number of rare species, CF, year 1-3, 35th abundance percentile, VST, osps
c_1_3_cf_osps_vst_rare_a3_n <- length(rownames(c_1_3_cf_osps_vst_rare_a3))
# number of rare species, CF, year 4-6, 35th abundance percentile, VST, osps
c_4_6_cf_osps_vst_rare_a3_n <- length(rownames(c_4_6_cf_osps_vst_rare_a3))
# merge information
c_cf_rare_osps_vst_a3_list <- c(c_0_cf_osps_vst_rare_a3_n, c_1_3_cf_osps_vst_rare_a3_n, c_4_6_cf_osps_vst_rare_a3_n)

# number of core species, Healthy, year 0, 15th abundance percentile, VST, osps
c_0_h_osps_vst_core_a1_n <- length(rownames(c_0_h_osps_vst_core_a1))
# number of core species, Healthy, year 1-3, 15th abundance percentile, VST, osps
c_1_3_h_osps_vst_core_a1_n <- length(rownames(c_1_3_h_osps_vst_core_a1))
# number of core species, Healthy, year 4-6, 15th abundance percentile, VST, osps
c_4_6_h_osps_vst_core_a1_n <- length(rownames(c_4_6_h_osps_vst_core_a1))
# merge information
c_h_core_osps_vst_a1_list <- c(c_0_h_osps_vst_core_a1_n, c_1_3_h_osps_vst_core_a1_n, c_4_6_h_osps_vst_core_a1_n)

# number of core species, Healthy, year 0, 25th abundance percentile, VST, osps
c_0_h_osps_vst_core_a2_n <- length(rownames(c_0_h_osps_vst_core_a2))
# number of core species, Healthy, year 1-3, 25th abundance percentile, VST, osps
c_1_3_h_osps_vst_core_a2_n <- length(rownames(c_1_3_h_osps_vst_core_a2))
# number of core species, Healthy, year 4-6, 25th abundance percentile, VST, osps
c_4_6_h_osps_vst_core_a2_n <- length(rownames(c_4_6_h_osps_vst_core_a2))
# merge information
c_h_core_osps_vst_a2_list <- c(c_0_h_osps_vst_core_a2_n, c_1_3_h_osps_vst_core_a2_n, c_4_6_h_osps_vst_core_a2_n)

# number of core species, Healthy, year 0, 35th abundance percentile, VST, osps
c_0_h_osps_vst_core_a3_n <- length(rownames(c_0_h_osps_vst_core_a3))
# number of core species, Healthy, year 1-3, 35th abundance percentile, VST, osps
c_1_3_h_osps_vst_core_a3_n <- length(rownames(c_1_3_h_osps_vst_core_a3))
# number of core species, Healthy, year 4-6, 35th abundance percentile, VST, osps
c_4_6_h_osps_vst_core_a3_n <- length(rownames(c_4_6_h_osps_vst_core_a3))
# merge information
c_h_core_osps_vst_a3_list <- c(c_0_h_osps_vst_core_a3_n, c_1_3_h_osps_vst_core_a3_n, c_4_6_h_osps_vst_core_a3_n)

# number of rare species, Healthy, year 0, 15th abundance percentile, VST, osps
c_0_h_osps_vst_rare_a1_n <- length(rownames(c_0_h_osps_vst_rare_a1))
# number of rare species, Healthy, year 0, 15th abundance percentile, VST, osps
c_1_3_h_osps_vst_rare_a1_n <- length(rownames(c_1_3_h_osps_vst_rare_a1))
# number of rare species, Healthy, year 0, 15th abundance percentile, VST, osps
c_4_6_h_osps_vst_rare_a1_n <- length(rownames(c_4_6_h_osps_vst_rare_a1))
# merge information
c_h_rare_osps_vst_a1_list <- c(c_0_h_osps_vst_rare_a1_n, c_1_3_h_osps_vst_rare_a1_n, c_4_6_h_osps_vst_rare_a1_n)

# number of rare species, Healthy, year 0, 25th abundance percentile, VST, osps
c_0_h_osps_vst_rare_a2_n <- length(rownames(c_0_h_osps_vst_rare_a2))
# number of rare species, Healthy, year 1-3, 25th abundance percentile, VST, osps
c_1_3_h_osps_vst_rare_a2_n <- length(rownames(c_1_3_h_osps_vst_rare_a2))
# number of rare species, Healthy, year 4-6, 25th abundance percentile, VST, osps
c_4_6_h_osps_vst_rare_a2_n <- length(rownames(c_4_6_h_osps_vst_rare_a2))
# merge information
c_h_rare_osps_vst_a2_list <- c(c_0_h_osps_vst_rare_a2_n, c_1_3_h_osps_vst_rare_a2_n, c_4_6_h_osps_vst_rare_a2_n)

# number of rare species, Healthy, year 0, 35th abundance percentile, VST, osps
c_0_h_osps_vst_rare_a3_n <- length(rownames(c_0_h_osps_vst_rare_a3))
# number of rare species, Healthy, year 1-3, 35th abundance percentile, VST, osps
c_1_3_h_osps_vst_rare_a3_n <- length(rownames(c_1_3_h_osps_vst_rare_a3))
# number of rare species, Healthy, year 4-6, 35th abundance percentile, VST, osps
c_4_6_h_osps_vst_rare_a3_n <- length(rownames(c_4_6_h_osps_vst_rare_a3))
# merge information
c_h_rare_osps_vst_a3_list <- c(c_0_h_osps_vst_rare_a3_n, c_1_3_h_osps_vst_rare_a3_n, c_4_6_h_osps_vst_rare_a3_n)


# merge all information 
# 15th abundance percentile, BCPHC 
age_osps_a1 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_osps_a1 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_osps_a1 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_osps_a1 <- c(c_cf_core_osps_a1_list, c_cf_rare_osps_a1_list, c_h_core_osps_a1_list, c_h_rare_osps_a1_list)

merge_data_venn_osps_a1 <- data.frame(cbind(age_osps_a1, state_osps_a1, species_type_osps_a1, c_osps_a1))
# add metadata
colnames(merge_data_venn_osps_a1) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_osps_a1$Threshold <- '15th percentile'
merge_data_venn_osps_a1$Normalisation <- 'BCPHC'
# convert age group to class factor and specify order
merge_data_venn_osps_a1$Age <- factor(merge_data_venn_osps_a1$Age, levels = c('0', '1-3', '4-6'))


# 25th abundance percentile, BCPHC 
age_osps_a2 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_osps_a2 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_osps_a2 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_osps_a2 <- c(c_cf_core_osps_a2_list, c_cf_rare_osps_a2_list, c_h_core_osps_a2_list, c_h_rare_osps_a2_list)

merge_data_venn_osps_a2 <- data.frame(cbind(age_osps_a2, state_osps_a2, species_type_osps_a2, c_osps_a2))
# add metadata
colnames(merge_data_venn_osps_a2) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_osps_a2$Threshold <- '25th percentile'
merge_data_venn_osps_a2$Normalisation <- 'BCPHC'
# convert age group to class factor and specify order
merge_data_venn_osps_a2$Age <- factor(merge_data_venn_osps_a2$Age, levels = c('0', '1-3', '4-6'))

# 35th abundance percentile, BCPHC 
age_osps_a3 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_osps_a3 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_osps_a3 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_osps_a3 <- c(c_cf_core_osps_a3_list, c_cf_rare_osps_a3_list, c_h_core_osps_a3_list, c_h_rare_osps_a3_list)
# add metadata 
merge_data_venn_osps_a3 <- data.frame(cbind(age_osps_a3, state_osps_a3, species_type_osps_a3, c_osps_a3))
colnames(merge_data_venn_osps_a3) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_osps_a3$Threshold <- '35th percentile'
merge_data_venn_osps_a3$Normalisation <- 'BCPHC'
# convert age group to class factor and specify order
merge_data_venn_osps_a3$Age <- factor(merge_data_venn_osps_a3$Age, levels = c('0', '1-3', '4-6'))


# 15th abundance percentile, RLE 
age_osps_rle_a1 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_osps_rle_a1 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_osps_rle_a1 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_osps_rle_a1 <- c(c_cf_core_osps_rle_a1_list, c_cf_rare_osps_rle_a1_list, c_h_core_osps_rle_a1_list, c_h_rare_osps_rle_a1_list)
# add metadata 
merge_data_venn_osps_rle_a1 <- data.frame(cbind(age_osps_rle_a1, state_osps_rle_a1, species_type_osps_rle_a1, c_osps_rle_a1))
colnames(merge_data_venn_osps_rle_a1) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_osps_rle_a1$Threshold <- '15th percentile'
merge_data_venn_osps_rle_a1$Normalisation <- 'RLE'
# convert age group to class factor and specify order
merge_data_venn_osps_rle_a1$Age <- factor(merge_data_venn_osps_rle_a1$Age, levels = c('0', '1-3', '4-6'))

# 25th abundance percentile, RLE 
age_osps_rle_a2 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_osps_rle_a2 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_osps_rle_a2 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_osps_rle_a2 <- c(c_cf_core_osps_rle_a2_list, c_cf_rare_osps_rle_a2_list, c_h_core_osps_rle_a2_list, c_h_rare_osps_rle_a2_list)
# add metadata 
merge_data_venn_osps_rle_a2 <- data.frame(cbind(age_osps_rle_a2, state_osps_rle_a2, species_type_osps_rle_a2, c_osps_rle_a2))
colnames(merge_data_venn_osps_rle_a2) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_osps_rle_a2$Threshold <- '25th percentile'
merge_data_venn_osps_rle_a2$Normalisation <- 'RLE'
# convert age group to class factor and specify order
merge_data_venn_osps_rle_a2$Age <- factor(merge_data_venn_osps_rle_a2$Age, levels = c('0', '1-3', '4-6'))

# 35th abundance percentile, RLE 
age_osps_rle_a3 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_osps_rle_a3 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_osps_rle_a3 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_osps_rle_a3 <- c(c_cf_core_osps_rle_a3_list, c_cf_rare_osps_rle_a3_list, c_h_core_osps_rle_a3_list, c_h_rare_osps_rle_a3_list)
# add metadata 
merge_data_venn_osps_rle_a3 <- data.frame(cbind(age_osps_rle_a3, state_osps_rle_a3, species_type_osps_rle_a3, c_osps_rle_a3))
colnames(merge_data_venn_osps_rle_a3) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_osps_rle_a3$Threshold <- '35th percentile'
merge_data_venn_osps_rle_a3$Normalisation <- 'RLE'
# convert age group to class factor and specify order
merge_data_venn_osps_rle_a3$Age <- factor(merge_data_venn_osps_rle_a3$Age, levels = c('0', '1-3', '4-6'))


# 15th abundance percentile, VST 
age_osps_vst_a1 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_osps_vst_a1 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_osps_vst_a1 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_osps_vst_a1 <- c(c_cf_core_osps_vst_a1_list, c_cf_rare_osps_vst_a1_list, c_h_core_osps_vst_a1_list, c_h_rare_osps_vst_a1_list)
# add metadata 
merge_data_venn_osps_vst_a1 <- data.frame(cbind(age_osps_vst_a1, state_osps_vst_a1, species_type_osps_vst_a1, c_osps_vst_a1))
colnames(merge_data_venn_osps_vst_a1) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_osps_vst_a1$Threshold <- '15th percentile'
merge_data_venn_osps_vst_a1$Normalisation <- 'VST'
# convert age group to class factor and specify order
merge_data_venn_osps_vst_a1$Age <- factor(merge_data_venn_osps_vst_a1$Age, levels = c('0', '1-3', '4-6'))

# 25th abundance percentile, VST 
age_osps_vst_a2 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_osps_vst_a2 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_osps_vst_a2 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_osps_vst_a2 <- c(c_cf_core_osps_vst_a2_list, c_cf_rare_osps_vst_a2_list, c_h_core_osps_vst_a2_list, c_h_rare_osps_vst_a2_list)
# add metadata 
merge_data_venn_osps_vst_a2 <- data.frame(cbind(age_osps_vst_a2, state_osps_vst_a2, species_type_osps_vst_a2, c_osps_vst_a2))
colnames(merge_data_venn_osps_vst_a2) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_osps_vst_a2$Threshold <- '25th percentile'
merge_data_venn_osps_vst_a2$Normalisation <- 'VST'
# convert age group to class factor and specify order
merge_data_venn_osps_vst_a2$Age <- factor(merge_data_venn_osps_vst_a2$Age, levels = c('0', '1-3', '4-6'))

# 35th abundance percentile, VST 
age_osps_vst_a3 <- c('0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6', '0', '1-3', '4-6')
state_osps_vst_a3 <- c('CF', 'CF', 'CF', 'CF', 'CF', 'CF', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy', 'Healthy')
species_type_osps_vst_a3 <- c('core', 'core', 'core', 'rare', 'rare', 'rare', 'core', 'core', 'core', 'rare', 'rare', 'rare')
c_osps_vst_a3 <- c(c_cf_core_osps_vst_a3_list, c_cf_rare_osps_vst_a3_list, c_h_core_osps_vst_a3_list, c_h_rare_osps_vst_a3_list)
# add metadata 
merge_data_venn_osps_vst_a3 <- data.frame(cbind(age_osps_vst_a3, state_osps_vst_a3, species_type_osps_vst_a3, c_osps_vst_a3))
colnames(merge_data_venn_osps_vst_a3) <- c('Age', 'State', 'Species_type', 'Species_count')
merge_data_venn_osps_vst_a3$Threshold <- '35th percentile'
merge_data_venn_osps_vst_a3$Normalisation <- 'VST'
# convert age group to class factor and specify order
merge_data_venn_osps_vst_a3$Age <- factor(merge_data_venn_osps_vst_a3$Age, levels = c('0', '1-3', '4-6'))

# merge all information in one data frame
merge_data_venn_osps <- data.frame(rbind(merge_data_venn_osps_a1, merge_data_venn_osps_a2, merge_data_venn_osps_a3,
                                         merge_data_venn_osps_rle_a1, merge_data_venn_osps_rle_a2, merge_data_venn_osps_rle_a3,
                                         merge_data_venn_osps_vst_a1, merge_data_venn_osps_vst_a2, merge_data_venn_osps_vst_a3))
# add column with osps information
merge_data_venn_osps$Database <- 'One-strain per species'

# bind pangenome and one strain per species database
merge_data_venn <- data.frame(rbind(merge_data_venn_pangenome, merge_data_venn_osps))

# add column with merged information
merge_data_venn$State_Threshold <- paste(merge_data_venn$State,'_',merge_data_venn$Threshold)
# convert species type to class factor
merge_data_venn$Species_type <- factor(merge_data_venn$Species_type, labels = c('Core species biosphere', 'Rare species biosphere'))
# remove empty spaces
merge_data_venn$State_Threshold <- str_replace_all(merge_data_venn$State_Threshold, ' _ ', '_')
# convert species number to class numeric
merge_data_venn$Species_count <- as.numeric(as.character(merge_data_venn$Species_count))

# subset dataset and get rare species
merge_data_venn_rare <- subset(merge_data_venn, Species_type == 'Rare species biosphere')
# subset dataset and get core species
merge_data_venn_core <- subset(merge_data_venn, Species_type == 'Core species biosphere')

# statistical comparison of species number and disease state, grouped by age, rare species
merge_data_venn_rare_wilcox <- compare_means(Species_count~State, data=merge_data_venn_rare, group.by = "Age")
merge_data_venn_rare_wilcox$Species_type <- "rare"
# statistical comparison of species number and disease state, grouped by age, core species
merge_data_venn_core_wilcox <- compare_means(Species_count~State, data=merge_data_venn_core, group.by = "Age")
merge_data_venn_core_wilcox$Species_type <- "core"

# generate a dataframe of statistics output
venn_core_rare_stats <- data.frame(rbind(merge_data_venn_rare_wilcox, merge_data_venn_core_wilcox))

# compute r effect size with confidence intervals, rare species and store in data frame
venn_rare_effsize <- merge_data_venn_rare %>% group_by(Age) %>% wilcox_effsize(Species_count~State, ci=TRUE)
venn_rare_effsize$Species_type <- "rare"

# compute r effect size with confidence intervals, core species and store in data frame
venn_core_effsize <- merge_data_venn_core %>% group_by(Age) %>% wilcox_effsize(Species_count~State, ci=TRUE)
venn_core_effsize$Species_type <- "core"

# make one large statistics table
venn_core_rare_stats_effsize <- data.frame(rbind(venn_rare_effsize, venn_core_effsize))
# round effect size to two decimal places
venn_core_rare_stats_effsize$effsize <- round(venn_core_rare_stats_effsize$effsize,2)

# generate plot for core species biosphere
merge_data_venn_core_osps <- subset(merge_data_venn_core, Database == "One-strain per species")
core_line_osps <- ggline(merge_data_venn_core_osps, x='Age', y='Species_count', color='State', add = c('mean_se'), size=0.2) +
  theme_bw(base_size=4) + xlab('Age (in years)') + ylab('Number of species\n') +
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        axis.text.x = element_text(size=11), axis.text.y = element_text(size=11), axis.title.x = element_text(size=11), axis.title.y = element_text(size=11),
        legend.text = element_text(size=11)) +
  scale_colour_manual(values=c('red', 'black')) +
  scale_x_discrete(labels=c('0', '1-3', '4-6')) + ylim(0,115) +
  geom_label(aes(x=2, y=108), label='One-strain-per-species database, \nCore species biosphere', size=3, bg='white')

merge_data_venn_rare_osps <- subset(merge_data_venn_rare, Database == "One-strain per species")
# generate plot for rare species biosphere
rare_line_osps <- ggline(merge_data_venn_rare_osps, x='Age', y='Species_count', color='State', add = c('mean_se'), size=0.2) +
  theme_bw(base_size=4) + xlab('Age (in years)') + ylab('\n ') +
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        axis.text.x = element_text(size=11), axis.text.y = element_text(size=11), axis.title.x = element_text(size=11), axis.title.y = element_text(size=11),
        legend.text = element_text(size=11)) + 
  scale_colour_manual(values=c('red', 'black')) +
  scale_x_discrete(labels=c('0', '1-3', '4-6')) + ylim(0,115) +
  geom_label(aes(x=2, y=108), label='One-strain-per-species database, \nRare species biosphere', size=3, bg='white')

merge_data_venn_core_pan<- subset(merge_data_venn_core, Database == "Pangenome")
core_line_pan <- ggline(merge_data_venn_core_pan, x='Age', y='Species_count', color='State', add = c('mean_se'), size=0.2) +
  theme_bw(base_size=4) + xlab('Age (in years)') + ylab('Number of species\n') +
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        axis.text.x = element_text(size=11), axis.text.y = element_text(size=11), axis.title.x = element_text(size=11), axis.title.y = element_text(size=11),
        legend.text = element_text(size=11)) +
  scale_colour_manual(values=c('red', 'black')) +
  scale_x_discrete(labels=c('0', '1-3', '4-6')) + ylim(0,115) +
  geom_label(aes(x=2, y=108), label='Pan-genome database, \nCore species biosphere', size=3, bg='white')

merge_data_venn_rare_pan <- subset(merge_data_venn_rare, Database == "Pangenome")
# generate plot for rare species biosphere
rare_line_pan<- ggline(merge_data_venn_rare_pan, x='Age', y='Species_count', color='State', add = c('mean_se'), size=0.2) +
  theme_bw(base_size=4) + xlab('Age (in years)') + ylab('\n ') +
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        axis.text.x = element_text(size=11), axis.text.y = element_text(size=11), axis.title.x = element_text(size=11), axis.title.y = element_text(size=11),
        legend.text = element_text(size=11)) + 
  scale_colour_manual(values=c('red', 'black')) +
  scale_x_discrete(labels=c('0', '1-3', '4-6')) + ylim(0,115) +
  geom_label(aes(x=2, y=108), label='Pan-genome database, \nRare species biosphere', size=3, bg='white')

# merge plot for core and rare species biosphere
line_plots <- ggarrange(core_line_osps, rare_line_osps, core_line_pan, rare_line_pan, labels=c('A', 'B', 'C', 'D'), common.legend = TRUE, font.label = list(size = 12, color = "black"))

#########################################################################################################################################
# Statistical comparison of count data
# Fisher's Exact Test for Count Data, pan-genome, 0 years, core species
merge_data_venn_core_pan_0 <- subset(merge_data_venn_core_pan, Age == "0")
table_it_0_core <- table(merge_data_venn_core_pan_0$Species_count,merge_data_venn_core_pan_0$State)
table_it_df_0_core <- data.frame(table_it_0_core)
table_it_df_h_0_core <- subset(table_it_df_0_core, Var2 == "Healthy")
table_it_df_cf_0_core <- subset(table_it_df_0_core, Var2 == "CF")
table_it_df_cf_0_core$Freq2 <- table_it_df_h_0_core$Freq
table_it_df_cf_0_core$Var1 <- NULL
table_it_df_cf_0_core$Var2 <- NULL
fisher.test(table_it_df_cf_0_core, conf.int = TRUE) 

# Fisher's Exact Test for Count Data, pan-genome, 1-3 years, core species
merge_data_venn_core_pan_1_3 <- subset(merge_data_venn_core_pan, Age == "1-3")
table_it_1_3_core <- table(merge_data_venn_core_pan_1_3$Species_count,merge_data_venn_core_pan_1_3$State)
table_it_df_1_3_core <- data.frame(table_it_1_3_core)
table_it_df_h_1_3_core <- subset(table_it_df_1_3_core, Var2 == "Healthy")
table_it_df_cf_1_3_core <- subset(table_it_df_1_3_core, Var2 == "CF")
table_it_df_cf_1_3_core$Freq2 <- table_it_df_h_1_3_core$Freq
table_it_df_cf_1_3_core$Var1 <- NULL
table_it_df_cf_1_3_core$Var2 <- NULL
fisher.test(table_it_df_cf_1_3_core, conf.int = TRUE)


# Fisher's Exact Test for Count Data, pan-genome, 1-3 years, core species
merge_data_venn_core_pan_4_6 <- subset(merge_data_venn_core_pan, Age == "4-6")
table_it_4_6_core <- table(merge_data_venn_core_pan_4_6$Species_count, merge_data_venn_core_pan_4_6$State)
table_it_df_4_6_core <- data.frame(table_it_4_6_core)
table_it_df_h_4_6_core <- subset(table_it_df_4_6_core, Var2 == "Healthy")
table_it_df_cf_4_6_core <- subset(table_it_df_4_6_core, Var2 == "CF")
table_it_df_cf_4_6_core$Freq2 <- table_it_df_h_4_6_core$Freq
table_it_df_cf_4_6_core$Var1 <- NULL
table_it_df_cf_4_6_core$Var2 <- NULL
fisher.test(table_it_df_cf_4_6_core, conf.int = TRUE) 


# Statistical comparison of count data
# Fisher's Exact Test for Count Data, pan-genome, 0 years, rare species
merge_data_venn_rare_pan_0 <- subset(merge_data_venn_rare_pan, Age == "0")
table_it_0_rare <- table(merge_data_venn_rare_pan_0$Species_count,merge_data_venn_rare_pan_0$State)
table_it_df_0_rare <- data.frame(table_it_0_rare)
table_it_df_h_0_rare <- subset(table_it_df_0_rare, Var2 == "Healthy")
table_it_df_cf_0_rare <- subset(table_it_df_0_rare, Var2 == "CF")
table_it_df_cf_0_rare$Freq2 <- table_it_df_h_0_rare$Freq
table_it_df_cf_0_rare$Var1 <- NULL
table_it_df_cf_0_rare$Var2 <- NULL
fisher.test(table_it_df_cf_0_rare, conf.int = TRUE) 

# Fisher's Exact Test for Count Data, pan-genome, 1-3 years, rare species
merge_data_venn_rare_pan_1_3 <- subset(merge_data_venn_rare_pan, Age == "1-3")
table_it_1_3_rare <- table(merge_data_venn_rare_pan_1_3$Species_count,merge_data_venn_rare_pan_1_3$State)
table_it_df_1_3_rare <- data.frame(table_it_1_3_rare)
table_it_df_h_1_3_rare <- subset(table_it_df_1_3_rare, Var2 == "Healthy")
table_it_df_cf_1_3_rare <- subset(table_it_df_1_3_rare, Var2 == "CF")
table_it_df_cf_1_3_rare$Freq2 <- table_it_df_h_1_3_rare$Freq
table_it_df_cf_1_3_rare$Var1 <- NULL
table_it_df_cf_1_3_rare$Var2 <- NULL
fisher.test(table_it_df_cf_1_3_rare, conf.int = TRUE)


# Fisher's Exact Test for Count Data, pan-genome, 4-6 years, rare species
merge_data_venn_rare_pan_4_6 <- subset(merge_data_venn_rare_pan, Age == "4-6")
table_it_4_6_rare <- table(merge_data_venn_rare_pan_4_6$Species_count, merge_data_venn_rare_pan_4_6$State)
table_it_df_4_6_rare <- data.frame(table_it_4_6_rare)
table_it_df_h_4_6_rare <- subset(table_it_df_4_6_rare, Var2 == "Healthy")
table_it_df_cf_4_6_rare <- subset(table_it_df_4_6_rare, Var2 == "CF")
table_it_df_cf_4_6_rare$Freq2 <- table_it_df_h_4_6_rare$Freq
table_it_df_cf_4_6_rare$Var1 <- NULL
table_it_df_cf_4_6_rare$Var2 <- NULL
fisher.test(table_it_df_cf_4_6_rare, conf.int = TRUE) 


# Statistical comparison of count data
# Fisher's Exact Test for Count Data, one-strain-per-species, 0 years, core species
merge_data_venn_core_osps_0 <- subset(merge_data_venn_core_osps, Age == "0")
table_it_0_core_osps <- table(merge_data_venn_core_osps_0$Species_count,merge_data_venn_core_osps_0$State)
table_it_df_0_core_osps <- data.frame(table_it_0_core_osps)
table_it_df_h_0_core_osps <- subset(table_it_df_0_core_osps, Var2 == "Healthy")
table_it_df_cf_0_core_osps <- subset(table_it_df_0_core_osps, Var2 == "CF")
table_it_df_cf_0_core_osps$Freq2 <- table_it_df_h_0_core_osps$Freq
table_it_df_cf_0_core_osps$Var1 <- NULL
table_it_df_cf_0_core_osps$Var2 <- NULL
fisher.test(table_it_df_cf_0_core_osps, conf.int = TRUE) 

# Fisher's Exact Test for Count Data, one-strain-per-species, 1-3 years, core species
merge_data_venn_core_osps_1_3 <- subset(merge_data_venn_core_osps, Age == "1-3")
table_it_1_3_core_osps <- table(merge_data_venn_core_osps_1_3$Species_count,merge_data_venn_core_osps_1_3$State)
table_it_df_1_3_core_osps <- data.frame(table_it_1_3_core_osps)
table_it_df_h_1_3_core_osps <- subset(table_it_df_1_3_core_osps, Var2 == "Healthy")
table_it_df_cf_1_3_core_osps <- subset(table_it_df_1_3_core_osps, Var2 == "CF")
table_it_df_cf_1_3_core_osps$Freq2 <- table_it_df_h_1_3_core_osps$Freq
table_it_df_cf_1_3_core_osps$Var1 <- NULL
table_it_df_cf_1_3_core_osps$Var2 <- NULL
fisher.test(table_it_df_cf_1_3_core, conf.int = TRUE) 


# Fisher's Exact Test for Count Data, one-strain-per-species, 1-3 years, core species
merge_data_venn_core_osps_4_6 <- subset(merge_data_venn_core_osps, Age == "4-6")
table_it_4_6_core_osps <- table(merge_data_venn_core_osps_4_6$Species_count, merge_data_venn_core_osps_4_6$State)
table_it_df_4_6_core_osps <- data.frame(table_it_4_6_core_osps)
table_it_df_h_4_6_core_osps <- subset(table_it_df_4_6_core_osps, Var2 == "Healthy")
table_it_df_cf_4_6_core_osps <- subset(table_it_df_4_6_core_osps, Var2 == "CF")
table_it_df_cf_4_6_core_osps$Freq2 <- table_it_df_h_4_6_core_osps$Freq
table_it_df_cf_4_6_core_osps$Var1 <- NULL
table_it_df_cf_4_6_core_osps$Var2 <- NULL
fisher.test(table_it_df_cf_4_6_core_osps, conf.int = TRUE) 


# Statistical comparison of count data
# Fisher's Exact Test for Count Data, one-strain-per-species, 0 years, rare species
merge_data_venn_rare_osps_0 <- subset(merge_data_venn_rare_osps, Age == "0")
table_it_0_rare_osps <- table(merge_data_venn_rare_osps_0$Species_count,merge_data_venn_rare_osps_0$State)
table_it_df_0_rare_osps <- data.frame(table_it_0_rare_osps)
table_it_df_h_0_rare_osps <- subset(table_it_df_0_rare_osps, Var2 == "Healthy")
table_it_df_cf_0_rare_osps <- subset(table_it_df_0_rare_osps, Var2 == "CF")
table_it_df_cf_0_rare_osps$Freq2 <- table_it_df_h_0_rare_osps$Freq
table_it_df_cf_0_rare_osps$Var1 <- NULL
table_it_df_cf_0_rare_osps$Var2 <- NULL
fisher.test(table_it_df_cf_0_rare_osps, conf.int = TRUE) 

# Fisher's Exact Test for Count Data, one-strain-per-species, 1-3 years, rare species
merge_data_venn_rare_osps_1_3 <- subset(merge_data_venn_rare_osps, Age == "1-3")
merge_data_venn_rare_osps
table_it_1_3_rare_osps <- table(merge_data_venn_rare_osps_1_3$Species_count,merge_data_venn_rare_osps_1_3$State)
table_it_df_1_3_rare_osps <- data.frame(table_it_1_3_rare_osps)
table_it_df_h_1_3_rare_osps <- subset(table_it_df_1_3_rare_osps, Var2 == "Healthy")
table_it_df_cf_1_3_rare_osps <- subset(table_it_df_1_3_rare_osps, Var2 == "CF")
table_it_df_cf_1_3_rare_osps$Freq2 <- table_it_df_h_1_3_rare_osps$Freq
table_it_df_cf_1_3_rare_osps$Var1 <- NULL
table_it_df_cf_1_3_rare_osps$Var2 <- NULL
fisher.test(table_it_df_cf_1_3_rare_osps, conf.int = TRUE) 


# Fisher's Exact Test for Count Data, one-strain-per-species, 1-3 years, rare species
merge_data_venn_rare_osps_4_6 <- subset(merge_data_venn_rare_osps, Age == "4-6")
table_it_4_6_rare_osps <- table(merge_data_venn_rare_osps_4_6$Species_count, merge_data_venn_rare_osps_4_6$State)
table_it_df_4_6_rare_osps <- data.frame(table_it_4_6_rare_osps)
table_it_df_h_4_6_rare_osps <- subset(table_it_df_4_6_rare_osps, Var2 == "Healthy")
table_it_df_cf_4_6_rare_osps <- subset(table_it_df_4_6_rare_osps, Var2 == "CF")
table_it_df_cf_4_6_rare_osps$Freq2 <- table_it_df_h_4_6_rare_osps$Freq
table_it_df_cf_4_6_rare_osps$Var1 <- NULL
table_it_df_cf_4_6_rare_osps$Var2 <- NULL
fisher.test(table_it_df_cf_4_6_rare_osps, conf.int = TRUE) 

############################################################################################################
# Preparation for random forest analysis
# Start with pangenome
# Obtain background core species (healty, cf), 15th abundance percentile, BCPHC, pangenome
f_core_species_pangenome_a1 <- c(rownames(background_core_h_pangenome_a1), rownames(background_core_cf_pangenome_a1))
# remove duplicate entries
f_core_species_pangenome_a1 <- f_core_species_pangenome_a1[!duplicated(f_core_species_pangenome_a1)]
# Obtain background rare species (healthy, cf), 15th abundance percentile, BCPHC, pangenome
f_rare_species_pangenome_a1 <- c(rownames(background_rare_h_pangenome_a1), rownames(background_rare_cf_pangenome_a1))
# remove duplicate entries
f_rare_species_pangenome_a1 <- f_rare_species_pangenome_a1[!duplicated(f_rare_species_pangenome_a1)]

# Obtain background core species (healty, cf), 25th abundance percentile, BCPHC, pangenome
f_core_species_pangenome_a2 <- c(rownames(background_core_h_pangenome_a2), rownames(background_core_cf_pangenome_a2))
# remove duplicate entries
f_core_species_pangenome_a2 <- f_core_species_pangenome_a2[!duplicated(f_core_species_pangenome_a2)]
# Obtain background rare species (healthy, cf), 25th abundance percentile, BCPHC, pangenome
f_rare_species_pangenome_a2 <- c(rownames(background_rare_h_pangenome_a2), rownames(background_rare_cf_pangenome_a2))
# remove duplicate entries
f_rare_species_pangenome_a2 <- f_rare_species_pangenome_a2[!duplicated(f_rare_species_pangenome_a2)]

# Obtain background core species (healty, cf), 35th abundance percentile, BCPHC, pangenome
f_core_species_pangenome_a3 <- c(rownames(background_core_h_pangenome_a3), rownames(background_core_cf_pangenome_a3))
# remove duplicate entries
f_core_species_pangenome_a3 <- f_core_species_pangenome_a3[!duplicated(f_core_species_pangenome_a3)]
# Obtain background rare species (healthy, cf), 35th abundance percentile, BCPHC, pangenome
f_rare_species_pangenome_a3 <- c(rownames(background_rare_h_pangenome_a3), rownames(background_rare_cf_pangenome_a3))
# remove duplicate entries
f_rare_species_pangenome_a3 <- f_rare_species_pangenome_a3[!duplicated(f_rare_species_pangenome_a3)]


# Obtain background core species (healty, cf), 15th abundance percentile, RLE, pangenome
f_core_species_pangenome_rle_a1 <- c(rownames(background_core_h_pangenome_rle_a1), rownames(background_core_cf_pangenome_rle_a1))
# remove duplicate entries
f_core_species_pangenome_rle_a1 <- f_core_species_pangenome_rle_a1[!duplicated(f_core_species_pangenome_rle_a1)]
# Obtain background rare species (healthy, cf), 15th abundance percentile, RLE, pangenome
f_rare_species_pangenome_rle_a1 <- c(rownames(background_rare_h_pangenome_rle_a1), rownames(background_rare_cf_pangenome_rle_a1))
# remove duplicate entries
f_rare_species_pangenome_rle_a1 <- f_rare_species_pangenome_rle_a1[!duplicated(f_rare_species_pangenome_rle_a1)]

# Obtain background core species (healty, cf), 25th abundance percentile, RLE, pangenome
f_core_species_pangenome_rle_a2 <- c(rownames(background_core_h_pangenome_rle_a2), rownames(background_core_cf_pangenome_rle_a2))
# remove duplicate entries
f_core_species_pangenome_rle_a2 <- f_core_species_pangenome_rle_a2[!duplicated(f_core_species_pangenome_rle_a2)]
# Obtain background rare species (healthy, cf), 25th abundance percentile, RLE, pangenome
f_rare_species_pangenome_rle_a2 <- c(rownames(background_rare_h_pangenome_rle_a2), rownames(background_rare_cf_pangenome_rle_a2))
# remove duplicate entries
f_rare_species_pangenome_rle_a2 <- f_rare_species_pangenome_rle_a2[!duplicated(f_rare_species_pangenome_rle_a2)]

# Obtain background core species (healty, cf), 35th abundance percentile, RLE, pangenome
f_core_species_pangenome_rle_a3 <- c(rownames(background_core_h_pangenome_rle_a3), rownames(background_core_cf_pangenome_rle_a3))
# remove duplicate entries
f_core_species_pangenome_rle_a3 <- f_core_species_pangenome_rle_a3[!duplicated(f_core_species_pangenome_rle_a3)]
# Obtain background rare species (healthy, cf), 35th abundance percentile, RLE, pangenome
f_rare_species_pangenome_rle_a3 <- c(rownames(background_rare_h_pangenome_rle_a3), rownames(background_rare_cf_pangenome_rle_a3))
# remove duplicate entries
f_rare_species_pangenome_rle_a3 <- f_rare_species_pangenome_rle_a3[!duplicated(f_rare_species_pangenome_rle_a3)]



# Obtain background core species (healthy, cf), 15th abundance percentile, VST, pangenome
f_core_species_pangenome_vst_a1 <- c(rownames(background_core_h_pangenome_vst_a1), rownames(background_core_cf_pangenome_vst_a1))
# remove duplicate entries
f_core_species_pangenome_vst_a1 <- f_core_species_pangenome_vst_a1[!duplicated(f_core_species_pangenome_vst_a1)]
# Obtain background rare species (healthy, cf), 15th abundance percentile, VST, pangenome
f_rare_species_pangenome_vst_a1 <- c(rownames(background_rare_h_pangenome_vst_a1), rownames(background_rare_cf_pangenome_vst_a1))
# remove duplicate entries
f_rare_species_pangenome_vst_a1 <- f_rare_species_pangenome_vst_a1[!duplicated(f_rare_species_pangenome_vst_a1)]

# Obtain background core species (healthy, cf), 25th abundance percentile, VST, pangenome
f_core_species_pangenome_vst_a2 <- c(rownames(background_core_h_pangenome_vst_a2), rownames(background_core_cf_pangenome_vst_a2))
# remove duplicate entries
f_core_species_pangenome_vst_a2 <- f_core_species_pangenome_vst_a2[!duplicated(f_core_species_pangenome_vst_a2)]
# Obtain background rare species (healthy, cf), 25th abundance percentile, VST, pangenome
f_rare_species_pangenome_vst_a2 <- c(rownames(background_rare_h_pangenome_vst_a2), rownames(background_rare_cf_pangenome_vst_a2))
# remove duplicate entries
f_rare_species_pangenome_vst_a2 <- f_rare_species_pangenome_vst_a2[!duplicated(f_rare_species_pangenome_vst_a2)]

# Obtain background core species (healthy, cf), 35th abundance percentile, VST, pangenome
f_core_species_pangenome_vst_a3 <- c(rownames(background_core_h_pangenome_vst_a3), rownames(background_core_cf_pangenome_vst_a3))
# remove duplicate entries
f_core_species_pangenome_vst_a3 <- f_core_species_pangenome_vst_a3[!duplicated(f_core_species_pangenome_vst_a3)]
# Obtain background rare species (healthy, cf), 35th abundance percentile, VST, pangenome
f_rare_species_pangenome_vst_a3 <- c(rownames(background_rare_h_pangenome_vst_a3), rownames(background_rare_cf_pangenome_vst_a3))
# remove duplicate entries
f_rare_species_pangenome_vst_a3 <- f_rare_species_pangenome_vst_a3[!duplicated(f_rare_species_pangenome_vst_a3)]


# Continue with one-strain per species database
# Obtain background core species (healty, cf), 15th abundance percentile, BCPHC, osps
f_core_species_osps_a1 <- c(rownames(background_core_h_osps_a1), rownames(background_core_cf_osps_a1))
# remove duplicate entries
f_core_species_osps_a1 <- f_core_species_osps_a1[!duplicated(f_core_species_osps_a1)]
# Obtain background rare species (healthy, cf), 15th abundance percentile, BCPHC, osps
f_rare_species_osps_a1 <- c(rownames(background_rare_h_osps_a1), rownames(background_rare_cf_osps_a1))
# remove duplicate entries
f_rare_species_osps_a1 <- f_rare_species_osps_a1[!duplicated(f_rare_species_osps_a1)]

# Obtain background core species (healty, cf), 25th abundance percentile, BCPHC, osps
f_core_species_osps_a2 <- c(rownames(background_core_h_osps_a2), rownames(background_core_cf_osps_a2))
# remove duplicate entries
f_core_species_osps_a2 <- f_core_species_osps_a2[!duplicated(f_core_species_osps_a2)]
# Obtain background rare species (healthy, cf), 25th abundance percentile, BCPHC, osps
f_rare_species_osps_a2 <- c(rownames(background_rare_h_osps_a2), rownames(background_rare_cf_osps_a2))
# remove duplicate entries
f_rare_species_osps_a2 <- f_rare_species_osps_a2[!duplicated(f_rare_species_osps_a2)]

# Obtain background core species (healty, cf), 35th abundance percentile, BCPHC, osps
f_core_species_osps_a3 <- c(rownames(background_core_h_osps_a3), rownames(background_core_cf_osps_a3))
# remove duplicate entries
f_core_species_osps_a3 <- f_core_species_osps_a3[!duplicated(f_core_species_osps_a3)]
# Obtain background rare species (healthy, cf), 35th abundance percentile, BCPHC, osps
f_rare_species_osps_a3 <- c(rownames(background_rare_h_osps_a3), rownames(background_rare_cf_osps_a3))
# remove duplicate entries
f_rare_species_osps_a3 <- f_rare_species_osps_a3[!duplicated(f_rare_species_osps_a3)]


# Obtain background core species (healty, cf), 15th abundance percentile, RLE, osps
f_core_species_osps_rle_a1 <- c(rownames(background_core_h_osps_rle_a1), rownames(background_core_cf_osps_rle_a1))
# remove duplicate entries
f_core_species_osps_rle_a1 <- f_core_species_osps_rle_a1[!duplicated(f_core_species_osps_rle_a1)]
# Obtain background rare species (healthy, cf), 15th abundance percentile, RLE, osps
f_rare_species_osps_rle_a1 <- c(rownames(background_rare_h_osps_rle_a1), rownames(background_rare_cf_osps_rle_a1))
# remove duplicate entries
f_rare_species_osps_rle_a1 <- f_rare_species_osps_rle_a1[!duplicated(f_rare_species_osps_rle_a1)]

# Obtain background core species (healty, cf), 25th abundance percentile, RLE, osps
f_core_species_osps_rle_a2 <- c(rownames(background_core_h_osps_rle_a2), rownames(background_core_cf_osps_rle_a2))
# remove duplicate entries
f_core_species_osps_rle_a2 <- f_core_species_osps_rle_a2[!duplicated(f_core_species_osps_rle_a2)]
# Obtain background rare species (healthy, cf), 25th abundance percentile, RLE, osps
f_rare_species_osps_rle_a2 <- c(rownames(background_rare_h_osps_rle_a2), rownames(background_rare_cf_osps_rle_a2))
# remove duplicate entries
f_rare_species_osps_rle_a2 <- f_rare_species_osps_rle_a2[!duplicated(f_rare_species_osps_rle_a2)]

# Obtain background core species (healty, cf), 35th abundance percentile, RLE, osps
f_core_species_osps_rle_a3 <- c(rownames(background_core_h_osps_rle_a3), rownames(background_core_cf_osps_rle_a3))
# remove duplicate entries
f_core_species_osps_rle_a3 <- f_core_species_osps_rle_a3[!duplicated(f_core_species_osps_rle_a3)]
# Obtain background rare species (healthy, cf), 35th abundance percentile, RLE, osps
f_rare_species_osps_rle_a3 <- c(rownames(background_rare_h_osps_rle_a3), rownames(background_rare_cf_osps_rle_a3))
# remove duplicate entries
f_rare_species_osps_rle_a3 <- f_rare_species_osps_rle_a3[!duplicated(f_rare_species_osps_rle_a3)]



# Obtain background core species (healthy, cf), 15th abundance percentile, VST, osps
f_core_species_osps_vst_a1 <- c(rownames(background_core_h_osps_vst_a1), rownames(background_core_cf_osps_vst_a1))
# remove duplicate entries
f_core_species_osps_vst_a1 <- f_core_species_osps_vst_a1[!duplicated(f_core_species_osps_vst_a1)]
# Obtain background rare species (healthy, cf), 15th abundance percentile, VST, osps
f_rare_species_osps_vst_a1 <- c(rownames(background_rare_h_osps_vst_a1), rownames(background_rare_cf_osps_vst_a1))
# remove duplicate entries
f_rare_species_osps_vst_a1 <- f_rare_species_osps_vst_a1[!duplicated(f_rare_species_osps_vst_a1)]

# Obtain background core species (healthy, cf), 25th abundance percentile, VST, osps
f_core_species_osps_vst_a2 <- c(rownames(background_core_h_osps_vst_a2), rownames(background_core_cf_osps_vst_a2))
# remove duplicate entries
f_core_species_osps_vst_a2 <- f_core_species_osps_vst_a2[!duplicated(f_core_species_osps_vst_a2)]
# Obtain background rare species (healthy, cf), 25th abundance percentile, VST, osps
f_rare_species_osps_vst_a2 <- c(rownames(background_rare_h_osps_vst_a2), rownames(background_rare_cf_osps_vst_a2))
# remove duplicate entries
f_rare_species_osps_vst_a2 <- f_rare_species_osps_vst_a2[!duplicated(f_rare_species_osps_vst_a2)]

# Obtain background core species (healthy, cf), 35th abundance percentile, VST, osps
f_core_species_osps_vst_a3 <- c(rownames(background_core_h_osps_vst_a3), rownames(background_core_cf_osps_vst_a3))
# remove duplicate entries
f_core_species_osps_vst_a3 <- f_core_species_osps_vst_a3[!duplicated(f_core_species_osps_vst_a3)]
# Obtain background rare species (healthy, cf), 35th abundance percentile, VST, osps
f_rare_species_osps_vst_a3 <- c(rownames(background_rare_h_osps_vst_a3), rownames(background_rare_cf_osps_vst_a3))
# remove duplicate entries
f_rare_species_osps_vst_a3 <- f_rare_species_osps_vst_a3[!duplicated(f_rare_species_osps_vst_a3)]


# random forest, pangenome, bcphc, 15th abundance percentile
# define the non-persisting rare species (only found in one or two age groups but not always)
fluctuating_rare_pangenome_a1 <- c(all_rare_a1_pangenome, f_rare_species_pangenome_a1)
# remove duplicate entries
fluctuating_rare_pangenome_a1 <- fluctuating_rare_pangenome_a1[!duplicated(fluctuating_rare_pangenome_a1)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_pangenome_a1 <- c(all_core_a1_pangenome, f_core_species_pangenome_a1)
# remove duplicate entries
fluctuating_core_pangenome_a1 <- fluctuating_core_pangenome_a1[!duplicated(fluctuating_core_pangenome_a1)]

# get all species found in children
uniform_species_pangenome_a1 = c(all_rare_a1_pangenome, all_core_a1_pangenome)
# extract genus information by splitting species name into two parts
uniform_genus_pangenome_a1 = unlist(map(strsplit(uniform_species_pangenome_a1, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_pangenome_a1 = data.frame(cbind(uniform_species_pangenome_a1, uniform_genus_pangenome_a1))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_pangenome_a1$Phylum <- 'X'
uniform_df_pangenome_a1$Class <- 'X'
uniform_df_pangenome_a1$Order <- 'X'
uniform_df_pangenome_a1$Family <- 'X'
# re-name columns
colnames(uniform_df_pangenome_a1) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_pangenome_a1 <- uniform_df_pangenome_a1[!duplicated(uniform_df_pangenome_a1), ]
# re-name row names
rownames(uniform_df_pangenome_a1) <- uniform_df_pangenome_a1$otu
# remove non-numeric information from data frame
uniform_df_pangenome_a1$otu <- NULL
# transpose data frame
uniform_df_pangenome_a1 <- data.frame(t(uniform_df_pangenome_a1))
# add columns with host-associated variables
uniform_df_pangenome_a1$Age <- 'Age'
uniform_df_pangenome_a1$BMI <- 'BMI'
uniform_df_pangenome_a1$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_pangenome_a1 <- data.frame(t(uniform_df_pangenome_a1))
# clean up the row names
rownames(uniform_df_pangenome_a1) <- str_replace(rownames(uniform_df_pangenome_a1), '\\.', ' ')

# subset pangenome database based on abundance estimations (15th abundance percentile)
ds_pangenome_all_a1 <- subset(ds_pangenome, rownames(ds_pangenome) %in% rownames(uniform_df_pangenome_a1))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_pangenome_all_a1 <- data.frame(t(ds_pangenome_all_a1))
# order pangenome data frame by patient id (alphabetically)
ds_pangenome_all_a1 <- ds_pangenome_all_a1[order(rownames(ds_pangenome_all_a1)),]

# add host-associated variables
ds_pangenome_all_a1$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_pangenome_all_a1$Age <- ifelse(ds_pangenome_all_a1$Age_1 == '0 years', 0, ifelse(ds_pangenome_all_a1$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_pangenome_all_a1$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_pangenome_all_a1$BMI <- md_sub$BMI
# remove age character column
ds_pangenome_all_a1$Age_1 <- NULL

# re-arrange data frames
ds_pangenome_all_a1_t <- data.frame(t(ds_pangenome_all_a1))
md_sub_2_pangenome_a1 <- data.frame(t(md_sub))
md_sub_3_pangenome_a1 <- data.frame(t(md_sub_2_pangenome_a1))

# convert data frames to matrix 
otu_mat_pangenome_a1 <- as.matrix(ds_pangenome_all_a1_t)
tax_mat_pangenome_a1 <- as.matrix(uniform_df_pangenome_a1)
# convert pangenome to otu_table (phyloseq object)
OTU_pangenome_a1 = otu_table(otu_mat_pangenome_a1, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_pangenome_a1) <- str_replace(rownames(OTU_pangenome_a1), '\\.', ' ')
# make tax table (phyloseq object)
TAX_pangenome_a1 = tax_table(tax_mat_pangenome_a1)
# make meta data table (phyloseq object)
samples_pangenome_a1 = sample_data(md_sub_3_pangenome_a1)
# generate final phyloseq object
erie.merge_pangenome_a1 <- phyloseq(OTU_pangenome_a1, TAX_pangenome_a1, samples_pangenome_a1)

# transpose data frame
predictors_pangenome_a1 <- t(otu_table(erie.merge_pangenome_a1))
# obtain response variable (CF vs healthy)
response_pangenome_a1 <- as.factor(sample_data(erie.merge_pangenome_a1)$State)
# add response variable to data frame
rf.data_pangenome_a1 <- data.frame(response_pangenome_a1, predictors_pangenome_a1)
# determine mtry for random forest
sample_val_pangenome_a1 = sqrt(ncol(rf.data_pangenome_a1))
# generate three empty variables
rf_pangenome_a1_pred_final = NULL
imp.sort.gini_pangenome_a1_final = NULL
error_rate_pangenome_a1_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_pangenome_a1 <- randomForest(response_pangenome_a1~., data = rf.data_pangenome_a1, 
                                             ntree = 80, 
                                             mtry=sample_val_pangenome_a1, 
                                             importance=TRUE,
                                             na.action = na.roughfix)
  
  # store error rate locally
  error_rate_pangenome_a1 <- erie.classify_pangenome_a1$err.rate
  # make data frame of error rate
  error_rate_pangenome_a1 <- as.data.frame(error_rate_pangenome_a1)
  # add meta data
  error_rate_pangenome_a1$Database <- 'Pangenome'
  error_rate_pangenome_a1$Normalisation <- 'BCPHC'
  error_rate_pangenome_a1$Threshold <- '15th percentile'
  error_rate_pangenome_a1$Seed <- i
  # re-name columns
  colnames(error_rate_pangenome_a1) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_pangenome_a1_final <- rbind(error_rate_pangenome_a1_final, error_rate_pangenome_a1)
  # obtain class probability
  rf_pangenome_a1_pred <- predict(erie.classify_pangenome_a1, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_pangenome_a1_pred <- as.data.frame(rf_pangenome_a1_pred)
  # add meta data
  rf_pangenome_a1_pred$State <- md$State
  rf_pangenome_a1_pred$Database <- 'Pangenome'
  rf_pangenome_a1_pred$Normalisation <- 'BCPHC'
  rf_pangenome_a1_pred$Threshold <- '15th percentile'
  rf_pangenome_a1_pred$Seed <- i
  rownames(rf_pangenome_a1_pred) <- NULL
  # transfer data frame to global environment
  rf_pangenome_a1_pred_final <- rbind(rf_pangenome_a1_pred_final, rf_pangenome_a1_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_pangenome_a1 <- Boruta(response_pangenome_a1~., data = rf.data_pangenome_a1, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_pangenome_a1_df <- as.data.frame(boruta_pangenome_a1$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_pangenome_a1 <- importance(erie.classify_pangenome_a1)
  imp_pangenome_a1 <- data.frame(predictors = rownames(imp_pangenome_a1), imp_pangenome_a1)
  # add meta data
  imp_pangenome_a1$Boruta_name <- rownames(boruta_pangenome_a1_df)
  imp_pangenome_a1$Boruta_predict <- boruta_pangenome_a1_df$`boruta_pangenome_a1$finalDecision`
  imp_pangenome_a1_sub <- subset(imp_pangenome_a1, MeanDecreaseAccuracy > 0.0)
  imp_pangenome_a1_sub$predictors <- str_replace(imp_pangenome_a1_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-assosciated
  imp_pangenome_a1_sub$species_type <- with(imp_pangenome_a1_sub,
                                            ifelse(predictors %in% f_rare_species_pangenome_a1, 
                                                   'Rare species', 
                                                   ifelse(predictors %in% f_core_species_pangenome_a1, 
                                                          'Core species', 
                                                          ifelse(predictors %in% fluctuating_rare_pangenome_a1, 
                                                                 'Rare species',
                                                                 ifelse(predictors %in% fluctuating_core_pangenome_a1, 
                                                                        'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_pangenome_a1_sub$species_type_2 <- with(imp_pangenome_a1_sub,
                                              ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                     ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                            ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                                   ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_pangenome_a1_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                      'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                      'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_pangenome_a1_sub <- subset(imp_pangenome_a1_sub, Boruta_predict != "Rejected")
  imp_pangenome_a1_short <- ddply(imp_pangenome_a1_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_pangenome_a1_short$Seed <- i
  imp_pangenome_a1_short$Threshold <- '15th percentile'
  imp_pangenome_a1_short$Normalisation <- 'BCPHC'
  # store short table globally
  imp.sort.gini_pangenome_a1_final <- rbind(imp.sort.gini_pangenome_a1_final, imp_pangenome_a1_short)
}



# random forest, pangenome, bcphc, 25th abundance percentile
# define the non-persisting rare species (only found in one or two age groups but not always)
fluctuating_rare_pangenome_a2 <- c(all_rare_a2_pangenome, f_rare_species_pangenome_a2)
# remove duplicate entries
fluctuating_rare_pangenome_a2 <- fluctuating_rare_pangenome_a2[!duplicated(fluctuating_rare_pangenome_a2)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_pangenome_a2 <- c(all_core_a2_pangenome, f_core_species_pangenome_a2)
# remove duplicate entries
fluctuating_core_pangenome_a2 <- fluctuating_core_pangenome_a2[!duplicated(fluctuating_core_pangenome_a2)]

# get all species found in children
uniform_species_pangenome_a2 = c(all_rare_a2_pangenome, all_core_a2_pangenome)
# extract genus information by splitting species name into two parts
uniform_genus_pangenome_a2 = unlist(map(strsplit(uniform_species_pangenome_a2, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_pangenome_a2 = data.frame(cbind(uniform_species_pangenome_a2, uniform_genus_pangenome_a2))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_pangenome_a2$Phylum <- 'X'
uniform_df_pangenome_a2$Class <- 'X'
uniform_df_pangenome_a2$Order <- 'X'
uniform_df_pangenome_a2$Family <- 'X'
# re-name columns
colnames(uniform_df_pangenome_a2) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_pangenome_a2 <- uniform_df_pangenome_a2[!duplicated(uniform_df_pangenome_a2), ]
# re-name row names
rownames(uniform_df_pangenome_a2) <- uniform_df_pangenome_a2$otu
# remove non-numeric information from data frame
uniform_df_pangenome_a2$otu <- NULL
# transpose data frame
uniform_df_pangenome_a2 <- data.frame(t(uniform_df_pangenome_a2))
# add columns with host-associated variables
uniform_df_pangenome_a2$Age <- 'Age'
uniform_df_pangenome_a2$BMI <- 'BMI'
uniform_df_pangenome_a2$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_pangenome_a2 <- data.frame(t(uniform_df_pangenome_a2))
# clean up the row names
rownames(uniform_df_pangenome_a2) <- str_replace(rownames(uniform_df_pangenome_a2), '\\.', ' ')

# subset pangenome database based on abundance estimations (25th abundance percentile)
ds_pangenome_all_a2 <- subset(ds_pangenome, rownames(ds_pangenome) %in% rownames(uniform_df_pangenome_a2))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_pangenome_all_a2 <- data.frame(t(ds_pangenome_all_a2))
# order pangenome data frame by patient id (alphabetically)
ds_pangenome_all_a2 <- ds_pangenome_all_a2[order(rownames(ds_pangenome_all_a2)),]

# add host-associated variables
ds_pangenome_all_a2$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_pangenome_all_a2$Age <- ifelse(ds_pangenome_all_a2$Age_1 == '0 years', 0, ifelse(ds_pangenome_all_a2$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_pangenome_all_a2$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_pangenome_all_a2$BMI <- md_sub$BMI
# remove age character column
ds_pangenome_all_a2$Age_1 <- NULL

# re-arrange data frames
ds_pangenome_all_a2_t <- data.frame(t(ds_pangenome_all_a2))
md_sub_2_pangenome_a2 <- data.frame(t(md_sub))
md_sub_3_pangenome_a2 <- data.frame(t(md_sub_2_pangenome_a2))

# convert data frames to matrix 
otu_mat_pangenome_a2 <- as.matrix(ds_pangenome_all_a2_t)
tax_mat_pangenome_a2 <- as.matrix(uniform_df_pangenome_a2)
# convert pangenome to otu_table (phyloseq object)
OTU_pangenome_a2 = otu_table(otu_mat_pangenome_a2, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_pangenome_a2) <- str_replace(rownames(OTU_pangenome_a2), '\\.', ' ')
# make tax table (phyloseq object)
TAX_pangenome_a2 = tax_table(tax_mat_pangenome_a2)
# make meta data table (phyloseq object)
samples_pangenome_a2 = sample_data(md_sub_3_pangenome_a2)
# generate final phyloseq object
erie.merge_pangenome_a2 <- phyloseq(OTU_pangenome_a2, TAX_pangenome_a2, samples_pangenome_a2)

# transpose data frame
predictors_pangenome_a2 <- t(otu_table(erie.merge_pangenome_a2))
# obtain response variable (CF vs healthy)
response_pangenome_a2 <- as.factor(sample_data(erie.merge_pangenome_a2)$State)
# add response variable to data frame
rf.data_pangenome_a2 <- data.frame(response_pangenome_a2, predictors_pangenome_a2)
# determine mtry for random forest
sample_val_pangenome_a2 = sqrt(ncol(rf.data_pangenome_a2))
# generate three empty variables
rf_pangenome_a2_pred_final = NULL
imp.sort.gini_pangenome_a2_final = NULL
error_rate_pangenome_a2_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_pangenome_a2 <- randomForest(response_pangenome_a2~., data = rf.data_pangenome_a2, 
                                             ntree = 80, 
                                             mtry=sample_val_pangenome_a2, 
                                             importance=TRUE,
                                             na.action = na.roughfix)
  
  # store error rate locally
  error_rate_pangenome_a2 <- erie.classify_pangenome_a2$err.rate
  # make data frame of error rate
  error_rate_pangenome_a2 <- as.data.frame(error_rate_pangenome_a2)
  # add meta data
  error_rate_pangenome_a2$Database <- 'Pangenome'
  error_rate_pangenome_a2$Normalisation <- 'BCPHC'
  error_rate_pangenome_a2$Threshold <- '25th percentile'
  error_rate_pangenome_a2$Seed <- i
  # re-name columns
  colnames(error_rate_pangenome_a2) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_pangenome_a2_final <- rbind(error_rate_pangenome_a2_final, error_rate_pangenome_a2)
  # obtain class probability
  rf_pangenome_a2_pred <- predict(erie.classify_pangenome_a2, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_pangenome_a2_pred <- as.data.frame(rf_pangenome_a2_pred)
  # add meta data
  rf_pangenome_a2_pred$State <- md$State
  rf_pangenome_a2_pred$Database <- 'Pangenome'
  rf_pangenome_a2_pred$Normalisation <- 'BCPHC'
  rf_pangenome_a2_pred$Threshold <- '25th percentile'
  rf_pangenome_a2_pred$Seed <- i
  rownames(rf_pangenome_a2_pred) <- NULL
  # transfer data frame to global environment
  rf_pangenome_a2_pred_final <- rbind(rf_pangenome_a2_pred_final, rf_pangenome_a2_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_pangenome_a2 <- Boruta(response_pangenome_a2~., data = rf.data_pangenome_a2, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_pangenome_a2_df <- as.data.frame(boruta_pangenome_a2$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_pangenome_a2 <- importance(erie.classify_pangenome_a2)
  imp_pangenome_a2 <- data.frame(predictors = rownames(imp_pangenome_a2), imp_pangenome_a2)
  # add meta data
  imp_pangenome_a2$Boruta_name <- rownames(boruta_pangenome_a2_df)
  imp_pangenome_a2$Boruta_predict <- boruta_pangenome_a2_df$`boruta_pangenome_a2$finalDecision`
  imp_pangenome_a2_sub <- subset(imp_pangenome_a2, MeanDecreaseAccuracy > 0.0)
  imp_pangenome_a2_sub$predictors <- str_replace(imp_pangenome_a2_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-assosciated
  imp_pangenome_a2_sub$species_type <- with(imp_pangenome_a2_sub,
                                            ifelse(predictors %in% f_rare_species_pangenome_a2, 
                                                   'Rare species', 
                                                   ifelse(predictors %in% f_core_species_pangenome_a2, 
                                                          'Core species', 
                                                          ifelse(predictors %in% fluctuating_rare_pangenome_a2, 
                                                                 'Rare species',
                                                                 ifelse(predictors %in% fluctuating_core_pangenome_a2, 
                                                                        'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_pangenome_a2_sub$species_type_2 <- with(imp_pangenome_a2_sub,
                                              ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                     ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                            ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                                   ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_pangenome_a2_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                      'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                      'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_pangenome_a2_sub <- subset(imp_pangenome_a2_sub, Boruta_predict != "Rejected")
  imp_pangenome_a2_short <- ddply(imp_pangenome_a2_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_pangenome_a2_short$Seed <- i
  imp_pangenome_a2_short$Threshold <- '25th percentile'
  imp_pangenome_a2_short$Normalisation <- 'BCPHC'
  # store short table globally
  imp.sort.gini_pangenome_a2_final <- rbind(imp.sort.gini_pangenome_a2_final, imp_pangenome_a2_short)
}

# random forest, pangenome, bcphc, 35th abundance percentile
# define the non-persisting rare species (only found in one or two age groups but not always)
fluctuating_rare_pangenome_a3 <- c(all_rare_a3_pangenome, f_rare_species_pangenome_a3)
# remove duplicate entries
fluctuating_rare_pangenome_a3 <- fluctuating_rare_pangenome_a3[!duplicated(fluctuating_rare_pangenome_a3)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_pangenome_a3 <- c(all_core_a3_pangenome, f_core_species_pangenome_a3)
# remove duplicate entries
fluctuating_core_pangenome_a3 <- fluctuating_core_pangenome_a3[!duplicated(fluctuating_core_pangenome_a3)]

# get all species found in children
uniform_species_pangenome_a3 = c(all_rare_a3_pangenome, all_core_a3_pangenome)
# extract genus information by splitting species name into two parts
uniform_genus_pangenome_a3 = unlist(map(strsplit(uniform_species_pangenome_a3, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_pangenome_a3 = data.frame(cbind(uniform_species_pangenome_a3, uniform_genus_pangenome_a3))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_pangenome_a3$Phylum <- 'X'
uniform_df_pangenome_a3$Class <- 'X'
uniform_df_pangenome_a3$Order <- 'X'
uniform_df_pangenome_a3$Family <- 'X'
# re-name columns
colnames(uniform_df_pangenome_a3) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_pangenome_a3 <- uniform_df_pangenome_a3[!duplicated(uniform_df_pangenome_a3), ]
# re-name row names
rownames(uniform_df_pangenome_a3) <- uniform_df_pangenome_a3$otu
# remove non-numeric information from data frame
uniform_df_pangenome_a3$otu <- NULL
# transpose data frame
uniform_df_pangenome_a3 <- data.frame(t(uniform_df_pangenome_a3))
# add columns with host-associated variables
uniform_df_pangenome_a3$Age <- 'Age'
uniform_df_pangenome_a3$BMI <- 'BMI'
uniform_df_pangenome_a3$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_pangenome_a3 <- data.frame(t(uniform_df_pangenome_a3))
# clean up the row names
rownames(uniform_df_pangenome_a3) <- str_replace(rownames(uniform_df_pangenome_a3), '\\.', ' ')

# subset pangenome database based on abundance estimations (35th abundance percentile)
ds_pangenome_all_a3 <- subset(ds_pangenome, rownames(ds_pangenome) %in% rownames(uniform_df_pangenome_a3))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_pangenome_all_a3 <- data.frame(t(ds_pangenome_all_a3))
# order pangenome data frame by patient id (alphabetically)
ds_pangenome_all_a3 <- ds_pangenome_all_a3[order(rownames(ds_pangenome_all_a3)),]

# add host-associated variables
ds_pangenome_all_a3$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_pangenome_all_a3$Age <- ifelse(ds_pangenome_all_a3$Age_1 == '0 years', 0, ifelse(ds_pangenome_all_a3$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_pangenome_all_a3$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_pangenome_all_a3$BMI <- md_sub$BMI
# remove age character column
ds_pangenome_all_a3$Age_1 <- NULL

# re-arrange data frames
ds_pangenome_all_a3_t <- data.frame(t(ds_pangenome_all_a3))
md_sub_2_pangenome_a3 <- data.frame(t(md_sub))
md_sub_3_pangenome_a3 <- data.frame(t(md_sub_2_pangenome_a3))

# convert data frames to matrix 
otu_mat_pangenome_a3 <- as.matrix(ds_pangenome_all_a3_t)
tax_mat_pangenome_a3 <- as.matrix(uniform_df_pangenome_a3)
# convert pangenome to otu_table (phyloseq object)
OTU_pangenome_a3 = otu_table(otu_mat_pangenome_a3, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_pangenome_a3) <- str_replace(rownames(OTU_pangenome_a3), '\\.', ' ')
# make tax table (phyloseq object)
TAX_pangenome_a3 = tax_table(tax_mat_pangenome_a3)
# make meta data table (phyloseq object)
samples_pangenome_a3 = sample_data(md_sub_3_pangenome_a3)
# generate final phyloseq object
erie.merge_pangenome_a3 <- phyloseq(OTU_pangenome_a3, TAX_pangenome_a3, samples_pangenome_a3)

# transpose data frame
predictors_pangenome_a3 <- t(otu_table(erie.merge_pangenome_a3))
# obtain response variable (CF vs healthy)
response_pangenome_a3 <- as.factor(sample_data(erie.merge_pangenome_a3)$State)
# add response variable to data frame
rf.data_pangenome_a3 <- data.frame(response_pangenome_a3, predictors_pangenome_a3)
# determine mtry for random forest
sample_val_pangenome_a3 = sqrt(ncol(rf.data_pangenome_a3))
# generate three empty variables
rf_pangenome_a3_pred_final = NULL
imp.sort.gini_pangenome_a3_final = NULL
error_rate_pangenome_a3_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_pangenome_a3 <- randomForest(response_pangenome_a3~., data = rf.data_pangenome_a3, 
                                             ntree = 80, 
                                             mtry=sample_val_pangenome_a3, 
                                             importance=TRUE,
                                             na.action = na.roughfix)
  
  # store error rate locally
  error_rate_pangenome_a3 <- erie.classify_pangenome_a3$err.rate
  # make data frame of error rate
  error_rate_pangenome_a3 <- as.data.frame(error_rate_pangenome_a3)
  # add meta data
  error_rate_pangenome_a3$Database <- 'Pangenome'
  error_rate_pangenome_a3$Normalisation <- 'BCPHC'
  error_rate_pangenome_a3$Threshold <- '35th percentile'
  error_rate_pangenome_a3$Seed <- i
  # re-name columns
  colnames(error_rate_pangenome_a3) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_pangenome_a3_final <- rbind(error_rate_pangenome_a3_final, error_rate_pangenome_a3)
  # obtain class probability
  rf_pangenome_a3_pred <- predict(erie.classify_pangenome_a3, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_pangenome_a3_pred <- as.data.frame(rf_pangenome_a3_pred)
  # add meta data
  rf_pangenome_a3_pred$State <- md$State
  rf_pangenome_a3_pred$Database <- 'Pangenome'
  rf_pangenome_a3_pred$Normalisation <- 'BCPHC'
  rf_pangenome_a3_pred$Threshold <- '35th percentile'
  rf_pangenome_a3_pred$Seed <- i
  rownames(rf_pangenome_a3_pred) <- NULL
  # transfer data frame to global environment
  rf_pangenome_a3_pred_final <- rbind(rf_pangenome_a3_pred_final, rf_pangenome_a3_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_pangenome_a3 <- Boruta(response_pangenome_a3~., data = rf.data_pangenome_a3, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_pangenome_a3_df <- as.data.frame(boruta_pangenome_a3$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_pangenome_a3 <- importance(erie.classify_pangenome_a3)
  imp_pangenome_a3 <- data.frame(predictors = rownames(imp_pangenome_a3), imp_pangenome_a3)
  # add meta data
  imp_pangenome_a3$Boruta_name <- rownames(boruta_pangenome_a3_df)
  imp_pangenome_a3$Boruta_predict <- boruta_pangenome_a3_df$`boruta_pangenome_a3$finalDecision`
  imp_pangenome_a3_sub <- subset(imp_pangenome_a3, MeanDecreaseAccuracy > 0.0)
  imp_pangenome_a3_sub$predictors <- str_replace(imp_pangenome_a3_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-assosciated
  imp_pangenome_a3_sub$species_type <- with(imp_pangenome_a3_sub,
                                            ifelse(predictors %in% f_rare_species_pangenome_a3, 
                                                   'Rare species', 
                                                   ifelse(predictors %in% f_core_species_pangenome_a3, 
                                                          'Core species', 
                                                          ifelse(predictors %in% fluctuating_rare_pangenome_a3, 
                                                                 'Rare species',
                                                                 ifelse(predictors %in% fluctuating_core_pangenome_a3, 
                                                                        'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_pangenome_a3_sub$species_type_2 <- with(imp_pangenome_a3_sub,
                                              ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                     ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                            ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                                   ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_pangenome_a3_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                      'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                      'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_pangenome_a3_sub <- subset(imp_pangenome_a3_sub, Boruta_predict != "Rejected")
  imp_pangenome_a3_short <-ddply(imp_pangenome_a3_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_pangenome_a3_short$Seed <- i
  imp_pangenome_a3_short$Threshold <- '35th percentile'
  imp_pangenome_a3_short$Normalisation <- 'BCPHC'
  # store short table globally
  imp.sort.gini_pangenome_a3_final <- rbind(imp.sort.gini_pangenome_a3_final, imp_pangenome_a3_short)
}


# merge all information, random forest, pangenome, BCPHC, 15th abundance percentile, 25th abundance percentile, 35th abundance percentile
rf.imp.sort.pangenome <- data.frame(rbind(imp.sort.gini_pangenome_a1_final,imp.sort.gini_pangenome_a2_final, imp.sort.gini_pangenome_a3_final))
rf.imp.sort.pangenome$Database <- 'Pangenome'
rf.pred.prob.pangenome <- data.frame(rbind(rf_pangenome_a1_pred_final,rf_pangenome_a2_pred_final, rf_pangenome_a3_pred_final))
rf.error.pangenome <- data.frame(rbind(error_rate_pangenome_a1_final,error_rate_pangenome_a2_final, error_rate_pangenome_a3_final))


# Continue with RLE-normalised data
# random forest for RLE normalised data, 15th abundance percentile
fluctuating_rare_pangenome_rle_a1 <- c(all_rare_a1_pangenome_rle, f_rare_species_pangenome_rle_a1)
# remove duplicate entries
fluctuating_rare_pangenome_rle_a1 <- fluctuating_rare_pangenome_rle_a1[!duplicated(fluctuating_rare_pangenome_rle_a1)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_pangenome_rle_a1 <- c(all_core_a1_pangenome_rle, f_core_species_pangenome_rle_a1)
# remove duplicate entries
fluctuating_core_pangenome_rle_a1 <- fluctuating_core_pangenome_rle_a1[!duplicated(fluctuating_core_pangenome_rle_a1)]

# get all species found in children
uniform_species_pangenome_rle_a1 = c(all_rare_a1_pangenome_rle, all_core_a1_pangenome_rle)
# extract genus information by splitting species name into two parts
uniform_genus_pangenome_rle_a1 = unlist(map(strsplit(uniform_species_pangenome_rle_a1, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_pangenome_rle_a1 = data.frame(cbind(uniform_species_pangenome_rle_a1, uniform_genus_pangenome_rle_a1))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_pangenome_rle_a1$Phylum <- 'X'
uniform_df_pangenome_rle_a1$Class <- 'X'
uniform_df_pangenome_rle_a1$Order <- 'X'
uniform_df_pangenome_rle_a1$Family <- 'X'
# re-name columns
colnames(uniform_df_pangenome_rle_a1) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_pangenome_rle_a1 <- uniform_df_pangenome_rle_a1[!duplicated(uniform_df_pangenome_rle_a1), ]
# re-name row names
rownames(uniform_df_pangenome_rle_a1) <- uniform_df_pangenome_rle_a1$otu
# remove non-numeric information from data frame
uniform_df_pangenome_rle_a1$otu <- NULL
# transpose data frame
uniform_df_pangenome_rle_a1 <- data.frame(t(uniform_df_pangenome_rle_a1))
# add columns with host-associated variables
uniform_df_pangenome_rle_a1$Age <- 'Age'
uniform_df_pangenome_rle_a1$BMI <- 'BMI'
uniform_df_pangenome_rle_a1$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_pangenome_rle_a1 <- data.frame(t(uniform_df_pangenome_rle_a1))
# clean up the row names
rownames(uniform_df_pangenome_rle_a1) <- str_replace(rownames(uniform_df_pangenome_rle_a1), '\\.', ' ')

# subset pangenome database based on abundance estimations (15th abundance percentile)
ds_pangenome_rle_all_a1 <- subset(ds_pangenome_rle, rownames(ds_pangenome_rle) %in% rownames(uniform_df_pangenome_rle_a1))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_pangenome_rle_all_a1 <- data.frame(t(ds_pangenome_rle_all_a1))
# order pangenome data frame by patient id (alphabetically)
ds_pangenome_rle_all_a1 <- ds_pangenome_rle_all_a1[order(rownames(ds_pangenome_rle_all_a1)),]

# add host-associated variables
ds_pangenome_rle_all_a1$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_pangenome_rle_all_a1$Age <- ifelse(ds_pangenome_rle_all_a1$Age_1 == '0 years', 0, ifelse(ds_pangenome_rle_all_a1$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_pangenome_rle_all_a1$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_pangenome_rle_all_a1$BMI <- md_sub$BMI
# remove age character column
ds_pangenome_rle_all_a1$Age_1 <- NULL

# re-arrange data frames
ds_pangenome_rle_all_a1_t <- data.frame(t(ds_pangenome_rle_all_a1))
md_sub_2_pangenome_rle_a1 <- data.frame(t(md_sub))
md_sub_3_pangenome_rle_a1 <- data.frame(t(md_sub_2_pangenome_rle_a1))

# convert data frames to matrix 
otu_mat_pangenome_rle_a1 <- as.matrix(ds_pangenome_rle_all_a1_t)
tax_mat_pangenome_rle_a1 <- as.matrix(uniform_df_pangenome_rle_a1)
# convert pangenome to otu_table (phyloseq object)
OTU_pangenome_rle_a1 = otu_table(otu_mat_pangenome_rle_a1, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_pangenome_rle_a1) <- str_replace(rownames(OTU_pangenome_rle_a1), '\\.', ' ')
# make tax table (phyloseq object)
TAX_pangenome_rle_a1 = tax_table(tax_mat_pangenome_rle_a1)
# make meta data table (phyloseq object)
samples_pangenome_rle_a1 = sample_data(md_sub_3_pangenome_rle_a1)
# generate final phyloseq object
erie.merge_pangenome_rle_a1 <- phyloseq(OTU_pangenome_rle_a1, TAX_pangenome_rle_a1, samples_pangenome_rle_a1)

# transpose data frame
predictors_pangenome_rle_a1 <- t(otu_table(erie.merge_pangenome_rle_a1))
# obtain response variable (CF vs healthy)
response_pangenome_rle_a1 <- as.factor(sample_data(erie.merge_pangenome_rle_a1)$State)
# add response variable to data frame
rf.data_pangenome_rle_a1 <- data.frame(response_pangenome_rle_a1, predictors_pangenome_rle_a1)
# determine mtry for random forest
sample_val_pangenome_rle_a1 = sqrt(ncol(rf.data_pangenome_rle_a1))
# generate three empty variables
rf_pangenome_rle_a1_pred_final = NULL
imp.sort.gini_pangenome_rle_a1_final = NULL
error_rate_pangenome_rle_a1_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_pangenome_rle_a1 <- randomForest(response_pangenome_rle_a1~., data = rf.data_pangenome_rle_a1, 
                                                 ntree = 80, 
                                                 mtry=sample_val_pangenome_rle_a1, 
                                                 importance=TRUE,
                                                 na.action = na.roughfix)
  
  # store error rate locally
  error_rate_pangenome_rle_a1 <- erie.classify_pangenome_rle_a1$err.rate
  # make data frame of error rate
  error_rate_pangenome_rle_a1 <- as.data.frame(error_rate_pangenome_rle_a1)
  # add meta data
  error_rate_pangenome_rle_a1$Database <- 'Pangenome'
  error_rate_pangenome_rle_a1$Normalisation <- 'RLE'
  error_rate_pangenome_rle_a1$Threshold <- '15th percentile'
  error_rate_pangenome_rle_a1$Seed <- i
  # re-name columns
  colnames(error_rate_pangenome_rle_a1) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_pangenome_rle_a1_final <- rbind(error_rate_pangenome_rle_a1_final, error_rate_pangenome_rle_a1)
  # obtain class probability
  rf_pangenome_rle_a1_pred <- predict(erie.classify_pangenome_rle_a1, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_pangenome_rle_a1_pred <- as.data.frame(rf_pangenome_rle_a1_pred)
  # add meta data
  rf_pangenome_rle_a1_pred$State <- md$State
  rf_pangenome_rle_a1_pred$Database <- 'Pangenome'
  rf_pangenome_rle_a1_pred$Normalisation <- 'RLE'
  rf_pangenome_rle_a1_pred$Threshold <- '15th percentile'
  rf_pangenome_rle_a1_pred$Seed <- i
  rownames(rf_pangenome_rle_a1_pred) <- NULL
  # transfer data frame to global environment
  rf_pangenome_rle_a1_pred_final <- rbind(rf_pangenome_rle_a1_pred_final, rf_pangenome_rle_a1_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_pangenome_rle_a1 <- Boruta(response_pangenome_rle_a1~., data = rf.data_pangenome_rle_a1, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_pangenome_rle_a1_df <- as.data.frame(boruta_pangenome_rle_a1$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_pangenome_rle_a1 <- importance(erie.classify_pangenome_rle_a1)
  imp_pangenome_rle_a1 <- data.frame(predictors = rownames(imp_pangenome_rle_a1), imp_pangenome_rle_a1)
  # add meta data
  imp_pangenome_rle_a1$Boruta_name <- rownames(boruta_pangenome_rle_a1_df)
  imp_pangenome_rle_a1$Boruta_predict <- boruta_pangenome_rle_a1_df$`boruta_pangenome_rle_a1$finalDecision`
  imp_pangenome_rle_a1_sub <- subset(imp_pangenome_rle_a1, MeanDecreaseAccuracy > 0.0)
  imp_pangenome_rle_a1_sub$predictors <- str_replace(imp_pangenome_rle_a1_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-associated
  imp_pangenome_rle_a1_sub$species_type <- with(imp_pangenome_rle_a1_sub,
                                                ifelse(predictors %in% f_rare_species_pangenome_rle_a1, 
                                                       'Rare species', 
                                                       ifelse(predictors %in% f_core_species_pangenome_rle_a1, 
                                                              'Core species', 
                                                              ifelse(predictors %in% fluctuating_rare_pangenome_rle_a1, 
                                                                     'Rare species',
                                                                     ifelse(predictors %in% fluctuating_core_pangenome_rle_a1, 
                                                                            'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_pangenome_rle_a1_sub$species_type_2 <- with(imp_pangenome_rle_a1_sub,
                                                  ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                         ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                                ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                                       ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_pangenome_rle_a1_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                          'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                          'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_pangenome_rle_a1_sub <- subset(imp_pangenome_rle_a1_sub, Boruta_predict != "Rejected")
  imp_pangenome_rle_a1_short <- ddply(imp_pangenome_rle_a1_sub, 'Species_type_2', numcolwise(sum))
  
  # add meta data
  imp_pangenome_rle_a1_short$Seed <- i
  imp_pangenome_rle_a1_short$Threshold <- '15th percentile'
  imp_pangenome_rle_a1_short$Normalisation <- 'RLE'
  # store short table globally
  imp.sort.gini_pangenome_rle_a1_final <- rbind(imp.sort.gini_pangenome_rle_a1_final, imp_pangenome_rle_a1_short)
}


# random forest, pangenome, RLE, 25th abundance percentile
# define the non-persisting rare species (only found in one or two age groups but not always)
fluctuating_rare_pangenome_rle_a2 <- c(all_rare_a2_pangenome_rle, f_rare_species_pangenome_rle_a2)
# remove duplicate entries
fluctuating_rare_pangenome_rle_a2 <- fluctuating_rare_pangenome_rle_a2[!duplicated(fluctuating_rare_pangenome_rle_a2)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_pangenome_rle_a2 <- c(all_core_a2_pangenome_rle, f_core_species_pangenome_rle_a2)
# remove duplicate entries
fluctuating_core_pangenome_rle_a2 <- fluctuating_core_pangenome_rle_a2[!duplicated(fluctuating_core_pangenome_rle_a2)]

# get all species found in children
uniform_species_pangenome_rle_a2 = c(all_rare_a2_pangenome_rle, all_core_a2_pangenome_rle)
# extract genus information by splitting species name into two parts
uniform_genus_pangenome_rle_a2 = unlist(map(strsplit(uniform_species_pangenome_rle_a2, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_pangenome_rle_a2 = data.frame(cbind(uniform_species_pangenome_rle_a2, uniform_genus_pangenome_rle_a2))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_pangenome_rle_a2$Phylum <- 'X'
uniform_df_pangenome_rle_a2$Class <- 'X'
uniform_df_pangenome_rle_a2$Order <- 'X'
uniform_df_pangenome_rle_a2$Family <- 'X'
# re-name columns
colnames(uniform_df_pangenome_rle_a2) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_pangenome_rle_a2 <- uniform_df_pangenome_rle_a2[!duplicated(uniform_df_pangenome_rle_a2), ]
# re-name row names
rownames(uniform_df_pangenome_rle_a2) <- uniform_df_pangenome_rle_a2$otu
# remove non-numeric information from data frame
uniform_df_pangenome_rle_a2$otu <- NULL
# transpose data frame
uniform_df_pangenome_rle_a2 <- data.frame(t(uniform_df_pangenome_rle_a2))
# add columns with host-associated variables
uniform_df_pangenome_rle_a2$Age <- 'Age'
uniform_df_pangenome_rle_a2$BMI <- 'BMI'
uniform_df_pangenome_rle_a2$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_pangenome_rle_a2 <- data.frame(t(uniform_df_pangenome_rle_a2))
# clean up the row names
rownames(uniform_df_pangenome_rle_a2) <- str_replace(rownames(uniform_df_pangenome_rle_a2), '\\.', ' ')

# subset pangenome database based on abundance estimations (25th abundance percentile)
ds_pangenome_rle_all_a2 <- subset(ds_pangenome_rle, rownames(ds_pangenome_rle) %in% rownames(uniform_df_pangenome_rle_a2))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_pangenome_rle_all_a2 <- data.frame(t(ds_pangenome_rle_all_a2))
# order pangenome data frame by patient id (alphabetically)
ds_pangenome_rle_all_a2 <- ds_pangenome_rle_all_a2[order(rownames(ds_pangenome_rle_all_a2)),]

# add host-associated variables
ds_pangenome_rle_all_a2$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_pangenome_rle_all_a2$Age <- ifelse(ds_pangenome_rle_all_a2$Age_1 == '0 years', 0, ifelse(ds_pangenome_rle_all_a2$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_pangenome_rle_all_a2$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_pangenome_rle_all_a2$BMI <- md_sub$BMI
# remove age character column
ds_pangenome_rle_all_a2$Age_1 <- NULL

# re-arrange data frames
ds_pangenome_rle_all_a2_t <- data.frame(t(ds_pangenome_rle_all_a2))
md_sub_2_pangenome_rle_a2 <- data.frame(t(md_sub))
md_sub_3_pangenome_rle_a2 <- data.frame(t(md_sub_2_pangenome_rle_a2))

# convert data frames to matrix 
otu_mat_pangenome_rle_a2 <- as.matrix(ds_pangenome_rle_all_a2_t)
tax_mat_pangenome_rle_a2 <- as.matrix(uniform_df_pangenome_rle_a2)
# convert pangenome_rle to otu_table (phyloseq object)
OTU_pangenome_rle_a2 = otu_table(otu_mat_pangenome_rle_a2, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_pangenome_rle_a2) <- str_replace(rownames(OTU_pangenome_rle_a2), '\\.', ' ')
# make tax table (phyloseq object)
TAX_pangenome_rle_a2 = tax_table(tax_mat_pangenome_rle_a2)
# make meta data table (phyloseq object)
samples_pangenome_rle_a2 = sample_data(md_sub_3_pangenome_rle_a2)
# generate final phyloseq object
erie.merge_pangenome_rle_a2 <- phyloseq(OTU_pangenome_rle_a2, TAX_pangenome_rle_a2, samples_pangenome_rle_a2)

# transpose data frame
predictors_pangenome_rle_a2 <- t(otu_table(erie.merge_pangenome_rle_a2))
# obtain response variable (CF vs healthy)
response_pangenome_rle_a2 <- as.factor(sample_data(erie.merge_pangenome_rle_a2)$State)
# add response variable to data frame
rf.data_pangenome_rle_a2 <- data.frame(response_pangenome_rle_a2, predictors_pangenome_rle_a2)
# determine mtry for random forest
sample_val_pangenome_rle_a2 = sqrt(ncol(rf.data_pangenome_rle_a2))
# generate three empty variables
rf_pangenome_rle_a2_pred_final = NULL
imp.sort.gini_pangenome_rle_a2_final = NULL
error_rate_pangenome_rle_a2_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_pangenome_rle_a2 <- randomForest(response_pangenome_rle_a2~., data = rf.data_pangenome_rle_a2, 
                                                 ntree = 80, 
                                                 mtry=sample_val_pangenome_rle_a2, 
                                                 importance=TRUE,
                                                 na.action = na.roughfix)
  
  # store error rate locally
  error_rate_pangenome_rle_a2 <- erie.classify_pangenome_rle_a2$err.rate
  # make data frame of error rate
  error_rate_pangenome_rle_a2 <- as.data.frame(error_rate_pangenome_rle_a2)
  # add meta data
  error_rate_pangenome_rle_a2$Database <- 'Pangenome'
  error_rate_pangenome_rle_a2$Normalisation <- 'RLE'
  error_rate_pangenome_rle_a2$Threshold <- '25th percentile'
  error_rate_pangenome_rle_a2$Seed <- i
  # re-name columns
  colnames(error_rate_pangenome_rle_a2) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_pangenome_rle_a2_final <- rbind(error_rate_pangenome_rle_a2_final, error_rate_pangenome_rle_a2)
  # obtain class probability
  rf_pangenome_rle_a2_pred <- predict(erie.classify_pangenome_rle_a2, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_pangenome_rle_a2_pred <- as.data.frame(rf_pangenome_rle_a2_pred)
  # add meta data
  rf_pangenome_rle_a2_pred$State <- md$State
  rf_pangenome_rle_a2_pred$Database <- 'Pangenome'
  rf_pangenome_rle_a2_pred$Normalisation <- 'RLE'
  rf_pangenome_rle_a2_pred$Threshold <- '25th percentile'
  rf_pangenome_rle_a2_pred$Seed <- i
  rownames(rf_pangenome_rle_a2_pred) <- NULL
  # transfer data frame to global environment
  rf_pangenome_rle_a2_pred_final <- rbind(rf_pangenome_rle_a2_pred_final, rf_pangenome_rle_a2_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_pangenome_rle_a2 <- Boruta(response_pangenome_rle_a2~., data = rf.data_pangenome_rle_a2, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_pangenome_rle_a2_df <- as.data.frame(boruta_pangenome_rle_a2$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_pangenome_rle_a2 <- importance(erie.classify_pangenome_rle_a2)
  imp_pangenome_rle_a2 <- data.frame(predictors = rownames(imp_pangenome_rle_a2), imp_pangenome_rle_a2)
  # add meta data
  imp_pangenome_rle_a2$Boruta_name <- rownames(boruta_pangenome_rle_a2_df)
  imp_pangenome_rle_a2$Boruta_predict <- boruta_pangenome_rle_a2_df$`boruta_pangenome_rle_a2$finalDecision`
  imp_pangenome_rle_a2_sub <- subset(imp_pangenome_rle_a2, MeanDecreaseAccuracy > 0.0)
  imp_pangenome_rle_a2_sub$predictors <- str_replace(imp_pangenome_rle_a2_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-assosciated
  imp_pangenome_rle_a2_sub$species_type <- with(imp_pangenome_rle_a2_sub,
                                                ifelse(predictors %in% f_rare_species_pangenome_rle_a2, 
                                                       'Rare species', 
                                                       ifelse(predictors %in% f_core_species_pangenome_rle_a2, 
                                                              'Core species', 
                                                              ifelse(predictors %in% fluctuating_rare_pangenome_rle_a2, 
                                                                     'Rare species',
                                                                     ifelse(predictors %in% fluctuating_core_pangenome_rle_a2, 
                                                                            'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_pangenome_rle_a2_sub$species_type_2 <- with(imp_pangenome_rle_a2_sub,
                                                  ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                         ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                                ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                                       ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_pangenome_rle_a2_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                          'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                          'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_pangenome_rle_a2_sub <- subset(imp_pangenome_rle_a2_sub, Boruta_predict != "Rejected")
  imp_pangenome_rle_a2_short <- ddply(imp_pangenome_rle_a2_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_pangenome_rle_a2_short$Seed <- i
  imp_pangenome_rle_a2_short$Threshold <- '25th percentile'
  imp_pangenome_rle_a2_short$Normalisation <- 'RLE'
  # store short table globally
  imp.sort.gini_pangenome_rle_a2_final <- rbind(imp.sort.gini_pangenome_rle_a2_final, imp_pangenome_rle_a2_short)
}


# random forest, pangenome, RLE, 35th abundance percentile
# define the non-persisting rare species (only found in one or two age groups but not always)
fluctuating_rare_pangenome_rle_a3 <- c(all_rare_a3_pangenome_rle, f_rare_species_pangenome_rle_a3)
# remove duplicate entries
fluctuating_rare_pangenome_rle_a3 <- fluctuating_rare_pangenome_rle_a3[!duplicated(fluctuating_rare_pangenome_rle_a3)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_pangenome_rle_a3 <- c(all_core_a3_pangenome_rle, f_core_species_pangenome_rle_a3)
# remove duplicate entries
fluctuating_core_pangenome_rle_a3 <- fluctuating_core_pangenome_rle_a3[!duplicated(fluctuating_core_pangenome_rle_a3)]

# get all species found in children
uniform_species_pangenome_rle_a3 = c(all_rare_a3_pangenome_rle, all_core_a3_pangenome_rle)
# extract genus information by splitting species name into two parts
uniform_genus_pangenome_rle_a3 = unlist(map(strsplit(uniform_species_pangenome_rle_a3, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_pangenome_rle_a3 = data.frame(cbind(uniform_species_pangenome_rle_a3, uniform_genus_pangenome_rle_a3))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_pangenome_rle_a3$Phylum <- 'X'
uniform_df_pangenome_rle_a3$Class <- 'X'
uniform_df_pangenome_rle_a3$Order <- 'X'
uniform_df_pangenome_rle_a3$Family <- 'X'
# re-name columns
colnames(uniform_df_pangenome_rle_a3) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_pangenome_rle_a3 <- uniform_df_pangenome_rle_a3[!duplicated(uniform_df_pangenome_rle_a3), ]
# re-name row names
rownames(uniform_df_pangenome_rle_a3) <- uniform_df_pangenome_rle_a3$otu
# remove non-numeric information from data frame
uniform_df_pangenome_rle_a3$otu <- NULL
# transpose data frame
uniform_df_pangenome_rle_a3 <- data.frame(t(uniform_df_pangenome_rle_a3))
# add columns with host-associated variables
uniform_df_pangenome_rle_a3$Age <- 'Age'
uniform_df_pangenome_rle_a3$BMI <- 'BMI'
uniform_df_pangenome_rle_a3$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_pangenome_rle_a3 <- data.frame(t(uniform_df_pangenome_rle_a3))
# clean up the row names
rownames(uniform_df_pangenome_rle_a3) <- str_replace(rownames(uniform_df_pangenome_rle_a3), '\\.', ' ')

# subset pangenome database based on abundance estimations (35th abundance percentile)
ds_pangenome_rle_all_a3 <- subset(ds_pangenome_rle, rownames(ds_pangenome_rle) %in% rownames(uniform_df_pangenome_rle_a3))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_pangenome_rle_all_a3 <- data.frame(t(ds_pangenome_rle_all_a3))
# order pangenome data frame by patient id (alphabetically)
ds_pangenome_rle_all_a3 <- ds_pangenome_rle_all_a3[order(rownames(ds_pangenome_rle_all_a3)),]

# add host-associated variables
ds_pangenome_rle_all_a3$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_pangenome_rle_all_a3$Age <- ifelse(ds_pangenome_rle_all_a3$Age_1 == '0 years', 0, ifelse(ds_pangenome_rle_all_a3$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_pangenome_rle_all_a3$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_pangenome_rle_all_a3$BMI <- md_sub$BMI
# remove age character column
ds_pangenome_rle_all_a3$Age_1 <- NULL

# re-arrange data frames
ds_pangenome_rle_all_a3_t <- data.frame(t(ds_pangenome_rle_all_a3))
md_sub_2_pangenome_rle_a3 <- data.frame(t(md_sub))
md_sub_3_pangenome_rle_a3 <- data.frame(t(md_sub_2_pangenome_rle_a3))

# convert data frames to matrix 
otu_mat_pangenome_rle_a3 <- as.matrix(ds_pangenome_rle_all_a3_t)
tax_mat_pangenome_rle_a3 <- as.matrix(uniform_df_pangenome_rle_a3)
# convert pangenome to otu_table (phyloseq object)
OTU_pangenome_rle_a3 = otu_table(otu_mat_pangenome_rle_a3, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_pangenome_rle_a3) <- str_replace(rownames(OTU_pangenome_rle_a3), '\\.', ' ')
# make tax table (phyloseq object)
TAX_pangenome_rle_a3 = tax_table(tax_mat_pangenome_rle_a3)
# make meta data table (phyloseq object)
samples_pangenome_rle_a3 = sample_data(md_sub_3_pangenome_rle_a3)
# generate final phyloseq object
erie.merge_pangenome_rle_a3 <- phyloseq(OTU_pangenome_rle_a3, TAX_pangenome_rle_a3, samples_pangenome_rle_a3)

# transpose data frame
predictors_pangenome_rle_a3 <- t(otu_table(erie.merge_pangenome_rle_a3))
# obtain response variable (CF vs healthy)
response_pangenome_rle_a3 <- as.factor(sample_data(erie.merge_pangenome_rle_a3)$State)
# add response variable to data frame
rf.data_pangenome_rle_a3 <- data.frame(response_pangenome_rle_a3, predictors_pangenome_rle_a3)
# determine mtry for random forest
sample_val_pangenome_rle_a3 = sqrt(ncol(rf.data_pangenome_rle_a3))
# generate three empty variables
rf_pangenome_rle_a3_pred_final = NULL
imp.sort.gini_pangenome_rle_a3_final = NULL
error_rate_pangenome_rle_a3_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_pangenome_rle_a3 <- randomForest(response_pangenome_rle_a3~., data = rf.data_pangenome_rle_a3, 
                                                 ntree = 80, 
                                                 mtry=sample_val_pangenome_rle_a3, 
                                                 importance=TRUE,
                                                 na.action = na.roughfix)
  
  # store error rate locally
  error_rate_pangenome_rle_a3 <- erie.classify_pangenome_rle_a3$err.rate
  # make data frame of error rate
  error_rate_pangenome_rle_a3 <- as.data.frame(error_rate_pangenome_rle_a3)
  # add meta data
  error_rate_pangenome_rle_a3$Database <- 'Pangenome'
  error_rate_pangenome_rle_a3$Normalisation <- 'RLE'
  error_rate_pangenome_rle_a3$Threshold <- '35th percentile'
  error_rate_pangenome_rle_a3$Seed <- i
  # re-name columns
  colnames(error_rate_pangenome_rle_a3) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_pangenome_rle_a3_final <- rbind(error_rate_pangenome_rle_a3_final, error_rate_pangenome_rle_a3)
  # obtain class probability
  rf_pangenome_rle_a3_pred <- predict(erie.classify_pangenome_rle_a3, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_pangenome_rle_a3_pred <- as.data.frame(rf_pangenome_rle_a3_pred)
  # add meta data
  rf_pangenome_rle_a3_pred$State <- md$State
  rf_pangenome_rle_a3_pred$Database <- 'Pangenome'
  rf_pangenome_rle_a3_pred$Normalisation <- 'RLE'
  rf_pangenome_rle_a3_pred$Threshold <- '35th percentile'
  rf_pangenome_rle_a3_pred$Seed <- i
  rownames(rf_pangenome_rle_a3_pred) <- NULL
  # transfer data frame to global environment
  rf_pangenome_rle_a3_pred_final <- rbind(rf_pangenome_rle_a3_pred_final, rf_pangenome_rle_a3_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_pangenome_rle_a3 <- Boruta(response_pangenome_rle_a3~., data = rf.data_pangenome_rle_a3, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_pangenome_rle_a3_df <- as.data.frame(boruta_pangenome_rle_a3$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_pangenome_rle_a3 <- importance(erie.classify_pangenome_rle_a3)
  imp_pangenome_rle_a3 <- data.frame(predictors = rownames(imp_pangenome_rle_a3), imp_pangenome_rle_a3)
  # add meta data
  imp_pangenome_rle_a3$Boruta_name <- rownames(boruta_pangenome_rle_a3_df)
  imp_pangenome_rle_a3$Boruta_predict <- boruta_pangenome_rle_a3_df$`boruta_pangenome_rle_a3$finalDecision`
  imp_pangenome_rle_a3_sub <- subset(imp_pangenome_rle_a3, MeanDecreaseAccuracy > 0.0)
  imp_pangenome_rle_a3_sub$predictors <- str_replace(imp_pangenome_rle_a3_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-assosciated
  imp_pangenome_rle_a3_sub$species_type <- with(imp_pangenome_rle_a3_sub,
                                                ifelse(predictors %in% f_rare_species_pangenome_rle_a3, 
                                                       'Rare species', 
                                                       ifelse(predictors %in% f_core_species_pangenome_rle_a3, 
                                                              'Core species', 
                                                              ifelse(predictors %in% fluctuating_rare_pangenome_rle_a3, 
                                                                     'Rare species',
                                                                     ifelse(predictors %in% fluctuating_core_pangenome_rle_a3, 
                                                                            'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_pangenome_rle_a3_sub$species_type_2 <- with(imp_pangenome_rle_a3_sub,
                                                  ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                         ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                                ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                                       ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_pangenome_rle_a3_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                          'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                          'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_pangenome_rle_a3_sub <- subset(imp_pangenome_rle_a3_sub, Boruta_predict != "Rejected")
  imp_pangenome_rle_a3_short <- ddply(imp_pangenome_rle_a3_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_pangenome_rle_a3_short$Seed <- i
  imp_pangenome_rle_a3_short$Threshold <- '35th percentile'
  imp_pangenome_rle_a3_short$Normalisation <- 'RLE'
  # store short table globally
  imp.sort.gini_pangenome_rle_a3_final <- rbind(imp.sort.gini_pangenome_rle_a3_final, imp_pangenome_rle_a3_short)
}


# merge all information, random forest, pangenome, RLE, 15th abundance percentile, 25th abundance percentile, 35th abundance percentile
rf.imp.sort.pangenome_rle <- data.frame(rbind(imp.sort.gini_pangenome_rle_a1_final,imp.sort.gini_pangenome_rle_a2_final, imp.sort.gini_pangenome_rle_a3_final))
rf.imp.sort.pangenome_rle$Database <- 'Pangenome'
rf.pred.prob.pangenome_rle <- data.frame(rbind(rf_pangenome_rle_a1_pred_final,rf_pangenome_rle_a2_pred_final, rf_pangenome_rle_a3_pred_final))
rf.error.pangenome_rle <- data.frame(rbind(error_rate_pangenome_rle_a1_final,error_rate_pangenome_rle_a2_final, error_rate_pangenome_rle_a3_final))


# Continue with VST-normalised data
# random forest for VST normalised data, 15th abundance percentile
fluctuating_rare_pangenome_vst_a1 <- c(all_rare_a1_pangenome_vst, f_rare_species_pangenome_vst_a1)
# remove duplicate entries
fluctuating_rare_pangenome_vst_a1 <- fluctuating_rare_pangenome_vst_a1[!duplicated(fluctuating_rare_pangenome_vst_a1)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_pangenome_vst_a1 <- c(all_core_a1_pangenome_vst, f_core_species_pangenome_vst_a1)
# remove duplicate entries
fluctuating_core_pangenome_vst_a1 <- fluctuating_core_pangenome_vst_a1[!duplicated(fluctuating_core_pangenome_vst_a1)]

# get all species found in children
uniform_species_pangenome_vst_a1 = c(all_rare_a1_pangenome_vst, all_core_a1_pangenome_vst)
# extract genus information by splitting species name into two parts
uniform_genus_pangenome_vst_a1 = unlist(map(strsplit(uniform_species_pangenome_vst_a1, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_pangenome_vst_a1 = data.frame(cbind(uniform_species_pangenome_vst_a1, uniform_genus_pangenome_vst_a1))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_pangenome_vst_a1$Phylum <- 'X'
uniform_df_pangenome_vst_a1$Class <- 'X'
uniform_df_pangenome_vst_a1$Order <- 'X'
uniform_df_pangenome_vst_a1$Family <- 'X'
# re-name columns
colnames(uniform_df_pangenome_vst_a1) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_pangenome_vst_a1 <- uniform_df_pangenome_vst_a1[!duplicated(uniform_df_pangenome_vst_a1), ]
# re-name row names
rownames(uniform_df_pangenome_vst_a1) <- uniform_df_pangenome_vst_a1$otu
# remove non-numeric information from data frame
uniform_df_pangenome_vst_a1$otu <- NULL
# transpose data frame
uniform_df_pangenome_vst_a1 <- data.frame(t(uniform_df_pangenome_vst_a1))
# add columns with host-associated variables
uniform_df_pangenome_vst_a1$Age <- 'Age'
uniform_df_pangenome_vst_a1$BMI <- 'BMI'
uniform_df_pangenome_vst_a1$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_pangenome_vst_a1 <- data.frame(t(uniform_df_pangenome_vst_a1))
# clean up the row names
rownames(uniform_df_pangenome_vst_a1) <- str_replace(rownames(uniform_df_pangenome_vst_a1), '\\.', ' ')

# subset pangenome database based on abundance estimations (15th abundance percentile)
ds_pangenome_vst_all_a1 <- subset(ds_pangenome_vst, rownames(ds_pangenome_vst) %in% rownames(uniform_df_pangenome_vst_a1))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_pangenome_vst_all_a1 <- data.frame(t(ds_pangenome_vst_all_a1))
# order pangenome data frame by patient id (alphabetically)
ds_pangenome_vst_all_a1 <- ds_pangenome_vst_all_a1[order(rownames(ds_pangenome_vst_all_a1)),]

# add host-associated variables
ds_pangenome_vst_all_a1$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_pangenome_vst_all_a1$Age <- ifelse(ds_pangenome_vst_all_a1$Age_1 == '0 years', 0, ifelse(ds_pangenome_vst_all_a1$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_pangenome_vst_all_a1$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_pangenome_vst_all_a1$BMI <- md_sub$BMI
# remove age character column
ds_pangenome_vst_all_a1$Age_1 <- NULL

# re-arrange data frames
ds_pangenome_vst_all_a1_t <- data.frame(t(ds_pangenome_vst_all_a1))
md_sub_2_pangenome_vst_a1 <- data.frame(t(md_sub))
md_sub_3_pangenome_vst_a1 <- data.frame(t(md_sub_2_pangenome_vst_a1))

# convert data frames to matrix 
otu_mat_pangenome_vst_a1 <- as.matrix(ds_pangenome_vst_all_a1_t)
tax_mat_pangenome_vst_a1 <- as.matrix(uniform_df_pangenome_vst_a1)
# convert pangenome to otu_table (phyloseq object)
OTU_pangenome_vst_a1 = otu_table(otu_mat_pangenome_vst_a1, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_pangenome_vst_a1) <- str_replace(rownames(OTU_pangenome_vst_a1), '\\.', ' ')
# make tax table (phyloseq object)
TAX_pangenome_vst_a1 = tax_table(tax_mat_pangenome_vst_a1)
# make meta data table (phyloseq object)
samples_pangenome_vst_a1 = sample_data(md_sub_3_pangenome_vst_a1)
# generate final phyloseq object
erie.merge_pangenome_vst_a1 <- phyloseq(OTU_pangenome_vst_a1, TAX_pangenome_vst_a1, samples_pangenome_vst_a1)

# transpose data frame
predictors_pangenome_vst_a1 <- t(otu_table(erie.merge_pangenome_vst_a1))
# obtain response variable (CF vs healthy)
response_pangenome_vst_a1 <- as.factor(sample_data(erie.merge_pangenome_vst_a1)$State)
# add response variable to data frame
rf.data_pangenome_vst_a1 <- data.frame(response_pangenome_vst_a1, predictors_pangenome_vst_a1)
# determine mtry for random forest
sample_val_pangenome_vst_a1 = sqrt(ncol(rf.data_pangenome_vst_a1))
# generate three empty variables
rf_pangenome_vst_a1_pred_final = NULL
imp.sort.gini_pangenome_vst_a1_final = NULL
error_rate_pangenome_vst_a1_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_pangenome_vst_a1 <- randomForest(response_pangenome_vst_a1~., data = rf.data_pangenome_vst_a1, 
                                                 ntree = 80, 
                                                 mtry=sample_val_pangenome_vst_a1, 
                                                 importance=TRUE,
                                                 na.action = na.roughfix)
  
  # store error rate locally
  error_rate_pangenome_vst_a1 <- erie.classify_pangenome_vst_a1$err.rate
  # make data frame of error rate
  error_rate_pangenome_vst_a1 <- as.data.frame(error_rate_pangenome_vst_a1)
  # add meta data
  error_rate_pangenome_vst_a1$Database <- 'Pangenome'
  error_rate_pangenome_vst_a1$Normalisation <- 'VST'
  error_rate_pangenome_vst_a1$Threshold <- '15th percentile'
  error_rate_pangenome_vst_a1$Seed <- i
  # re-name columns
  colnames(error_rate_pangenome_vst_a1) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_pangenome_vst_a1_final <- rbind(error_rate_pangenome_vst_a1_final, error_rate_pangenome_vst_a1)
  # obtain class probability
  rf_pangenome_vst_a1_pred <- predict(erie.classify_pangenome_vst_a1, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_pangenome_vst_a1_pred <- as.data.frame(rf_pangenome_vst_a1_pred)
  # add meta data
  rf_pangenome_vst_a1_pred$State <- md$State
  rf_pangenome_vst_a1_pred$Database <- 'Pangenome'
  rf_pangenome_vst_a1_pred$Normalisation <- 'VST'
  rf_pangenome_vst_a1_pred$Threshold <- '15th percentile'
  rf_pangenome_vst_a1_pred$Seed <- i
  rownames(rf_pangenome_vst_a1_pred) <- NULL
  # transfer data frame to global environment
  rf_pangenome_vst_a1_pred_final <- rbind(rf_pangenome_vst_a1_pred_final, rf_pangenome_vst_a1_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_pangenome_vst_a1 <- Boruta(response_pangenome_vst_a1~., data = rf.data_pangenome_vst_a1, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_pangenome_vst_a1_df <- as.data.frame(boruta_pangenome_vst_a1$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_pangenome_vst_a1 <- importance(erie.classify_pangenome_vst_a1)
  imp_pangenome_vst_a1 <- data.frame(predictors = rownames(imp_pangenome_vst_a1), imp_pangenome_vst_a1)
  # add meta data
  imp_pangenome_vst_a1$Boruta_name <- rownames(boruta_pangenome_vst_a1_df)
  imp_pangenome_vst_a1$Boruta_predict <- boruta_pangenome_vst_a1_df$`boruta_pangenome_vst_a1$finalDecision`
  imp_pangenome_vst_a1_sub <- subset(imp_pangenome_vst_a1, MeanDecreaseAccuracy > 0.0)
  imp_pangenome_vst_a1_sub$predictors <- str_replace(imp_pangenome_vst_a1_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-associated
  imp_pangenome_vst_a1_sub$species_type <- with(imp_pangenome_vst_a1_sub,
                                                ifelse(predictors %in% f_rare_species_pangenome_vst_a1, 
                                                       'Rare species', 
                                                       ifelse(predictors %in% f_core_species_pangenome_vst_a1, 
                                                              'Core species', 
                                                              ifelse(predictors %in% fluctuating_rare_pangenome_vst_a1, 
                                                                     'Rare species',
                                                                     ifelse(predictors %in% fluctuating_core_pangenome_vst_a1, 
                                                                            'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_pangenome_vst_a1_sub$species_type_2 <- with(imp_pangenome_vst_a1_sub,
                                                  ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                         ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                                ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                                       ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_pangenome_vst_a1_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                          'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                          'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_pangenome_vst_a1_sub <- subset(imp_pangenome_vst_a1_sub, Boruta_predict != "Rejected")
  imp_pangenome_vst_a1_short <-ddply(imp_pangenome_vst_a1_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_pangenome_vst_a1_short$Seed <- i
  imp_pangenome_vst_a1_short$Threshold <- '15th percentile'
  imp_pangenome_vst_a1_short$Normalisation <- 'VST'
  # store short table globally
  imp.sort.gini_pangenome_vst_a1_final <- rbind(imp.sort.gini_pangenome_vst_a1_final, imp_pangenome_vst_a1_short)
}


# random forest, pangenome, VST, 25th abundance percentile
# define the non-persisting rare species (only found in one or two age groups but not always)
fluctuating_rare_pangenome_vst_a2 <- c(all_rare_a2_pangenome_vst, f_rare_species_pangenome_vst_a2)
# remove duplicate entries
fluctuating_rare_pangenome_vst_a2 <- fluctuating_rare_pangenome_vst_a2[!duplicated(fluctuating_rare_pangenome_vst_a2)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_pangenome_vst_a2 <- c(all_core_a2_pangenome_vst, f_core_species_pangenome_vst_a2)
# remove duplicate entries
fluctuating_core_pangenome_vst_a2 <- fluctuating_core_pangenome_vst_a2[!duplicated(fluctuating_core_pangenome_vst_a2)]

# get all species found in children
uniform_species_pangenome_vst_a2 = c(all_rare_a2_pangenome_vst, all_core_a2_pangenome_vst)
# extract genus information by splitting species name into two parts
uniform_genus_pangenome_vst_a2 = unlist(map(strsplit(uniform_species_pangenome_vst_a2, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_pangenome_vst_a2 = data.frame(cbind(uniform_species_pangenome_vst_a2, uniform_genus_pangenome_vst_a2))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_pangenome_vst_a2$Phylum <- 'X'
uniform_df_pangenome_vst_a2$Class <- 'X'
uniform_df_pangenome_vst_a2$Order <- 'X'
uniform_df_pangenome_vst_a2$Family <- 'X'
# re-name columns
colnames(uniform_df_pangenome_vst_a2) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_pangenome_vst_a2 <- uniform_df_pangenome_vst_a2[!duplicated(uniform_df_pangenome_vst_a2), ]
# re-name row names
rownames(uniform_df_pangenome_vst_a2) <- uniform_df_pangenome_vst_a2$otu
# remove non-numeric information from data frame
uniform_df_pangenome_vst_a2$otu <- NULL
# transpose data frame
uniform_df_pangenome_vst_a2 <- data.frame(t(uniform_df_pangenome_vst_a2))
# add columns with host-associated variables
uniform_df_pangenome_vst_a2$Age <- 'Age'
uniform_df_pangenome_vst_a2$BMI <- 'BMI'
uniform_df_pangenome_vst_a2$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_pangenome_vst_a2 <- data.frame(t(uniform_df_pangenome_vst_a2))
# clean up the row names
rownames(uniform_df_pangenome_vst_a2) <- str_replace(rownames(uniform_df_pangenome_vst_a2), '\\.', ' ')

# subset pangenome database based on abundance estimations (25th abundance percentile)
ds_pangenome_vst_all_a2 <- subset(ds_pangenome_vst, rownames(ds_pangenome_vst) %in% rownames(uniform_df_pangenome_vst_a2))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_pangenome_vst_all_a2 <- data.frame(t(ds_pangenome_vst_all_a2))
# order pangenome data frame by patient id (alphabetically)
ds_pangenome_vst_all_a2 <- ds_pangenome_vst_all_a2[order(rownames(ds_pangenome_vst_all_a2)),]

# add host-associated variables
ds_pangenome_vst_all_a2$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_pangenome_vst_all_a2$Age <- ifelse(ds_pangenome_vst_all_a2$Age_1 == '0 years', 0, ifelse(ds_pangenome_vst_all_a2$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_pangenome_vst_all_a2$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_pangenome_vst_all_a2$BMI <- md_sub$BMI
# remove age character column
ds_pangenome_vst_all_a2$Age_1 <- NULL

# re-arrange data frames
ds_pangenome_vst_all_a2_t <- data.frame(t(ds_pangenome_vst_all_a2))
md_sub_2_pangenome_vst_a2 <- data.frame(t(md_sub))
md_sub_3_pangenome_vst_a2 <- data.frame(t(md_sub_2_pangenome_vst_a2))

# convert data frames to matrix 
otu_mat_pangenome_vst_a2 <- as.matrix(ds_pangenome_vst_all_a2_t)
tax_mat_pangenome_vst_a2 <- as.matrix(uniform_df_pangenome_vst_a2)
# convert pangenome_vst to otu_table (phyloseq object)
OTU_pangenome_vst_a2 = otu_table(otu_mat_pangenome_vst_a2, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_pangenome_vst_a2) <- str_replace(rownames(OTU_pangenome_vst_a2), '\\.', ' ')
# make tax table (phyloseq object)
TAX_pangenome_vst_a2 = tax_table(tax_mat_pangenome_vst_a2)
# make meta data table (phyloseq object)
samples_pangenome_vst_a2 = sample_data(md_sub_3_pangenome_vst_a2)
# generate final phyloseq object
erie.merge_pangenome_vst_a2 <- phyloseq(OTU_pangenome_vst_a2, TAX_pangenome_vst_a2, samples_pangenome_vst_a2)

# transpose data frame
predictors_pangenome_vst_a2 <- t(otu_table(erie.merge_pangenome_vst_a2))
# obtain response variable (CF vs healthy)
response_pangenome_vst_a2 <- as.factor(sample_data(erie.merge_pangenome_vst_a2)$State)
# add response variable to data frame
rf.data_pangenome_vst_a2 <- data.frame(response_pangenome_vst_a2, predictors_pangenome_vst_a2)
# determine mtry for random forest
sample_val_pangenome_vst_a2 = sqrt(ncol(rf.data_pangenome_vst_a2))
# generate three empty variables
rf_pangenome_vst_a2_pred_final = NULL
imp.sort.gini_pangenome_vst_a2_final = NULL
error_rate_pangenome_vst_a2_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_pangenome_vst_a2 <- randomForest(response_pangenome_vst_a2~., data = rf.data_pangenome_vst_a2, 
                                                 ntree = 80, 
                                                 mtry=sample_val_pangenome_vst_a2, 
                                                 importance=TRUE,
                                                 na.action = na.roughfix)
  
  # store error rate locally
  error_rate_pangenome_vst_a2 <- erie.classify_pangenome_vst_a2$err.rate
  # make data frame of error rate
  error_rate_pangenome_vst_a2 <- as.data.frame(error_rate_pangenome_vst_a2)
  # add meta data
  error_rate_pangenome_vst_a2$Database <- 'Pangenome'
  error_rate_pangenome_vst_a2$Normalisation <- 'VST'
  error_rate_pangenome_vst_a2$Threshold <- '25th percentile'
  error_rate_pangenome_vst_a2$Seed <- i
  # re-name columns
  colnames(error_rate_pangenome_vst_a2) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_pangenome_vst_a2_final <- rbind(error_rate_pangenome_vst_a2_final, error_rate_pangenome_vst_a2)
  # obtain class probability
  rf_pangenome_vst_a2_pred <- predict(erie.classify_pangenome_vst_a2, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_pangenome_vst_a2_pred <- as.data.frame(rf_pangenome_vst_a2_pred)
  # add meta data
  rf_pangenome_vst_a2_pred$State <- md$State
  rf_pangenome_vst_a2_pred$Database <- 'Pangenome'
  rf_pangenome_vst_a2_pred$Normalisation <- 'VST'
  rf_pangenome_vst_a2_pred$Threshold <- '25th percentile'
  rf_pangenome_vst_a2_pred$Seed <- i
  rownames(rf_pangenome_vst_a2_pred) <- NULL
  # transfer data frame to global environment
  rf_pangenome_vst_a2_pred_final <- rbind(rf_pangenome_vst_a2_pred_final, rf_pangenome_vst_a2_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_pangenome_vst_a2 <- Boruta(response_pangenome_vst_a2~., data = rf.data_pangenome_vst_a2, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_pangenome_vst_a2_df <- as.data.frame(boruta_pangenome_vst_a2$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_pangenome_vst_a2 <- importance(erie.classify_pangenome_vst_a2)
  imp_pangenome_vst_a2 <- data.frame(predictors = rownames(imp_pangenome_vst_a2), imp_pangenome_vst_a2)
  # add meta data
  imp_pangenome_vst_a2$Boruta_name <- rownames(boruta_pangenome_vst_a2_df)
  imp_pangenome_vst_a2$Boruta_predict <- boruta_pangenome_vst_a2_df$`boruta_pangenome_vst_a2$finalDecision`
  imp_pangenome_vst_a2_sub <- subset(imp_pangenome_vst_a2, MeanDecreaseAccuracy > 0.0)
  imp_pangenome_vst_a2_sub$predictors <- str_replace(imp_pangenome_vst_a2_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-assosciated
  imp_pangenome_vst_a2_sub$species_type <- with(imp_pangenome_vst_a2_sub,
                                                ifelse(predictors %in% f_rare_species_pangenome_vst_a2, 
                                                       'Rare species', 
                                                       ifelse(predictors %in% f_core_species_pangenome_vst_a2, 
                                                              'Core species', 
                                                              ifelse(predictors %in% fluctuating_rare_pangenome_vst_a2, 
                                                                     'Rare species',
                                                                     ifelse(predictors %in% fluctuating_core_pangenome_vst_a2, 
                                                                            'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_pangenome_vst_a2_sub$species_type_2 <- with(imp_pangenome_vst_a2_sub,
                                                  ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                         ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                                ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                                       ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_pangenome_vst_a2_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                          'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                          'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_pangenome_vst_a2_sub <- subset(imp_pangenome_vst_a2_sub, Boruta_predict != "Rejected")
  imp_pangenome_vst_a2_short <- ddply(imp_pangenome_vst_a2_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_pangenome_vst_a2_short$Seed <- i
  imp_pangenome_vst_a2_short$Threshold <- '25th percentile'
  imp_pangenome_vst_a2_short$Normalisation <- 'VST'
  # store short table globally
  imp.sort.gini_pangenome_vst_a2_final <- rbind(imp.sort.gini_pangenome_vst_a2_final, imp_pangenome_vst_a2_short)
}


# random forest, pangenome, VST, 35th abundance percentile
# define the non-persisting rare species (only found in one or two age groups but not always)
fluctuating_rare_pangenome_vst_a3 <- c(all_rare_a3_pangenome_vst, f_rare_species_pangenome_vst_a3)
# remove duplicate entries
fluctuating_rare_pangenome_vst_a3 <- fluctuating_rare_pangenome_vst_a3[!duplicated(fluctuating_rare_pangenome_vst_a3)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_pangenome_vst_a3 <- c(all_core_a3_pangenome_vst, f_core_species_pangenome_vst_a3)
# remove duplicate entries
fluctuating_core_pangenome_vst_a3 <- fluctuating_core_pangenome_vst_a3[!duplicated(fluctuating_core_pangenome_vst_a3)]

# get all species found in children
uniform_species_pangenome_vst_a3 = c(all_rare_a3_pangenome_vst, all_core_a3_pangenome_vst)
# extract genus information by splitting species name into two parts
uniform_genus_pangenome_vst_a3 = unlist(map(strsplit(uniform_species_pangenome_vst_a3, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_pangenome_vst_a3 = data.frame(cbind(uniform_species_pangenome_vst_a3, uniform_genus_pangenome_vst_a3))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_pangenome_vst_a3$Phylum <- 'X'
uniform_df_pangenome_vst_a3$Class <- 'X'
uniform_df_pangenome_vst_a3$Order <- 'X'
uniform_df_pangenome_vst_a3$Family <- 'X'
# re-name columns
colnames(uniform_df_pangenome_vst_a3) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_pangenome_vst_a3 <- uniform_df_pangenome_vst_a3[!duplicated(uniform_df_pangenome_vst_a3), ]
# re-name row names
rownames(uniform_df_pangenome_vst_a3) <- uniform_df_pangenome_vst_a3$otu
# remove non-numeric information from data frame
uniform_df_pangenome_vst_a3$otu <- NULL
# transpose data frame
uniform_df_pangenome_vst_a3 <- data.frame(t(uniform_df_pangenome_vst_a3))
# add columns with host-associated variables
uniform_df_pangenome_vst_a3$Age <- 'Age'
uniform_df_pangenome_vst_a3$BMI <- 'BMI'
uniform_df_pangenome_vst_a3$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_pangenome_vst_a3 <- data.frame(t(uniform_df_pangenome_vst_a3))
# clean up the row names
rownames(uniform_df_pangenome_vst_a3) <- str_replace(rownames(uniform_df_pangenome_vst_a3), '\\.', ' ')

# subset pangenome database based on abundance estimations (35th abundance percentile)
ds_pangenome_vst_all_a3 <- subset(ds_pangenome_vst, rownames(ds_pangenome_vst) %in% rownames(uniform_df_pangenome_vst_a3))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_pangenome_vst_all_a3 <- data.frame(t(ds_pangenome_vst_all_a3))
# order pangenome data frame by patient id (alphabetically)
ds_pangenome_vst_all_a3 <- ds_pangenome_vst_all_a3[order(rownames(ds_pangenome_vst_all_a3)),]

# add host-associated variables
ds_pangenome_vst_all_a3$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_pangenome_vst_all_a3$Age <- ifelse(ds_pangenome_vst_all_a3$Age_1 == '0 years', 0, ifelse(ds_pangenome_vst_all_a3$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_pangenome_vst_all_a3$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_pangenome_vst_all_a3$BMI <- md_sub$BMI
# remove age character column
ds_pangenome_vst_all_a3$Age_1 <- NULL

# re-arrange data frames
ds_pangenome_vst_all_a3_t <- data.frame(t(ds_pangenome_vst_all_a3))
md_sub_2_pangenome_vst_a3 <- data.frame(t(md_sub))
md_sub_3_pangenome_vst_a3 <- data.frame(t(md_sub_2_pangenome_vst_a3))

# convert data frames to matrix 
otu_mat_pangenome_vst_a3 <- as.matrix(ds_pangenome_vst_all_a3_t)
tax_mat_pangenome_vst_a3 <- as.matrix(uniform_df_pangenome_vst_a3)
# convert pangenome to otu_table (phyloseq object)
OTU_pangenome_vst_a3 = otu_table(otu_mat_pangenome_vst_a3, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_pangenome_vst_a3) <- str_replace(rownames(OTU_pangenome_vst_a3), '\\.', ' ')
# make tax table (phyloseq object)
TAX_pangenome_vst_a3 = tax_table(tax_mat_pangenome_vst_a3)
# make meta data table (phyloseq object)
samples_pangenome_vst_a3 = sample_data(md_sub_3_pangenome_vst_a3)
# generate final phyloseq object
erie.merge_pangenome_vst_a3 <- phyloseq(OTU_pangenome_vst_a3, TAX_pangenome_vst_a3, samples_pangenome_vst_a3)

# transpose data frame
predictors_pangenome_vst_a3 <- t(otu_table(erie.merge_pangenome_vst_a3))
# obtain response variable (CF vs healthy)
response_pangenome_vst_a3 <- as.factor(sample_data(erie.merge_pangenome_vst_a3)$State)
# add response variable to data frame
rf.data_pangenome_vst_a3 <- data.frame(response_pangenome_vst_a3, predictors_pangenome_vst_a3)
# determine mtry for random forest
sample_val_pangenome_vst_a3 = sqrt(ncol(rf.data_pangenome_vst_a3))
# generate three empty variables
rf_pangenome_vst_a3_pred_final = NULL
imp.sort.gini_pangenome_vst_a3_final = NULL
error_rate_pangenome_vst_a3_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_pangenome_vst_a3 <- randomForest(response_pangenome_vst_a3~., data = rf.data_pangenome_vst_a3, 
                                                 ntree = 80, 
                                                 mtry=sample_val_pangenome_vst_a3, 
                                                 importance=TRUE,
                                                 na.action = na.roughfix)
  
  # store error rate locally
  error_rate_pangenome_vst_a3 <- erie.classify_pangenome_vst_a3$err.rate
  # make data frame of error rate
  error_rate_pangenome_vst_a3 <- as.data.frame(error_rate_pangenome_vst_a3)
  # add meta data
  error_rate_pangenome_vst_a3$Database <- 'Pangenome'
  error_rate_pangenome_vst_a3$Normalisation <- 'VST'
  error_rate_pangenome_vst_a3$Threshold <- '35th percentile'
  error_rate_pangenome_vst_a3$Seed <- i
  # re-name columns
  colnames(error_rate_pangenome_vst_a3) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_pangenome_vst_a3_final <- rbind(error_rate_pangenome_vst_a3_final, error_rate_pangenome_vst_a3)
  # obtain class probability
  rf_pangenome_vst_a3_pred <- predict(erie.classify_pangenome_vst_a3, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_pangenome_vst_a3_pred <- as.data.frame(rf_pangenome_vst_a3_pred)
  # add meta data
  rf_pangenome_vst_a3_pred$State <- md$State
  rf_pangenome_vst_a3_pred$Database <- 'Pangenome'
  rf_pangenome_vst_a3_pred$Normalisation <- 'VST'
  rf_pangenome_vst_a3_pred$Threshold <- '35th percentile'
  rf_pangenome_vst_a3_pred$Seed <- i
  rownames(rf_pangenome_vst_a3_pred) <- NULL
  # transfer data frame to global environment
  rf_pangenome_vst_a3_pred_final <- rbind(rf_pangenome_vst_a3_pred_final, rf_pangenome_vst_a3_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_pangenome_vst_a3 <- Boruta(response_pangenome_vst_a3~., data = rf.data_pangenome_vst_a3, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_pangenome_vst_a3_df <- as.data.frame(boruta_pangenome_vst_a3$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_pangenome_vst_a3 <- importance(erie.classify_pangenome_vst_a3)
  imp_pangenome_vst_a3 <- data.frame(predictors = rownames(imp_pangenome_vst_a3), imp_pangenome_vst_a3)
  # add meta data
  imp_pangenome_vst_a3$Boruta_name <- rownames(boruta_pangenome_vst_a3_df)
  imp_pangenome_vst_a3$Boruta_predict <- boruta_pangenome_vst_a3_df$`boruta_pangenome_vst_a3$finalDecision`
  imp_pangenome_vst_a3_sub <- subset(imp_pangenome_vst_a3, MeanDecreaseAccuracy > 0.0)
  imp_pangenome_vst_a3_sub$predictors <- str_replace(imp_pangenome_vst_a3_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-assosciated
  imp_pangenome_vst_a3_sub$species_type <- with(imp_pangenome_vst_a3_sub,
                                                ifelse(predictors %in% f_rare_species_pangenome_vst_a3, 
                                                       'Rare species', 
                                                       ifelse(predictors %in% f_core_species_pangenome_vst_a3, 
                                                              'Core species', 
                                                              ifelse(predictors %in% fluctuating_rare_pangenome_vst_a3, 
                                                                     'Rare species',
                                                                     ifelse(predictors %in% fluctuating_core_pangenome_vst_a3, 
                                                                            'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_pangenome_vst_a3_sub$species_type_2 <- with(imp_pangenome_vst_a3_sub,
                                                  ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                         ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                                ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                                       ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_pangenome_vst_a3_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                          'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                          'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_pangenome_vst_a3_sub <- subset(imp_pangenome_vst_a3_sub, Boruta_predict != "Rejected")
  imp_pangenome_vst_a3_short <- ddply(imp_pangenome_vst_a3_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_pangenome_vst_a3_short$Seed <- i
  imp_pangenome_vst_a3_short$Threshold <- '35th percentile'
  imp_pangenome_vst_a3_short$Normalisation <- 'VST'
  # store short table globally
  imp.sort.gini_pangenome_vst_a3_final <- rbind(imp.sort.gini_pangenome_vst_a3_final, imp_pangenome_vst_a3_short)
}


# merge all information, random forest, pangenome, VST, 15th abundance percentile, 25th abundance percentile, 35th abundance percentile
rf.imp.sort.pangenome_vst <- data.frame(rbind(imp.sort.gini_pangenome_vst_a1_final,imp.sort.gini_pangenome_vst_a2_final, imp.sort.gini_pangenome_vst_a3_final))
rf.imp.sort.pangenome_vst$Database <- 'Pangenome'
rf.pred.prob.pangenome_vst <- data.frame(rbind(rf_pangenome_vst_a1_pred_final,rf_pangenome_vst_a2_pred_final, rf_pangenome_vst_a3_pred_final))
rf.error.pangenome_vst <- data.frame(rbind(error_rate_pangenome_vst_a1_final,error_rate_pangenome_vst_a2_final, error_rate_pangenome_vst_a3_final))

# merge all information, random forest, pangenome, VST, BCPHC, RLE, 15th abundance, 25th abundance, 35th abundance
rf.imp.sort.pangenome.all <- data.frame(rbind(rf.imp.sort.pangenome, rf.imp.sort.pangenome_rle, rf.imp.sort.pangenome_vst))
rf.pred.prob.pangenome.all <- data.frame(rbind(rf.pred.prob.pangenome, rf.pred.prob.pangenome_rle, rf.pred.prob.pangenome_vst))
rf.error.pangenome.all <- data.frame(rbind(rf.error.pangenome,rf.error.pangenome_rle, rf.error.pangenome_vst))



# Continue random forest with one-strain-per species reference database data
# Random forest, osps, bcphc, 15th abundance percentile
# define the non-persisting rare species (only found in one or two age groups but not always)
fluctuating_rare_osps_a1 <- c(all_rare_a1_osps, f_rare_species_osps_a1)
# remove duplicate entries
fluctuating_rare_osps_a1 <- fluctuating_rare_osps_a1[!duplicated(fluctuating_rare_osps_a1)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_osps_a1 <- c(all_core_a1_osps, f_core_species_osps_a1)
# remove duplicate entries
fluctuating_core_osps_a1 <- fluctuating_core_osps_a1[!duplicated(fluctuating_core_osps_a1)]

# get all species found in children
uniform_species_osps_a1 = c(all_rare_a1_osps, all_core_a1_osps)
# extract genus information by splitting species name into two parts
uniform_genus_osps_a1 = unlist(map(strsplit(uniform_species_osps_a1, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_osps_a1 = data.frame(cbind(uniform_species_osps_a1, uniform_genus_osps_a1))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_osps_a1$Phylum <- 'X'
uniform_df_osps_a1$Class <- 'X'
uniform_df_osps_a1$Order <- 'X'
uniform_df_osps_a1$Family <- 'X'
# re-name columns
colnames(uniform_df_osps_a1) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_osps_a1 <- uniform_df_osps_a1[!duplicated(uniform_df_osps_a1), ]
# re-name row names
rownames(uniform_df_osps_a1) <- uniform_df_osps_a1$otu
# remove non-numeric information from data frame
uniform_df_osps_a1$otu <- NULL
# transpose data frame
uniform_df_osps_a1 <- data.frame(t(uniform_df_osps_a1))
# add columns with host-associated variables
uniform_df_osps_a1$Age <- 'Age'
uniform_df_osps_a1$BMI <- 'BMI'
uniform_df_osps_a1$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_osps_a1 <- data.frame(t(uniform_df_osps_a1))
# clean up the row names
rownames(uniform_df_osps_a1) <- str_replace(rownames(uniform_df_osps_a1), '\\.', ' ')

# subset osps database based on abundance estimations (15th abundance percentile)
ds_osps_all_a1 <- subset(ds_osps, rownames(ds_osps) %in% rownames(uniform_df_osps_a1))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_osps_all_a1 <- data.frame(t(ds_osps_all_a1))
# order osps data frame by patient id (alphabetically)
ds_osps_all_a1 <- ds_osps_all_a1[order(rownames(ds_osps_all_a1)),]

# add host-associated variables
ds_osps_all_a1$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_osps_all_a1$Age <- ifelse(ds_osps_all_a1$Age_1 == '0 years', 0, ifelse(ds_osps_all_a1$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_osps_all_a1$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_osps_all_a1$BMI <- md_sub$BMI
# remove age character column
ds_osps_all_a1$Age_1 <- NULL

# re-arrange data frames
ds_osps_all_a1_t <- data.frame(t(ds_osps_all_a1))
md_sub_2_osps_a1 <- data.frame(t(md_sub))
md_sub_3_osps_a1 <- data.frame(t(md_sub_2_osps_a1))

# convert data frames to matrix 
otu_mat_osps_a1 <- as.matrix(ds_osps_all_a1_t)
tax_mat_osps_a1 <- as.matrix(uniform_df_osps_a1)
# convert osps to otu_table (phyloseq object)
OTU_osps_a1 = otu_table(otu_mat_osps_a1, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_osps_a1) <- str_replace(rownames(OTU_osps_a1), '\\.', ' ')
# make tax table (phyloseq object)
TAX_osps_a1 = tax_table(tax_mat_osps_a1)
# make meta data table (phyloseq object)
samples_osps_a1 = sample_data(md_sub_3_osps_a1)
# generate final phyloseq object
erie.merge_osps_a1 <- phyloseq(OTU_osps_a1, TAX_osps_a1, samples_osps_a1)

# transpose data frame
predictors_osps_a1 <- t(otu_table(erie.merge_osps_a1))
# obtain response variable (CF vs healthy)
response_osps_a1 <- as.factor(sample_data(erie.merge_osps_a1)$State)
# add response variable to data frame
rf.data_osps_a1 <- data.frame(response_osps_a1, predictors_osps_a1)
# determine mtry for random forest
sample_val_osps_a1 = sqrt(ncol(rf.data_osps_a1))
# generate three empty variables
rf_osps_a1_pred_final = NULL
imp.sort.gini_osps_a1_final = NULL
error_rate_osps_a1_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_osps_a1 <- randomForest(response_osps_a1~., data = rf.data_osps_a1, 
                                        ntree = 80, 
                                        mtry=sample_val_osps_a1, 
                                        importance=TRUE,
                                        na.action = na.roughfix)
  
  # store error rate locally
  error_rate_osps_a1 <- erie.classify_osps_a1$err.rate
  # make data frame of error rate
  error_rate_osps_a1 <- as.data.frame(error_rate_osps_a1)
  # add meta data
  error_rate_osps_a1$Database <- 'One strain per species'
  error_rate_osps_a1$Normalisation <- 'BCPHC'
  error_rate_osps_a1$Threshold <- '15th percentile'
  error_rate_osps_a1$Seed <- i
  # re-name columns
  colnames(error_rate_osps_a1) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_osps_a1_final <- rbind(error_rate_osps_a1_final, error_rate_osps_a1)
  # obtain class probability
  rf_osps_a1_pred <- predict(erie.classify_osps_a1, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_osps_a1_pred <- as.data.frame(rf_osps_a1_pred)
  # add meta data
  rf_osps_a1_pred$State <- md$State
  rf_osps_a1_pred$Database <- 'One strain per species'
  rf_osps_a1_pred$Normalisation <- 'BCPHC'
  rf_osps_a1_pred$Threshold <- '15th percentile'
  rf_osps_a1_pred$Seed <- i
  rownames(rf_osps_a1_pred) <- NULL
  # transfer data frame to global environment
  rf_osps_a1_pred_final <- rbind(rf_osps_a1_pred_final, rf_osps_a1_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_osps_a1 <- Boruta(response_osps_a1~., data = rf.data_osps_a1, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_osps_a1_df <- as.data.frame(boruta_osps_a1$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_osps_a1 <- importance(erie.classify_osps_a1)
  imp_osps_a1 <- data.frame(predictors = rownames(imp_osps_a1), imp_osps_a1)
  # add meta data
  imp_osps_a1$Boruta_name <- rownames(boruta_osps_a1_df)
  imp_osps_a1$Boruta_predict <- boruta_osps_a1_df$`boruta_osps_a1$finalDecision`
  imp_osps_a1_sub <- subset(imp_osps_a1, MeanDecreaseAccuracy > 0.0)
  imp_osps_a1_sub$predictors <- str_replace(imp_osps_a1_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-assosciated
  imp_osps_a1_sub$species_type <- with(imp_osps_a1_sub,
                                       ifelse(predictors %in% f_rare_species_osps_a1, 
                                              'Rare species', 
                                              ifelse(predictors %in% f_core_species_osps_a1, 
                                                     'Core species', 
                                                     ifelse(predictors %in% fluctuating_rare_osps_a1, 
                                                            'Rare species',
                                                            ifelse(predictors %in% fluctuating_core_osps_a1, 
                                                                   'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_osps_a1_sub$species_type_2 <- with(imp_osps_a1_sub,
                                         ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                       ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                              ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_osps_a1_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                 'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                 'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_osps_a1_sub <- subset(imp_osps_a1_sub, Boruta_predict != "Rejected")
  imp_osps_a1_short <- ddply(imp_osps_a1_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_osps_a1_short$Seed <- i
  imp_osps_a1_short$Threshold <- '15th percentile'
  imp_osps_a1_short$Normalisation <- 'BCPHC'
  # store short table globally
  imp.sort.gini_osps_a1_final <- rbind(imp.sort.gini_osps_a1_final, imp_osps_a1_short)
}


# random forest, osps, bcphc, 25th abundance percentile
# define the non-persisting rare species (only found in one or two age groups but not always)
fluctuating_rare_osps_a2 <- c(all_rare_a2_osps, f_rare_species_osps_a2)
# remove duplicate entries
fluctuating_rare_osps_a2 <- fluctuating_rare_osps_a2[!duplicated(fluctuating_rare_osps_a2)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_osps_a2 <- c(all_core_a2_osps, f_core_species_osps_a2)
# remove duplicate entries
fluctuating_core_osps_a2 <- fluctuating_core_osps_a2[!duplicated(fluctuating_core_osps_a2)]

# get all species found in children
uniform_species_osps_a2 = c(all_rare_a2_osps, all_core_a2_osps)
# extract genus information by splitting species name into two parts
uniform_genus_osps_a2 = unlist(map(strsplit(uniform_species_osps_a2, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_osps_a2 = data.frame(cbind(uniform_species_osps_a2, uniform_genus_osps_a2))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_osps_a2$Phylum <- 'X'
uniform_df_osps_a2$Class <- 'X'
uniform_df_osps_a2$Order <- 'X'
uniform_df_osps_a2$Family <- 'X'
# re-name columns
colnames(uniform_df_osps_a2) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_osps_a2 <- uniform_df_osps_a2[!duplicated(uniform_df_osps_a2), ]
# re-name row names
rownames(uniform_df_osps_a2) <- uniform_df_osps_a2$otu
# remove non-numeric information from data frame
uniform_df_osps_a2$otu <- NULL
# transpose data frame
uniform_df_osps_a2 <- data.frame(t(uniform_df_osps_a2))
# add columns with host-associated variables
uniform_df_osps_a2$Age <- 'Age'
uniform_df_osps_a2$BMI <- 'BMI'
uniform_df_osps_a2$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_osps_a2 <- data.frame(t(uniform_df_osps_a2))
# clean up the row names
rownames(uniform_df_osps_a2) <- str_replace(rownames(uniform_df_osps_a2), '\\.', ' ')

# subset osps database based on abundance estimations (25th abundance percentile)
ds_osps_all_a2 <- subset(ds_osps, rownames(ds_osps) %in% rownames(uniform_df_osps_a2))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_osps_all_a2 <- data.frame(t(ds_osps_all_a2))
# order osps data frame by patient id (alphabetically)
ds_osps_all_a2 <- ds_osps_all_a2[order(rownames(ds_osps_all_a2)),]

# add host-associated variables
ds_osps_all_a2$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_osps_all_a2$Age <- ifelse(ds_osps_all_a2$Age_1 == '0 years', 0, ifelse(ds_osps_all_a2$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_osps_all_a2$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_osps_all_a2$BMI <- md_sub$BMI
# remove age character column
ds_osps_all_a2$Age_1 <- NULL

# re-arrange data frames
ds_osps_all_a2_t <- data.frame(t(ds_osps_all_a2))
md_sub_2_osps_a2 <- data.frame(t(md_sub))
md_sub_3_osps_a2 <- data.frame(t(md_sub_2_osps_a2))

# convert data frames to matrix 
otu_mat_osps_a2 <- as.matrix(ds_osps_all_a2_t)
tax_mat_osps_a2 <- as.matrix(uniform_df_osps_a2)
# convert osps to otu_table (phyloseq object)
OTU_osps_a2 = otu_table(otu_mat_osps_a2, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_osps_a2) <- str_replace(rownames(OTU_osps_a2), '\\.', ' ')
# make tax table (phyloseq object)
TAX_osps_a2 = tax_table(tax_mat_osps_a2)
# make meta data table (phyloseq object)
samples_osps_a2 = sample_data(md_sub_3_osps_a2)
# generate final phyloseq object
erie.merge_osps_a2 <- phyloseq(OTU_osps_a2, TAX_osps_a2, samples_osps_a2)

# transpose data frame
predictors_osps_a2 <- t(otu_table(erie.merge_osps_a2))
# obtain response variable (CF vs healthy)
response_osps_a2 <- as.factor(sample_data(erie.merge_osps_a2)$State)
# add response variable to data frame
rf.data_osps_a2 <- data.frame(response_osps_a2, predictors_osps_a2)
# determine mtry for random forest
sample_val_osps_a2 = sqrt(ncol(rf.data_osps_a2))
# generate three empty variables
rf_osps_a2_pred_final = NULL
imp.sort.gini_osps_a2_final = NULL
error_rate_osps_a2_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_osps_a2 <- randomForest(response_osps_a2~., data = rf.data_osps_a2, 
                                        ntree = 80, 
                                        mtry=sample_val_osps_a2, 
                                        importance=TRUE,
                                        na.action = na.roughfix)
  
  # store error rate locally
  error_rate_osps_a2 <- erie.classify_osps_a2$err.rate
  # make data frame of error rate
  error_rate_osps_a2 <- as.data.frame(error_rate_osps_a2)
  # add meta data
  error_rate_osps_a2$Database <- 'One strain per species'
  error_rate_osps_a2$Normalisation <- 'BCPHC'
  error_rate_osps_a2$Threshold <- '25th percentile'
  error_rate_osps_a2$Seed <- i
  # re-name columns
  colnames(error_rate_osps_a2) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_osps_a2_final <- rbind(error_rate_osps_a2_final, error_rate_osps_a2)
  # obtain class probability
  rf_osps_a2_pred <- predict(erie.classify_osps_a2, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_osps_a2_pred <- as.data.frame(rf_osps_a2_pred)
  # add meta data
  rf_osps_a2_pred$State <- md$State
  rf_osps_a2_pred$Database <- 'One strain per species'
  rf_osps_a2_pred$Normalisation <- 'BCPHC'
  rf_osps_a2_pred$Threshold <- '25th percentile'
  rf_osps_a2_pred$Seed <- i
  rownames(rf_osps_a2_pred) <- NULL
  # transfer data frame to global environment
  rf_osps_a2_pred_final <- rbind(rf_osps_a2_pred_final, rf_osps_a2_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_osps_a2 <- Boruta(response_osps_a2~., data = rf.data_osps_a2, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_osps_a2_df <- as.data.frame(boruta_osps_a2$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_osps_a2 <- importance(erie.classify_osps_a2)
  imp_osps_a2 <- data.frame(predictors = rownames(imp_osps_a2), imp_osps_a2)
  # add meta data
  imp_osps_a2$Boruta_name <- rownames(boruta_osps_a2_df)
  imp_osps_a2$Boruta_predict <- boruta_osps_a2_df$`boruta_osps_a2$finalDecision`
  imp_osps_a2_sub <- subset(imp_osps_a2, MeanDecreaseAccuracy > 0.0)
  imp_osps_a2_sub$predictors <- str_replace(imp_osps_a2_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-assosciated
  imp_osps_a2_sub$species_type <- with(imp_osps_a2_sub,
                                       ifelse(predictors %in% f_rare_species_osps_a2, 
                                              'Rare species', 
                                              ifelse(predictors %in% f_core_species_osps_a2, 
                                                     'Core species', 
                                                     ifelse(predictors %in% fluctuating_rare_osps_a2, 
                                                            'Rare species',
                                                            ifelse(predictors %in% fluctuating_core_osps_a2, 
                                                                   'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_osps_a2_sub$species_type_2 <- with(imp_osps_a2_sub,
                                         ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                       ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                              ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_osps_a2_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                 'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                 'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_osps_a2_sub <- subset(imp_osps_a2_sub, Boruta_predict != "Rejected")
  imp_osps_a2_short <- ddply(imp_osps_a2_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_osps_a2_short$Seed <- i
  imp_osps_a2_short$Threshold <- '25th percentile'
  imp_osps_a2_short$Normalisation <- 'BCPHC'
  # store short table globally
  imp.sort.gini_osps_a2_final <- rbind(imp.sort.gini_osps_a2_final, imp_osps_a2_short)
}

# random forest, osps, bcphc, 35th abundance percentile
# define the non-persisting rare species (only found in one or two age groups but not always)
fluctuating_rare_osps_a3 <- c(all_rare_a3_osps, f_rare_species_osps_a3)
# remove duplicate entries
fluctuating_rare_osps_a3 <- fluctuating_rare_osps_a3[!duplicated(fluctuating_rare_osps_a3)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_osps_a3 <- c(all_core_a3_osps, f_core_species_osps_a3)
# remove duplicate entries
fluctuating_core_osps_a3 <- fluctuating_core_osps_a3[!duplicated(fluctuating_core_osps_a3)]

# get all species found in children
uniform_species_osps_a3 = c(all_rare_a3_osps, all_core_a3_osps)
# extract genus information by splitting species name into two parts
uniform_genus_osps_a3 = unlist(map(strsplit(uniform_species_osps_a3, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_osps_a3 = data.frame(cbind(uniform_species_osps_a3, uniform_genus_osps_a3))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_osps_a3$Phylum <- 'X'
uniform_df_osps_a3$Class <- 'X'
uniform_df_osps_a3$Order <- 'X'
uniform_df_osps_a3$Family <- 'X'
# re-name columns
colnames(uniform_df_osps_a3) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_osps_a3 <- uniform_df_osps_a3[!duplicated(uniform_df_osps_a3), ]
# re-name row names
rownames(uniform_df_osps_a3) <- uniform_df_osps_a3$otu
# remove non-numeric information from data frame
uniform_df_osps_a3$otu <- NULL
# transpose data frame
uniform_df_osps_a3 <- data.frame(t(uniform_df_osps_a3))
# add columns with host-associated variables
uniform_df_osps_a3$Age <- 'Age'
uniform_df_osps_a3$BMI <- 'BMI'
uniform_df_osps_a3$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_osps_a3 <- data.frame(t(uniform_df_osps_a3))
# clean up the row names
rownames(uniform_df_osps_a3) <- str_replace(rownames(uniform_df_osps_a3), '\\.', ' ')

# subset osps database based on abundance estimations (35th abundance percentile)
ds_osps_all_a3 <- subset(ds_osps, rownames(ds_osps) %in% rownames(uniform_df_osps_a3))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_osps_all_a3 <- data.frame(t(ds_osps_all_a3))
# order osps data frame by patient id (alphabetically)
ds_osps_all_a3 <- ds_osps_all_a3[order(rownames(ds_osps_all_a3)),]

# add host-associated variables
ds_osps_all_a3$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_osps_all_a3$Age <- ifelse(ds_osps_all_a3$Age_1 == '0 years', 0, ifelse(ds_osps_all_a3$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_osps_all_a3$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_osps_all_a3$BMI <- md_sub$BMI
# remove age character column
ds_osps_all_a3$Age_1 <- NULL

# re-arrange data frames
ds_osps_all_a3_t <- data.frame(t(ds_osps_all_a3))
md_sub_2_osps_a3 <- data.frame(t(md_sub))
md_sub_3_osps_a3 <- data.frame(t(md_sub_2_osps_a3))

# convert data frames to matrix 
otu_mat_osps_a3 <- as.matrix(ds_osps_all_a3_t)
tax_mat_osps_a3 <- as.matrix(uniform_df_osps_a3)
# convert osps to otu_table (phyloseq object)
OTU_osps_a3 = otu_table(otu_mat_osps_a3, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_osps_a3) <- str_replace(rownames(OTU_osps_a3), '\\.', ' ')
# make tax table (phyloseq object)
TAX_osps_a3 = tax_table(tax_mat_osps_a3)
# make meta data table (phyloseq object)
samples_osps_a3 = sample_data(md_sub_3_osps_a3)
# generate final phyloseq object
erie.merge_osps_a3 <- phyloseq(OTU_osps_a3, TAX_osps_a3, samples_osps_a3)

# transpose data frame
predictors_osps_a3 <- t(otu_table(erie.merge_osps_a3))
# obtain response variable (CF vs healthy)
response_osps_a3 <- as.factor(sample_data(erie.merge_osps_a3)$State)
# add response variable to data frame
rf.data_osps_a3 <- data.frame(response_osps_a3, predictors_osps_a3)
# determine mtry for random forest
sample_val_osps_a3 = sqrt(ncol(rf.data_osps_a3))
# generate three empty variables
rf_osps_a3_pred_final = NULL
imp.sort.gini_osps_a3_final = NULL
error_rate_osps_a3_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_osps_a3 <- randomForest(response_osps_a3~., data = rf.data_osps_a3, 
                                        ntree = 80, 
                                        mtry=sample_val_osps_a3, 
                                        importance=TRUE,
                                        na.action = na.roughfix)
  
  # store error rate locally
  error_rate_osps_a3 <- erie.classify_osps_a3$err.rate
  # make data frame of error rate
  error_rate_osps_a3 <- as.data.frame(error_rate_osps_a3)
  # add meta data
  error_rate_osps_a3$Database <- 'One strain per species'
  error_rate_osps_a3$Normalisation <- 'BCPHC'
  error_rate_osps_a3$Threshold <- '35th percentile'
  error_rate_osps_a3$Seed <- i
  # re-name columns
  colnames(error_rate_osps_a3) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_osps_a3_final <- rbind(error_rate_osps_a3_final, error_rate_osps_a3)
  # obtain class probability
  rf_osps_a3_pred <- predict(erie.classify_osps_a3, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_osps_a3_pred <- as.data.frame(rf_osps_a3_pred)
  # add meta data
  rf_osps_a3_pred$State <- md$State
  rf_osps_a3_pred$Database <- 'One strain per species'
  rf_osps_a3_pred$Normalisation <- 'BCPHC'
  rf_osps_a3_pred$Threshold <- '35th percentile'
  rf_osps_a3_pred$Seed <- i
  rownames(rf_osps_a3_pred) <- NULL
  # transfer data frame to global environment
  rf_osps_a3_pred_final <- rbind(rf_osps_a3_pred_final, rf_osps_a3_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_osps_a3 <- Boruta(response_osps_a3~., data = rf.data_osps_a3, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_osps_a3_df <- as.data.frame(boruta_osps_a3$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_osps_a3 <- importance(erie.classify_osps_a3)
  imp_osps_a3 <- data.frame(predictors = rownames(imp_osps_a3), imp_osps_a3)
  # add meta data
  imp_osps_a3$Boruta_name <- rownames(boruta_osps_a3_df)
  imp_osps_a3$Boruta_predict <- boruta_osps_a3_df$`boruta_osps_a3$finalDecision`
  imp_osps_a3_sub <- subset(imp_osps_a3, MeanDecreaseAccuracy > 0.0)
  imp_osps_a3_sub$predictors <- str_replace(imp_osps_a3_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-assosciated
  imp_osps_a3_sub$species_type <- with(imp_osps_a3_sub,
                                       ifelse(predictors %in% f_rare_species_osps_a3, 
                                              'Rare species', 
                                              ifelse(predictors %in% f_core_species_osps_a3, 
                                                     'Core species', 
                                                     ifelse(predictors %in% fluctuating_rare_osps_a3, 
                                                            'Rare species',
                                                            ifelse(predictors %in% fluctuating_core_osps_a3, 
                                                                   'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_osps_a3_sub$species_type_2 <- with(imp_osps_a3_sub,
                                         ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                       ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                              ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_osps_a3_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                 'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                 'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_osps_a3_sub <- subset(imp_osps_a3_sub, Boruta_predict != "Rejected")
  imp_osps_a3_short <- ddply(imp_osps_a3_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_osps_a3_short$Seed <- i
  imp_osps_a3_short$Threshold <- '35th percentile'
  imp_osps_a3_short$Normalisation <- 'BCPHC'
  # store short table globally
  imp.sort.gini_osps_a3_final <- rbind(imp.sort.gini_osps_a3_final, imp_osps_a3_short)
}


# merge all information, random forest, osps, BCPHC, 15th abundance percentile, 25th abundance percentile, 35th abundance percentile
rf.imp.sort.osps <- data.frame(rbind(imp.sort.gini_osps_a1_final,imp.sort.gini_osps_a2_final, imp.sort.gini_osps_a3_final))
rf.imp.sort.osps$Database <- 'One strain per species'
rf.pred.prob.osps <- data.frame(rbind(rf_osps_a1_pred_final,rf_osps_a2_pred_final, rf_osps_a3_pred_final))
rf.error.osps <- data.frame(rbind(error_rate_osps_a1_final,error_rate_osps_a2_final, error_rate_osps_a3_final))


# Continue with RLE-normalised data
# random forest for RLE normalised data, 15th abundance percentile
fluctuating_rare_osps_rle_a1 <- c(all_rare_a1_osps_rle, f_rare_species_osps_rle_a1)
# remove duplicate entries
fluctuating_rare_osps_rle_a1 <- fluctuating_rare_osps_rle_a1[!duplicated(fluctuating_rare_osps_rle_a1)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_osps_rle_a1 <- c(all_core_a1_osps_rle, f_core_species_osps_rle_a1)
# remove duplicate entries
fluctuating_core_osps_rle_a1 <- fluctuating_core_osps_rle_a1[!duplicated(fluctuating_core_osps_rle_a1)]

# get all species found in children
uniform_species_osps_rle_a1 = c(all_rare_a1_osps_rle, all_core_a1_osps_rle)
# extract genus information by splitting species name into two parts
uniform_genus_osps_rle_a1 = unlist(map(strsplit(uniform_species_osps_rle_a1, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_osps_rle_a1 = data.frame(cbind(uniform_species_osps_rle_a1, uniform_genus_osps_rle_a1))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_osps_rle_a1$Phylum <- 'X'
uniform_df_osps_rle_a1$Class <- 'X'
uniform_df_osps_rle_a1$Order <- 'X'
uniform_df_osps_rle_a1$Family <- 'X'
# re-name columns
colnames(uniform_df_osps_rle_a1) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_osps_rle_a1 <- uniform_df_osps_rle_a1[!duplicated(uniform_df_osps_rle_a1), ]
# re-name row names
rownames(uniform_df_osps_rle_a1) <- uniform_df_osps_rle_a1$otu
# remove non-numeric information from data frame
uniform_df_osps_rle_a1$otu <- NULL
# transpose data frame
uniform_df_osps_rle_a1 <- data.frame(t(uniform_df_osps_rle_a1))
# add columns with host-associated variables
uniform_df_osps_rle_a1$Age <- 'Age'
uniform_df_osps_rle_a1$BMI <- 'BMI'
uniform_df_osps_rle_a1$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_osps_rle_a1 <- data.frame(t(uniform_df_osps_rle_a1))
# clean up the row names
rownames(uniform_df_osps_rle_a1) <- str_replace(rownames(uniform_df_osps_rle_a1), '\\.', ' ')

# subset osps database based on abundance estimations (15th abundance percentile)
ds_osps_rle_all_a1 <- subset(ds_osps_rle, rownames(ds_osps_rle) %in% rownames(uniform_df_osps_rle_a1))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_osps_rle_all_a1 <- data.frame(t(ds_osps_rle_all_a1))
# order osps data frame by patient id (alphabetically)
ds_osps_rle_all_a1 <- ds_osps_rle_all_a1[order(rownames(ds_osps_rle_all_a1)),]

# add host-associated variables
ds_osps_rle_all_a1$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_osps_rle_all_a1$Age <- ifelse(ds_osps_rle_all_a1$Age_1 == '0 years', 0, ifelse(ds_osps_rle_all_a1$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_osps_rle_all_a1$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_osps_rle_all_a1$BMI <- md_sub$BMI
# remove age character column
ds_osps_rle_all_a1$Age_1 <- NULL

# re-arrange data frames
ds_osps_rle_all_a1_t <- data.frame(t(ds_osps_rle_all_a1))
md_sub_2_osps_rle_a1 <- data.frame(t(md_sub))
md_sub_3_osps_rle_a1 <- data.frame(t(md_sub_2_osps_rle_a1))

# convert data frames to matrix 
otu_mat_osps_rle_a1 <- as.matrix(ds_osps_rle_all_a1_t)
tax_mat_osps_rle_a1 <- as.matrix(uniform_df_osps_rle_a1)
# convert osps to otu_table (phyloseq object)
OTU_osps_rle_a1 = otu_table(otu_mat_osps_rle_a1, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_osps_rle_a1) <- str_replace(rownames(OTU_osps_rle_a1), '\\.', ' ')
# make tax table (phyloseq object)
TAX_osps_rle_a1 = tax_table(tax_mat_osps_rle_a1)
# make meta data table (phyloseq object)
samples_osps_rle_a1 = sample_data(md_sub_3_osps_rle_a1)
# generate final phyloseq object
erie.merge_osps_rle_a1 <- phyloseq(OTU_osps_rle_a1, TAX_osps_rle_a1, samples_osps_rle_a1)

# transpose data frame
predictors_osps_rle_a1 <- t(otu_table(erie.merge_osps_rle_a1))
# obtain response variable (CF vs healthy)
response_osps_rle_a1 <- as.factor(sample_data(erie.merge_osps_rle_a1)$State)
# add response variable to data frame
rf.data_osps_rle_a1 <- data.frame(response_osps_rle_a1, predictors_osps_rle_a1)
# determine mtry for random forest
sample_val_osps_rle_a1 = sqrt(ncol(rf.data_osps_rle_a1))
# generate three empty variables
rf_osps_rle_a1_pred_final = NULL
imp.sort.gini_osps_rle_a1_final = NULL
error_rate_osps_rle_a1_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_osps_rle_a1 <- randomForest(response_osps_rle_a1~., data = rf.data_osps_rle_a1, 
                                            ntree = 80, 
                                            mtry=sample_val_osps_rle_a1, 
                                            importance=TRUE,
                                            na.action = na.roughfix)
  
  # store error rate locally
  error_rate_osps_rle_a1 <- erie.classify_osps_rle_a1$err.rate
  # make data frame of error rate
  error_rate_osps_rle_a1 <- as.data.frame(error_rate_osps_rle_a1)
  # add meta data
  error_rate_osps_rle_a1$Database <- 'One strain per species'
  error_rate_osps_rle_a1$Normalisation <- 'RLE'
  error_rate_osps_rle_a1$Threshold <- '15th percentile'
  error_rate_osps_rle_a1$Seed <- i
  # re-name columns
  colnames(error_rate_osps_rle_a1) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_osps_rle_a1_final <- rbind(error_rate_osps_rle_a1_final, error_rate_osps_rle_a1)
  # obtain class probability
  rf_osps_rle_a1_pred <- predict(erie.classify_osps_rle_a1, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_osps_rle_a1_pred <- as.data.frame(rf_osps_rle_a1_pred)
  # add meta data
  rf_osps_rle_a1_pred$State <- md$State
  rf_osps_rle_a1_pred$Database <- 'One strain per species'
  rf_osps_rle_a1_pred$Normalisation <- 'RLE'
  rf_osps_rle_a1_pred$Threshold <- '15th percentile'
  rf_osps_rle_a1_pred$Seed <- i
  rownames(rf_osps_rle_a1_pred) <- NULL
  # transfer data frame to global environment
  rf_osps_rle_a1_pred_final <- rbind(rf_osps_rle_a1_pred_final, rf_osps_rle_a1_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_osps_rle_a1 <- Boruta(response_osps_rle_a1~., data = rf.data_osps_rle_a1, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_osps_rle_a1_df <- as.data.frame(boruta_osps_rle_a1$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_osps_rle_a1 <- importance(erie.classify_osps_rle_a1)
  imp_osps_rle_a1 <- data.frame(predictors = rownames(imp_osps_rle_a1), imp_osps_rle_a1)
  # add meta data
  imp_osps_rle_a1$Boruta_name <- rownames(boruta_osps_rle_a1_df)
  imp_osps_rle_a1$Boruta_predict <- boruta_osps_rle_a1_df$`boruta_osps_rle_a1$finalDecision`
  imp_osps_rle_a1_sub <- subset(imp_osps_rle_a1, MeanDecreaseAccuracy > 0.0)
  imp_osps_rle_a1_sub$predictors <- str_replace(imp_osps_rle_a1_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-associated
  imp_osps_rle_a1_sub$species_type <- with(imp_osps_rle_a1_sub,
                                           ifelse(predictors %in% f_rare_species_osps_rle_a1, 
                                                  'Rare species', 
                                                  ifelse(predictors %in% f_core_species_osps_rle_a1, 
                                                         'Core species', 
                                                         ifelse(predictors %in% fluctuating_rare_osps_rle_a1, 
                                                                'Rare species',
                                                                ifelse(predictors %in% fluctuating_core_osps_rle_a1, 
                                                                       'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_osps_rle_a1_sub$species_type_2 <- with(imp_osps_rle_a1_sub,
                                             ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                    ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                           ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                                  ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_osps_rle_a1_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                     'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                     'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_osps_rle_a1_sub <- subset(imp_osps_rle_a1_sub, Boruta_predict != "Rejected")
  imp_osps_rle_a1_short <- ddply(imp_osps_rle_a1_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_osps_rle_a1_short$Seed <- i
  imp_osps_rle_a1_short$Threshold <- '15th percentile'
  imp_osps_rle_a1_short$Normalisation <- 'RLE'
  # store short table globally
  imp.sort.gini_osps_rle_a1_final <- rbind(imp.sort.gini_osps_rle_a1_final, imp_osps_rle_a1_short)
}


# random forest, osps, RLE, 25th abundance percentile
# define the non-persisting rare species (only found in one or two age groups but not always)
fluctuating_rare_osps_rle_a2 <- c(all_rare_a2_osps_rle, f_rare_species_osps_rle_a2)
# remove duplicate entries
fluctuating_rare_osps_rle_a2 <- fluctuating_rare_osps_rle_a2[!duplicated(fluctuating_rare_osps_rle_a2)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_osps_rle_a2 <- c(all_core_a2_osps_rle, f_core_species_osps_rle_a2)
# remove duplicate entries
fluctuating_core_osps_rle_a2 <- fluctuating_core_osps_rle_a2[!duplicated(fluctuating_core_osps_rle_a2)]

# get all species found in children
uniform_species_osps_rle_a2 = c(all_rare_a2_osps_rle, all_core_a2_osps_rle)
# extract genus information by splitting species name into two parts
uniform_genus_osps_rle_a2 = unlist(map(strsplit(uniform_species_osps_rle_a2, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_osps_rle_a2 = data.frame(cbind(uniform_species_osps_rle_a2, uniform_genus_osps_rle_a2))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_osps_rle_a2$Phylum <- 'X'
uniform_df_osps_rle_a2$Class <- 'X'
uniform_df_osps_rle_a2$Order <- 'X'
uniform_df_osps_rle_a2$Family <- 'X'
# re-name columns
colnames(uniform_df_osps_rle_a2) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_osps_rle_a2 <- uniform_df_osps_rle_a2[!duplicated(uniform_df_osps_rle_a2), ]
# re-name row names
rownames(uniform_df_osps_rle_a2) <- uniform_df_osps_rle_a2$otu
# remove non-numeric information from data frame
uniform_df_osps_rle_a2$otu <- NULL
# transpose data frame
uniform_df_osps_rle_a2 <- data.frame(t(uniform_df_osps_rle_a2))
# add columns with host-associated variables
uniform_df_osps_rle_a2$Age <- 'Age'
uniform_df_osps_rle_a2$BMI <- 'BMI'
uniform_df_osps_rle_a2$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_osps_rle_a2 <- data.frame(t(uniform_df_osps_rle_a2))
# clean up the row names
rownames(uniform_df_osps_rle_a2) <- str_replace(rownames(uniform_df_osps_rle_a2), '\\.', ' ')

# subset osps database based on abundance estimations (25th abundance percentile)
ds_osps_rle_all_a2 <- subset(ds_osps_rle, rownames(ds_osps_rle) %in% rownames(uniform_df_osps_rle_a2))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_osps_rle_all_a2 <- data.frame(t(ds_osps_rle_all_a2))
# order osps data frame by patient id (alphabetically)
ds_osps_rle_all_a2 <- ds_osps_rle_all_a2[order(rownames(ds_osps_rle_all_a2)),]

# add host-associated variables
ds_osps_rle_all_a2$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_osps_rle_all_a2$Age <- ifelse(ds_osps_rle_all_a2$Age_1 == '0 years', 0, ifelse(ds_osps_rle_all_a2$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_osps_rle_all_a2$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_osps_rle_all_a2$BMI <- md_sub$BMI
# remove age character column
ds_osps_rle_all_a2$Age_1 <- NULL

# re-arrange data frames
ds_osps_rle_all_a2_t <- data.frame(t(ds_osps_rle_all_a2))
md_sub_2_osps_rle_a2 <- data.frame(t(md_sub))
md_sub_3_osps_rle_a2 <- data.frame(t(md_sub_2_osps_rle_a2))

# convert data frames to matrix 
otu_mat_osps_rle_a2 <- as.matrix(ds_osps_rle_all_a2_t)
tax_mat_osps_rle_a2 <- as.matrix(uniform_df_osps_rle_a2)
# convert osps_rle to otu_table (phyloseq object)
OTU_osps_rle_a2 = otu_table(otu_mat_osps_rle_a2, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_osps_rle_a2) <- str_replace(rownames(OTU_osps_rle_a2), '\\.', ' ')
# make tax table (phyloseq object)
TAX_osps_rle_a2 = tax_table(tax_mat_osps_rle_a2)
# make meta data table (phyloseq object)
samples_osps_rle_a2 = sample_data(md_sub_3_osps_rle_a2)
# generate final phyloseq object
erie.merge_osps_rle_a2 <- phyloseq(OTU_osps_rle_a2, TAX_osps_rle_a2, samples_osps_rle_a2)

# transpose data frame
predictors_osps_rle_a2 <- t(otu_table(erie.merge_osps_rle_a2))
# obtain response variable (CF vs healthy)
response_osps_rle_a2 <- as.factor(sample_data(erie.merge_osps_rle_a2)$State)
# add response variable to data frame
rf.data_osps_rle_a2 <- data.frame(response_osps_rle_a2, predictors_osps_rle_a2)
# determine mtry for random forest
sample_val_osps_rle_a2 = sqrt(ncol(rf.data_osps_rle_a2))
# generate three empty variables
rf_osps_rle_a2_pred_final = NULL
imp.sort.gini_osps_rle_a2_final = NULL
error_rate_osps_rle_a2_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_osps_rle_a2 <- randomForest(response_osps_rle_a2~., data = rf.data_osps_rle_a2, 
                                            ntree = 80, 
                                            mtry=sample_val_osps_rle_a2, 
                                            importance=TRUE,
                                            na.action = na.roughfix)
  
  # store error rate locally
  error_rate_osps_rle_a2 <- erie.classify_osps_rle_a2$err.rate
  # make data frame of error rate
  error_rate_osps_rle_a2 <- as.data.frame(error_rate_osps_rle_a2)
  # add meta data
  error_rate_osps_rle_a2$Database <- 'One strain per species'
  error_rate_osps_rle_a2$Normalisation <- 'RLE'
  error_rate_osps_rle_a2$Threshold <- '25th percentile'
  error_rate_osps_rle_a2$Seed <- i
  # re-name columns
  colnames(error_rate_osps_rle_a2) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_osps_rle_a2_final <- rbind(error_rate_osps_rle_a2_final, error_rate_osps_rle_a2)
  # obtain class probability
  rf_osps_rle_a2_pred <- predict(erie.classify_osps_rle_a2, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_osps_rle_a2_pred <- as.data.frame(rf_osps_rle_a2_pred)
  # add meta data
  rf_osps_rle_a2_pred$State <- md$State
  rf_osps_rle_a2_pred$Database <- 'One strain per species'
  rf_osps_rle_a2_pred$Normalisation <- 'RLE'
  rf_osps_rle_a2_pred$Threshold <- '25th percentile'
  rf_osps_rle_a2_pred$Seed <- i
  rownames(rf_osps_rle_a2_pred) <- NULL
  # transfer data frame to global environment
  rf_osps_rle_a2_pred_final <- rbind(rf_osps_rle_a2_pred_final, rf_osps_rle_a2_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_osps_rle_a2 <- Boruta(response_osps_rle_a2~., data = rf.data_osps_rle_a2, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_osps_rle_a2_df <- as.data.frame(boruta_osps_rle_a2$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_osps_rle_a2 <- importance(erie.classify_osps_rle_a2)
  imp_osps_rle_a2 <- data.frame(predictors = rownames(imp_osps_rle_a2), imp_osps_rle_a2)
  # add meta data
  imp_osps_rle_a2$Boruta_name <- rownames(boruta_osps_rle_a2_df)
  imp_osps_rle_a2$Boruta_predict <- boruta_osps_rle_a2_df$`boruta_osps_rle_a2$finalDecision`
  imp_osps_rle_a2_sub <- subset(imp_osps_rle_a2, MeanDecreaseAccuracy > 0.0)
  imp_osps_rle_a2_sub$predictors <- str_replace(imp_osps_rle_a2_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-assosciated
  imp_osps_rle_a2_sub$species_type <- with(imp_osps_rle_a2_sub,
                                           ifelse(predictors %in% f_rare_species_osps_rle_a2, 
                                                  'Rare species', 
                                                  ifelse(predictors %in% f_core_species_osps_rle_a2, 
                                                         'Core species', 
                                                         ifelse(predictors %in% fluctuating_rare_osps_rle_a2, 
                                                                'Rare species',
                                                                ifelse(predictors %in% fluctuating_core_osps_rle_a2, 
                                                                       'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_osps_rle_a2_sub$species_type_2 <- with(imp_osps_rle_a2_sub,
                                             ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                    ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                           ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                                  ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_osps_rle_a2_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                     'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                     'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_osps_rle_a2_sub <- subset(imp_osps_rle_a2_sub, Boruta_predict != "Rejected")
  imp_osps_rle_a2_short <- ddply(imp_osps_rle_a2_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_osps_rle_a2_short$Seed <- i
  imp_osps_rle_a2_short$Threshold <- '25th percentile'
  imp_osps_rle_a2_short$Normalisation <- 'RLE'
  # store short table globally
  imp.sort.gini_osps_rle_a2_final <- rbind(imp.sort.gini_osps_rle_a2_final, imp_osps_rle_a2_short)
}


# random forest, osps, RLE, 35th abundance percentile
# define the non-persisting rare species (only found in one or two age groups but not always)
fluctuating_rare_osps_rle_a3 <- c(all_rare_a3_osps_rle, f_rare_species_osps_rle_a3)
# remove duplicate entries
fluctuating_rare_osps_rle_a3 <- fluctuating_rare_osps_rle_a3[!duplicated(fluctuating_rare_osps_rle_a3)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_osps_rle_a3 <- c(all_core_a3_osps_rle, f_core_species_osps_rle_a3)
# remove duplicate entries
fluctuating_core_osps_rle_a3 <- fluctuating_core_osps_rle_a3[!duplicated(fluctuating_core_osps_rle_a3)]

# get all species found in children
uniform_species_osps_rle_a3 = c(all_rare_a3_osps_rle, all_core_a3_osps_rle)
# extract genus information by splitting species name into two parts
uniform_genus_osps_rle_a3 = unlist(map(strsplit(uniform_species_osps_rle_a3, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_osps_rle_a3 = data.frame(cbind(uniform_species_osps_rle_a3, uniform_genus_osps_rle_a3))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_osps_rle_a3$Phylum <- 'X'
uniform_df_osps_rle_a3$Class <- 'X'
uniform_df_osps_rle_a3$Order <- 'X'
uniform_df_osps_rle_a3$Family <- 'X'
# re-name columns
colnames(uniform_df_osps_rle_a3) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_osps_rle_a3 <- uniform_df_osps_rle_a3[!duplicated(uniform_df_osps_rle_a3), ]
# re-name row names
rownames(uniform_df_osps_rle_a3) <- uniform_df_osps_rle_a3$otu
# remove non-numeric information from data frame
uniform_df_osps_rle_a3$otu <- NULL
# transpose data frame
uniform_df_osps_rle_a3 <- data.frame(t(uniform_df_osps_rle_a3))
# add columns with host-associated variables
uniform_df_osps_rle_a3$Age <- 'Age'
uniform_df_osps_rle_a3$BMI <- 'BMI'
uniform_df_osps_rle_a3$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_osps_rle_a3 <- data.frame(t(uniform_df_osps_rle_a3))
# clean up the row names
rownames(uniform_df_osps_rle_a3) <- str_replace(rownames(uniform_df_osps_rle_a3), '\\.', ' ')

# subset osps database based on abundance estimations (35th abundance percentile)
ds_osps_rle_all_a3 <- subset(ds_osps_rle, rownames(ds_osps_rle) %in% rownames(uniform_df_osps_rle_a3))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_osps_rle_all_a3 <- data.frame(t(ds_osps_rle_all_a3))
# order osps data frame by patient id (alphabetically)
ds_osps_rle_all_a3 <- ds_osps_rle_all_a3[order(rownames(ds_osps_rle_all_a3)),]

# add host-associated variables
ds_osps_rle_all_a3$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_osps_rle_all_a3$Age <- ifelse(ds_osps_rle_all_a3$Age_1 == '0 years', 0, ifelse(ds_osps_rle_all_a3$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_osps_rle_all_a3$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_osps_rle_all_a3$BMI <- md_sub$BMI
# remove age character column
ds_osps_rle_all_a3$Age_1 <- NULL

# re-arrange data frames
ds_osps_rle_all_a3_t <- data.frame(t(ds_osps_rle_all_a3))
md_sub_2_osps_rle_a3 <- data.frame(t(md_sub))
md_sub_3_osps_rle_a3 <- data.frame(t(md_sub_2_osps_rle_a3))

# convert data frames to matrix 
otu_mat_osps_rle_a3 <- as.matrix(ds_osps_rle_all_a3_t)
tax_mat_osps_rle_a3 <- as.matrix(uniform_df_osps_rle_a3)
# convert osps to otu_table (phyloseq object)
OTU_osps_rle_a3 = otu_table(otu_mat_osps_rle_a3, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_osps_rle_a3) <- str_replace(rownames(OTU_osps_rle_a3), '\\.', ' ')
# make tax table (phyloseq object)
TAX_osps_rle_a3 = tax_table(tax_mat_osps_rle_a3)
# make meta data table (phyloseq object)
samples_osps_rle_a3 = sample_data(md_sub_3_osps_rle_a3)
# generate final phyloseq object
erie.merge_osps_rle_a3 <- phyloseq(OTU_osps_rle_a3, TAX_osps_rle_a3, samples_osps_rle_a3)

# transpose data frame
predictors_osps_rle_a3 <- t(otu_table(erie.merge_osps_rle_a3))
# obtain response variable (CF vs healthy)
response_osps_rle_a3 <- as.factor(sample_data(erie.merge_osps_rle_a3)$State)
# add response variable to data frame
rf.data_osps_rle_a3 <- data.frame(response_osps_rle_a3, predictors_osps_rle_a3)
# determine mtry for random forest
sample_val_osps_rle_a3 = sqrt(ncol(rf.data_osps_rle_a3))
# generate three empty variables
rf_osps_rle_a3_pred_final = NULL
imp.sort.gini_osps_rle_a3_final = NULL
error_rate_osps_rle_a3_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_osps_rle_a3 <- randomForest(response_osps_rle_a3~., data = rf.data_osps_rle_a3, 
                                            ntree = 80, 
                                            mtry=sample_val_osps_rle_a3, 
                                            importance=TRUE,
                                            na.action = na.roughfix)
  
  # store error rate locally
  error_rate_osps_rle_a3 <- erie.classify_osps_rle_a3$err.rate
  # make data frame of error rate
  error_rate_osps_rle_a3 <- as.data.frame(error_rate_osps_rle_a3)
  # add meta data
  error_rate_osps_rle_a3$Database <- 'One strain per species'
  error_rate_osps_rle_a3$Normalisation <- 'RLE'
  error_rate_osps_rle_a3$Threshold <- '35th percentile'
  error_rate_osps_rle_a3$Seed <- i
  # re-name columns
  colnames(error_rate_osps_rle_a3) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_osps_rle_a3_final <- rbind(error_rate_osps_rle_a3_final, error_rate_osps_rle_a3)
  # obtain class probability
  rf_osps_rle_a3_pred <- predict(erie.classify_osps_rle_a3, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_osps_rle_a3_pred <- as.data.frame(rf_osps_rle_a3_pred)
  # add meta data
  rf_osps_rle_a3_pred$State <- md$State
  rf_osps_rle_a3_pred$Database <- 'One strain per species'
  rf_osps_rle_a3_pred$Normalisation <- 'RLE'
  rf_osps_rle_a3_pred$Threshold <- '35th percentile'
  rf_osps_rle_a3_pred$Seed <- i
  rownames(rf_osps_rle_a3_pred) <- NULL
  # transfer data frame to global environment
  rf_osps_rle_a3_pred_final <- rbind(rf_osps_rle_a3_pred_final, rf_osps_rle_a3_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_osps_rle_a3 <- Boruta(response_osps_rle_a3~., data = rf.data_osps_rle_a3, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_osps_rle_a3_df <- as.data.frame(boruta_osps_rle_a3$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_osps_rle_a3 <- importance(erie.classify_osps_rle_a3)
  imp_osps_rle_a3 <- data.frame(predictors = rownames(imp_osps_rle_a3), imp_osps_rle_a3)
  # add meta data
  imp_osps_rle_a3$Boruta_name <- rownames(boruta_osps_rle_a3_df)
  imp_osps_rle_a3$Boruta_predict <- boruta_osps_rle_a3_df$`boruta_osps_rle_a3$finalDecision`
  imp_osps_rle_a3_sub <- subset(imp_osps_rle_a3, MeanDecreaseAccuracy > 0.0)
  imp_osps_rle_a3_sub$predictors <- str_replace(imp_osps_rle_a3_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-assosciated
  imp_osps_rle_a3_sub$species_type <- with(imp_osps_rle_a3_sub,
                                           ifelse(predictors %in% f_rare_species_osps_rle_a3, 
                                                  'Rare species', 
                                                  ifelse(predictors %in% f_core_species_osps_rle_a3, 
                                                         'Core species', 
                                                         ifelse(predictors %in% fluctuating_rare_osps_rle_a3, 
                                                                'Rare species',
                                                                ifelse(predictors %in% fluctuating_core_osps_rle_a3, 
                                                                       'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_osps_rle_a3_sub$species_type_2 <- with(imp_osps_rle_a3_sub,
                                             ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                    ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                           ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                                  ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_osps_rle_a3_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                     'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                     'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_osps_rle_a3_sub <- subset(imp_osps_rle_a3_sub, Boruta_predict != "Rejected")
  imp_osps_rle_a3_short <- ddply(imp_osps_rle_a3_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_osps_rle_a3_short$Seed <- i
  imp_osps_rle_a3_short$Threshold <- '35th percentile'
  imp_osps_rle_a3_short$Normalisation <- 'RLE'
  # store short table globally
  imp.sort.gini_osps_rle_a3_final <- rbind(imp.sort.gini_osps_rle_a3_final, imp_osps_rle_a3_short)
}


# merge all information, random forest, osps, RLE, 15th abundance percentile, 25th abundance percentile, 35th abundance percentile
rf.imp.sort.osps_rle <- data.frame(rbind(imp.sort.gini_osps_rle_a1_final,imp.sort.gini_osps_rle_a2_final, imp.sort.gini_osps_rle_a3_final))
rf.imp.sort.osps_rle$Database <- 'One strain per species'
rf.pred.prob.osps_rle <- data.frame(rbind(rf_osps_rle_a1_pred_final,rf_osps_rle_a2_pred_final, rf_osps_rle_a3_pred_final))
rf.error.osps_rle <- data.frame(rbind(error_rate_osps_rle_a1_final,error_rate_osps_rle_a2_final, error_rate_osps_rle_a3_final))


# Continue with VST-normalised data
# random forest for VST normalised data, 15th abundance percentile
fluctuating_rare_osps_vst_a1 <- c(all_rare_a1_osps_vst, f_rare_species_osps_vst_a1)
# remove duplicate entries
fluctuating_rare_osps_vst_a1 <- fluctuating_rare_osps_vst_a1[!duplicated(fluctuating_rare_osps_vst_a1)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_osps_vst_a1 <- c(all_core_a1_osps_vst, f_core_species_osps_vst_a1)
# remove duplicate entries
fluctuating_core_osps_vst_a1 <- fluctuating_core_osps_vst_a1[!duplicated(fluctuating_core_osps_vst_a1)]

# get all species found in children
uniform_species_osps_vst_a1 = c(all_rare_a1_osps_vst, all_core_a1_osps_vst)
# extract genus information by splitting species name into two parts
uniform_genus_osps_vst_a1 = unlist(map(strsplit(uniform_species_osps_vst_a1, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_osps_vst_a1 = data.frame(cbind(uniform_species_osps_vst_a1, uniform_genus_osps_vst_a1))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_osps_vst_a1$Phylum <- 'X'
uniform_df_osps_vst_a1$Class <- 'X'
uniform_df_osps_vst_a1$Order <- 'X'
uniform_df_osps_vst_a1$Family <- 'X'
# re-name columns
colnames(uniform_df_osps_vst_a1) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_osps_vst_a1 <- uniform_df_osps_vst_a1[!duplicated(uniform_df_osps_vst_a1), ]
# re-name row names
rownames(uniform_df_osps_vst_a1) <- uniform_df_osps_vst_a1$otu
# remove non-numeric information from data frame
uniform_df_osps_vst_a1$otu <- NULL
# transpose data frame
uniform_df_osps_vst_a1 <- data.frame(t(uniform_df_osps_vst_a1))
# add columns with host-associated variables
uniform_df_osps_vst_a1$Age <- 'Age'
uniform_df_osps_vst_a1$BMI <- 'BMI'
uniform_df_osps_vst_a1$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_osps_vst_a1 <- data.frame(t(uniform_df_osps_vst_a1))
# clean up the row names
rownames(uniform_df_osps_vst_a1) <- str_replace(rownames(uniform_df_osps_vst_a1), '\\.', ' ')

# subset osps database based on abundance estimations (15th abundance percentile)
ds_osps_vst_all_a1 <- subset(ds_osps_vst, rownames(ds_osps_vst) %in% rownames(uniform_df_osps_vst_a1))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_osps_vst_all_a1 <- data.frame(t(ds_osps_vst_all_a1))
# order osps data frame by patient id (alphabetically)
ds_osps_vst_all_a1 <- ds_osps_vst_all_a1[order(rownames(ds_osps_vst_all_a1)),]

# add host-associated variables
ds_osps_vst_all_a1$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_osps_vst_all_a1$Age <- ifelse(ds_osps_vst_all_a1$Age_1 == '0 years', 0, ifelse(ds_osps_vst_all_a1$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_osps_vst_all_a1$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_osps_vst_all_a1$BMI <- md_sub$BMI
# remove age character column
ds_osps_vst_all_a1$Age_1 <- NULL

# re-arrange data frames
ds_osps_vst_all_a1_t <- data.frame(t(ds_osps_vst_all_a1))
md_sub_2_osps_vst_a1 <- data.frame(t(md_sub))
md_sub_3_osps_vst_a1 <- data.frame(t(md_sub_2_osps_vst_a1))

# convert data frames to matrix 
otu_mat_osps_vst_a1 <- as.matrix(ds_osps_vst_all_a1_t)
tax_mat_osps_vst_a1 <- as.matrix(uniform_df_osps_vst_a1)
# convert osps to otu_table (phyloseq object)
OTU_osps_vst_a1 = otu_table(otu_mat_osps_vst_a1, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_osps_vst_a1) <- str_replace(rownames(OTU_osps_vst_a1), '\\.', ' ')
# make tax table (phyloseq object)
TAX_osps_vst_a1 = tax_table(tax_mat_osps_vst_a1)
# make meta data table (phyloseq object)
samples_osps_vst_a1 = sample_data(md_sub_3_osps_vst_a1)
# generate final phyloseq object
erie.merge_osps_vst_a1 <- phyloseq(OTU_osps_vst_a1, TAX_osps_vst_a1, samples_osps_vst_a1)

# transpose data frame
predictors_osps_vst_a1 <- t(otu_table(erie.merge_osps_vst_a1))
# obtain response variable (CF vs healthy)
response_osps_vst_a1 <- as.factor(sample_data(erie.merge_osps_vst_a1)$State)
# add response variable to data frame
rf.data_osps_vst_a1 <- data.frame(response_osps_vst_a1, predictors_osps_vst_a1)
# determine mtry for random forest
sample_val_osps_vst_a1 = sqrt(ncol(rf.data_osps_vst_a1))
# generate three empty variables
rf_osps_vst_a1_pred_final = NULL
imp.sort.gini_osps_vst_a1_final = NULL
error_rate_osps_vst_a1_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_osps_vst_a1 <- randomForest(response_osps_vst_a1~., data = rf.data_osps_vst_a1, 
                                            ntree = 80, 
                                            mtry=sample_val_osps_vst_a1, 
                                            importance=TRUE,
                                            na.action = na.roughfix)
  
  # store error rate locally
  error_rate_osps_vst_a1 <- erie.classify_osps_vst_a1$err.rate
  # make data frame of error rate
  error_rate_osps_vst_a1 <- as.data.frame(error_rate_osps_vst_a1)
  # add meta data
  error_rate_osps_vst_a1$Database <- 'One strain per species'
  error_rate_osps_vst_a1$Normalisation <- 'VST'
  error_rate_osps_vst_a1$Threshold <- '15th percentile'
  error_rate_osps_vst_a1$Seed <- i
  # re-name columns
  colnames(error_rate_osps_vst_a1) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_osps_vst_a1_final <- rbind(error_rate_osps_vst_a1_final, error_rate_osps_vst_a1)
  # obtain class probability
  rf_osps_vst_a1_pred <- predict(erie.classify_osps_vst_a1, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_osps_vst_a1_pred <- as.data.frame(rf_osps_vst_a1_pred)
  # add meta data
  rf_osps_vst_a1_pred$State <- md$State
  rf_osps_vst_a1_pred$Database <- 'One strain per species'
  rf_osps_vst_a1_pred$Normalisation <- 'VST'
  rf_osps_vst_a1_pred$Threshold <- '15th percentile'
  rf_osps_vst_a1_pred$Seed <- i
  rownames(rf_osps_vst_a1_pred) <- NULL
  # transfer data frame to global environment
  rf_osps_vst_a1_pred_final <- rbind(rf_osps_vst_a1_pred_final, rf_osps_vst_a1_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_osps_vst_a1 <- Boruta(response_osps_vst_a1~., data = rf.data_osps_vst_a1, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_osps_vst_a1_df <- as.data.frame(boruta_osps_vst_a1$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_osps_vst_a1 <- importance(erie.classify_osps_vst_a1)
  imp_osps_vst_a1 <- data.frame(predictors = rownames(imp_osps_vst_a1), imp_osps_vst_a1)
  # add meta data
  imp_osps_vst_a1$Boruta_name <- rownames(boruta_osps_vst_a1_df)
  imp_osps_vst_a1$Boruta_predict <- boruta_osps_vst_a1_df$`boruta_osps_vst_a1$finalDecision`
  imp_osps_vst_a1_sub <- subset(imp_osps_vst_a1, MeanDecreaseAccuracy > 0.0)
  imp_osps_vst_a1_sub$predictors <- str_replace(imp_osps_vst_a1_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-associated
  imp_osps_vst_a1_sub$species_type <- with(imp_osps_vst_a1_sub,
                                           ifelse(predictors %in% f_rare_species_osps_vst_a1, 
                                                  'Rare species', 
                                                  ifelse(predictors %in% f_core_species_osps_vst_a1, 
                                                         'Core species', 
                                                         ifelse(predictors %in% fluctuating_rare_osps_vst_a1, 
                                                                'Rare species',
                                                                ifelse(predictors %in% fluctuating_core_osps_vst_a1, 
                                                                       'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_osps_vst_a1_sub$species_type_2 <- with(imp_osps_vst_a1_sub,
                                             ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                    ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                           ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                                  ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_osps_vst_a1_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                     'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                     'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_osps_vst_a1_sub <- subset(imp_osps_vst_a1_sub, Boruta_predict != "Rejected")
  imp_osps_vst_a1_short <- ddply(imp_osps_vst_a1_sub, 'Species_type_2', numcolwise(sum))
  
  # add meta data
  imp_osps_vst_a1_short$Seed <- i
  imp_osps_vst_a1_short$Threshold <- '15th percentile'
  imp_osps_vst_a1_short$Normalisation <- 'VST'
  # store short table globally
  imp.sort.gini_osps_vst_a1_final <- rbind(imp.sort.gini_osps_vst_a1_final, imp_osps_vst_a1_short)
}


# random forest, osps, VST, 25th abundance percentile
# define the non-persisting rare species (only found in one or two age groups but not always)
fluctuating_rare_osps_vst_a2 <- c(all_rare_a2_osps_vst, f_rare_species_osps_vst_a2)
# remove duplicate entries
fluctuating_rare_osps_vst_a2 <- fluctuating_rare_osps_vst_a2[!duplicated(fluctuating_rare_osps_vst_a2)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_osps_vst_a2 <- c(all_core_a2_osps_vst, f_core_species_osps_vst_a2)
# remove duplicate entries
fluctuating_core_osps_vst_a2 <- fluctuating_core_osps_vst_a2[!duplicated(fluctuating_core_osps_vst_a2)]

# get all species found in children
uniform_species_osps_vst_a2 = c(all_rare_a2_osps_vst, all_core_a2_osps_vst)
# extract genus information by splitting species name into two parts
uniform_genus_osps_vst_a2 = unlist(map(strsplit(uniform_species_osps_vst_a2, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_osps_vst_a2 = data.frame(cbind(uniform_species_osps_vst_a2, uniform_genus_osps_vst_a2))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_osps_vst_a2$Phylum <- 'X'
uniform_df_osps_vst_a2$Class <- 'X'
uniform_df_osps_vst_a2$Order <- 'X'
uniform_df_osps_vst_a2$Family <- 'X'
# re-name columns
colnames(uniform_df_osps_vst_a2) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_osps_vst_a2 <- uniform_df_osps_vst_a2[!duplicated(uniform_df_osps_vst_a2), ]
# re-name row names
rownames(uniform_df_osps_vst_a2) <- uniform_df_osps_vst_a2$otu
# remove non-numeric information from data frame
uniform_df_osps_vst_a2$otu <- NULL
# transpose data frame
uniform_df_osps_vst_a2 <- data.frame(t(uniform_df_osps_vst_a2))
# add columns with host-associated variables
uniform_df_osps_vst_a2$Age <- 'Age'
uniform_df_osps_vst_a2$BMI <- 'BMI'
uniform_df_osps_vst_a2$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_osps_vst_a2 <- data.frame(t(uniform_df_osps_vst_a2))
# clean up the row names
rownames(uniform_df_osps_vst_a2) <- str_replace(rownames(uniform_df_osps_vst_a2), '\\.', ' ')

# subset osps database based on abundance estimations (25th abundance percentile)
ds_osps_vst_all_a2 <- subset(ds_osps_vst, rownames(ds_osps_vst) %in% rownames(uniform_df_osps_vst_a2))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_osps_vst_all_a2 <- data.frame(t(ds_osps_vst_all_a2))
# order osps data frame by patient id (alphabetically)
ds_osps_vst_all_a2 <- ds_osps_vst_all_a2[order(rownames(ds_osps_vst_all_a2)),]

# add host-associated variables
ds_osps_vst_all_a2$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_osps_vst_all_a2$Age <- ifelse(ds_osps_vst_all_a2$Age_1 == '0 years', 0, ifelse(ds_osps_vst_all_a2$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_osps_vst_all_a2$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_osps_vst_all_a2$BMI <- md_sub$BMI
# remove age character column
ds_osps_vst_all_a2$Age_1 <- NULL

# re-arrange data frames
ds_osps_vst_all_a2_t <- data.frame(t(ds_osps_vst_all_a2))
md_sub_2_osps_vst_a2 <- data.frame(t(md_sub))
md_sub_3_osps_vst_a2 <- data.frame(t(md_sub_2_osps_vst_a2))

# convert data frames to matrix 
otu_mat_osps_vst_a2 <- as.matrix(ds_osps_vst_all_a2_t)
tax_mat_osps_vst_a2 <- as.matrix(uniform_df_osps_vst_a2)
# convert osps_vst to otu_table (phyloseq object)
OTU_osps_vst_a2 = otu_table(otu_mat_osps_vst_a2, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_osps_vst_a2) <- str_replace(rownames(OTU_osps_vst_a2), '\\.', ' ')
# make tax table (phyloseq object)
TAX_osps_vst_a2 = tax_table(tax_mat_osps_vst_a2)
# make meta data table (phyloseq object)
samples_osps_vst_a2 = sample_data(md_sub_3_osps_vst_a2)
# generate final phyloseq object
erie.merge_osps_vst_a2 <- phyloseq(OTU_osps_vst_a2, TAX_osps_vst_a2, samples_osps_vst_a2)

# transpose data frame
predictors_osps_vst_a2 <- t(otu_table(erie.merge_osps_vst_a2))
# obtain response variable (CF vs healthy)
response_osps_vst_a2 <- as.factor(sample_data(erie.merge_osps_vst_a2)$State)
# add response variable to data frame
rf.data_osps_vst_a2 <- data.frame(response_osps_vst_a2, predictors_osps_vst_a2)
# determine mtry for random forest
sample_val_osps_vst_a2 = sqrt(ncol(rf.data_osps_vst_a2))
# generate three empty variables
rf_osps_vst_a2_pred_final = NULL
imp.sort.gini_osps_vst_a2_final = NULL
error_rate_osps_vst_a2_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_osps_vst_a2 <- randomForest(response_osps_vst_a2~., data = rf.data_osps_vst_a2, 
                                            ntree = 80, 
                                            mtry=sample_val_osps_vst_a2, 
                                            importance=TRUE,
                                            na.action = na.roughfix)
  
  # store error rate locally
  error_rate_osps_vst_a2 <- erie.classify_osps_vst_a2$err.rate
  # make data frame of error rate
  error_rate_osps_vst_a2 <- as.data.frame(error_rate_osps_vst_a2)
  # add meta data
  error_rate_osps_vst_a2$Database <- 'One strain per species'
  error_rate_osps_vst_a2$Normalisation <- 'VST'
  error_rate_osps_vst_a2$Threshold <- '25th percentile'
  error_rate_osps_vst_a2$Seed <- i
  # re-name columns
  colnames(error_rate_osps_vst_a2) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_osps_vst_a2_final <- rbind(error_rate_osps_vst_a2_final, error_rate_osps_vst_a2)
  # obtain class probability
  rf_osps_vst_a2_pred <- predict(erie.classify_osps_vst_a2, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_osps_vst_a2_pred <- as.data.frame(rf_osps_vst_a2_pred)
  # add meta data
  rf_osps_vst_a2_pred$State <- md$State
  rf_osps_vst_a2_pred$Database <- 'One strain per species'
  rf_osps_vst_a2_pred$Normalisation <- 'VST'
  rf_osps_vst_a2_pred$Threshold <- '25th percentile'
  rf_osps_vst_a2_pred$Seed <- i
  rownames(rf_osps_vst_a2_pred) <- NULL
  # transfer data frame to global environment
  rf_osps_vst_a2_pred_final <- rbind(rf_osps_vst_a2_pred_final, rf_osps_vst_a2_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_osps_vst_a2 <- Boruta(response_osps_vst_a2~., data = rf.data_osps_vst_a2, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_osps_vst_a2_df <- as.data.frame(boruta_osps_vst_a2$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_osps_vst_a2 <- importance(erie.classify_osps_vst_a2)
  imp_osps_vst_a2 <- data.frame(predictors = rownames(imp_osps_vst_a2), imp_osps_vst_a2)
  # add meta data
  imp_osps_vst_a2$Boruta_name <- rownames(boruta_osps_vst_a2_df)
  imp_osps_vst_a2$Boruta_predict <- boruta_osps_vst_a2_df$`boruta_osps_vst_a2$finalDecision`
  imp_osps_vst_a2_sub <- subset(imp_osps_vst_a2, MeanDecreaseAccuracy > 0.0)
  imp_osps_vst_a2_sub$predictors <- str_replace(imp_osps_vst_a2_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-assosciated
  imp_osps_vst_a2_sub$species_type <- with(imp_osps_vst_a2_sub,
                                           ifelse(predictors %in% f_rare_species_osps_vst_a2, 
                                                  'Rare species', 
                                                  ifelse(predictors %in% f_core_species_osps_vst_a2, 
                                                         'Core species', 
                                                         ifelse(predictors %in% fluctuating_rare_osps_vst_a2, 
                                                                'Rare species',
                                                                ifelse(predictors %in% fluctuating_core_osps_vst_a2, 
                                                                       'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_osps_vst_a2_sub$species_type_2 <- with(imp_osps_vst_a2_sub,
                                             ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                    ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                           ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                                  ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_osps_vst_a2_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                     'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                     'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_osps_vst_a2_sub <- subset(imp_osps_vst_a2_sub, Boruta_predict != "Rejected")
  imp_osps_vst_a2_short <- ddply(imp_osps_vst_a2_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_osps_vst_a2_short$Seed <- i
  imp_osps_vst_a2_short$Threshold <- '25th percentile'
  imp_osps_vst_a2_short$Normalisation <- 'VST'
  # store short table globally
  imp.sort.gini_osps_vst_a2_final <- rbind(imp.sort.gini_osps_vst_a2_final, imp_osps_vst_a2_short)
}


# random forest, osps, VST, 35th abundance percentile
# define the non-persisting rare species (only found in one or two age groups but not always)
fluctuating_rare_osps_vst_a3 <- c(all_rare_a3_osps_vst, f_rare_species_osps_vst_a3)
# remove duplicate entries
fluctuating_rare_osps_vst_a3 <- fluctuating_rare_osps_vst_a3[!duplicated(fluctuating_rare_osps_vst_a3)]
# define the non-persisting core species (only found in one or two age groups but not always)
fluctuating_core_osps_vst_a3 <- c(all_core_a3_osps_vst, f_core_species_osps_vst_a3)
# remove duplicate entries
fluctuating_core_osps_vst_a3 <- fluctuating_core_osps_vst_a3[!duplicated(fluctuating_core_osps_vst_a3)]

# get all species found in children
uniform_species_osps_vst_a3 = c(all_rare_a3_osps_vst, all_core_a3_osps_vst)
# extract genus information by splitting species name into two parts
uniform_genus_osps_vst_a3 = unlist(map(strsplit(uniform_species_osps_vst_a3, split = ' '), 1))
# generate a dataframe of species and genus information
uniform_df_osps_vst_a3 = data.frame(cbind(uniform_species_osps_vst_a3, uniform_genus_osps_vst_a3))
# add pseudo information (taxa levels), which is required to make phyloseq object but not used for random forest analysis
uniform_df_osps_vst_a3$Phylum <- 'X'
uniform_df_osps_vst_a3$Class <- 'X'
uniform_df_osps_vst_a3$Order <- 'X'
uniform_df_osps_vst_a3$Family <- 'X'
# re-name columns
colnames(uniform_df_osps_vst_a3) <- c('otu', 'Genus', 'Phylum','Class', 'Order', 'Family')
# remove duplicate entries
uniform_df_osps_vst_a3 <- uniform_df_osps_vst_a3[!duplicated(uniform_df_osps_vst_a3), ]
# re-name row names
rownames(uniform_df_osps_vst_a3) <- uniform_df_osps_vst_a3$otu
# remove non-numeric information from data frame
uniform_df_osps_vst_a3$otu <- NULL
# transpose data frame
uniform_df_osps_vst_a3 <- data.frame(t(uniform_df_osps_vst_a3))
# add columns with host-associated variables
uniform_df_osps_vst_a3$Age <- 'Age'
uniform_df_osps_vst_a3$BMI <- 'BMI'
uniform_df_osps_vst_a3$Gender <- 'Gender'
# transpose data frame back to original space
uniform_df_osps_vst_a3 <- data.frame(t(uniform_df_osps_vst_a3))
# clean up the row names
rownames(uniform_df_osps_vst_a3) <- str_replace(rownames(uniform_df_osps_vst_a3), '\\.', ' ')

# subset osps database based on abundance estimations (35th abundance percentile)
ds_osps_vst_all_a3 <- subset(ds_osps_vst, rownames(ds_osps_vst) %in% rownames(uniform_df_osps_vst_a3))
# re-name metadata table, so that the original table is not altered
md_sub <- md
# order meta data frame by patient id (alphabetically)
md_sub <- md_sub[order(rownames(md_sub)),]
ds_osps_vst_all_a3 <- data.frame(t(ds_osps_vst_all_a3))
# order osps data frame by patient id (alphabetically)
ds_osps_vst_all_a3 <- ds_osps_vst_all_a3[order(rownames(ds_osps_vst_all_a3)),]

# add host-associated variables
ds_osps_vst_all_a3$Age_1 <- md_sub$AgeGroup
# convert age group information from character to factor
ds_osps_vst_all_a3$Age <- ifelse(ds_osps_vst_all_a3$Age_1 == '0 years', 0, ifelse(ds_osps_vst_all_a3$Age_1 == '1-3 years', 1, 2))
# convert gender from character to factor
ds_osps_vst_all_a3$Gender <- ifelse(md_sub$Gender == 'm', 0, 1)
# add BMI
ds_osps_vst_all_a3$BMI <- md_sub$BMI
# remove age character column
ds_osps_vst_all_a3$Age_1 <- NULL

# re-arrange data frames
ds_osps_vst_all_a3_t <- data.frame(t(ds_osps_vst_all_a3))
md_sub_2_osps_vst_a3 <- data.frame(t(md_sub))
md_sub_3_osps_vst_a3 <- data.frame(t(md_sub_2_osps_vst_a3))

# convert data frames to matrix 
otu_mat_osps_vst_a3 <- as.matrix(ds_osps_vst_all_a3_t)
tax_mat_osps_vst_a3 <- as.matrix(uniform_df_osps_vst_a3)
# convert osps to otu_table (phyloseq object)
OTU_osps_vst_a3 = otu_table(otu_mat_osps_vst_a3, taxa_are_rows = TRUE)
# clean row names
rownames(OTU_osps_vst_a3) <- str_replace(rownames(OTU_osps_vst_a3), '\\.', ' ')
# make tax table (phyloseq object)
TAX_osps_vst_a3 = tax_table(tax_mat_osps_vst_a3)
# make meta data table (phyloseq object)
samples_osps_vst_a3 = sample_data(md_sub_3_osps_vst_a3)
# generate final phyloseq object
erie.merge_osps_vst_a3 <- phyloseq(OTU_osps_vst_a3, TAX_osps_vst_a3, samples_osps_vst_a3)

# transpose data frame
predictors_osps_vst_a3 <- t(otu_table(erie.merge_osps_vst_a3))
# obtain response variable (CF vs healthy)
response_osps_vst_a3 <- as.factor(sample_data(erie.merge_osps_vst_a3)$State)
# add response variable to data frame
rf.data_osps_vst_a3 <- data.frame(response_osps_vst_a3, predictors_osps_vst_a3)
# determine mtry for random forest
sample_val_osps_vst_a3 = sqrt(ncol(rf.data_osps_vst_a3))
# generate three empty variables
rf_osps_vst_a3_pred_final = NULL
imp.sort.gini_osps_vst_a3_final = NULL
error_rate_osps_vst_a3_final = NULL
# loop over random forest to generate results for different seeds set
for (i in random_seeds){
  set.seed(i)
  # run random forest
  erie.classify_osps_vst_a3 <- randomForest(response_osps_vst_a3~., data = rf.data_osps_vst_a3, 
                                            ntree = 80, 
                                            mtry=sample_val_osps_vst_a3, 
                                            importance=TRUE,
                                            na.action = na.roughfix)
  
  # store error rate locally
  error_rate_osps_vst_a3 <- erie.classify_osps_vst_a3$err.rate
  # make data frame of error rate
  error_rate_osps_vst_a3 <- as.data.frame(error_rate_osps_vst_a3)
  # add meta data
  error_rate_osps_vst_a3$Database <- 'One strain per species'
  error_rate_osps_vst_a3$Normalisation <- 'VST'
  error_rate_osps_vst_a3$Threshold <- '35th percentile'
  error_rate_osps_vst_a3$Seed <- i
  # re-name columns
  colnames(error_rate_osps_vst_a3) <- c('OOB_all', 'Error_CF', 'Error_H', 'Database', 'Normalisation', 'Threshold', 'Seed')
  # transfer error rate data frame to global environment
  error_rate_osps_vst_a3_final <- rbind(error_rate_osps_vst_a3_final, error_rate_osps_vst_a3)
  # obtain class probability
  rf_osps_vst_a3_pred <- predict(erie.classify_osps_vst_a3, type='prob', norm.votes = TRUE, predict.all = TRUE)
  # store in data frame locally
  rf_osps_vst_a3_pred <- as.data.frame(rf_osps_vst_a3_pred)
  # add meta data
  rf_osps_vst_a3_pred$State <- md$State
  rf_osps_vst_a3_pred$Database <- 'One strain per species'
  rf_osps_vst_a3_pred$Normalisation <- 'VST'
  rf_osps_vst_a3_pred$Threshold <- '35th percentile'
  rf_osps_vst_a3_pred$Seed <- i
  rownames(rf_osps_vst_a3_pred) <- NULL
  # transfer data frame to global environment
  rf_osps_vst_a3_pred_final <- rbind(rf_osps_vst_a3_pred_final, rf_osps_vst_a3_pred)
  
  # confirm random forest with boruta
  # run boruta
  boruta_osps_vst_a3 <- Boruta(response_osps_vst_a3~., data = rf.data_osps_vst_a3, pValue = 0.05, mcAdj=TRUE)
  # store final decision in data frame
  boruta_osps_vst_a3_df <- as.data.frame(boruta_osps_vst_a3$finalDecision)
  
  # Make a data frame with predictor names and their importance
  imp_osps_vst_a3 <- importance(erie.classify_osps_vst_a3)
  imp_osps_vst_a3 <- data.frame(predictors = rownames(imp_osps_vst_a3), imp_osps_vst_a3)
  # add meta data
  imp_osps_vst_a3$Boruta_name <- rownames(boruta_osps_vst_a3_df)
  imp_osps_vst_a3$Boruta_predict <- boruta_osps_vst_a3_df$`boruta_osps_vst_a3$finalDecision`
  imp_osps_vst_a3_sub <- subset(imp_osps_vst_a3, MeanDecreaseAccuracy > 0.0)
  imp_osps_vst_a3_sub$predictors <- str_replace(imp_osps_vst_a3_sub$predictors, '\\.', ' ') 
  # get information which feature is core/rare species or host-associated
  imp_osps_vst_a3_sub$species_type <- with(imp_osps_vst_a3_sub,
                                           ifelse(predictors %in% f_rare_species_osps_vst_a3, 
                                                  'Rare species', 
                                                  ifelse(predictors %in% f_core_species_osps_vst_a3, 
                                                         'Core species', 
                                                         ifelse(predictors %in% fluctuating_rare_osps_vst_a3, 
                                                                'Rare species',
                                                                ifelse(predictors %in% fluctuating_core_osps_vst_a3, 
                                                                       'Core species', 'Host-associated')))))
  # add Boruta information to data frame
  imp_osps_vst_a3_sub$species_type_2 <- with(imp_osps_vst_a3_sub,
                                             ifelse(Boruta_predict == 'Rejected' & species_type == "Rare species", 'Random_rare',
                                                    ifelse(Boruta_predict == 'Rejected' & species_type == "Core species", 'Random_core',
                                                           ifelse(Boruta_predict == 'Rejected' & species_type == "Host-associated", 'Random_host',
                                                                  ifelse(Boruta_predict == 'Tentative', 'Tenative', species_type)))))
  # add column names
  colnames(imp_osps_vst_a3_sub) <- c('Predictors', 'CF', 'Healthy', 'MeanDecreaseAccuracy', 
                                     'MeanDecreaseGini', 'Boruta_name', 'Boruta_predict', 
                                     'Species_type', 'Species_type_2')
  # generate short table based on species type, concatenate and sum duplicate rows
  imp_osps_vst_a3_sub <- subset(imp_osps_vst_a3_sub, Boruta_predict != "Rejected")
  imp_osps_vst_a3_short <- ddply(imp_osps_vst_a3_sub, 'Species_type_2', numcolwise(sum))
  # add meta data
  imp_osps_vst_a3_short$Seed <- i
  imp_osps_vst_a3_short$Threshold <- '35th percentile'
  imp_osps_vst_a3_short$Normalisation <- 'VST'
  # store short table globally
  imp.sort.gini_osps_vst_a3_final <- rbind(imp.sort.gini_osps_vst_a3_final, imp_osps_vst_a3_short)
}

# merge all information, random forest, osps, VST, 15th abundance percentile, 25th abundance percentile, 35th abundance percentile
rf.imp.sort.osps_vst <- data.frame(rbind(imp.sort.gini_osps_vst_a1_final,imp.sort.gini_osps_vst_a2_final, imp.sort.gini_osps_vst_a3_final))
rf.imp.sort.osps_vst$Database <- 'One strain per species'
rf.pred.prob.osps_vst <- data.frame(rbind(rf_osps_vst_a1_pred_final,rf_osps_vst_a2_pred_final, rf_osps_vst_a3_pred_final))
rf.error.osps_vst <- data.frame(rbind(error_rate_osps_vst_a1_final,error_rate_osps_vst_a2_final, error_rate_osps_vst_a3_final))

############################################################################################################
# merge all information, random forest, osps, VST, BCPHC, RLE, 15th abundance, 25th abundance, 35th abundance
rf.imp.sort.osps.all <- data.frame(rbind(rf.imp.sort.osps, rf.imp.sort.osps_rle, rf.imp.sort.osps_vst))
rf.pred.prob.osps.all <- data.frame(rbind(rf.pred.prob.osps, rf.pred.prob.osps_rle, rf.pred.prob.osps_vst))
rf.error.osps.all <- data.frame(rbind(rf.error.osps,rf.error.osps_rle, rf.error.osps_vst))

# merge information of pangenome and one strain per species database outputs
rf.imp.sort.all  <- data.frame(rbind(rf.imp.sort.pangenome.all, rf.imp.sort.osps.all))
rf.imp.sort.all.nr <- subset(rf.imp.sort.all, Species_type_2 != "Tenative")

# convert classification information to class factor
rf.imp.sort.all.nr$Species_type_2 <- factor(rf.imp.sort.all.nr$Species_type_2, levels = c('Host-associated', 'Core species', 'Rare species'))
# rename normalisation if necessary
rf.imp.sort.all.nr$Normalisation_2 <- with(rf.imp.sort.all.nr, ifelse(Normalisation == 'BCPHC', 'BCPHC', ifelse(Normalisation == 'RLE', 'RLE', 'VST')))

# generate plot of mean decrease accuracy with statistics
rf.imp.sort.all.nr$merged_all <- paste(rf.imp.sort.all.nr$Normalisation_2, ";", rf.imp.sort.all.nr$Database, ";", rf.imp.sort.all.nr$Threshold)

# re-name column entries
rf.imp.sort.all.nr$merged_all <- str_replace_all(rf.imp.sort.all.nr$merged_all, " ;", ";")
rf.imp.sort.all.nr$merged_all <- str_replace_all(rf.imp.sort.all.nr$merged_all, "One strain per species", "DB1")
rf.imp.sort.all.nr$merged_all <- str_replace_all(rf.imp.sort.all.nr$merged_all, "Pangenome", "DB2")
rf.imp.sort.all.nr$merged_all <- str_replace_all(rf.imp.sort.all.nr$merged_all, "percentile", "-PCTL")
rf.imp.sort.all.nr$merged_all <- str_replace_all(rf.imp.sort.all.nr$merged_all, "15th -PCTL", "15th-PCTL")
rf.imp.sort.all.nr$merged_all <- str_replace_all(rf.imp.sort.all.nr$merged_all, "25th -PCTL", "25th-PCTL")
rf.imp.sort.all.nr$merged_all <- str_replace_all(rf.imp.sort.all.nr$merged_all, "35th -PCTL", "35th-PCTL")
rf.imp.sort.all.nr$Species_type_2 <- str_replace_all(rf.imp.sort.all.nr$Species_type_2, "Core species", "Core spp.")
rf.imp.sort.all.nr$Species_type_2 <- str_replace_all(rf.imp.sort.all.nr$Species_type_2, "Rare species", "Rare spp.")

# generate plot of mean decrease accuracy with statistics
plot_accuracy <- ggplot(rf.imp.sort.all.nr, aes(x=merged_all, y=MeanDecreaseAccuracy)) +
  geom_point(aes(colour=Species_type_2),position = position_jitterdodge(dodge.width = 0.5), alpha=0.8, size=0.6) +
  scale_color_manual(values = c('darkorange', 'black', 'blue')) +
  theme_bw(base_size = 12) +  xlab(' ') + ylab('Mean Decrease Accuracy') +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,110)) +
  theme(panel.grid = element_blank(), 
        legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size=12),
        axis.title = element_text(size=12), axis.text.x = element_text(size=12), axis.text.y = element_text(size=12)) + 
  coord_flip() + guides(colour=guide_legend(ncol=3, override.aes = list(size=3, alpha=1)))


# generate plot of mean decrease gini with statistics
plot_gini <- ggplot(rf.imp.sort.all.nr, aes(x=merged_all, y=MeanDecreaseGini)) +
  geom_point(aes(colour=Species_type_2),position = position_jitterdodge(dodge.width = 0.5), alpha=0.8, size=0.6) +
  scale_color_manual(values = c('darkorange', 'black', 'blue')) +
  theme_bw(base_size = 12) +  xlab(' ') + ylab('Mean Decrease Gini') +
  scale_y_continuous(breaks = c(0,20,40,60,80,100), limits = c(0,110)) +
  theme(panel.grid = element_blank(), 
        legend.title = element_blank(), legend.position = "bottom", legend.text = element_text(size=12),
        axis.title = element_text(size=12), axis.text = element_text(size=12), axis.text.y = element_blank()) + 
  coord_flip() + guides(colour=guide_legend(ncol=3, override.aes = list(size=2, alpha=1)))

# merge both plots
plot_rf <- ggarrange(plot_accuracy, plot_gini, common.legend = TRUE, legend='bottom', labels = c('A', 'B'), nrow=1, widths = c(1,0.65), font.label = list(size = 14, color = "black"))

############################################################################################################
# statistical comparison of mean decrease accuracy with variables
acc_compare <- compare_means(MeanDecreaseAccuracy ~ merged_all, data = rf.imp.sort.all.nr, group.by=c("Normalisation_2","Threshold", "Database"))
# store in table
acc_compare_table <- data.frame(acc_compare)

# statistical comparison of mean decrease gini with variables
gini_compare <- compare_means(MeanDecreaseGini ~ merged_all, data = rf.imp.sort.all.nr, group.by=c("Normalisation_2","Threshold", "Database"))
# store in table
gini_compare_table <- data.frame(gini_compare)

# merge statistics to one data table, which will be exported later on
compare_rf_table <- data.frame(rbind(acc_compare_table, gini_compare_table))
compare_rf_table$p.format <- NULL
compare_rf_table$p <- NULL
compare_rf_table$method <- NULL

# get effect size information, mean decrease gini per species type
rf_table_effsize_gini <- rf.imp.sort.all.nr %>% group_by(Normalisation_2, Threshold, Database) %>% wilcox_effsize(MeanDecreaseGini~merged_all, ci=TRUE)

# get effect size information, mean decrease accuracy per species type
rf_table_effsize_acc <- rf.imp.sort.all.nr %>% group_by(Normalisation_2,Threshold, Database) %>% wilcox_effsize(MeanDecreaseAccuracy~merged_all, ci=TRUE)

# store in large table
compare_rf_table_effsize <- data.frame(rbind(rf_table_effsize_gini, rf_table_effsize_acc))
# round effect size to two decimal places
compare_rf_table_effsize$effsize <- round(compare_rf_table_effsize$effsize, 2)
compare_rf_table_effsize$n2 <- NULL
compare_rf_table_effsize$n1 <- NULL
compare_rf_table_effsize$conf.low <- round(compare_rf_table_effsize$conf.low,2)
# merge data table
rf_statistcs <- merge(compare_rf_table, compare_rf_table_effsize, by=c(".y.", "group1", "group2", "Threshold", "Database", "Normalisation_2"))
# re-name variables
rf_statistcs$group1 <- ifelse(grepl("Core species", rf_statistcs$group1), "Core species", ifelse(grepl("Rare species", rf_statistcs$group1), "Rare species", "Host-associated"))
rf_statistcs$group2 <- ifelse(grepl("Core species", rf_statistcs$group2), "Core species", ifelse(grepl("Rare species", rf_statistcs$group2), "Rare species", "Host-associated"))

############################################################################################################
# get information on class errors
# bind pangenome and one-strain per species database
rf.error.all <- data.frame(rbind(rf.error.osps.all, rf.error.pangenome.all))
rf.error.all <- subset(rf.error.all, OOB_all > 0)

# mean OOB estimate of error (total)
mean_OOB_all <- mean(rf.error.all$OOB_all)
sd_OOB_all <- sd(rf.error.all$OOB_all)
mean_OOB_all # 0.11 
sd_OOB_all # 0.06

# mean class error (CF)
mean_OOB_CF <- mean(rf.error.all$Error_CF)
sd_OOB_CF <- sd(rf.error.all$Error_CF)
mean_OOB_CF # 0.09
sd_OOB_CF # 0.08

# mean class error (Healthy)
mean_OOB_H <- mean(rf.error.all$Error_H)
sd_OOB_H <- sd(rf.error.all$Error_H)
mean_OOB_H # 0.12
sd_OOB_H #0.07

############################################################################################################
# Final output, figures
ggsave(filename="output_figures/Supplementary_Figure_S1.pdf", line_plots, device="pdf", height=18, width=18, units="cm")
ggsave(filename="output_figures/Figure_01.pdf", vennDiagramsPlot_pangenome, device="pdf", height=16, width=16, units="cm")
ggsave(filename="output_figures/Figure_02.pdf", plot_rf, height=15, width=23, units="cm", device="pdf")

############################################################################################################
# Final output, tables
write.table(rf_statistcs, file='output_figures/rf_statistics.csv', sep=';', row.names = FALSE, col.names = TRUE)
write.table(venn_core_rare_stats, file='output_figures/venn_statistics.csv', sep=';', row.names = FALSE, col.names = TRUE)
write.table(venn_core_rare_stats_effsize, file='output_figures/venn_statistics_effsize.csv', sep=';', row.names = FALSE, col.names = TRUE)
write.table(background_table, file = 'input_files/taxonomic_data/background_species.csv', row.names = FALSE, col.names = TRUE, sep=';')

############################################################################################################
