# title: "Taxonomic analysis of the human airway rare species"
# author: "Marie-Madlen Pust"
# date: "17 04 2021"

# clean global environment
rm(list=ls())

# set working directory
setwd("C:/Users/marie/Desktop/Phd studies3/e_rare_biosphere/R")

# load packages
library('readr')
library('viridis')
library('dplyr')
library('stringr')
library('tidyr')
library('ggplot2')
library('factoextra')
library('gmodels')
library('ggpubr')
library('plyr')
library('purrr')
library('Hmisc')
library('reshape')
library('vegan')
library('magrittr')
library('scales')
library('grid')
library('reshape2')
library('phyloseq')
library('rcompanion')
library('SimilarityMeasures')
library('viridis')
library('randomForest')
library('knitr')
library('ggrepel')
library('forcats')
# remotes::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
library('pheatmap')
#library('devtools')
#install_github('knapply/homophily')
#install_github('dyerlab/popgraph')
#install_github('tomwenseleers/export')
library('homophily')
library('export')
library('popgraph')
library('qgraph')
library('sna')
library('tidygraph')
library('expm')
library('CINNA')
library('igraph')
library('NetSwan')
library('CINNA')
library('graphkernels')

# set global variables 
abund_spec = 95.0

# global functions
#function which will bootstrap the standard error of the median
bootstrap_median <- function(data, num) {
  resamples <- lapply(1:num, function(i) sample(data, replace=T))
  r.median <- sapply(resamples, median)
  std.err <- sqrt(var(r.median))
  list(std.err=std.err, resamples=resamples, medians=r.median)}

# Network statistics
all_indices <- function(g){
  res <- matrix(0,vcount(g),4)
  res[,1] <- igraph::degree(g)
  res[,2] <- igraph::betweenness(g)
  res[,3] <- igraph::closeness(g)
  res[,4] <- igraph::hub_score(g)$vector
  apply(res,2,function(x) round(x,8))}

# import data
# normalised count table
dataset <- read_delim(
  "data_input/merged_project3_rare_2.csv", ";", 
  escape_double = FALSE, trim_ws = TRUE, skip_empty_rows = TRUE)
dataset <- data.frame(dataset)
# rename columns
colnames(dataset) <- sub("_S.*","",colnames(dataset)) 
# ID is row name
rownames(dataset) <- dataset$Id
dataset$Id <- NULL
# keep only newborns/infants count data
keep_newborns <-
  c("EMCF01", "EMCF03", "EMCF04", "EMCF17", "EMCF19", "EMCF21", "EMCF22", "EMCF25", "EMCF26", 
    "EMCF30", "EMCF32", "EMCF33", "EMCF34", "EMCF35", "EMCF36", "EMCF37", "KGCF01", "KGCF05", 
    "KGCF19", "KGCF21", "KGCF25", "KGCF26", "KGCF27", "KGCF28", "KGCF29", "KGCF30", "KGCF33", 
    "KGCF35", "KGCF38", "KGCF40", "EMCF02", "EMCF05", "EMCF06", "EMCF09", "EMCF10", "EMCF13", 
    "EMCF14", "EMCF16", "EMCF18", "EMCF27", "EMCF28", "EMCF29", "EMCF31", "KGCF03", "KGCF13", 
    "KGCF15", "KGCF22", "KGCF24", "KGCF37", "KGCF50", "MCF04s1", "MCF05s1", "MCF06s1", "MCF09s1", 
    "MCF10s1", "MCF11s1", "MCF12s1", "KGCF02", "KGCF04", "KGCF07", "KGCF14", "KGCF16", "KGCF17", 
    "KGCF18", "KGCF32", "KGCF36", "KGCF39", "KGCF41", "KGCF42", "KGCF44", "KGCF45", "KGCF46", 
    "KGCF47", "KGCF48", "KGCF49", "KGCF51", "KGCF52", "KGCF53","KGCF55", "KGCF56", "KGCF57", 
    "KGCF58", "MCF01s1", "MCF02s1", "MCF07s1", "MCF08s1", "MCF13s1")
dataset_sub <- select(dataset, keep_newborns)

# import metadata
metadata <- read_delim(
  "data_input/spatial_metaData_2020_12.csv",  ";", 
  escape_double = FALSE, trim_ws = TRUE, skip_empty_rows = TRUE)
metadata <- data.frame(metadata)
rownames(metadata) <- metadata$Sample 
metadata$Sample <- NULL
metadata_sub <- subset(metadata, rownames(metadata) %in% keep_newborns)

# raspir tables
healthy_0_raspir <- read_csv("data_input/healthy_0_raspir.csv", col_names = FALSE)
healthy_0_raspir <- healthy_0_raspir$X1

healthy_A_raspir <- read_csv("data_input/healthy_1_3_raspir.csv", col_names = FALSE)
healthy_A_raspir <- healthy_A_raspir$X1

healthy_C_raspir <- read_csv("data_input/healthy_4_6_raspir.csv", col_names = FALSE)
healthy_C_raspir <- healthy_C_raspir$X1

cf_0_raspir <- read_csv("data_input/cf_0_raspir.csv", col_names = FALSE)
cf_0_raspir <- cf_0_raspir$X1

cf_A_raspir <- read_csv("data_input/cf_1_3_raspir.csv", col_names = FALSE)
cf_A_raspir <- cf_A_raspir$X1

cf_C_raspir <- read_csv("data_input/cf_4_6_raspir.csv", col_names = FALSE)
cf_C_raspir <- cf_C_raspir$X1

# Age group (0 years)
metaData_0 <- subset(metadata_sub, ageGroup == "0 years")
# Healthy
metaData_0_h <- subset(metaData_0, state == "Healthy")
count_0_h <- select(dataset_sub, rownames(metaData_0_h))
count_0_h <- subset(count_0_h, rownames(count_0_h) %in% healthy_0_raspir)
count_0_h$species <- rownames(count_0_h)
count_0_h$species <- substring(count_0_h$species, 13, nchar(count_0_h$species))
count_0_h$species <- str_replace(count_0_h$species, "1_", "")
count_0_h$species <- str_replace(count_0_h$species, "2_", "")
count_0_h$species2 <- map(strsplit(count_0_h$species, split = "_"), 1)
count_0_h$species3 <- map(strsplit(count_0_h$species, split = "_"), 2)
count_0_h$species4 <- paste(count_0_h$species2, count_0_h$species3, sep= "_")
count_0_h$species <- NULL
count_0_h$species2 <- NULL
count_0_h$species3 <- NULL
count_0_h <- ddply(count_0_h,"species4",numcolwise(sum))
rownames(count_0_h) <- count_0_h$species4
count_0_h$species4 <- NULL
sum_healthy_0 = colSums(count_0_h)
sum_all_healthy_0 = sum(sum_healthy_0)
count_0_h$abundance <- (rowSums(
  count_0_h[,1:ncol(count_0_h)]) / sum_all_healthy_0) * 100
# sort abundance decreasing
count_0_h <- count_0_h[with(count_0_h, order(-abundance)), ]
count_0_h$cumsum <- cumsum(count_0_h$abundance)
# table with abundant species
count_0_h_core <- subset(count_0_h, cumsum <= abund_spec)
count_0_h_core$cumsum <- NULL
count_0_h_core$abundance <- NULL
# table with rare species
count_0_h_rare <- subset(count_0_h, cumsum > abund_spec)
count_0_h_rare$cumsum <- NULL
count_0_h_rare$abundance <- NULL

#CF
metaData_0_cf <- subset(metaData_0, state == "CF")
count_0_cf <- select(dataset_sub, rownames(metaData_0_cf))
count_0_cf <- subset(count_0_cf, rownames(count_0_cf) %in% cf_0_raspir)
count_0_cf$species <- rownames(count_0_cf)
count_0_cf$species <- substring(count_0_cf$species, 13, nchar(count_0_cf$species))
count_0_cf$species <- str_replace(count_0_cf$species, "1_", "")
count_0_cf$species <- str_replace(count_0_cf$species, "2_", "")
count_0_cf$species2 <- map(strsplit(count_0_cf$species, split = "_"), 1)
count_0_cf$species3 <- map(strsplit(count_0_cf$species, split = "_"), 2)
count_0_cf$species4 <- paste(count_0_cf$species2, count_0_cf$species3, sep= "_")
count_0_cf$species <- NULL
count_0_cf$species2 <- NULL
count_0_cf$species3 <- NULL
count_0_cf <- ddply(count_0_cf,"species4",numcolwise(sum))
rownames(count_0_cf) <- count_0_cf$species4
count_0_cf$species4 <- NULL
sum_cf_0 = colSums(count_0_cf)
sum_all_cf_0 = sum(sum_cf_0)
count_0_cf$abundance <- (rowSums(count_0_cf[,1:ncol(count_0_cf)]) / sum_all_cf_0) * 100
# sort abundance decreasing
count_0_cf <- count_0_cf[with(count_0_cf, order(-abundance)), ]
count_0_cf$cumsum <- cumsum(count_0_cf$abundance)
# table with abundant species
count_0_cf_core <- subset(count_0_cf, cumsum <= abund_spec)
count_0_cf_core$cumsum <- NULL
count_0_cf_core$abundance <- NULL
# table with rare species
count_0_cf_rare <- subset(count_0_cf, cumsum > abund_spec)
count_0_cf_rare$cumsum <- NULL
count_0_cf_rare$abundance <- NULL

# Age group (1-3 years)
metaData_1_3 <- subset(metadata_sub, ageGroup == "1-3 years")
# Healthy
metaData_1_3_h <- subset(metaData_1_3, state == "Healthy")
count_1_3_h <- select(dataset_sub, rownames(metaData_1_3_h))
count_1_3_h <- subset(count_1_3_h, rownames(count_1_3_h) %in% healthy_A_raspir)
count_1_3_h$species <- rownames(count_1_3_h)
count_1_3_h$species <- substring(count_1_3_h$species, 13, nchar(count_1_3_h$species))
count_1_3_h$species <- str_replace(count_1_3_h$species, "1_", "")
count_1_3_h$species <- str_replace(count_1_3_h$species, "2_", "")
count_1_3_h$species2 <- map(strsplit(count_1_3_h$species, split = "_"), 1)
count_1_3_h$species3 <- map(strsplit(count_1_3_h$species, split = "_"), 2)
count_1_3_h$species4 <- paste(count_1_3_h$species2, count_1_3_h$species3, sep= "_")
count_1_3_h$species <- NULL
count_1_3_h$species2 <- NULL
count_1_3_h$species3 <- NULL
count_1_3_h <- ddply(count_1_3_h, "species4", numcolwise(sum))
rownames(count_1_3_h) <- count_1_3_h$species4
count_1_3_h$species4 <- NULL
sum_healthy_1_3 = colSums(count_1_3_h)
sum_all_healthy_1_3 = sum(sum_healthy_1_3)
count_1_3_h$abundance <- (rowSums(
  count_1_3_h[,1:ncol(count_1_3_h)]) / sum_all_healthy_1_3) * 100
# sort abundance decreasing
count_1_3_h <- count_1_3_h[with(count_1_3_h, order(-abundance)), ]
count_1_3_h$cumsum <- cumsum(count_1_3_h$abundance)
# table with abundant species
count_1_3_h_core <- subset(count_1_3_h, cumsum <= abund_spec)
count_1_3_h_core$cumsum <- NULL
count_1_3_h_core$abundance <- NULL
# table with rare species
count_1_3_h_rare <- subset(count_1_3_h, cumsum > abund_spec)
count_1_3_h_rare$cumsum <- NULL
count_1_3_h_rare$abundance <- NULL

# CF
metaData_1_3_cf <- subset(metaData_1_3, state == "CF")
count_1_3_cf <- select(dataset_sub, rownames(metaData_1_3_cf))
count_1_3_cf <- subset(count_1_3_cf, rownames(count_1_3_cf) %in% cf_A_raspir)
count_1_3_cf$species <- rownames(count_1_3_cf)
count_1_3_cf$species <- substring(count_1_3_cf$species, 13, nchar(count_1_3_cf$species))
count_1_3_cf$species <- str_replace(count_1_3_cf$species, "1_", "")
count_1_3_cf$species <- str_replace(count_1_3_cf$species, "2_", "")
count_1_3_cf$species2 <- map(strsplit(count_1_3_cf$species, split = "_"), 1)
count_1_3_cf$species3 <- map(strsplit(count_1_3_cf$species, split = "_"), 2)
count_1_3_cf$species4 <- paste(count_1_3_cf$species2, count_1_3_cf$species3, sep= "_")
count_1_3_cf$species <- NULL
count_1_3_cf$species2 <- NULL
count_1_3_cf$species3 <- NULL
count_1_3_cf <- ddply(count_1_3_cf, "species4", numcolwise(sum))
rownames(count_1_3_cf) <- count_1_3_cf$species4
count_1_3_cf$species4 <- NULL
sum_cf_1_3 = colSums(count_1_3_cf)
sum_all_cf_1_3 = sum(sum_cf_1_3)
count_1_3_cf$abundance <- (rowSums(
  count_1_3_cf[,1:ncol(count_1_3_cf)]) / sum_all_cf_1_3) * 100
# sort abundance decreasing
count_1_3_cf <- count_1_3_cf[with(count_1_3_cf, order(-abundance)), ]
count_1_3_cf$cumsum <- cumsum(count_1_3_cf$abundance)
# table with abundant species
count_1_3_cf_core <- subset(count_1_3_cf, cumsum <= abund_spec)
count_1_3_cf_core$cumsum <- NULL
count_1_3_cf_core$abundance <- NULL
# table with rare species
count_1_3_cf_rare <- subset(count_1_3_cf, cumsum > abund_spec)
count_1_3_cf_rare$cumsum <- NULL
count_1_3_cf_rare$abundance <- NULL

# Age group (4-6 years)
metaData_4_6 <- subset(metadata_sub, ageGroup == "4-6 years")
# Healthy
metaData_4_6_h <- subset(metaData_4_6, state == "Healthy")
count_4_6_h <- select(dataset_sub, rownames(metaData_4_6_h))
count_4_6_h <- subset(count_4_6_h, rownames(count_4_6_h) %in% healthy_C_raspir)
count_4_6_h$species <- rownames(count_4_6_h)
count_4_6_h$species <- substring(count_4_6_h$species, 13, nchar(count_4_6_h$species))
count_4_6_h$species <- str_replace(count_4_6_h$species, "1_", "")
count_4_6_h$species <- str_replace(count_4_6_h$species, "2_", "")
count_4_6_h$species2 <- map(strsplit(count_4_6_h$species, split = "_"), 1)
count_4_6_h$species3 <- map(strsplit(count_4_6_h$species, split = "_"), 2)
count_4_6_h$species4 <- paste(count_4_6_h$species2, count_4_6_h$species3, sep= "_")
count_4_6_h$species <- NULL
count_4_6_h$species2 <- NULL
count_4_6_h$species3 <- NULL
count_4_6_h <- ddply(count_4_6_h, "species4", numcolwise(sum))
rownames(count_4_6_h) <- count_4_6_h$species4
count_4_6_h$species4 <- NULL
sum_healthy_4_6 = colSums(count_4_6_h)
sum_all_healthy_4_6 = sum(sum_healthy_4_6)
count_4_6_h$abundance <- (rowSums(
  count_4_6_h[,1:ncol(count_4_6_h)]) / sum_all_healthy_4_6) * 100
# sort abundance decreasing
count_4_6_h <- count_4_6_h[with(count_4_6_h, order(-abundance)), ]
count_4_6_h$cumsum <- cumsum(count_4_6_h$abundance)
# table with abundant species
count_4_6_h_core <- subset(count_4_6_h, cumsum <= abund_spec)
count_4_6_h_core$cumsum <- NULL
count_4_6_h_core$abundance <- NULL
# table with rare species
count_4_6_h_rare <- subset(count_4_6_h, cumsum > abund_spec)
count_4_6_h_rare$cumsum <- NULL
count_4_6_h_rare$abundance <- NULL

# CF
metaData_4_6_cf <- subset(metaData_4_6, state == "CF")
count_4_6_cf <- select(dataset_sub, rownames(metaData_4_6_cf))
count_4_6_cf <- subset(count_4_6_cf, rownames(count_4_6_cf) %in% cf_C_raspir)
count_4_6_cf$species <- rownames(count_4_6_cf)
count_4_6_cf$species <- substring(count_4_6_cf$species, 13, nchar(count_4_6_cf$species))
count_4_6_cf$species <- str_replace(count_4_6_cf$species, "1_", "")
count_4_6_cf$species <- str_replace(count_4_6_cf$species, "2_", "")
count_4_6_cf$species2 <- map(strsplit(count_4_6_cf$species, split = "_"), 1)
count_4_6_cf$species3 <- map(strsplit(count_4_6_cf$species, split = "_"), 2)
count_4_6_cf$species4 <- paste(count_4_6_cf$species2, count_4_6_cf$species3, sep= "_")
count_4_6_cf$species <- NULL
count_4_6_cf$species2 <- NULL
count_4_6_cf$species3 <- NULL
count_4_6_cf <- ddply(count_4_6_cf, "species4", numcolwise(sum))
rownames(count_4_6_cf) <- count_4_6_cf$species4
count_4_6_cf$species4 <- NULL
sum_cf_4_6 = colSums(count_4_6_cf)
sum_all_cf_4_6 = sum(sum_cf_4_6)
count_4_6_cf$abundance <- (rowSums(
  count_4_6_cf[,1:ncol(count_4_6_cf)]) / sum_all_cf_4_6) * 100
# sort abundance decreasing
count_4_6_cf <- count_4_6_cf[with(count_4_6_cf, order(-abundance)), ]
count_4_6_cf$cumsum <- cumsum(count_4_6_cf$abundance)
# table with abundant species
count_4_6_cf_core <- subset(count_4_6_cf, cumsum <= abund_spec)
count_4_6_cf_core$cumsum <- NULL
count_4_6_cf_core$abundance <- NULL
# table with rare species
count_4_6_cf_rare <- subset(count_4_6_cf, cumsum > abund_spec)
count_4_6_cf_rare$cumsum <- NULL
count_4_6_cf_rare$abundance <- NULL

dataset_sub$species <- rownames(dataset_sub)
dataset_sub$species <- substring(dataset_sub$species, 13, nchar(dataset_sub$species))
dataset_sub$species <- str_replace(dataset_sub$species, "1_", "")
dataset_sub$species <- str_replace(dataset_sub$species, "2_", "")
dataset_sub$species2 <- map(strsplit(dataset_sub$species, split = "_"), 1)
dataset_sub$species3 <- map(strsplit(dataset_sub$species, split = "_"), 2)
dataset_sub$species4 <- paste(dataset_sub$species2, dataset_sub$species3, sep= "_")
dataset_sub$species <- NULL
dataset_sub$species2 <- NULL
dataset_sub$species3 <- NULL
dataset_sub <- ddply(dataset_sub, "species4", numcolwise(sum))
rownames(dataset_sub) <- dataset_sub$species4
dataset_sub$species4 <- NULL

healthy_rare <- c(rownames(count_0_h_rare), rownames(count_1_3_h_rare), rownames(count_4_6_h_rare))
healthy_rare <- healthy_rare[!duplicated(healthy_rare)]
healthy_core <- c(rownames(count_0_h_core), rownames(count_1_3_h_core), rownames(count_4_6_h_core))
healthy_core <- healthy_core[!duplicated(healthy_core)]

cf_rare <- c(rownames(count_0_cf_rare), rownames(count_1_3_cf_rare), rownames(count_4_6_cf_rare))
cf_rare <- cf_rare[!duplicated(cf_rare)]
cf_core <- c(rownames(count_0_cf_core), rownames(count_1_3_cf_core), rownames(count_4_6_cf_core))
cf_core <- cf_core[!duplicated(cf_core)]

all_rare <- c(healthy_rare, cf_rare)
all_core <- c(healthy_core, cf_core)

# obtain diversity information
dataset_species_rare <- subset(dataset_sub, rownames(dataset_sub) %in% all_rare)
dataset_species_rare <- dataset_species_rare[order(rownames(dataset_species_rare), decreasing = FALSE),]
dataset_species_rare <- data.frame(t(dataset_species_rare))
dataset_species_rare <- dataset_species_rare[order(rownames(dataset_species_rare), decreasing = FALSE),]

dataset_species_core <- subset(dataset_sub, rownames(dataset_sub) %in% all_core)
dataset_species_core <- dataset_species_core[order(rownames(dataset_species_core), decreasing = FALSE),]
dataset_species_core <- data.frame(t(dataset_species_core))
dataset_species_core <- dataset_species_core[order(rownames(dataset_species_core), decreasing = FALSE),]
metadata_sub <- metadata_sub[order(rownames(metadata_sub), decreasing = FALSE),]

metadata_sub$shannon_core <- vegan::diversity(dataset_species_core, index = "shannon")
metadata_sub$shannon_rare <- vegan::diversity(dataset_species_rare, index = "shannon")
metadata_sub$load_core <- rowSums(dataset_species_core)
metadata_sub$load_rare <- rowSums(dataset_species_rare)

metadata_sub$shannon_rare_cato <- ifelse(metadata_sub$shannon_rare > median(metadata_sub$shannon_core), "high", "low")
metadata_sub$shannon_core_cato <- ifelse(metadata_sub$shannon_core > median(metadata_sub$shannon_rare), "high", "low")
metadata_sub$load_rare_cato <- ifelse(metadata_sub$load_rare > median(metadata_sub$load_rare), "high", "low")
metadata_sub$load_core_cato <- ifelse(metadata_sub$load_core > median(metadata_sub$load_core), "high", "low")

# Venn diagrams
venn_all_rare_species_h <- data.frame(dataset_sub)
venn_all_rare_species_h$h_0_summary <- ifelse(rownames(venn_all_rare_species_h) %in% rownames(count_0_h_rare), 1, 0)
venn_all_rare_species_h$h_1_3_summary <- ifelse(rownames(venn_all_rare_species_h) %in% rownames(count_1_3_h_rare), 1, 0)
venn_all_rare_species_h$h_4_6_summary <- ifelse(rownames(venn_all_rare_species_h) %in% rownames(count_4_6_h_rare), 1, 0)
venn_all_rare_species_h <- venn_all_rare_species_h %>% select(h_0_summary, h_1_3_summary, h_4_6_summary)
venn_all_rare_species_h <- venn_all_rare_species_h[rowSums(venn_all_rare_species_h) > 0,]
background_rare_h <- venn_all_rare_species_h[rowSums(venn_all_rare_species_h)>2,]
non_persistent_rare_h <- venn_all_rare_species_h[rowSums(venn_all_rare_species_h)<3,]

venn_all_core_species_h <- data.frame(dataset_sub)
venn_all_core_species_h$h_0_summary <- ifelse(rownames(venn_all_core_species_h) %in% rownames(count_0_h_core), 1, 0)
venn_all_core_species_h$h_1_3_summary <- ifelse(rownames(venn_all_core_species_h) %in% rownames(count_1_3_h_core), 1, 0)
venn_all_core_species_h$h_4_6_summary <- ifelse(rownames(venn_all_core_species_h) %in% rownames(count_4_6_h_core), 1, 0)
venn_all_core_species_h <- venn_all_core_species_h %>% select(h_0_summary, h_1_3_summary, h_4_6_summary)
venn_all_core_species_h <- venn_all_core_species_h[rowSums(venn_all_core_species_h) > 0,]
background_core_h <- venn_all_core_species_h[rowSums(venn_all_core_species_h)>2,]
non_persistent_core_h <- venn_all_core_species_h[rowSums(venn_all_core_species_h)<3,]

venn_all_rare_species_cf <- data.frame(dataset_sub)
venn_all_rare_species_cf$cf_0_summary <- ifelse(rownames(venn_all_rare_species_cf) %in% rownames(count_0_cf_rare), 1, 0)
venn_all_rare_species_cf$cf_1_3_summary <- ifelse(rownames(venn_all_rare_species_cf) %in% rownames(count_1_3_cf_rare), 1, 0)
venn_all_rare_species_cf$cf_4_6_summary <- ifelse(rownames(venn_all_rare_species_cf) %in% rownames(count_4_6_cf_rare), 1, 0)
venn_all_rare_species_cf <- venn_all_rare_species_cf %>% select(cf_0_summary, cf_1_3_summary, cf_4_6_summary)
venn_all_rare_species_cf <- venn_all_rare_species_cf[rowSums(venn_all_rare_species_cf) > 0,]
background_rare_cf <- venn_all_rare_species_cf[rowSums(venn_all_rare_species_cf)>2,]
non_persistent_rare_cf <- venn_all_rare_species_cf[rowSums(venn_all_rare_species_cf)<3,]

venn_all_core_species_cf <- data.frame(dataset_sub)
venn_all_core_species_cf$cf_0_summary <- ifelse(rownames(venn_all_core_species_cf) %in% rownames(count_0_cf_core), 1, 0)
venn_all_core_species_cf$cf_1_3_summary <- ifelse(rownames(venn_all_core_species_cf) %in% rownames(count_1_3_cf_core), 1, 0)
venn_all_core_species_cf$cf_4_6_summary <- ifelse(rownames(venn_all_core_species_cf) %in% rownames(count_4_6_cf_core), 1, 0)
venn_all_core_species_cf <- venn_all_core_species_cf %>% select(cf_0_summary, cf_1_3_summary, cf_4_6_summary)
venn_all_core_species_cf <- venn_all_core_species_cf[rowSums(venn_all_core_species_cf) > 0,]
background_core_cf <- venn_all_core_species_cf[rowSums(venn_all_core_species_cf)>2,]
non_persistent_core_cf <- venn_all_core_species_cf[rowSums(venn_all_core_species_cf)<3,]

# venn diagram (rare species)
h_x_rare <- list(
  A = rownames(count_0_h_rare), 
  B = rownames(count_1_3_h_rare), 
  C = rownames(count_4_6_h_rare))

cf_x_rare <- list(
  A = rownames(count_0_cf_rare), 
  B = rownames(count_1_3_cf_rare), 
  C = rownames(count_4_6_cf_rare))

# venn diagram (core species)
h_x_core <- list(
  A = rownames(count_0_h_core), 
  B = rownames(count_1_3_h_core), 
  C = rownames(count_4_6_h_core))
cf_x_core <- list(
  A = rownames(count_0_cf_core), 
  B = rownames(count_1_3_cf_core), 
  C = rownames(count_4_6_cf_core))


# Venn diagram (rare species) ####
cf_rare_venn <- 
  ggVennDiagram(cf_x_rare, size=1.5, category.names = c(
    "0", "1-3", "4-6"), label = "both", percent_digits = 0) +
  theme_pubr(border=TRUE, legend="none") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  geom_label(aes(x=-2, y=9), label="CF, rare species")


h_rare_venn <-
  ggVennDiagram(h_x_rare, size=1.5, category.names = c(
    "0", "1-3", "4-6"), label = "both", percent_digits = 0) +
  theme_pubr(border=TRUE, legend="none") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  geom_label(aes(x=-1.5, y=9), label="Healthy, rare species")


cf_core_venn <- 
  ggVennDiagram(cf_x_core, size=1.5, category.names = c(
    "0", "1-3", "4-6"), label = "both", percent_digits = 0) +
  theme_pubr(border=TRUE, legend="none") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  geom_label(aes(x=-2, y=9), label="CF, core species")

h_core_venn <-
  ggVennDiagram(h_x_core, size=1.5, category.names = c(
    "0", "1-3", "4-6"), label = "both", percent_digits = 0) +
  theme_pubr(border=TRUE, legend="none") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  geom_label(aes(x=-1.3, y=9), label="Healthy, core species")


vennDiagramsPlot <-
  ggarrange(h_rare_venn, h_core_venn,
            cf_rare_venn, cf_core_venn, 
            labels = c("a", "b", "c", "d"), 
            nrow=2, ncol=2, legend = FALSE)


# venn diagram (from rare to core)
h_x_rare_to_core <- list(
  A = c(rownames(count_0_h_rare), rownames(count_1_3_h_rare)),
  B = rownames(count_4_6_h_core))
cf_x_rare_to_core <- list(
  A =  c(rownames(count_0_cf_rare), rownames(count_1_3_cf_rare)),
  B = rownames(count_4_6_cf_core))

# venn diagram (from core to rare)
h_x_core_to_rare <- list(
  A = c(rownames(count_0_h_core), rownames(count_1_3_h_core)),
  B = rownames(count_4_6_h_rare))
cf_x_core_to_rare <- list(
  A =  c(rownames(count_0_cf_core), rownames(count_1_3_cf_core)),
  B = rownames(count_4_6_cf_rare))


cf_rare_to_core_venn <- 
  ggVennDiagram(cf_x_rare_to_core, category.names = c(
    "0-3", "4-6"), label = "both", percent_digits = 0) +
  theme_pubr(border=TRUE, legend="none") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  geom_label(aes(x=-1.3, y=9), label="CF, rare to core") 

h_rare_to_core_venn <-
  ggVennDiagram(h_x_rare_to_core, category.names = c(
    "0-3", "4-6"), label = "both", percent_digits = 0) +
  theme_pubr(border=TRUE, legend="none") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  geom_label(aes(x=-1.3, y=9), label="Healthy, rare to core") 


cf_core_to_rare_venn <- 
  ggVennDiagram(cf_x_core_to_rare, category.names = c(
    "0-3", "4-6"), label = "both", percent_digits = 0) +
  theme_pubr(border=TRUE, legend="none") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  geom_label(aes(x=-1.3, y=9), label="CF, core to rare") 


h_core_to_rare_venn <-
  ggVennDiagram(h_x_core_to_rare, category.names = c(
    "0-3", "4-6"), label = "both", percent_digits = 0) +
  theme_pubr(border=TRUE, legend="none") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  geom_label(aes(x=-1.3, y=9), label="Healthy, core to rare") 

from_where_to_where <- ggarrange(h_rare_to_core_venn, cf_rare_to_core_venn,
                                 h_core_to_rare_venn, cf_core_to_rare_venn,
                                 labels = c("a", "b", "c", "d"), nrow=2, ncol=2,
                                 legend = FALSE)

# Obtain background core species (healty, cf)
f_core_species <- c(rownames(background_core_h), rownames(background_core_cf))
f_core_species <- f_core_species[!duplicated(f_core_species)]

# Obtain background rare species (healthy, cf)
f_rare_species <- c(rownames(background_rare_h), rownames(background_rare_cf))
f_rare_species <- f_rare_species[!duplicated(f_rare_species)]

heatmap_plot <- rbind(background_core_h, background_rare_h)
heatmap_plot$cf_0_summary <- ifelse(rownames(heatmap_plot) %in% rownames(background_rare_cf) |
                                      rownames(heatmap_plot) %in% rownames(background_core_cf), 1, 0)
heatmap_plot$cf_1_3_summary <- ifelse(rownames(heatmap_plot) %in% rownames(background_rare_cf) |
                                        rownames(heatmap_plot) %in% rownames(background_core_cf), 1, 0)
heatmap_plot$cf_4_6_summary <- ifelse(rownames(heatmap_plot) %in% rownames(background_rare_cf) |
                                        rownames(heatmap_plot) %in% rownames(background_core_cf), 1, 0)
heatmap_plot <- heatmap_plot[order(rownames(heatmap_plot)),]
heatmap_plot$blank <- 0
heatmap_plot$water <- 0
rownames(heatmap_plot) <- str_replace(rownames(heatmap_plot), "_", " ")

# print rownames italic
newnames <- lapply(rownames(heatmap_plot),function(x) bquote(italic(.(x))))
plot_pheatmap <- pheatmap(heatmap_plot, cluster_rows = FALSE, cluster_cols = FALSE,
                          color = c("azure", "cyan4"), legend = FALSE, cellwidth = 50,
                          cellheight = 15, angle_col = 0, gaps_col = c(1,2,3,4,5,6,7),
                          labels_col = c("H-0", "H-1-3", "H-4-6", "CF-0", "CF-1-3", "CF-4-6", "Blank", "Water"),
                          labels_row = as.expression(newnames))

# Where do the rare species come from?
cf_y_rare_spec_A <- c(
  "Upper respiratory tract",                                    # Neisseria subflava
  "Upper respiratory tract",                                    # Neisseria lactamica
  "Upper respiratory tract",                                    # Neisseria polysaccharea
  "Upper respiratory tract",                                    # Streptococcus koreensis
  "Upper respiratory tract",                                    # Gemella sanguinis
  "Various human body sites",                                   # Neisseria flavescens
  "Various human body sites",                                   # Veillonella dispar
  "Various human body sites",                                   # Haemophilus pittmaniae
  "Upper respiratory tract",                                    # Neisseria sicca
  "Environment, animal or food-associated")                     # Rodentibacter heylii

cf_y_rare_spec_B <- c(
  "Upper respiratory tract",                                    # Neisseria_lactamica
  "Various human body sites",                                   # Neisseria_flavescens
  "Upper respiratory tract",                                    # Neisseria_sicca
  "Environment, animal or food-associated",                     # Streptococcus_equinus
  "Upper respiratory tract",                                    # Capnocytophaga_gingivalis
  "Various human body sites",                                   # Fusobacterium_periodonticum
  "Upper respiratory tract",                                    # Leptotrichia_wadei
  "Upper respiratory tract",                                    # Mogibacterium_diversum
  "Various human body sites",                                   # Veillonella_parvula
  "Upper respiratory tract",                                    # Leptotrichia_shahii
  "Upper respiratory tract",                                    # Streptococcus_sanguinis
  "Various human body sites",                                   # Abiotrophia_defectiva
  "Various human body sites",                                   # Gemella_morbillorum
  "Upper respiratory tract",                                    # Streptococcus_vestibularis
  "Upper respiratory tract",                                    # Neisseria_polysaccharea
  "Upper respiratory tract",                                    # Capnocytophaga_sputigena
  "Upper respiratory tract",                                    # Capnocytophaga_leadbetteri
  "Upper respiratory tract",                                    # Neisseria_elongata
  "Upper respiratory tract",                                    # Prevotella_oris
  "Various human body sites",                                   # Haemophilus_pittmaniae
  "Various human body sites",                                   # Fusobacterium_nucleatum
  "Upper respiratory tract",                                    # Aggregatibacter_segnis
  "Upper respiratory tract",                                    # Cardiobacterium_hominis
  "Various human body sites",                                   # Schaalia_meyeri
  "Upper respiratory tract",                                    # Kingella_oralis
  "Upper respiratory tract",                                    # Leptotrichia_hongkongensis
  "Upper respiratory tract",                                    # Haemophilus_influenzae
  "Various human body sites",                                   # Pseudoleptotrichia_goodfellowii
  "Upper respiratory tract",                                    # Prevotella_enoeca
  "Environment, animal or food-associated",                     # Veillonella_rodentium
  "Upper respiratory tract",                                    # Leptotrichia_hofstadii
  "Upper respiratory tract",                                    # Actinomyces_viscosus
  "Upper respiratory tract",                                    # Leptotrichia_buccalis
  "Various human body sites",                                   # Sneathia_amnii
  "Various human body sites",                                   # Porphyromonas_asaccharolytica
  "Environment, animal or food-associated",                     # Porphyromonas_cangingivalis
  "Environment, animal or food-associated",                     # Bacteroides_heparinolyticus
  "Environment, animal or food-associated",                     # Prevotella_ruminicola
  "Environment, animal or food-associated",                     # Barnesiella_viscericola
  "Various human body sites",                                   # Paraprevotella_xylaniphila
  "Various human body sites",                                   # Fusobacterium_mortiferum
  "Environment, animal or food-associated",                     # Phocaeicola_salanitronis
  "Various human body sites")                                   # Bacteroides_uniformis

cf_y_rare_spec_C <- c(
  "Upper respiratory tract",                                    # Fusobacterium_pseudoperiodonticum
  "Upper respiratory tract",                                    # Streptococcus_koreensis
  "Environment, animal or food-associated",                     # Streptococcus_equinus
  "Upper respiratory tract",                                    # Streptococcus_oralis
  "Upper respiratory tract",                                    # Streptococcus_pseudopneumoniae
  "Various human body sites",                                   # Atopobium_parvulum
  "Upper respiratory tract",                                    # Neisseria_sicca
  "Various human body sites",                                   # Abiotrophia_defectiva
  "Various human body sites",                                   # Fusobacterium_periodonticum
  "Upper respiratory tract",                                    # Prevotella_fusca
  "Upper respiratory tract",                                    # Prevotella_enoeca
  "Upper respiratory tract",                                    # Rothia_aeria
  "Various human body sites",                                   # Campylobacter_rectus
  "Various human body sites",                                   # Lachnospiraceae_bacterium
  "Upper respiratory tract",                                    # Filifactor_alocis
  "Upper respiratory tract",                                    # Prevotella_denticola
  "Upper respiratory tract",                                    # Corynebacterium_matruchotii
  "Upper respiratory tract",                                    # Streptococcus_sanguinis
  "Various human body sites",                                   # Schaalia_cardiffensis
  "Various human body sites",                                   # Streptococcus_milleri
  "Upper respiratory tract",                                    # Rothia_dentocariosa
  "Upper respiratory tract",                                    # Treponema_denticola
  "Upper respiratory tract")                                    # Treponema_putidum


h_y_rare_spec_A <- c(
  "Various human body sites",                                    # Veillonella_parvula
  "Upper respiratory tract",                                     # Neisseria_cinerea
  "Environment, animal or food-associated",                      # Streptococcus_equinus
  "Upper respiratory tract",                                     # Streptococcus_pseudopneumoniae
  "Various human body sites",                                    # Neisseria_flavescens
  "Upper respiratory tract",                                     # Fusobacterium_pseudoperiodonticum
  "Upper respiratory tract",                                     # Streptococcus_vestibularis
  "Various human body sites",                                    # Schaalia_meyeri
  "Upper respiratory tract",                                     # Streptococcus_pneumoniae
  "Upper respiratory tract",                                     # Haemophilus_haemolyticus
  "Upper respiratory tract",                                     # Leptotrichia_shahii
  "Environment, animal or food-associated",                      # Veillonella_rodentium
  "Upper respiratory tract",                                     # Streptococcus_gordonii
  "Upper respiratory tract",                                     # Neisseria_lactamica
  "Upper respiratory tract",                                     # Prevotella_oris
  "Upper respiratory tract",                                     # Neisseria_sicca
  "Upper respiratory tract",                                     # Prevotella_fusca
  "Various human body sites",                                    # Leptotrichia_trevisanii
  "Upper respiratory tract",                                     # Capnocytophaga_sputigena
  "Various human body sites",                                    # Gemella_morbillorum
  "Upper respiratory tract",                                     # Neisseria_polysaccharea
  "Environment, animal or food-associated",                      # Streptococcus_himalayensis
  "Various human body sites",                                    # Pseudoleptotrichia_goodfellowii
  "Environment, animal or food-associated",                      # Rothia_nasimurium
  "Upper respiratory tract",                                     # Leptotrichia_hofstadii
  "Upper respiratory tract",                                     # Leptotrichia_buccalis
  "Upper respiratory tract",                                     # Capnocytophaga_endodontalis
  "Various human body sites",                                    # Fusobacterium_nucleatum
  "Upper respiratory tract",                                     # Leptotrichia_hongkongensis
  "Upper respiratory tract",                                     # Prevotella_enoeca
  "Upper respiratory tract",                                     # Prevotella_denticola
  "Environment, animal or food-associated",                      # Bacteroides_heparinolyticus
  "Various human body sites",                                    # Actinomyces_israelii
  "Environment, animal or food-associated",                      # Bergeyella_cardium
  "Various human body sites",                                    # Porphyromonas_asaccharolytica
  "Upper respiratory tract",                                     # Actinomyces_slackii
  "Environment, animal or food-associated",                      # Porphyromonas_cangingivalis
  "Environment, animal or food-associated")                      # Prevotella_ruminicola

h_y_rare_spec_B <- c(
  "Upper respiratory tract",                                    # Haemophilus_haemolyticus
  "Upper respiratory tract",                                    # Streptococcus_gwangjuense
  "Various human body sites",                                   # Neisseria_flavescens
  "Upper respiratory tract",                                    # Neisseria_cinerea
  "Environment, animal or food-associated",                     # Streptococcus_equinus
  "Various human body sites",                                   # Veillonella_parvula
  "Upper respiratory tract",                                    # Streptococcus_pneumoniae
  "Upper respiratory tract",                                    # Neisseria_lactamica
  "Upper respiratory tract",                                    # Neisseria_sicca
  "Various human body sites",                                   # Fusobacterium_periodonticum
  "Various human body sites",                                   # Haemophilus_pittmaniae
  "Various human body sites",                                   # Schaalia_meyeri
  "Upper respiratory tract",                                    # Neisseria_polysaccharea
  "Various human body sites",                                   # Gemella_morbillorum
  "Upper respiratory tract",                                    # Neisseria_elongata
  "Upper respiratory tract",                                    # Streptococcus_vestibularis
  "Upper respiratory tract",                                    # Prevotella_oris
  "Upper respiratory tract",                                    # Aggregatibacter_segnis
  "Upper respiratory tract",                                    # Prevotella_fusca
  "Upper respiratory tract",                                    # Capnocytophaga_gingivalis
  "Various human body sites",                                   # Fusobacterium_nucleatum
  "Upper respiratory tract",                                    # Streptococcus_gordonii
  "Upper respiratory tract",                                    # Leptotrichia_hofstadii
  "Environment, animal or food-associated",                     # Veillonella_rodentium
  "Upper respiratory tract",                                    # Streptococcus_cristatus
  "Various human body sites",                                   # Pseudoleptotrichia_goodfellowii
  "Upper respiratory tract",                                    # Leptotrichia_buccalis
  "Environment, animal or food-associated",                     # Bergeyella_cardium
  "Upper respiratory tract",                                    # Leptotrichia_hongkongensis
  "Upper respiratory tract",                                    # Prevotella_enoeca
  "Various human body sites",                                   # Porphyromonas_asaccharolytica
  "Environment, animal or food-associated",                     # Actinobacillus_porcitonsillarum
  "Upper respiratory tract",                                    # Streptococcus_sanguinis
  "Upper respiratory tract",                                    # Capnocytophaga_endodontalis
  "Upper respiratory tract",                                    # Filifactor_alocis
  "Upper respiratory tract",                                    # Capnocytophaga_sputigena
  "Environment, animal or food-associated",                     # Streptococcus_himalayensis
  "Environment, animal or food-associated",                     # Porphyromonas_cangingivalis
  "Various human body sites",                                   # Eikenella_corrodens
  "Various human body sites",                                   # Lachnospiraceae_bacterium
  "Environment, animal or food-associated",                     # Streptococcus_marmotae
  "Upper respiratory tract",                                    # Haemophilus_influenzae
  "Various human body sites",                                   # Actinomyces_israelii
  "Environment, animal or food-associated",                     # Streptococcus_respiraculi
  "Environment, animal or food-associated",                     # Barnesiella_viscericola
  "Environment, animal or food-associated",                     # Otariodibacter_oris
  "Environment, animal or food-associated",                     # Actinobacillus_indolicus
  "Environment, animal or food-associated",                     # Jeotgalibaca_arthritidis
  "Environment, animal or food-associated",                     # Streptococcus_merionis
  "Environment, animal or food-associated",                     # Prevotella_ruminicola
  "Environment, animal or food-associated",                     # Rodentibacter_heylii
  "Environment, animal or food-associated",                     # Mannheimia_granulomatis
  "Environment, animal or food-associated")                     # Vagococcus_carniphilus

h_y_rare_spec_C <-c(
  "Upper respiratory tract",                                    # Capnocytophaga_sputigena
  "Upper respiratory tract",                                    # Haemophilus_haemolyticus
  "Upper respiratory tract",                                    # Prevotella_oris 
  "Various human body sites",                                   # Veillonella_parvula
  "Upper respiratory tract",                                    # Streptococcus_koreensis
  "Various human body sites",                                   # Schaalia_meyeri
  "Upper respiratory tract",                                    # Leptotrichia_shahii
  "Upper respiratory tract",                                    # Leptotrichia_buccalis
  "Various human body sites",                                   # Haemophilus_pittmaniae
  "Various human body sites",                                   # Leptotrichia_trevisanii
  "Upper respiratory tract",                                    # Capnocytophaga_endodontalis
  "Environment, animal or food-associated",                     # Streptococcus_equinus
  "Upper respiratory tract",                                    # Streptococcus_pneumoniae
  "Upper respiratory tract",                                    # Neisseria_sicca
  "Upper respiratory tract",                                    # Leptotrichia_hofstadii
  "Various human body sites",                                   # Eikenella_corrodens
  "Various human body sites",                                   # Fusobacterium_nucleatum
  "Upper respiratory tract",                                    # Selenomonas_sputigena
  "Various human body sites",                                   # Prevotella_intermedia
  "Upper respiratory tract",                                    # Neisseria_elongata
  "Upper respiratory tract",                                    # Streptococcus_sanguinis
  "Various human body sites",                                   # Abiotrophia_defectiva
  "Upper respiratory tract",                                    # Prevotella_denticola
  "Upper respiratory tract",                                    # Prevotella_enoeca
  "Upper respiratory tract",                                    # Filifactor_alocis
  "Environment, animal or food-associated",                     # Veillonella_rodentium
  "Upper respiratory tract",                                    # Streptococcus_vestibularis
  "Various human body sites",                                   # Simonsiella_muelleri
  "Various human body sites",                                   # Gemella_morbillorum
  "Upper respiratory tract",                                    # Prevotella_dentalis
  "Upper respiratory tract",                                    # Leptotrichia_hongkongensis
  "Various human body sites",                                   # Pseudoleptotrichia_goodfellowii
  "Various human body sites",                                   # Lachnospiraceae_bacterium
  "Upper respiratory tract",                                    # Neisseria_polysaccharea
  "Upper respiratory tract",                                    # Rothia_aeria
  "Upper respiratory tract",                                    # Cardiobacterium_hominis
  "Environment, animal or food-associated",                     # Bacteroides_heparinolyticus
  "Upper respiratory tract",                                    # Streptococcus_gordonii
  "Various human body sites",                                   # Schaalia_cardiffensis
  "Various human body sites",                                   # Porphyromonas_asaccharolytica
  "Upper respiratory tract",                                    # Streptococcus_cristatus
  "Environment, animal or food-associated",                     # Streptococcus_himalayensis
  "Various human body sites",                                   # Actinomyces_israelii
  "Environment, animal or food-associated",                     # Prevotella_ruminicola
  "Environment, animal or food-associated",                     # Barnesiella_viscericola
  "Various human body sites",                                   # Paraprevotella_xylaniphila
  "Upper respiratory tract",                                    # Actinomyces_viscosus
  "Upper respiratory tract",                                    # Actinomyces_radicidentis
  "Environment, animal or food-associated",                     # Muribaculaceae_bacterium
  "Environment, animal or food-associated",                     # Bacteroides_zoogleoformans
  "Upper respiratory tract",                                    # Actinomyces_slackii
  "Various human body sites")                                   # Bacteroides_uniformis


my_samples <- c("cf_0", "cf_1_3", "cf_4_6", "h_0", "h_1_3", "h_4_6")
env_info <- c(
  table(cf_y_rare_spec_A)[1],table(cf_y_rare_spec_B)[1],table(cf_y_rare_spec_C)[1],
  table(h_y_rare_spec_A)[1],table(h_y_rare_spec_B)[1],table(h_y_rare_spec_C)[1])
upper_info <- c(
  table(cf_y_rare_spec_A)[2],table(cf_y_rare_spec_B)[2],table(cf_y_rare_spec_C)[2],
  table(h_y_rare_spec_A)[2],table(h_y_rare_spec_B)[2],table(h_y_rare_spec_C)[2])
micro_info <- c(
  table(cf_y_rare_spec_A)[3],table(cf_y_rare_spec_B)[3],table(cf_y_rare_spec_C)[3],
  table(h_y_rare_spec_A)[3],table(h_y_rare_spec_B)[3],table(h_y_rare_spec_C)[3])

merge_loc_df <- data.frame(rbind(env_info, upper_info, micro_info))
colnames(merge_loc_df) <- my_samples
merge_loc_df <- data.frame(t(merge_loc_df))
merge_loc_df$TOTAL <- rowSums(merge_loc_df[,1:3])
merge_loc_df_2 <- merge_loc_df[,1:(ncol(merge_loc_df)-1)] / merge_loc_df$TOTAL
merge_loc_df_2$TOTAL <- NULL
merge_loc_df_2 <- merge_loc_df_2 * 100
merge_loc_df_2 <- round(merge_loc_df_2, 1)
merge_loc_df_2$id <- rownames(merge_loc_df_2)
rownames(merge_loc_df_2) <- NULL
merge_loc_df_2$state <- c("CF", "CF", "CF", "Healthy", "Healthy", "Healthy")
merge_loc_df_2$age <- c("0", "1-3", "4-6", "0", "1-3", "4-6")

merge_loc_df_3 <- gather(merge_loc_df_2, key="where", value="perc", -c(id, state, age))
merge_loc_df_3$where2 <- with(merge_loc_df_3, 
                              ifelse(where == "env_info", "Environment, animal or food-associated",
                                     ifelse(where == "micro_info", 
                                            "Various human body sites",
                                            ifelse(where == "upper_info", "Upper respiratory tract", NA))))

where_rare_plot <- ggplot(merge_loc_df_3) + geom_bar(aes(fill=where2, y=perc, x=age), position="fill", stat="identity") +
  facet_wrap(~state) + xlab("Age group") + ylab(" \n") + scale_fill_viridis(discrete = T) +
  theme_pubr(border=TRUE, legend="bottom") + theme(legend.title = element_blank()) + guides(fill=guide_legend(nrow=1,byrow=TRUE))

# Where do the core species come from?
cf_y_core_spec_A <- c(
  "Upper respiratory tract",                                    # Streptococcus_mitis
  "Upper respiratory tract",                                    # Neisseria_cinerea
  "Various human body sites",                                   # Gemella_haemolysans
  "Upper respiratory tract",                                    # Streptococcus_pseudopneumoniae
  "Upper respiratory tract",                                    # Haemophilus_haemolyticus
  "Upper respiratory tract",                                    # Rothia_mucilaginosa
  "Upper respiratory tract",                                    # Neisseria_mucosa
  "Upper respiratory tract",                                    # Streptococcus_oralis
  "Upper respiratory tract",                                    # Streptococcus_pneumoniae
  "Upper respiratory tract",                                    # Streptococcus_australis
  "Upper respiratory tract")                                    # Haemophilus_parainfluenzae

cf_y_core_spec_B <- c(
  "Upper respiratory tract",                                    # Rothia_mucilaginosa
  "Upper respiratory tract",                                    # Streptococcus_mitis
  "Various human body sites",                                   # Gemella_haemolysans
  "Upper respiratory tract",                                    # Gemella_sanguinis
  "Upper respiratory tract",                                    # Haemophilus_parainfluenzae
  "Various human body sites",                                   # Veillonella_dispar
  "Various human body sites",                                   # Veillonella_atypica
  "Upper respiratory tract",                                    # Neisseria_subflava
  "Upper respiratory tract",                                    # Neisseria_mucosa
  "Upper respiratory tract",                                    # Schaalia_odontolytica
  "Various human body sites",                                   # Haemophilus_parahaemolyticus
  "Upper respiratory tract",                                    # Streptococcus_pseudopneumoniae
  "Upper respiratory tract",                                    # Streptococcus_australis
  "Upper respiratory tract",                                    # Neisseria_cinerea
  "Upper respiratory tract",                                    # Streptococcus_parasanguinis
  "Various human body sites",                                   # Prevotella_melaninogenica
  "Upper respiratory tract",                                    # Lautropia_mirabilis
  "Upper respiratory tract",                                    # Streptococcus_salivarius
  "Upper respiratory tract",                                    # Streptococcus_oralis
  "Upper respiratory tract",                                    # Streptococcus_gwangjuense
  "Upper respiratory tract",                                    # Fusobacterium_pseudoperiodonticum
  "Upper respiratory tract",                                    # Haemophilus_haemolyticus
  "Various human body sites",                                   # Prevotella_jejuni
  "Upper respiratory tract",                                    # Streptococcus_koreensis
  "Upper respiratory tract",                                    # Streptococcus_pneumoniae
  "Various human body sites")                                   # Campylobacter_concisus


cf_y_core_spec_C <- c(
  "Upper respiratory tract",                                    # Haemophilus_parainfluenzae
  "Various human body sites",                                   # Veillonella_dispar
  "Various human body sites",                                   # Veillonella_atypica
  "Upper respiratory tract",                                    # Lautropia_mirabilis
  "Upper respiratory tract",                                    # Rothia_mucilaginosa
  "Upper respiratory tract",                                    # Mogibacterium_diversum
  "Various human body sites",                                   # Prevotella_melaninogenica
  "Various human body sites",                                   # Gemella_haemolysans
  "Upper respiratory tract",                                    # Streptococcus_mitis
  "Upper respiratory tract",                                    # Capnocytophaga_gingivalis
  "Upper respiratory tract",                                    # Neisseria_mucosa
  "Various human body sites",                                   # Prevotella_jejuni
  "Upper respiratory tract",                                    # Streptococcus_parasanguinis
  "Upper respiratory tract",                                    # Schaalia_odontolytica
  "Upper respiratory tract",                                    # Gemella_sanguinis
  "Upper respiratory tract",                                    # Streptococcus_australis
  "Upper respiratory tract",                                    # Eubacterium_sulci
  "Various human body sites",                                   # Campylobacter_concisus
  "Upper respiratory tract",                                    # Streptococcus_salivarius
  "Various human body sites",                                   # Haemophilus_parahaemolyticus
  "Upper respiratory tract",                                    # Leptotrichia_wadei
  "Upper respiratory tract",                                    # Capnocytophaga_leadbetteri
  "Various human body sites",                                   # Veillonella_parvula
  "Upper respiratory tract",                                    # Streptococcus_gwangjuense
  "Upper respiratory tract")                                    # Prevotella_oris

h_y_core_spec_A <- c(
  "Various human body sites",                                   # Veillonella_dispar
  "Various human body sites",                                   # Veillonella_atypica
  "Upper respiratory tract",                                    # Rothia_mucilaginosa
  "Upper respiratory tract",                                    # Schaalia_odontolytica
  "Upper respiratory tract",                                    # Neisseria_mucosa
  "Upper respiratory tract",                                    # Streptococcus_mitis
  "Upper respiratory tract",                                    # Streptococcus_parasanguinis
  "Upper respiratory tract",                                    # Actinomyces_pacaensis
  "Upper respiratory tract",                                    # Gemella_sanguinis
  "Various human body sites",                                   # Prevotella_melaninogenica
  "Various human body sites",                                   # Gemella_haemolysans
  "Upper respiratory tract",                                    # Neisseria_subflava
  "Upper respiratory tract",                                    # Streptococcus_australis
  "Upper respiratory tract",                                    # Streptococcus_oralis
  "Various human body sites",                                   # Prevotella_jejuni
  "Upper respiratory tract",                                    # Streptococcus_salivarius
  "Upper respiratory tract",                                    # Haemophilus_parainfluenzae
  "Upper respiratory tract",                                    # Leptotrichia_wadei
  "Upper respiratory tract",                                    # Atopobium_parvulum
  "Various human body sites",                                   # Campylobacter_concisus
  "Upper respiratory tract",                                    # Streptococcus_koreensis
  "Upper respiratory tract")                                    # Streptococcus_gwangjuense


h_y_core_spec_B <- c(
  "Upper respiratory tract",                                    # Rothia_mucilaginosa
  "Upper respiratory tract",                                    # Haemophilus_parainfluenzae
  "Various human body sites",                                   # Prevotella_jejuni
  "Upper respiratory tract",                                    # Gemella_sanguinis
  "Upper respiratory tract",                                    # Neisseria_mucosa
  "Various human body sites",                                   # Veillonella_dispar
  "Various human body sites",                                   # Veillonella_atypica
  "Upper respiratory tract",                                    # Streptococcus_mitis
  "Various human body sites",                                   # Prevotella_melaninogenica
  "Various human body sites",                                   # Gemella_haemolysans
  "Upper respiratory tract",                                    # Schaalia_odontolytica
  "Upper respiratory tract",                                    # Streptococcus_australis
  "Upper respiratory tract",                                    # Streptococcus_parasanguinis
  "Upper respiratory tract",                                    # Neisseria_subflava
  "Various human body sites",                                   # Haemophilus_parahaemolyticus
  "Upper respiratory tract",                                    # Capnocytophaga_leadbetteri
  "Upper respiratory tract",                                    # Streptococcus_salivarius
  "Upper respiratory tract",                                    # Leptotrichia_wadei
  "Upper respiratory tract",                                    # Streptococcus_oralis
  "Upper respiratory tract",                                    # Fusobacterium_pseudoperiodonticum
  "Upper respiratory tract",                                    # Streptococcus_koreensis
  "Upper respiratory tract",                                    # Streptococcus_pseudopneumoniae
  "Various human body sites",                                   # Simonsiella_muelleri
  "Various human body sites")                                   # Campylobacter_concisus

h_y_core_spec_C <- c(
  "Various human body sites",                                   # Prevotella_jejuni
  "Upper respiratory tract",                                    # Leptotrichia_wadei
  "Various human body sites",                                   # Veillonella_dispar
  "Various human body sites",                                   # Prevotella_melaninogenica
  "Upper respiratory tract",                                    # Rothia_mucilaginosa
  "Various human body sites",                                   # Veillonella_atypica
  "Upper respiratory tract",                                    # Neisseria_mucosa
  "Upper respiratory tract",                                    # Haemophilus_parainfluenzae
  "Various human body sites",                                   # Haemophilus_parahaemolyticus
  "Upper respiratory tract",                                    # Schaalia_odontolytica
  "Upper respiratory tract",                                    # Eubacterium_sulci
  "Upper respiratory tract",                                    # Mogibacterium_diversum
  "Upper respiratory tract",                                    # Neisseria_subflava
  "Upper respiratory tract",                                    # Gemella_sanguinis
  "Upper respiratory tract",                                    # Streptococcus_mitis
  "Various human body sites",                                   # Campylobacter_concisus
  "Upper respiratory tract",                                    # Capnocytophaga_leadbetteri
  "Various human body sites",                                   # Gemella_haemolysans
  "Upper respiratory tract",                                    # Actinomyces_pacaensis
  "Upper respiratory tract",                                    # Streptococcus_parasanguinis
  "Upper respiratory tract",                                    # Fusobacterium_pseudoperiodonticum
  "Upper respiratory tract",                                    # Capnocytophaga_gingivalis
  "Upper respiratory tract",                                    # Atopobium_parvulum
  "Upper respiratory tract",                                    # Streptococcus_salivarius
  "Upper respiratory tract",                                    # Streptococcus_australis
  "Upper respiratory tract",                                    # Lautropia_mirabilis
  "Upper respiratory tract",                                    # Streptococcus_oralis
  "Upper respiratory tract",                                    # Streptococcus_pseudopneumoniae
  "Upper respiratory tract",                                    # Neisseria_cinerea
  "Various human body sites",                                   # Neisseria_flavescens
  "Various human body sites",                                   # Fusobacterium_periodonticum
  "Various human body sites",                                   # Campylobacter_showae
  "Upper respiratory tract",                                    # Streptococcus_gwangjuense
  "Various human body sites")                                   # Prevotella_fusca


env_info_core <- c(0,0,0,0,0,0)
upper_info_core <- 
  c(table(cf_y_core_spec_A)[1], table(cf_y_core_spec_B)[1],table(cf_y_core_spec_C)[1],
    table(h_y_core_spec_A)[1], table(h_y_core_spec_B)[1], table(h_y_core_spec_C)[1])

micro_info_core <- c(
  table(cf_y_core_spec_A)[2], table(cf_y_core_spec_B)[2], table(cf_y_core_spec_C)[2],
  table(h_y_core_spec_A)[2], table(h_y_core_spec_B)[2], table(h_y_core_spec_C)[2])

merge_loc_df_core <- data.frame(rbind(env_info_core, upper_info_core, micro_info_core))
colnames(merge_loc_df_core) <- my_samples
merge_loc_df_core <- data.frame(t(merge_loc_df_core))
merge_loc_df_core$TOTAL <- rowSums(merge_loc_df_core[,1:3])
merge_loc_df_core_2 <- merge_loc_df_core[,1:(ncol(merge_loc_df_core)-1)] / merge_loc_df_core$TOTAL
merge_loc_df_core_2$TOTAL <- NULL
merge_loc_df_core_2 <- merge_loc_df_core_2 * 100
merge_loc_df_core_2 <- round(merge_loc_df_core_2, 1)
merge_loc_df_core_2$id <- rownames(merge_loc_df_core_2)
rownames(merge_loc_df_core_2) <- NULL
merge_loc_df_core_2$state <- c("CF", "CF", "CF", "Healthy", "Healthy", "Healthy")
merge_loc_df_core_2$age <- c("0", "1-3", "4-6", "0", "1-3", "4-6")

merge_loc_df_core_3 <- gather(merge_loc_df_core_2, key="where", value="perc", -c(id, state, age))
merge_loc_df_core_3$where2 <- with(merge_loc_df_core_3, 
                                   ifelse(where == "env_info_core", "Environment, animal or food-associated",
                                          ifelse(where == "gastro_info_core", "Gastrointestinal tract",
                                                 ifelse(where == "micro_info_core", 
                                                        "Isolated from various human body sites",
                                                        ifelse(where == "upper_info_core", "Upper respiratory tract", NA)))))

where_core_plot <- ggplot(merge_loc_df_core_3) + geom_bar(aes(fill=where2, y=perc, x=age), position="fill", stat="identity") +
  facet_wrap(~state) + xlab("Age group") + ylab("Proportion of bacterial species\n") + scale_fill_viridis(discrete = T) +
  theme_pubr(border=TRUE, legend="bottom") + theme(legend.title = element_blank()) + guides(fill=guide_legend(nrow=1,byrow=TRUE))

where_species <- ggarrange(where_core_plot, where_rare_plot,common.legend = TRUE, legend = "bottom", labels=c("a", "b"))

# n per age group
metadata_sub_h <- subset(metadata_sub, state == "Healthy")
table(metadata_sub_h$ageGroup)

metadata_sub_cf <- subset(metadata_sub, state == "CF")
table(metadata_sub_cf$ageGroup)

# make taxonomy data table
uniform_species = c(all_rare, all_core)
uniform_genus = unlist(map(strsplit(uniform_species, split = "_"), 1))

uniform_df = data.frame(cbind(uniform_species, uniform_genus))
uniform_df$Phylum <- with(uniform_df, ifelse(
  uniform_genus == "Rothia", "Actinobacteria",
  ifelse(uniform_genus == "Streptococcus", "Firmicutes",
         ifelse(uniform_genus == "Veillonella", "Firmicutes",
                ifelse(uniform_genus == "Prevotella", "Bacteroidetes",
                       ifelse( uniform_genus == "Neisseria", "Proteobacteria", 
                               ifelse(uniform_genus == "Gemella", "Firmicutes",
                                      ifelse(uniform_genus == "Haemophilus", "Proteobacteria",
                                             ifelse(uniform_genus == "Campylobacter", "Proteobacteria",
                                                    ifelse(uniform_genus == "Schaalia", "Actinobacteria",
                                                           ifelse(uniform_genus == "Lautropia", "Proteobacteria",
                                                                  ifelse(uniform_genus == "Actinomyces", "Actinobacteria",
                                                                         ifelse(uniform_genus == "Capnocytophaga", "Bacteroidetes",
                                                                                ifelse(uniform_genus == "Eubacterium", "Firmicutes", 
                                                                                       ifelse(uniform_genus == "Fusobacterium", "Fusobacteria",
                                                                                              ifelse(uniform_genus == "Mogibacterium", "Firmicutes",
                                                                                                     ifelse(uniform_genus == "Atopobium", "Actinobacteria",
                                                                                                            ifelse(uniform_genus == "Abiotrophia", "Firmicutes",
                                                                                                                   ifelse(uniform_genus == "Bacteroides", "Bacteroidetes",
                                                                                                                          ifelse(uniform_genus == "Porphyromonas", "Bacteroidetes",
                                                                                                                                 ifelse(uniform_genus == "Aggregatibacter", "Proteobacteria",
                                                                                                                                        ifelse(uniform_genus == "Pseudoleptotrichia", "Fusobacteria",
                                                                                                                                               ifelse(uniform_genus == "Bergeyella", "Bacteroidetes",
                                                                                                                                                      ifelse(uniform_genus == "Filifactor", "Firmicutes",
                                                                                                                                                             ifelse(uniform_genus == "Paraprevotella", "Bacteroidetes",
                                                                                                                                                                    ifelse(uniform_genus == "Barnesiella","Bacteroidetes",
                                                                                                                                                                           ifelse(uniform_genus == "Cardiobacterium", "Proteobacteria",
                                                                                                                                                                                  ifelse(uniform_genus == "Lachnospiraceae", "Firmicutes",
                                                                                                                                                                                         ifelse(uniform_genus == "Eikenella", "Proteobacteria",
                                                                                                                                                                                                ifelse(uniform_genus == "Rodentibacter", "Proteobacteria",
                                                                                                                                                                                                       ifelse(uniform_genus == "Leptotrichia", "Fusobacteria", 
                                                                                                                                                                                                              NA)))))))))))))))))))))))))))))))


uniform_df$Class <- with(uniform_df, ifelse(
  uniform_genus == "Rothia", "Actinobacteria",
  ifelse(uniform_genus == "Streptococcus", "Bacilli",
         ifelse(uniform_genus == "Veillonella", "Negativicutes",
                ifelse(uniform_genus == "Prevotella", "Bacteroidia",
                       ifelse(uniform_genus == "Neisseria", "Betaproteobacteria", 
                              ifelse(uniform_genus == "Gemella", "Bacilli",
                                     ifelse(uniform_genus == "Haemophilus", "Gammaproteobacteria",
                                            ifelse(uniform_genus == "Campylobacter", "Epsilonproteobacteria",
                                                   ifelse(uniform_genus == "Schaalia", "Actinobacteria",
                                                          ifelse(uniform_genus == "Lautropia", "Betaproteobacteria",
                                                                 ifelse(uniform_genus == "Actinomyces", "Actinobacteria",
                                                                        ifelse(uniform_genus == "Capnocytophaga", "Flavobacteriia",
                                                                               ifelse(uniform_genus == "Eubacterium", "Clostridia", 
                                                                                      ifelse(uniform_genus == "Fusobacterium", "Fusobacteriia",
                                                                                             ifelse(uniform_genus == "Mogibacterium", "Clostridia",
                                                                                                    ifelse(uniform_genus == "Atopobium", "Actinobacteria",
                                                                                                           ifelse(uniform_genus == "Abiotrophia", "Bacilli",
                                                                                                                  ifelse(uniform_genus == "Bacteroides", "Bacteroidia",
                                                                                                                         ifelse(uniform_genus == "Porphyromonas", "Bacteroidia",
                                                                                                                                ifelse(uniform_genus == "Aggregatibacter", "Gammaproteobacteria",
                                                                                                                                       ifelse(uniform_genus == "Pseudoleptotrichia", "Fusobacteriia",
                                                                                                                                              ifelse(uniform_genus == "Bergeyella", "Flavobacteriia",
                                                                                                                                                     ifelse(uniform_genus == "Filifactor", "Clostridia",
                                                                                                                                                            ifelse(uniform_genus == "Paraprevotella", "Bacteroidia",
                                                                                                                                                                   ifelse(uniform_genus == "Barnesiella", "Bacteroidia",
                                                                                                                                                                          ifelse(uniform_genus == "Cardiobacterium", "Gammaproteobacteria",
                                                                                                                                                                                 ifelse(uniform_genus == "Lachnospiraceae", "Clostridia",
                                                                                                                                                                                        ifelse(uniform_genus == "Eikenella", "Betaproteobacteria",
                                                                                                                                                                                               ifelse(uniform_genus == "Rodentibacter", "Gammaproteobacteria",
                                                                                                                                                                                                      ifelse(uniform_genus == "Leptotrichia","Fusobacteria", 
                                                                                                                                                                                                             NA)))))))))))))))))))))))))))))))

uniform_df$Order <- with(uniform_df, ifelse(
  uniform_genus == "Rothia", "Actinomycetales",
  ifelse(uniform_genus == "Streptococcus", "Lactobacillales",
         ifelse(uniform_genus == "Veillonella", "Vellionellales",
                ifelse(uniform_genus == "Prevotella", "Bacteroidales",
                       ifelse(uniform_genus == "Neisseria", "Neisseriales", 
                              ifelse(uniform_genus == "Gemella", "Bacillales",
                                     ifelse(uniform_genus == "Haemophilus", "Pasteurellales",
                                            ifelse(uniform_genus == "Campylobacter", "Campylobacterales",
                                                   ifelse(uniform_genus == "Schaalia", "Actinomycetales",
                                                          ifelse(uniform_genus == "Lautropia", "Burkholderiales",
                                                                 ifelse(uniform_genus == "Actinomyces", "Actinomycetales",
                                                                        ifelse(uniform_genus == "Capnocytophaga", "Flavobacteriales",
                                                                               ifelse(uniform_genus == "Eubacterium", "Clostridiales", 
                                                                                      ifelse(uniform_genus == "Fusobacterium", "Fusobacteriales",
                                                                                             ifelse(uniform_genus == "Mogibacterium", "Clostridiales",
                                                                                                    ifelse(uniform_genus == "Atopobium", "Coriobacteriales",
                                                                                                           ifelse(uniform_genus == "Abiotrophia", "Lactobacillales",
                                                                                                                  ifelse(uniform_genus == "Bacteroides", "Bacteroidales",
                                                                                                                         ifelse(uniform_genus == "Porphyromonas", "Bacteroidales",
                                                                                                                                ifelse(uniform_genus == "Aggregatibacter", "Pasteurellales",
                                                                                                                                       ifelse(uniform_genus == "Pseudoleptotrichia", "Fusobacteriales",
                                                                                                                                              ifelse(uniform_genus == "Bergeyella", "Flavobacteriales",
                                                                                                                                                     ifelse(uniform_genus == "Filifactor", "Clostridiales",
                                                                                                                                                            ifelse(uniform_genus == "Paraprevotella", "Bacteroidales",
                                                                                                                                                                   ifelse(uniform_genus == "Barnesiella", "Bacteroidales",
                                                                                                                                                                          ifelse(uniform_genus == "Cardiobacterium", "Cardiobacteriales",
                                                                                                                                                                                 ifelse(uniform_genus == "Lachnospiraceae", "Lachnospirales",
                                                                                                                                                                                        ifelse(uniform_genus == "Eikenella", "Neisseriales",
                                                                                                                                                                                               ifelse(uniform_genus == "Rodentibacter", "Pasteurellales",
                                                                                                                                                                                                      ifelse(uniform_genus == "Leptotrichia", "Fusobacterales", 
                                                                                                                                                                                                             NA)))))))))))))))))))))))))))))))

uniform_df$Family <- with(
  uniform_df, ifelse(uniform_genus == "Rothia", "Micrococcaceae",
                     ifelse(uniform_genus == "Streptococcus", "Streptococcaceae",
                            ifelse(uniform_genus == "Veillonella", "Veillonellaceae",
                                   ifelse(uniform_genus == "Prevotella", "Prevotellaceae",
                                          ifelse(uniform_genus == "Neisseria", "Neisseriaceae", 
                                                 ifelse(uniform_genus == "Gemella", "Streptococcaceae",
                                                        ifelse(uniform_genus == "Haemophilus", "Pasteurellaceae",
                                                               ifelse(uniform_genus == "Campylobacter", "Campylobacteraceae",
                                                                      ifelse(uniform_genus == "Schaalia", "Actinomycetaceae",
                                                                             ifelse(uniform_genus == "Lautropia", "Burkholderiaceae",
                                                                                    ifelse(uniform_genus == "Actinomyces", "Actinomycetaceae",
                                                                                           ifelse(uniform_genus == "Capnocytophaga", "Flavobacteriaceae",
                                                                                                  ifelse(uniform_genus == "Eubacterium", "Eubacteriaceae", 
                                                                                                         ifelse(uniform_genus == "Fusobacterium", "Fusobacteriaceae",
                                                                                                                ifelse(uniform_genus == "Mogibacterium", "Eubacteriaceae",
                                                                                                                       ifelse(uniform_genus == "Atopobium", "Coriobacteriaceae",
                                                                                                                              ifelse(uniform_genus == "Abiotrophia", "Aerococcaceae",
                                                                                                                                     ifelse(uniform_genus == "Bacteroides", "Bacteroidaceae",
                                                                                                                                            ifelse(uniform_genus == "Porphyromonas", "Porphyromonadaceae",
                                                                                                                                                   ifelse(uniform_genus == "Aggregatibacter", "Pasteurellaceae",
                                                                                                                                                          ifelse(uniform_genus == "Pseudoleptotrichia", "Leptotrichiaceae",
                                                                                                                                                                 ifelse(uniform_genus == "Bergeyella", "Flavobacteriaceae",
                                                                                                                                                                        ifelse(uniform_genus == "Filifactor", "Peptostreptococcaceae",
                                                                                                                                                                               ifelse(uniform_genus == "Paraprevotella", "Prevotellaceae",
                                                                                                                                                                                      ifelse(uniform_genus == "Barnesiella", "Porphyromonadaceae",
                                                                                                                                                                                             ifelse(uniform_genus == "Cardiobacterium", "Cardiobacteriaceae",
                                                                                                                                                                                                    ifelse(uniform_genus == "Lachnospiraceae", "Lachnospiraceae",
                                                                                                                                                                                                           ifelse(uniform_genus == "Eikenella", "Neisseriaceae",
                                                                                                                                                                                                                  ifelse(uniform_genus == "Rodentibacter", "Pasteurellaceae",
                                                                                                                                                                                                                         ifelse(uniform_genus == "Leptotrichia", "Leptotrichiaceae", 
                                                                                                                                                                                                                                NA)))))))))))))))))))))))))))))))

colnames(uniform_df) <- c("otu", "Genus", "Phylum","Class", "Order", "Family")
uniform_df <- uniform_df[!duplicated(uniform_df), ]
rownames(uniform_df) <- uniform_df$otu
uniform_df$otu <- NULL
uniform_df <- data.frame(t(uniform_df))
uniform_df$Age <- "Age"
uniform_df$BMI <- "BMI"
uniform_df$LCI <- "LCI"
uniform_df$PI_PS_2 <- "PI_PS"
uniform_df <- data.frame(t(uniform_df))
dataset_sub_all <- subset(dataset_sub, rownames(dataset_sub) %in% rownames(uniform_df))
metadata_sub$sample <- rownames(metadata_sub)

dataset_sub_all <- data.frame(t(dataset_sub_all))
dataset_sub_all$Age_1 <- metadata_sub$ageGroup
dataset_sub_all$Age <- ifelse(dataset_sub_all$Age_1 == "0 years", 0,
                              ifelse(dataset_sub_all$Age_1 == "1-3 years", 1, 2))
dataset_sub_all$LCI <- metadata_sub$LCI
dataset_sub_all$PI_PS <- metadata_sub$PI_PS
dataset_sub_all$PI_PS_2 <- ifelse(dataset_sub_all$PI_PS == "PI", 2, 
                                  ifelse(dataset_sub_all$PI_PS == "Healthy", 0, 0))

dataset_sub_all$BMI <- metadata_sub$BMI
dataset_sub_all$PI_PS <- NULL
dataset_sub_all$Age_1 <- NULL

dataset_sub_all <- as.data.frame(scale(dataset_sub_all))

dataset_sub_all_t <- data.frame(t(dataset_sub_all))
metadata_sub_2 <- data.frame(t(metadata_sub))
metadata_sub_3 <- data.frame(t(metadata_sub_2))

otu_mat <- as.matrix(dataset_sub_all_t)
tax_mat <- as.matrix(uniform_df)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)

TAX = tax_table(tax_mat)
samples = sample_data(metadata_sub_3)
erie.merge <- phyloseq(OTU, TAX, samples)

ntaxa(erie.merge)
predictors <- t(otu_table(erie.merge))
dim(predictors)
# Make one column for our outcome/response variable 
response <- as.factor(sample_data(erie.merge)$state)
# Combine them into 1 data frame
rf.data <- data.frame(response, predictors)

sample_val = ncol(rf.data) / 3
set.seed(112)
erie.classify <- randomForest(response~., data = rf.data, 
                              ntree = 19, 
                              mtry=sample_val, 
                              importance=TRUE,
                              na.action = na.roughfix)

print(erie.classify)
plot(erie.classify)

# Make a data frame with predictor names and their importance
imp <- importance(erie.classify)
imp <- data.frame(predictors = rownames(imp), imp)
# Order the predictor levels by importance
imp.sort.gini <- imp[order(imp$MeanDecreaseGini, decreasing = TRUE),]
imp.sort.gini$predictors <- factor(imp.sort.gini$predictors, levels = imp.sort.gini$predictors)
imp.sort.gini <- subset(imp.sort.gini, MeanDecreaseAccuracy > 0)

fluctuating_rare <- c(all_rare, f_rare_species)
fluctuating_rare <- fluctuating_rare[!duplicated(fluctuating_rare)]
fluctuating_core <- c(all_core, f_core_species)
fluctuating_core <- fluctuating_core[!duplicated(fluctuating_core)]

imp.sort.gini$species_type <- with(imp.sort.gini,
                                   ifelse(predictors %in% f_rare_species, "Rare species (background)", 
                                          ifelse(predictors %in% f_core_species, "Core species (background)", 
                                                 ifelse(predictors %in% fluctuating_rare, "Rare species (non-persistent)",
                                                        ifelse(predictors %in% fluctuating_core, "Core species (non-persistent)", 
                                                               "Host-associated variables")))))

imp.sort.gini$predictors <- str_replace(imp.sort.gini$predictors, "_", " ")

# most important predictors/features for classifying infant samples into CF and healthy
randomForest_plot_Gini <- ggplot(imp.sort.gini, aes(x = reorder(predictors, MeanDecreaseGini), 
                                                    y = MeanDecreaseGini)) +
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + xlab("") +
  theme(legend.title = element_blank(), text = element_text(size = 12), 
        axis.text.y = element_text(face = "italic"), legend.position = "none")

imp.sort.gini.short <-ddply(imp.sort.gini, "species_type", numcolwise(sum))
aa <- c("Host-associated variables", 0, 0, 0, 0)
imp.sort.gini.short <- rbind(imp.sort.gini.short, aa)
imp.sort.gini.short$MeanDecreaseAccuracy <- as.numeric(as.character(imp.sort.gini.short$MeanDecreaseAccuracy))
imp.sort.gini.short$MeanDecreaseGini <- as.numeric(as.character(imp.sort.gini.short$MeanDecreaseGini))
imp.sort.gini.short$CF <- as.numeric(as.character(imp.sort.gini.short$CF))
imp.sort.gini.short$Healthy <- as.numeric(as.character(imp.sort.gini.short$Healthy))

randomForest_plot_Gini1 <- ggplot(imp.sort.gini.short, aes(x = reorder(species_type, MeanDecreaseGini), 
                                                           y = MeanDecreaseGini, fill=species_type)) +
  geom_bar(stat = "identity") + coord_flip() + theme_pubr(border=TRUE) + xlab("") +
  scale_fill_manual(values = c("orange", "orange", "black", "lightblue", "lightblue")) +
  theme(legend.title = element_blank(), text = element_text(size = 12), legend.position = "none", axis.text.y = element_blank())

randomForest_plot_Acc1 <- ggplot(imp.sort.gini.short, aes(x = reorder(species_type, MeanDecreaseAccuracy), 
                                                          y = MeanDecreaseAccuracy, fill=species_type)) +
  geom_bar(stat = "identity") + coord_flip() + theme_pubr(border=TRUE) + xlab("") +
  scale_fill_manual(values = c("orange", "orange", "black", "lightblue", "lightblue")) +
  theme(legend.title = element_blank(), text = element_text(size = 12), legend.position = "none") 

sum_all_rf <- sum(table(imp.sort.gini$species_type))
sum_rare_rf <- table(imp.sort.gini$species_type)[3] + table(imp.sort.gini$species_type)[4]
sum_core_rf <- table(imp.sort.gini$species_type)[1] + table(imp.sort.gini$species_type)[2]
percentage_rare_species = (sum_rare_rf / sum_all_rf) * 100
percentage_core_species = (sum_core_rf / sum_all_rf) * 100

plot_pie <- data.frame(value1=c(percentage_rare_species, percentage_core_species),
                       group1=c("Rare species",  "Core species"))

bp <- ggplot(plot_pie, aes(x="", y=value1, fill=fct_inorder(group1))) + geom_bar(width = 1, stat = "identity") +
  coord_polar("y") + theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
                           panel.grid=element_blank(), panel.background=element_blank(), 
                           legend.direction = "vertical", legend.position = "bottom", legend.title = element_blank()) +
  scale_fill_manual(values = c("lightblue", "orange")) +
  geom_label(aes(label = paste(round(value1), "%")), size=4, show.legend = F, nudge_x = 0) 

randomForest_plot_final <- ggarrange(randomForest_plot_Acc1, randomForest_plot_Gini1, bp, 
                                     nrow=1, labels = c("a", "b", "c"), widths=c(2, 1.28, 1))
randomForest_plot_final


# Principal component analysis ####
pca_patternRec <- rf.data
pca_patternRec <- subset(pca_patternRec, rownames(pca_patternRec) != "KGCF17")
pca_patternRec <- subset(pca_patternRec, rownames(pca_patternRec) != "EMCF09")
pca_patternRec <- subset(pca_patternRec, rownames(pca_patternRec) != "EMCF28")

pca_patternRec <- na.omit(pca_patternRec)
state_plotting <- pca_patternRec$response
pca_patternRec$response <- NULL
pca_patternRec <- prcomp(pca_patternRec, center = FALSE, scale. = FALSE)

pcaInd <- fviz_pca_ind(pca_patternRec, geom.ind = c("arrow"), pointshape = 21, pointsize = 3, 
                       palette = c("#00AFBB", "#E7B800"), fill.ind = state_plotting, col.ind = state_plotting, 
                       repel = FALSE, legend.title = "") + theme_pubr(border = TRUE, base_size = 12, legend = "right") +
  theme(title = element_blank(), legend.title = element_blank())


pcaVar <- fviz_pca_var(pca_patternRec, select.var = list(contrib = 40), col.var = "contrib", geom = c("arrow"),
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) +
  theme_pubr(border=TRUE, base_size = 12, legend = "right") + theme(title = element_blank(), legend.title = element_blank())

# Extract variable contributions of PCA
varContribution <- data.frame(pcaVar$data)
# rename columns
colnames(varContribution) <- c("variable", "Dim1", "Dim2", "coord", "cos2", "contrib")
# remove rownames
rownames(varContribution) <- NULL
# order table
varContribution <- varContribution[order(-varContribution$contrib),]
# first 40 variables explain most of the variance 
sum(varContribution$contrib[1:40])
# extract 40 variables
varContribution$species_type <- with(varContribution,
                                     ifelse(variable %in% f_rare_species, "Rare species (background)", 
                                            ifelse(variable %in% f_core_species, "Core species (background)", 
                                                   ifelse(variable %in% fluctuating_rare, "Rare species (non-persistent)",
                                                          ifelse(variable %in% fluctuating_core, "Core species (non-persistent)", 
                                                                 "Host-associated variables")))))
varContribution40_PCA_randomForest <- varContribution[1:40,]
sum_pca <- sum(table(varContribution40_PCA_randomForest$species_type))
table(varContribution40_PCA_randomForest$species_type)
sum_rare_pca <- table(varContribution40_PCA_randomForest$species_type)[3] + table(varContribution40_PCA_randomForest$species_type)[4]
sum_core_pca <- table(varContribution40_PCA_randomForest$species_type)[1] + table(varContribution40_PCA_randomForest$species_type)[2]
percentage_rare_species2 = (sum_rare_pca / sum_pca) * 100
percentage_core_species2 = (sum_core_pca / sum_pca) * 100

plot_pie2 <- data.frame(value1=c(percentage_rare_species2, percentage_core_species2),
                        group1=c("Rare species", "Core species"))

bp2 <- ggplot(plot_pie2, aes(x="", y=value1, fill=fct_inorder(group1))) +
  geom_bar(width = 1, stat = "identity") + coord_polar("y") +
  theme(axis.text=element_blank(), axis.title=element_blank(), panel.grid=element_blank(), axis.ticks = element_blank(),
        panel.background=element_blank(), legend.direction = "vertical", legend.position = "right",
        legend.title = element_blank()) + scale_fill_manual(values = c("lightblue", "orange")) +
  geom_label(aes(label = paste(round(value1), "%")), size=4, show.legend = F, nudge_x = 0) 

pca_plot <- ggarrange(pcaInd, ggarrange(pcaVar, bp2, nrow = 2, labels = c("b", "c")), nrow = 1, labels = "a")


ggsave("image_output/Fig_1.tif", vennDiagramsPlot, device="tiff",
       width=8.93, height=8.02, dpi=600)

ggsave("image_output/Fig_2.tif", randomForest_plot_final, dpi=600,
       width=12.1, height=5.73, device="tiff")

ggsave("image_output/Supplementary_Fig_1.tif", plot_pheatmap, device="tiff",
       width=8.75, height=9.29, dpi=300)

ggsave("image_output/Supplementary_Fig_3.tif", where_species, dpi=300,
       device="tiff", width=9.79, height=6.9)

ggsave("image_output/Supplementary_Fig_4.tif", from_where_to_where, device="tiff",
       width=8.93, height=8.02, dpi=300)

ggsave("image_output/Supplementary_Fig_5.tif", pca_plot, dpi=300,
       width=9.58, height=7.98, device="tiff")

write.csv2(varContribution40_PCA_randomForest, "PCA_randomForest_table.csv", row.names = FALSE)
