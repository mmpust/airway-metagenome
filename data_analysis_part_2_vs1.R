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
vennDiagramsPlot


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
uniform_df <- data.frame(t(uniform_df))
dataset_sub_all <- subset(dataset_sub, rownames(dataset_sub) %in% rownames(uniform_df))
metadata_sub$sample <- rownames(metadata_sub)
dataset_sub_all <- data.frame(t(dataset_sub_all))
dataset_sub_all <- as.data.frame(scale(dataset_sub_all))
dataset_sub_all$Age <- metadata_sub$Age_in_months
dataset_sub_all$BMI <- metadata_sub$BMI
dataset_sub_all$LCI <- metadata_sub$LCI

# effect of adding mutation class 
#Type of random forest: classification
#                    Number of trees: 28
#No. of variables tried at each split: 39

#       OOB estimate of  error rate: 30.59%
#Confusion matrix:
#       CF Healthy class.error
#CF      27      14   0.3414634
#Healthy 12      32   0.2727273
# Age	Age	0.0000000	-1.01835015	-1.0183502	0.23735225
# BMI	BMI	-1.0183502	-0.48302500	-0.8057138	0.63232213
# LCI	LCI	0.0000000	0.00000000	0.0000000	0.19331998
# genotype2	genotype2	0.0000000	0.00000000	0.0000000


# Random forest analysis ####
dataset_sub_all <- data.frame(t(dataset_sub_all))

otu_mat <- as.matrix(dataset_sub_all)
tax_mat <- as.matrix(uniform_df)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)

TAX = tax_table(tax_mat)
samples = sample_data(metadata_sub)
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
                              ntree = 28, 
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
  theme(legend.title = element_blank(), text = element_text(size = 12), legend.position = "none")

randomForest_plot_Acc1 <- ggplot(imp.sort.gini.short, aes(x = reorder(species_type, MeanDecreaseAccuracy), 
                                                          y = MeanDecreaseAccuracy, fill=species_type)) +
  geom_bar(stat = "identity") + coord_flip() + theme_pubr(border=TRUE) + xlab("") +
  scale_fill_manual(values = c("orange", "orange", "black", "lightblue", "lightblue")) +
  theme(legend.title = element_blank(), text = element_text(size = 12),
        legend.position = "none", axis.text.y = element_blank()) 

table(imp.sort.gini$species_type)
percentage_rare_species = (42 / 50) * 100
percentage_core_species = (8 / 50) * 100

plot_pie <- data.frame(value1=c(percentage_rare_species, percentage_core_species),
                       group1=c("Rare species",  "Core species"))

bp <- ggplot(plot_pie, aes(x="", y=value1, fill=fct_inorder(group1))) + geom_bar(width = 1, stat = "identity") +
  coord_polar("y") + theme(axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank(),
                           panel.grid=element_blank(), panel.background=element_blank(), 
                           legend.direction = "vertical", legend.position = "bottom", legend.title = element_blank()) +
  scale_fill_manual(values = c("lightblue", "orange")) +
  geom_label(aes(label = paste(round(value1), "%")), size=4, show.legend = F, nudge_x = 0) 

randomForest_plot_final <- ggarrange(randomForest_plot_Gini1, randomForest_plot_Acc1, bp, 
                                     nrow=1, labels = c("a", "b", "c"), widths=c(2, 1.28, 1))

# Principal component analysis ####
pca_patternRec <- rf.data
pca_patternRec <- subset(pca_patternRec, rownames(pca_patternRec) != "KGCF17")
pca_patternRec <- subset(pca_patternRec, rownames(pca_patternRec) != "EMCF09")
pca_patternRec <- subset(pca_patternRec, rownames(pca_patternRec) != "EMCF28")
state_plotting <- pca_patternRec$response
pca_patternRec$response <- NULL
pca_patternRec$LCI <- NULL
pca_patternRec$BMI <- NULL

pca_patternRec <- prcomp(pca_patternRec, center = FALSE, scale. = TRUE)

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
table(varContribution40_PCA_randomForest$species_type)

percentage_rare_species2 = (33 / 40) * 100
percentage_core_species2 = (7 / 40) * 100

plot_pie2 <- data.frame(value1=c(percentage_rare_species2, percentage_core_species2),
                        group1=c("Rare species", "Core species"))

bp2 <- ggplot(plot_pie2, aes(x="", y=value1, fill=fct_inorder(group1))) +
  geom_bar(width = 1, stat = "identity") + coord_polar("y") +
  theme(axis.text=element_blank(), axis.title=element_blank(), panel.grid=element_blank(),
        panel.background=element_blank(), legend.direction = "vertical", legend.position = "right",
        legend.title = element_blank()) + scale_fill_manual(values = c("lightblue", "orange")) +
  geom_label(aes(label = paste(round(value1), "%")), size=4, show.legend = F, nudge_x = 0) 

pca_plot <- ggarrange(pcaInd, ggarrange(pcaVar, bp2, nrow = 2, labels = c("b", "c")), nrow = 1, labels = "a")

# Generate network input ####
# set global variables 
weight_val = 0.60 
weight_val2 = -1.0
sig_level = 0.01

# correlation between rare and core species
# generate correlation matrix
# healthy, background rare
df_h_rare <- subset(dataset_sub,
                    rownames(dataset_sub) %in% rownames(background_core_h) |
                      rownames(dataset_sub) %in% rownames(background_rare_h))
df_h_rare <- data.frame(t(df_h_rare))
df_h_rare$state <- metadata_sub$state
df_h_rare <- subset(df_h_rare, state == "Healthy")
df_h_rare$state <- NULL
df_h_rare <- rcorr(as.matrix(df_h_rare), type = 'spearman')
# create node and edge lists
# healthy
df_h_rare_Pval <- df_h_rare$P
df_h_rare_COR <- df_h_rare$r
# generate edgelist
df_h_rare_COR_edges <- melt(df_h_rare_Pval)
df_h_rare_coredges <- melt(df_h_rare_COR)
df_h_rare_COR_edges$COR <- df_h_rare_coredges$value
df_h_rare_COR_edges$value <- round(df_h_rare_COR_edges$value, 5)
df_h_rare_COR_edges$Label <- row.names(df_h_rare_coredges)
colnames(df_h_rare_COR_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_h_rare_COR_edges) <- NULL
df_h_rare_COR_edges <- df_h_rare_COR_edges[complete.cases(df_h_rare_COR_edges), ]
df_h_rare_COR_edges_short <- subset(df_h_rare_COR_edges, pValue < sig_level & Weight >= weight_val | 
                                      pValue < sig_level & Weight <= weight_val2)
df_h_rare_COR_edges_short$cortype <- ifelse(df_h_rare_COR_edges_short$Weight > 0,  "pos", "neg")
df_h_rare_nodes <- select(df_h_rare_COR_edges_short, c("Source", "cortype"))
df_h_rare_nodes = df_h_rare_nodes[!duplicated(df_h_rare_nodes$Source),]
df_h_rare_nodes <- data.frame(df_h_rare_nodes)
rownames(df_h_rare_nodes) <- NULL
colnames(df_h_rare_nodes) <- c("Id", "CorType")
df_h_rare_nodes$species_type <- ifelse(df_h_rare_nodes$Id %in% rownames(background_rare_h), "rare", "core")
df_h_rare_nodes$Genus <- df_h_rare_nodes$Id
df_h_rare_nodes$Genus <- sapply(strsplit(as.character(df_h_rare_nodes$Genus),"\\."), `[`, 1)

# # CF, background rare
df_cf_rare <- subset(dataset_sub, 
                     rownames(dataset_sub) %in% rownames(background_core_cf) |
                       rownames(dataset_sub) %in% rownames(background_rare_cf))
df_cf_rare <- data.frame(t(df_cf_rare))
df_cf_rare$state <- metadata_sub$state
df_cf_rare <- subset(df_cf_rare, state == "CF")
df_cf_rare$state <- NULL
df_cf_rare <- rcorr(as.matrix(df_cf_rare), type = 'spearman')
# create node and edge lists
df_cf_rare_Pval <- df_cf_rare$P
df_cf_rare_COR <- df_cf_rare$r
# generate edgelist
df_cf_rare_COR_edges <- melt(df_cf_rare_Pval)
df_cf_rare_coredges <- melt(df_cf_rare_COR)
df_cf_rare_COR_edges$COR <- df_cf_rare_coredges$value
df_cf_rare_COR_edges$value <- round(df_cf_rare_COR_edges$value, 5)
df_cf_rare_COR_edges$Label <- row.names(df_cf_rare_coredges)
colnames(df_cf_rare_COR_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
rownames(df_cf_rare_COR_edges) <- NULL
df_cf_rare_COR_edges <- df_cf_rare_COR_edges[complete.cases(df_cf_rare_COR_edges), ]
df_cf_rare_COR_edges_short <- subset(df_cf_rare_COR_edges, 
                                     pValue < sig_level & Weight >= weight_val | 
                                       pValue < sig_level & Weight <= weight_val2)
df_cf_rare_COR_edges_short$cortype <- ifelse(df_cf_rare_COR_edges_short$Weight > 0,  "pos", "neg")
df_cf_rare_nodes <- select(df_cf_rare_COR_edges_short, c("Source", "cortype"))
df_cf_rare_nodes = df_cf_rare_nodes[!duplicated(df_cf_rare_nodes$Source),]
df_cf_rare_nodes <- data.frame(df_cf_rare_nodes)
rownames(df_cf_rare_nodes) <- NULL
colnames(df_cf_rare_nodes) <- c("Id", "CorType")
df_cf_rare_nodes$species_type <- ifelse(df_cf_rare_nodes$Id %in% rownames(background_rare_cf), "rare", "core")
df_cf_rare_nodes$Genus <- df_cf_rare_nodes$Id
df_cf_rare_nodes$Genus <- sapply(strsplit(as.character(df_cf_rare_nodes$Genus),"\\."), `[`, 1)

# Background networks ####
directed_network = FALSE
# healthy, background core and rare species ####
df_h_rare_nodes$species_type <- as.factor(as.character(df_h_rare_nodes$species_type))
h_rare_net <- graph_from_data_frame(d=df_h_rare_COR_edges_short, vertices=df_h_rare_nodes, directed = directed_network)
h_rare_net <- simplify(h_rare_net, remove.multiple = F, remove.loops = T)
V(h_rare_net)$color <- ifelse(V(h_rare_net)$species_type == "rare", "lightsteelblue2", "orange2")

h.net.plot <- layout.fruchterman.reingold(h_rare_net)
h.net.plot.df <- as.data.frame(h.net.plot)
h.net.plot.df$Id <- df_h_rare_nodes$Id
h.net.plot.df$species_type <- df_h_rare_nodes$species_type

h.net.plot.edges <- get.data.frame(h_rare_net)
h.net.plot.edges$from.x <- h.net.plot.df$V1[match(h.net.plot.edges$from, h.net.plot.df$Id)]  
h.net.plot.edges$from.y <- h.net.plot.df$V2[match(h.net.plot.edges$from, h.net.plot.df$Id)]
h.net.plot.edges$to.x <- h.net.plot.df$V1[match(h.net.plot.edges$to, h.net.plot.df$Id)]  
h.net.plot.edges$to.y <- h.net.plot.df$V2[match(h.net.plot.edges$to, h.net.plot.df$Id)]
h.net.plot.edges$species_type <- ifelse(h.net.plot.edges$from %in% f_rare_species, "rare", "core")

# cf, background core and rare ##################
cf_rare_edges <- df_cf_rare_COR_edges_short
cf_rare_nodes <- df_cf_rare_nodes
cf_rare_nodes$species_type <- as.factor(as.character(cf_rare_nodes$species_type))

cf_rare_net <- graph_from_data_frame(d=cf_rare_edges, vertices=cf_rare_nodes, directed = directed_network)
cf_rare_net <- simplify(cf_rare_net, remove.multiple = F, remove.loops = T)
V(cf_rare_net)$color <- ifelse(V(cf_rare_net)$species_type == "rare", "lightsteelblue2", "orange2")

cf.net.plot <- layout.fruchterman.reingold(cf_rare_net)
cf.net.plot.df <- as.data.frame(cf.net.plot)
cf.net.plot.df$Id <- cf_rare_nodes$Id
cf.net.plot.df$species_type <- cf_rare_nodes$species_type
cf.net.plot.df$node_size <- V(cf_rare_net)$size

cf.net.plot.edges <- get.data.frame(cf_rare_net)
cf.net.plot.edges$from.x <- cf.net.plot.df$V1[match(cf.net.plot.edges$from, cf.net.plot.df$Id)]  
cf.net.plot.edges$from.y <- cf.net.plot.df$V2[match(cf.net.plot.edges$from, cf.net.plot.df$Id)]
cf.net.plot.edges$to.x <- cf.net.plot.df$V1[match(cf.net.plot.edges$to, cf.net.plot.df$Id)]  
cf.net.plot.edges$to.y <- cf.net.plot.df$V2[match(cf.net.plot.edges$to, cf.net.plot.df$Id)]
cf.net.plot.edges$species_type <- ifelse(cf.net.plot.edges$from %in% f_rare_species, "rare", "core")

h_net_plot <- ggplot() + geom_segment(data = h.net.plot.edges, 
                                      aes(x = from.x, xend = to.x, y = from.y, yend = to.y, 
                                          colour=species_type), alpha=0.4) +
  geom_point(data = h.net.plot.df, aes(x=V1,y=V2,colour=species_type,size=2)) +  
  theme_pubr(border=TRUE, legend = "none")  +
  scale_colour_manual(values=c("orange2", "darkblue")) +
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(), 
        axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

cf_net_plot <- ggplot() + geom_segment(data = cf.net.plot.edges, 
                                       aes(x = from.x, xend = to.x, y = from.y, yend = to.y,
                                           colour=species_type), alpha=0.4) +
  geom_point(data = cf.net.plot.df, aes(x=V1,y=V2,colour=species_type,size=2)) +  
  theme_pubr(border=TRUE, legend = "none")  +
  scale_colour_manual(values=c("orange2", "darkblue")) +
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(), 
        axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

# "Distance Congruence" and is based upon a non-parametric correlation of path distances.
# test_congruence(h_rare_net, h_rare_net, method="distance")
test_congruence(h_rare_net, cf_rare_net, method="distance")

# Network robustness and vulnaribiltiy ####
# healthy, background core and rare species
h_rare_net_attack <-swan_combinatory(h_rare_net,30)
h_rare_net_attack <- data.frame(h_rare_net_attack)
colnames(h_rare_net_attack) <- c("fraction_removed", "betweeness", "degree", "cascading", "random")
h_rare_net_attack$fraction_removed <- round(h_rare_net_attack$fraction_removed,1)
h_rare_net_attack <- ddply(h_rare_net_attack,"fraction_removed",numcolwise(median))
h_rare_net_attack[1, ] <- c(0.0, 0.0, 0.0, 0.0, 0.0)

# cf, background core and rare species
cf_rare_net_attack <-swan_combinatory(cf_rare_net,30)
cf_rare_net_attack <- data.frame(cf_rare_net_attack)
colnames(cf_rare_net_attack) <- c("fraction_removed", "betweeness", "degree", "cascading", "random")
cf_rare_net_attack$fraction_removed <- round(cf_rare_net_attack$fraction_removed,1)
cf_rare_net_attack <- ddply(cf_rare_net_attack,"fraction_removed",numcolwise(median))
cf_rare_net_attack <- rbind(h_rare_net_attack[1,], cf_rare_net_attack)

# similarity measures
select_items <- c("fraction_removed", "degree")
cf_rare_net_matrix_degree <- select(cf_rare_net_attack, select_items)
cf_rare_net_matrix_degree <- as.matrix(cf_rare_net_matrix_degree)
h_rare_net_matrix_degree <- select(h_rare_net_attack, select_items)
h_rare_net_matrix_degree <- h_rare_net_matrix_degree[-c(3,6,9),]
h_rare_net_matrix_degree <- as.matrix(h_rare_net_matrix_degree)
original_frechet_h_cf_degree <- Frechet(h_rare_net_matrix_degree, cf_rare_net_matrix_degree)
# 0.25 
original_frechet_h_cf_degree
select_items2 <- c("fraction_removed", "random")
cf_rare_net_matrix_random <- select(cf_rare_net_attack, select_items2)
cf_rare_net_matrix_random <- as.matrix(cf_rare_net_matrix_random)
h_rare_net_matrix_random <- select(h_rare_net_attack, select_items2)
h_rare_net_matrix_random <- h_rare_net_matrix_random[-c(3,6,9),]
h_rare_net_matrix_random <- as.matrix(h_rare_net_matrix_random)
original_frechet_h_cf_random <- Frechet(h_rare_net_matrix_random, cf_rare_net_matrix_random)
# 0.12

attack_degree <- ggplot() + 
  geom_line(data = cf_rare_net_attack, aes(x = fraction_removed, y = degree), color = "red", size=1) + 
  geom_line(data = h_rare_net_attack, aes(x = fraction_removed, y = degree), color = "black", size=1) +
  xlab("Fraction of nodes removed") + geom_label(aes(x=0.9, y=0.05), label="Degree") +
  ylim(0.0, 1.0) + ylab("Connectivity loss") + theme_pubr(border=TRUE)

attack_cascase <- ggplot() +
  geom_line(data = cf_rare_net_attack, aes(x = fraction_removed, y = cascading), color = "red", size=1) +
  geom_line(data = h_rare_net_attack, aes(x = fraction_removed, y = cascading), color = "black", size=1) +
  xlab(" ") + ylab("Connectivity loss") + geom_label(aes(x=0.85, y=0.05), label="Betweeness") +
  ylim(0.0, 1.0) + theme_pubr(border=TRUE)

attack_random <- ggplot() +
  geom_line(data = cf_rare_net_attack, aes(x = fraction_removed, y = random), color = "red", size=1) +
  geom_line(data = h_rare_net_attack, aes(x = fraction_removed, y = random), color = "black", size=1) +
  xlab("Fraction of nodes removed") + ylab(" ") + geom_label(aes(x=0.9, y=0.05), label="Random") +
  ylim(0.0, 1.0) + theme_pubr(border=TRUE)

ab1 <- ggarrange(h_net_plot, cf_net_plot, nrow = 1, ncol = 2, labels = c("a", "b"), label.x = 0.025, label.y = 0.98) 
ab2 <- ggarrange(attack_degree, attack_random, nrow = 1, ncol = 2, labels = c("c", "d"), label.x = 0.12, label.y = 0.98)
first_network_plot <- ggarrange(ab1, ab2, nrow=2)

# Healthy statistics (background species)
h_rare_statistics <- all_indices(h_rare_net)
h_rare_statistics <- data.frame(h_rare_statistics)
h_rare_statistics$Id <- df_h_rare_nodes$Id
h_rare_statistics$type <- df_h_rare_nodes$species_type
colnames(h_rare_statistics) <- c("Degree", "Betweeness", "Closeness", 
                                 "Authority","Id", "Species_type")
h_rare_statistics$state <- "Healthy"
h_rare_statistics$Closeness <- h_rare_statistics$Closeness * 100
h_rare_statistics$Closeness <- round(h_rare_statistics$Closeness,1)

# Cf statistics (background species)
cf_rare_statistics <- all_indices(cf_rare_net)
cf_rare_statistics <- data.frame(cf_rare_statistics)
cf_rare_statistics$Id <- cf_rare_nodes$Id
cf_rare_statistics$type <- cf_rare_nodes$species_type
colnames(cf_rare_statistics) <- c("Degree", "Betweeness", "Closeness", 
                                  "Authority", "Id", "Species_type")
cf_rare_statistics$state <- "CF"
cf_rare_statistics$Closeness <- round(cf_rare_statistics$Closeness,1)

rare_statistics <- rbind(h_rare_statistics, cf_rare_statistics)
rare_statistics_L <- gather(rare_statistics, key="statistics", value="value", 
                            -c("Id", "Species_type", "state"))
rare_statistics_L <- subset(rare_statistics_L, state == "Healthy")
rare_statistics_L_1 <- subset(rare_statistics_L, statistics != "Authority")

stats_rare_core <- ggplot(rare_statistics_L_1, aes(x=Species_type,y=value)) +
  geom_boxplot() + geom_jitter(aes(colour=Species_type), width = 0.02, size=3) +
  facet_wrap(~statistics, scales = "free_y", nrow=1) + stat_compare_means(label="p.signif", 
                                                                          hide.ns = TRUE, size = 5,
                                                                          label.x.npc = c("center"),
                                                                          label.y.npc = c("top")) +
  theme_pubr(border=TRUE) + xlab("") + ylab("Centrality value") + scale_colour_manual(values=c("orange", "darkblue")) +
  theme(legend.position = "bottom", legend.title = element_blank())

wilcox.test(Betweeness ~ Species_type, data=rare_statistics)
wilcoxonR(rare_statistics$Betweeness, g=rare_statistics$Species_type, ci=TRUE)

wilcox.test(Closeness ~ Species_type, data=rare_statistics)
wilcoxonR(rare_statistics$Closeness, g=rare_statistics$Species_type, ci=TRUE)

wilcox.test(Degree ~ Species_type, data=rare_statistics)
wilcoxonR(rare_statistics$Degree, g=rare_statistics$Species_type, ci=TRUE)


# Network modulation ####
# Generate random integers between 1 and 200
random_seeds <- sample(1:200, 100) #20
add_nodes <- c(1:18)
frechet_attack_similarity = NULL
network_congruency = NULL
network_statistics = NULL
netwerk_kernel_randomWalk = NULL
network_kernel_shortestPath = NULL
network_kernel_WL = NULL
selected_rare_species = NULL
selected_core_species = NULL
selected_bothr_species = NULL
selected_bothc_species = NULL
sample_with_replacement = FALSE

for (i in random_seeds){
  for (n_nodes in add_nodes){
    set.seed(i)
    # CF, rescue healthy rare
    sample_h_rare <- sample(rownames(background_rare_h), n_nodes, replace = sample_with_replacement)
    selected_rare_species = rbind(selected_rare_species, data.frame(i, n_nodes, sample_h_rare))
    df_cf_rare_res <- subset(dataset_sub, 
                             rownames(dataset_sub) %in% rownames(background_core_cf) |
                               rownames(dataset_sub) %in% sample_h_rare)
    df_cf_rare_res <- data.frame(t(df_cf_rare_res))
    df_cf_rare_res$state <- metadata_sub$state
    df_cf_rare_res <- subset(df_cf_rare_res, state == "CF")
    df_cf_rare_res$state <- NULL
    df_cf_rare_res <- rcorr(as.matrix(df_cf_rare_res), type = 'spearman')
    
    # create node and edge lists
    df_cf_rare_Pval_res <- df_cf_rare_res$P
    df_cf_rare_COR_res <- df_cf_rare_res$r
    # generate edgelist
    df_cf_rare_COR_edges_res <- melt(df_cf_rare_Pval_res)
    df_cf_rare_coredges_res <- melt(df_cf_rare_COR_res)
    df_cf_rare_COR_edges_res$COR <- df_cf_rare_coredges_res$value
    df_cf_rare_COR_edges_res$value <- round(df_cf_rare_COR_edges_res$value, 5)
    df_cf_rare_COR_edges_res$Label <- row.names(df_cf_rare_coredges_res)
    colnames(df_cf_rare_COR_edges_res) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
    rownames(df_cf_rare_COR_edges_res) <- NULL
    df_cf_rare_COR_edges_res <- df_cf_rare_COR_edges_res[complete.cases(df_cf_rare_COR_edges_res), ]
    df_cf_rare_COR_edges_short_res <- subset(
      df_cf_rare_COR_edges_res, pValue < sig_level & Weight >= weight_val | pValue < sig_level & Weight <= weight_val2)
    df_cf_rare_COR_edges_short_res$cortype <- ifelse(df_cf_rare_COR_edges_short_res$Weight > 0,  "pos", "neg")
    df_cf_rare_nodes_res <- select(df_cf_rare_COR_edges_short_res, c("Source", "cortype"))
    df_cf_rare_nodes_res = df_cf_rare_nodes_res[!duplicated(df_cf_rare_nodes_res$Source),]
    df_cf_rare_nodes_res <- data.frame(df_cf_rare_nodes_res)
    rownames(df_cf_rare_nodes_res) <- NULL
    colnames(df_cf_rare_nodes_res) <- c("Id", "CorType")
    df_cf_rare_nodes_res$species_type <- ifelse(df_cf_rare_nodes_res$Id %in% rownames(background_rare_h), "rare", "core")
    df_cf_rare_nodes_res$Genus <- df_cf_rare_nodes_res$Id
    df_cf_rare_nodes_res$Genus <- sapply(strsplit(as.character(
      df_cf_rare_nodes_res$Genus),"\\."), `[`, 1)
    
    # rescue cf, healthy background rare species ###
    cf_rare_edges_res <- df_cf_rare_COR_edges_short_res
    cf_rare_nodes_res <- df_cf_rare_nodes_res
    cf_rare_nodes_res$species_type <- as.factor(as.character(cf_rare_nodes_res$species_type))
    
    cf_rare_net_res <- graph_from_data_frame(d=cf_rare_edges_res, vertices=cf_rare_nodes_res, directed = directed_network)
    cf_rare_net_res <- simplify(cf_rare_net_res, remove.multiple = F, remove.loops = T)
    V(cf_rare_net_res)$color <- ifelse(V(cf_rare_net_res)$species_type == "rare", "lightsteelblue2", "orange2")
    
    cf.rare.rescue.net.plot <- layout.fruchterman.reingold(cf_rare_net_res)
    cf.rare.rescue.net.plot.df <- as.data.frame(cf.rare.rescue.net.plot)
    cf.rare.rescue.net.plot.df$Id <- cf_rare_nodes_res$Id
    cf.rare.rescue.net.plot.df$species_type <- cf_rare_nodes_res$species_type
    cf.rare.rescue.net.plot.df$node_size <- V(cf_rare_net_res)$size
    
    cf.rare.rescue.net.plot.edges <- get.data.frame(cf_rare_net_res)
    cf.rare.rescue.net.plot.edges$from.x <- cf.rare.rescue.net.plot.df$V1[
      match(cf.rare.rescue.net.plot.edges$from, cf.rare.rescue.net.plot.df$Id)]  
    cf.rare.rescue.net.plot.edges$from.y <- cf.rare.rescue.net.plot.df$V2[
      match(cf.rare.rescue.net.plot.edges$from, cf.rare.rescue.net.plot.df$Id)]
    cf.rare.rescue.net.plot.edges$to.x <- cf.rare.rescue.net.plot.df$V1[
      match(cf.rare.rescue.net.plot.edges$to, cf.rare.rescue.net.plot.df$Id)]  
    cf.rare.rescue.net.plot.edges$to.y <- cf.rare.rescue.net.plot.df$V2[
      match(cf.rare.rescue.net.plot.edges$to, cf.rare.rescue.net.plot.df$Id)]
    cf.rare.rescue.net.plot.edges$species_type <- ifelse(
      cf.rare.rescue.net.plot.edges$from %in% f_rare_species, "rare", "core")
    
    # CF, rescue healthy core
    sample_h_core <- sample(rownames(background_core_h), n_nodes, replace = sample_with_replacement)
    selected_core_species = rbind(selected_core_species, data.frame(i, n_nodes, sample_h_core))
    df_cf_core_rare_res <- subset(dataset_sub, 
                                  rownames(dataset_sub) %in% sample_h_core |
                                    rownames(dataset_sub) %in% rownames(background_core_cf) |
                                    rownames(dataset_sub) %in% rownames(background_rare_cf))
    df_cf_core_rare_res <- data.frame(t(df_cf_core_rare_res))
    df_cf_core_rare_res$state <- metadata_sub$state
    df_cf_core_rare_res <- subset(df_cf_core_rare_res, state == "CF")
    df_cf_core_rare_res$state <- NULL
    df_cf_core_rare_res <- rcorr(as.matrix(df_cf_core_rare_res), type = 'spearman')
    # create node and edge lists
    df_cf_core_rare_Pval_res <- df_cf_core_rare_res$P
    df_cf_core_rare_COR_res <- df_cf_core_rare_res$r
    # generate edgelist
    df_cf_core_rare_COR_edges_res <- melt(df_cf_core_rare_Pval_res)
    df_cf_core_rare_coredges_res <- melt(df_cf_core_rare_COR_res)
    df_cf_core_rare_COR_edges_res$COR <- df_cf_core_rare_coredges_res$value
    df_cf_core_rare_COR_edges_res$value <- round(df_cf_core_rare_COR_edges_res$value, 5)
    df_cf_core_rare_COR_edges_res$Label <- row.names(df_cf_core_rare_coredges_res)
    colnames(df_cf_core_rare_COR_edges_res) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
    rownames(df_cf_core_rare_COR_edges_res) <- NULL
    df_cf_core_rare_COR_edges_res <- df_cf_core_rare_COR_edges_res[
      complete.cases(df_cf_core_rare_COR_edges_res), ]
    df_cf_core_rare_COR_edges_short_res <- subset(df_cf_core_rare_COR_edges_res, 
                                                  pValue < sig_level & Weight >= weight_val | 
                                                    pValue < sig_level & Weight <= weight_val2)
    df_cf_core_rare_COR_edges_short_res$cortype <- ifelse(
      df_cf_core_rare_COR_edges_short_res$Weight > 0,  "pos", "neg")
    df_cf_core_rare_nodes_res <- select(
      df_cf_core_rare_COR_edges_short_res, c("Source", "cortype"))
    df_cf_core_rare_nodes_res = df_cf_core_rare_nodes_res[!duplicated(df_cf_core_rare_nodes_res$Source),]
    df_cf_core_rare_nodes_res <- data.frame(df_cf_core_rare_nodes_res)
    rownames(df_cf_core_rare_nodes_res) <- NULL
    colnames(df_cf_core_rare_nodes_res) <- c("Id", "CorType")
    df_cf_core_rare_nodes_res$species_type <-
      ifelse(df_cf_core_rare_nodes_res$Id %in% rownames(background_rare_cf), "rare", "core")
    df_cf_core_rare_nodes_res$Genus <- df_cf_core_rare_nodes_res$Id
    df_cf_core_rare_nodes_res$Genus <- sapply(strsplit(as.character(
      df_cf_core_rare_nodes_res$Genus),"\\."), `[`, 1)
    
    # rescue cf, healthy background core species ###
    cf_core_rare_edges_res <- df_cf_core_rare_COR_edges_short_res
    cf_core_rare_nodes_res <- df_cf_core_rare_nodes_res
    cf_core_rare_nodes_res$species_type <- as.factor(as.character(cf_core_rare_nodes_res$species_type))
    
    cf_core_rare_net_res <- graph_from_data_frame(d=cf_core_rare_edges_res, vertices=cf_core_rare_nodes_res, 
                                                  directed = directed_network)
    cf_core_rare_net_res <- simplify(cf_core_rare_net_res, remove.multiple = F, remove.loops = T)
    # Generate colors based on media type:
    V(cf_core_rare_net_res)$color <- ifelse(V(cf_core_rare_net_res)$species_type == "rare", "lightsteelblue2", "orange2")
    cf.core.rescue.net.plot <- layout.fruchterman.reingold(cf_core_rare_net_res)
    cf.core.rescue.net.plot.df <- as.data.frame(cf.core.rescue.net.plot)
    cf.core.rescue.net.plot.df$Id <- cf_core_rare_nodes_res$Id
    cf.core.rescue.net.plot.df$species_type <- cf_core_rare_nodes_res$species_type
    cf.core.rescue.net.plot.df$node_size <- V(cf_core_rare_net_res)$size
    
    cf.core.rescue.net.plot.edges <- get.data.frame(cf_core_rare_net_res)
    cf.core.rescue.net.plot.edges$from.x <- cf.core.rescue.net.plot.df$V1[
      match(cf.core.rescue.net.plot.edges$from, cf.core.rescue.net.plot.df$Id)]  
    cf.core.rescue.net.plot.edges$from.y <- cf.core.rescue.net.plot.df$V2[
      match(cf.core.rescue.net.plot.edges$from, cf.core.rescue.net.plot.df$Id)]
    cf.core.rescue.net.plot.edges$to.x <- cf.core.rescue.net.plot.df$V1[
      match(cf.core.rescue.net.plot.edges$to, cf.core.rescue.net.plot.df$Id)]  
    cf.core.rescue.net.plot.edges$to.y <- cf.core.rescue.net.plot.df$V2[
      match(cf.core.rescue.net.plot.edges$to, cf.core.rescue.net.plot.df$Id)]
    cf.core.rescue.net.plot.edges$species_type <- ifelse(
      cf.core.rescue.net.plot.edges$from %in% f_rare_species, "rare", "core")
    
    
    # # CF, rescue half healthy core and half rare
    set.seed(i)
    n_nodes_half0 <- n_nodes / 2
    n_nodes_half1 <- ceiling(n_nodes_half0)
    sample_core <- sample(rownames(background_core_h), n_nodes_half1, replace = sample_with_replacement)
    sample_rare <- sample(rownames(background_rare_h), n_nodes_half1, replace = sample_with_replacement)
    selected_bothr_species = rbind(selected_bothr_species, data.frame(i, n_nodes, sample_rare))
    selected_bothc_species = rbind(selected_bothc_species, data.frame(i, n_nodes, sample_core))
    df_cf_res_both <- subset(dataset_sub, 
                             rownames(dataset_sub) %in% sample_core |
                               rownames(dataset_sub) %in% sample_rare |
                               rownames(dataset_sub) %in% rownames(background_core_cf) |
                               rownames(dataset_sub) %in% rownames(background_rare_cf))
    df_cf_res_both <- data.frame(t(df_cf_res_both))
    df_cf_res_both$state <- metadata_sub$state
    df_cf_res_both <- subset(df_cf_res_both, state == "CF")
    df_cf_res_both$state <- NULL
    df_cf_res_both <- rcorr(as.matrix(df_cf_res_both), type = 'spearman')
    # create node and edge lists
    df_cf_res_both_Pval_res <- df_cf_res_both$P
    df_cf_res_both_COR_res <- df_cf_res_both$r
    # generate edgelist
    df_cf_res_both_COR_edges_res <- melt(df_cf_res_both_Pval_res)
    df_cf_res_both_coredges_res <- melt(df_cf_res_both_COR_res)
    df_cf_res_both_COR_edges_res$COR <- df_cf_res_both_coredges_res$value
    df_cf_res_both_COR_edges_res$value <- round(df_cf_res_both_COR_edges_res$value, 5)
    df_cf_res_both_COR_edges_res$Label <- row.names(df_cf_res_both_coredges_res)
    colnames(df_cf_res_both_COR_edges_res) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
    rownames(df_cf_res_both_COR_edges_res) <- NULL
    df_cf_res_both_COR_edges_res <- df_cf_res_both_COR_edges_res[complete.cases(df_cf_res_both_COR_edges_res), ]
    df_cf_res_both_COR_edges_short_res <- subset(df_cf_res_both_COR_edges_res, 
                                                 pValue < sig_level & Weight >= weight_val | 
                                                   pValue < sig_level & Weight <= weight_val2)
    df_cf_res_both_COR_edges_short_res$cortype <- ifelse(
      df_cf_res_both_COR_edges_short_res$Weight > 0,  "pos", "neg")
    df_cf_res_both_nodes_res <- select(df_cf_res_both_COR_edges_short_res, c("Source", "cortype"))
    df_cf_res_both_nodes_res = df_cf_res_both_nodes_res[!duplicated(df_cf_res_both_nodes_res$Source),]
    df_cf_res_both_nodes_res <- data.frame(df_cf_res_both_nodes_res)
    rownames(df_cf_res_both_nodes_res) <- NULL
    colnames(df_cf_res_both_nodes_res) <- c("Id", "CorType")
    df_cf_res_both_nodes_res$species_type <- ifelse(df_cf_res_both_nodes_res$Id %in% rownames(background_rare_h), "rare", "core")
    df_cf_res_both_nodes_res$Genus <- df_cf_res_both_nodes_res$Id
    df_cf_res_both_nodes_res$Genus <- sapply(strsplit(as.character(
      df_cf_res_both_nodes_res$Genus),"\\."), `[`, 1)
    
    # rescue cf, healthy background core AND rare species ###
    df_cf_res_both_nodes_res$species_type <- as.factor(as.character(df_cf_res_both_nodes_res$species_type))
    cf_both_net_res <- graph_from_data_frame(d=df_cf_res_both_COR_edges_short_res, 
                                             vertices=df_cf_res_both_nodes_res, directed = directed_network)
    cf_both_net_res <- simplify(cf_both_net_res, remove.multiple = F, remove.loops = T)
    # Generate colors based on media type:
    V(cf_both_net_res)$color <- ifelse(V(cf_both_net_res)$species_type == "rare", "lightsteelblue2", "orange2")
    cf.both.rescue.net.plot <- layout.fruchterman.reingold(cf_both_net_res)
    cf.both.rescue.net.plot.df <- as.data.frame(cf.both.rescue.net.plot)
    cf.both.rescue.net.plot.df$Id <- df_cf_res_both_nodes_res$Id
    cf.both.rescue.net.plot.df$species_type <- df_cf_res_both_nodes_res$species_type
    cf.both.rescue.net.plot.df$node_size <- V(cf_both_net_res)$size
    
    cf.both.rescue.net.plot.edges <- get.data.frame(cf_both_net_res)
    cf.both.rescue.net.plot.edges$from.x <- cf.both.rescue.net.plot.df$V1[
      match(cf.both.rescue.net.plot.edges$from, cf.both.rescue.net.plot.df$Id)]  
    cf.both.rescue.net.plot.edges$from.y <- cf.both.rescue.net.plot.df$V2[
      match(cf.both.rescue.net.plot.edges$from, cf.both.rescue.net.plot.df$Id)]
    cf.both.rescue.net.plot.edges$to.x <- cf.both.rescue.net.plot.df$V1[
      match(cf.both.rescue.net.plot.edges$to, cf.both.rescue.net.plot.df$Id)]  
    cf.both.rescue.net.plot.edges$to.y <- cf.both.rescue.net.plot.df$V2[
      match(cf.both.rescue.net.plot.edges$to, cf.both.rescue.net.plot.df$Id)]
    cf.both.rescue.net.plot.edges$species_type <- ifelse(
      cf.both.rescue.net.plot.edges$from %in% f_rare_species, "rare", "core")
    
    
    # Attack statistics (healthy rare)
    cf_rescue_net_attack <- swan_combinatory(cf_rare_net_res,10)
    cf_rescue_net_attack <- data.frame(cf_rescue_net_attack)
    colnames(cf_rescue_net_attack) <- c("fraction_removed", "betweeness", "degree", "cascading", "random")
    cf_rescue_net_attack$fraction_removed <- round(cf_rescue_net_attack$fraction_removed,1)
    cf_rescue_net_attack <- ddply(cf_rescue_net_attack,"fraction_removed",numcolwise(median))
    cf_rescue_net_attack <- rbind(h_rare_net_attack[1,], cf_rescue_net_attack)
    
    cf_rescue_core_rare_net_attack <- swan_combinatory(cf_core_rare_net_res,10)
    cf_rescue_core_rare_net_attack <- data.frame(cf_rescue_core_rare_net_attack)
    colnames(cf_rescue_core_rare_net_attack) <- c("fraction_removed", "betweeness", "degree", "cascading", "random")
    cf_rescue_core_rare_net_attack$fraction_removed <- round(cf_rescue_core_rare_net_attack$fraction_removed,1)
    cf_rescue_core_rare_net_attack <- ddply(cf_rescue_core_rare_net_attack,"fraction_removed",numcolwise(median))
    cf_rescue_core_rare_net_attack <- rbind(h_rare_net_attack[1,], cf_rescue_core_rare_net_attack)
    
    cf_rescue_both_net_attack <- swan_combinatory(cf_both_net_res,10)
    cf_rescue_both_net_attack <- data.frame(cf_rescue_both_net_attack)
    colnames(cf_rescue_both_net_attack) <- c("fraction_removed", "betweeness", "degree", "cascading", "random")
    cf_rescue_both_net_attack$fraction_removed <- round(cf_rescue_both_net_attack$fraction_removed,1)
    cf_rescue_both_net_attack <- ddply(cf_rescue_both_net_attack,"fraction_removed",numcolwise(median))
    cf_rescue_both_net_attack <- rbind(cf_both_net_res[1,], cf_rescue_both_net_attack)
    
    # rescue CF (healthy rare), degree
    cf_rare_res_matrix_degree <- select(cf_rescue_net_attack, select_items)
    cf_rare_res_matrix_degree <- as.matrix(cf_rare_res_matrix_degree)
    h_rare_net_matrix_degree <- select(h_rare_net_attack, select_items)
    h_rare_net_matrix_degree <- as.matrix(h_rare_net_matrix_degree)
    
    # rescue CF (healthy core), degree
    cf_core_res_matrix_degree <- select(cf_rescue_core_rare_net_attack, select_items)
    cf_core_res_matrix_degree <- as.matrix(cf_core_res_matrix_degree)
    h_rare_net_matrix_degree <- select(h_rare_net_attack, select_items)
    h_rare_net_matrix_degree <- as.matrix(h_rare_net_matrix_degree)
    
    # rescue CF (both), degree
    cf_both_res_matrix_degree <- select(cf_rescue_both_net_attack, select_items)
    cf_both_res_matrix_degree <- as.matrix(cf_both_res_matrix_degree)
    cf_both_res_matrix_degree <- cf_both_res_matrix_degree[-1,]
    h_rare_net_matrix_degree <- select(h_rare_net_attack, select_items)
    h_rare_net_matrix_degree <- as.matrix(h_rare_net_matrix_degree)
    
    # rescue CF (healthy rare), random
    cf_rare_res_matrix_random <- select(cf_rescue_net_attack, select_items2)
    cf_rare_res_matrix_random <- as.matrix(cf_rare_res_matrix_random)
    h_rare_net_matrix_random <- select(h_rare_net_attack, select_items2)
    h_rare_net_matrix_random <- as.matrix(h_rare_net_matrix_random)
    
    # rescue CF (healthy core), random
    cf_core_res_matrix_random <- select(cf_rescue_core_rare_net_attack, select_items2)
    cf_core_res_matrix_random <- as.matrix(cf_core_res_matrix_random)
    h_rare_net_matrix_random <- select(h_rare_net_attack, select_items2)
    h_rare_net_matrix_random <- as.matrix(h_rare_net_matrix_random)
    
    # rescue CF (both), random
    cf_both_res_matrix_random <- select(cf_rescue_both_net_attack, select_items2)
    cf_both_res_matrix_random <- as.matrix(cf_both_res_matrix_random)
    cf_both_res_matrix_random <- cf_both_res_matrix_random[-1,]
    h_rare_net_matrix_random <- select(h_rare_net_attack, select_items2)
    h_rare_net_matrix_random <- as.matrix(h_rare_net_matrix_random)
    
    my_list1 <- list(h_rare_net, cf_rare_net, cf_rare_net_res, cf_core_rare_net_res, cf_both_net_res)
    my_K <- CalculateKStepRandomWalkKernel(my_list1, rep(1,2))
    my_K_df <- data.frame(my_K)
    my_K_df <- my_K_df / my_K_df[1,1]
    colnames(my_K_df) <- c(
      "Healthy", "CF", "CF + Healthy (rare)", "CF + Healthy (core)", "CF + Healthy (both)")
    rownames(my_K_df) <- c(
      "Healthy", "CF", "CF + Healthy (rare)", "CF + Healthy (core)", "CF + Healthy (both)")
    my_K_df <- my_K_df[1,]
    
    my_K2 <- CalculateShortestPathKernel(my_list1)
    my_K_df2 <- data.frame(my_K2)
    my_K_df2 <- my_K_df2 / my_K_df2[1,1]
    colnames(my_K_df2) <- c(
      "Healthy", "CF", "CF + Healthy (rare)", "CF + Healthy (core)", "CF + Healthy (both)")
    rownames(my_K_df2) <- c(
      "Healthy", "CF", "CF + Healthy (rare)", "CF + Healthy (core)", "CF + Healthy (both)")
    my_K_df2 <- my_K_df2[1,]
    
    my_K3 <-CalculateWLKernel(my_list1, 5)
    my_K_df3 <- data.frame(my_K3)
    my_K_df3 <- my_K_df3 / 	my_K_df3[1,1]
    colnames(my_K_df3) <- c(
      "Healthy", "CF", "CF + Healthy (rare)", "CF + Healthy (core)", "CF + Healthy (both)")
    rownames(my_K_df3) <- c(
      "Healthy", "CF", "CF + Healthy (rare)", "CF + Healthy (core)", "CF + Healthy (both)")
    my_K_df3 <- my_K_df3[1,]
    
    #my_K <- data.frame(rbind(my_K_df, my_K_df2, my_K_df3))
    #"Random Walk Kernel", "Shortest-Path Kernel", "Weisfeiler-Lehman Graph Kernel")
    netwerk_kernel_randomWalk <- rbind(netwerk_kernel_randomWalk, data.frame(i, n_nodes, my_K_df))
    network_kernel_shortestPath <- rbind(network_kernel_shortestPath, data.frame(i, n_nodes, my_K_df2))
    network_kernel_WL <- rbind(network_kernel_WL, data.frame(i, n_nodes, my_K_df3))
    
    cf_hrare_deg_frechet <- Frechet(h_rare_net_matrix_degree,cf_rare_res_matrix_degree)
    cf_hcore_deg_frechet <- Frechet(h_rare_net_matrix_degree,cf_core_res_matrix_degree)
    cf_both_deg_frechet <- Frechet(h_rare_net_matrix_degree,cf_both_res_matrix_degree)
    cf_hrare_rand_frechet <-Frechet(h_rare_net_matrix_random,cf_rare_res_matrix_random)
    cf_hcore_rand_frechet <-Frechet(h_rare_net_matrix_random,cf_core_res_matrix_random)
    cf_both_rand_frechet <-Frechet(h_rare_net_matrix_random,cf_both_res_matrix_random)
    frechet_attack_similarity = rbind(frechet_attack_similarity, 
                                      data.frame(i, n_nodes, cf_hrare_deg_frechet, cf_hcore_deg_frechet,
                                                 cf_both_deg_frechet, cf_hrare_rand_frechet,
                                                 cf_hcore_rand_frechet, cf_both_rand_frechet))
    
    cf_hrare_cong <- test_congruence(h_rare_net, cf_rare_net_res, method="distance")
    cf_hrare_cong_p <- cf_hrare_cong$p.value
    cf_hrare_cong_ci <- cf_hrare_cong$conf.int
    cf_hrare_cong_estimate <- cf_hrare_cong$estimate
    cf_hcore_cong <- test_congruence(cf_core_rare_net_res, h_rare_net,  method="distance")
    cf_hcore_cong_p <- cf_hcore_cong$p.value
    cf_hcore_cong_ci <- cf_hcore_cong$conf.int
    cf_hcore_cong_estimate <- cf_hcore_cong$estimate
    cf_both_cong <- test_congruence(h_rare_net, cf_both_net_res, method="distance")
    cf_both_cong_p <- cf_both_cong$p.value
    cf_both_cong_ci <- cf_both_cong$conf.int
    cf_both_cong_estimate <- cf_both_cong$estimate
    network_congruency <- rbind(network_congruency,
                                data.frame(i, n_nodes,
                                           cbind(cf_hrare_cong_p, cf_hrare_cong_estimate, cf_hrare_cong_ci),
                                           cbind(cf_hcore_cong_p, cf_hcore_cong_estimate, cf_hcore_cong_ci),
                                           cbind(cf_both_cong_p, cf_both_cong_estimate, cf_both_cong_ci)))
    
    # degree
    degree_cf_hrare <- sum(degree(cf_rare_net_res))
    degree_cf_hcore <- sum(degree(cf_core_rare_net_res))
    degree_cf_both <- sum(degree(cf_both_net_res))
    
    # number of nodes
    nNodes_cf_hrare <- gorder(cf_rare_net_res) 
    nNodes_cf_hcore <- gorder(cf_core_rare_net_res) 
    nNodes_cf_both <- gorder(cf_both_net_res)
    
    # number of edges
    nEdges_cf_hrare <- gsize(cf_rare_net_res) 
    nEdges_cf_hcore <- gsize(cf_core_rare_net_res) 
    nEdges_cf_both <- gsize(cf_both_net_res) 
    
    # network diameter
    diameter_cf_hrare <- diameter(cf_rare_net_res) 
    diameter_cf_hcore <- diameter(cf_core_rare_net_res) 
    diameter_cf_both <- diameter(cf_both_net_res) 
    
    # edge density
    density_cf_hrare <- edge_density(cf_rare_net_res) 
    density_cf_hcore <- edge_density(cf_core_rare_net_res) 
    density_cf_both <- edge_density(cf_both_net_res) 
    
    # transitivity
    trans_cf_hrare <- transitivity(cf_rare_net_res) 
    trans_cf_hcore <- transitivity(cf_core_rare_net_res) 
    trans_cf_both <- transitivity(cf_both_net_res)
    
    # eigencentrality
    eigen_cf_hrare <- centr_eigen(cf_rare_net_res)$centralization 
    eigen_cf_hcore <- centr_eigen(cf_core_rare_net_res)$centralization 
    eigen_cf_both <- centr_eigen(cf_both_net_res)$centralization
    
    network_statistics <- rbind(network_statistics, data.frame(
      i, n_nodes, degree_cf_hrare, degree_cf_hcore, degree_cf_both,
      nNodes_cf_hrare, nNodes_cf_hcore, nNodes_cf_both,
      nEdges_cf_hrare, nEdges_cf_hcore, nEdges_cf_both,
      diameter_cf_hrare, diameter_cf_hcore, diameter_cf_both,
      density_cf_hrare, density_cf_hcore, density_cf_both,
      trans_cf_hrare, trans_cf_hcore, trans_cf_both,
      eigen_cf_hrare, eigen_cf_hcore, eigen_cf_both))
  }
}

frechet_attack_similarity_median <-
  ddply(frechet_attack_similarity,"n_nodes",numcolwise(median))
colnames(frechet_attack_similarity_median) <- c(
  "n_nodes_med", "i_med","cf_hrare_deg_frechet_med", "cf_hcore_deg_frechet_med","cf_both_deg_frechet_med", 
  "cf_hrare_rand_frechet_med", "cf_hcore_rand_frechet_med", "cf_both_rand_frechet_med")
frechet_attack_similarity_median$cf_vs_healthy_degree <- original_frechet_h_cf_degree
frechet_attack_similarity_median$cf_vs_healthy_random <- original_frechet_h_cf_random
frechet_attack_similarity_minimum <-
  ddply(frechet_attack_similarity,"n_nodes",numcolwise(min))
colnames(frechet_attack_similarity_minimum) <- c(
  "n_nodes_min", "i_min","cf_hrare_deg_frechet_min", "cf_hcore_deg_frechet_min","cf_both_deg_frechet_min", 
  "cf_hrare_rand_frechet_min", "cf_hcore_rand_frechet_min", "cf_both_rand_frechet_min")

frechet_attack_similarity_maximum <-
  ddply(frechet_attack_similarity,"n_nodes",numcolwise(max))
colnames(frechet_attack_similarity_maximum) <- c(
  "n_nodes_max", "i_max","cf_hrare_deg_frechet_max", "cf_hcore_deg_frechet_max","cf_both_deg_frechet_max", 
  "cf_hrare_rand_frechet_max", "cf_hcore_rand_frechet_max", "cf_both_rand_frechet_max")

frechet_attack_similartiy_degree <- data.frame(cbind(
  frechet_attack_similarity_median, frechet_attack_similarity_minimum,frechet_attack_similarity_maximum))

frechet_TargetedAttack_similarity_plot <-ggplot(data=frechet_attack_similartiy_degree) +
  geom_pointrange(aes(x=n_nodes_med, y=cf_hrare_deg_frechet_med, ymin=cf_hrare_deg_frechet_min,
                      ymax=cf_hrare_deg_frechet_max), color="darkblue") +
  geom_line(aes(x=n_nodes_med, y=cf_hrare_deg_frechet_med), color="darkblue", size=1, alpha=0.5) +
  geom_pointrange(aes(x=n_nodes_med, y=cf_hcore_deg_frechet_med, ymin=cf_hcore_deg_frechet_min,
                      ymax=cf_hcore_deg_frechet_max), color="orange2") +
  geom_line(aes(x=n_nodes_med, y=cf_hcore_deg_frechet_med), color="orange2", size=1, alpha=0.5) +
  geom_pointrange(aes(x=n_nodes_med, y=cf_both_deg_frechet_med, ymin=cf_both_deg_frechet_min,
                      ymax=cf_both_deg_frechet_max), color="darkgreen") +
  geom_line(aes(x=n_nodes_med, y=cf_both_deg_frechet_med), color="darkgreen", size=1, alpha=0.5) +
  geom_line(aes(x=n_nodes_med, y=cf_vs_healthy_degree), color="black", size=1, alpha=1, linetype="dashed") +
  xlab("Number of added nodes") + ylim(0.0,0.8) + ylab("Frechet distance\n") + theme_pubr(border=TRUE) 

colors <- c("CF + Healthy (rare)" = "darkblue", "CF + Healthy (core)" = "orange2", 
            "CF + Healthy (both)" = "darkgreen", "CF" = "black")
frechet_RandomAttack_similarity_plot <- ggplot(data=frechet_attack_similartiy_degree) +
  geom_pointrange(aes(x=n_nodes_med, y=cf_hrare_rand_frechet_med, ymin=cf_hrare_rand_frechet_min,
                      ymax=cf_hrare_rand_frechet_max, color = "CF + Healthy (rare)")) +
  geom_line(aes(x=n_nodes_med, y=cf_hrare_rand_frechet_med, color = "CF + Healthy (rare)"), size=1, alpha=0.5) +
  geom_pointrange(aes(x=n_nodes_med, y=cf_hcore_rand_frechet_med, ymin=cf_hcore_rand_frechet_min,
                      ymax=cf_hcore_rand_frechet_max, color="CF + Healthy (core)")) +
  geom_line(aes(x=n_nodes_med, y=cf_hcore_rand_frechet_med, color="CF + Healthy (core)"), size=1, alpha=0.5) +
  geom_pointrange(aes(x=n_nodes_med, y=cf_both_rand_frechet_med, ymin=cf_both_rand_frechet_min,
                      ymax=cf_both_rand_frechet_max, color="CF + Healthy (both)")) +
  geom_line(aes(x=n_nodes_med, y=cf_both_rand_frechet_med, color="CF + Healthy (both)"), size=1, alpha=0.5) +
  geom_line(aes(x=n_nodes_med, y=cf_vs_healthy_random, color="CF"), size=1, alpha=1, linetype="dashed") +
  xlab("Number of added nodes") + ylim(0.0,0.8) + ylab("Frechet distance\n") + theme_pubr(border=TRUE)  +
  scale_color_manual(values = colors) + theme(legend.position = c(0.7, 0.87), legend.title = element_blank())


# Network kernel (Random Walk)
kernel_randomWalk_median <- ddply(netwerk_kernel_randomWalk,"n_nodes",numcolwise(median))
colnames(kernel_randomWalk_median) <- c("n_nodes_med", "i_med","Healthy_med", "CF_med", 
                                        "CF_Hrare_med", "CF_Hcore_med", "CF_Hboth_med")

kernel_randomWalk_min <- ddply(netwerk_kernel_randomWalk,"n_nodes",numcolwise(min))
colnames(kernel_randomWalk_min) <- c("n_nodes_min", "i_min","Healthy_min", "CF_min",
                                     "CF_Hrare_min", "CF_Hcore_min", "CF_Hboth_min")

kernel_randomWalk_max <- ddply(netwerk_kernel_randomWalk,"n_nodes",numcolwise(max))
colnames(kernel_randomWalk_max) <- c("n_nodes_max", "i_max","Healthy_max", "CF_max",
                                     "CF_Hrare_max", "CF_Hcore_max", "CF_Hboth_max")

kernel_randomWalk_df <- data.frame(cbind(kernel_randomWalk_median, kernel_randomWalk_min,kernel_randomWalk_max))

kernel_randomWalk_plot <- ggplot(data=kernel_randomWalk_df) + 
  geom_pointrange(aes(x=n_nodes_med, y=CF_Hrare_med, ymin=CF_Hrare_min, ymax=CF_Hrare_max), color="darkblue") +
  geom_line(aes(x=n_nodes_med, y=CF_Hrare_med), color="darkblue", size=1, alpha=0.5) +
  geom_pointrange(aes(x=n_nodes_med, y=CF_Hcore_med, ymin=CF_Hcore_min,ymax=CF_Hcore_max), color="orange2") +
  geom_line(aes(x=n_nodes_med, y=CF_Hcore_med), color="orange2", size=1, alpha=0.5) +
  geom_pointrange(aes(x=n_nodes_med, y=CF_Hboth_med, ymin=CF_Hboth_min, ymax=CF_Hboth_max), color="darkgreen") +
  geom_line(aes(x=n_nodes_med, y=CF_Hboth_med), color="darkgreen", size=1, alpha=0.5) +
  geom_line(aes(x=n_nodes_med, y=CF_med), color="black", size=1, alpha=1, linetype="dashed") +
  xlab(" ") + ylim(0.0,0.6) + ylab("Similarity random walk kernel\n") + theme_pubr(border=TRUE) 


# shortest pathway kernel
kernel_shortestPahway_median <- ddply(network_kernel_shortestPath,"n_nodes",numcolwise(median))
colnames(kernel_shortestPahway_median) <- c("n_nodes_med", "i_med","Healthy_med", 
                                            "CF_med","CF_Hrare_med", "CF_Hcore_med", "CF_Hboth_med")

kernel_shortestPahway_min <- ddply(network_kernel_shortestPath,"n_nodes",numcolwise(min))
colnames(kernel_shortestPahway_min) <- c("n_nodes_min", "i_min","Healthy_min", "CF_min",
                                         "CF_Hrare_min", "CF_Hcore_min", "CF_Hboth_min")

kernel_shortestPahway_max <- ddply(network_kernel_shortestPath,"n_nodes",numcolwise(max))
colnames(kernel_shortestPahway_max) <- c("n_nodes_max", "i_max","Healthy_max", "CF_max",
                                         "CF_Hrare_max", "CF_Hcore_max", "CF_Hboth_max")

kernel_shortestPahway_df <- data.frame(cbind(kernel_shortestPahway_median, kernel_shortestPahway_min,kernel_shortestPahway_max))

kernel_shortestPahway_plot <- ggplot(data=kernel_shortestPahway_df) +
  geom_pointrange(aes(x=n_nodes_med, y=CF_Hrare_med, ymin=CF_Hrare_min, ymax=CF_Hrare_max), color="darkblue") +
  geom_line(aes(x=n_nodes_med, y=CF_Hrare_med), color="darkblue", size=1, alpha=0.5) +
  geom_pointrange(aes(x=n_nodes_med, y=CF_Hcore_med, ymin=CF_Hcore_min, ymax=CF_Hcore_max), color="orange2") +
  geom_line(aes(x=n_nodes_med, y=CF_Hcore_med), color="orange2", size=1, alpha=0.5) +
  geom_pointrange(aes(x=n_nodes_med, y=CF_Hboth_med, ymin=CF_Hboth_min, ymax=CF_Hboth_max), color="darkgreen") +
  geom_line(aes(x=n_nodes_med, y=CF_Hboth_med), color="darkgreen", size=1, alpha=0.5) +
  geom_line(aes(x=n_nodes_med, y=CF_med), color="black", size=1, alpha=1, linetype="dashed") +
  xlab(" ") + ylab("Similarity shortest pathway kernel\n") + ylim(0.0,0.6) + theme_pubr(border=TRUE) 

# Weisfeiler-Lehman Graph Kernels
kernel_WL_median <- ddply(network_kernel_WL,"n_nodes",numcolwise(median))
colnames(kernel_WL_median) <- c("n_nodes_med", "i_med","Healthy_med", 
                                "CF_med","CF_Hrare_med", "CF_Hcore_med", "CF_Hboth_med")

kernel_WL_min <- ddply(network_kernel_WL,"n_nodes",numcolwise(min))
colnames(kernel_WL_min) <- c( "n_nodes_min", "i_min","Healthy_min", 
                              "CF_min","CF_Hrare_min", "CF_Hcore_min", "CF_Hboth_min")

kernel_WL_max <- ddply(network_kernel_WL,"n_nodes",numcolwise(max))
colnames(kernel_WL_max) <- c("n_nodes_max", "i_max","Healthy_max", 
                             "CF_max","CF_Hrare_max", "CF_Hcore_max", "CF_Hboth_max")

kernel_WL_df <- data.frame(cbind(kernel_WL_median, kernel_WL_min,kernel_WL_max))

kernel_WL_plot <- ggplot(data=kernel_WL_df) + 
  geom_pointrange(aes(x=n_nodes_med,  y=CF_Hrare_med, ymin=CF_Hrare_min, ymax=CF_Hrare_max), color="darkblue") +
  geom_line(aes(x=n_nodes_med, y=CF_Hrare_med), color="darkblue", size=1, alpha=0.5) +
  geom_pointrange(aes(x=n_nodes_med, y=CF_Hcore_med, ymin=CF_Hcore_min, ymax=CF_Hcore_max), color="orange2") +
  geom_line(aes(x=n_nodes_med, y=CF_Hcore_med), color="orange2", size=1, alpha=0.5) +
  geom_pointrange(aes(x=n_nodes_med, y=CF_Hboth_med, ymin=CF_Hboth_min, ymax=CF_Hboth_max), color="darkgreen") +
  geom_line(aes(x=n_nodes_med, y=CF_Hboth_med), color="darkgreen", size=1, alpha=0.5) +
  geom_line(aes(x=n_nodes_med, y=CF_med), color="black", size=1, alpha=1, linetype="dashed") +
  ylim(0.0,0.6) + xlab(" ") + ylab("Similarity Weisfeiler-Lehman kernel\n") + theme_pubr(border=TRUE) 

modulated_network_similarities <- ggarrange(kernel_shortestPahway_plot, kernel_WL_plot,
                                            frechet_TargetedAttack_similarity_plot,
                                            frechet_RandomAttack_similarity_plot, labels=c("a", "b", "c", "d"),
                                            common.legend = TRUE, legend = "bottom")

rownames(frechet_attack_similarity) <- NULL
colnames(frechet_attack_similarity) <- c("Seed", "N", "cf_hrare_deg", "cf_hcore_deg", "cf_hboth_deg",
                                         "cf_hrare_rand", "CF_hcore_rand", "cf_hboth_rand")

# Both (core and rare species) choice on network performance (degree)
cf_hboth <- data.frame(cbind(selected_bothc_species, selected_bothr_species$sample_rare))
colnames(cf_hboth) <- c("i", "n_nodes", "sample_core", "sample_rare")

cf_hboth_n4 <- subset(cf_hboth, n_nodes == 4)
frechet_hboth_deg_n4 <- subset(frechet_attack_similarity, N == 4)
frechet_hboth_deg_n4 <- frechet_hboth_deg_n4[rep(seq_len(nrow(frechet_hboth_deg_n4)), each = 2), ]
cf_hboth_n4_final <- data.frame(cbind(cf_hboth_n4, frechet_hboth_deg_n4$cf_hboth_deg))
colnames(cf_hboth_n4_final) <- c("i", "n_nodes", "sample_core", "sample_rare", "frechet_deg")
cf_hboth_n4_final$cat <- ifelse(cf_hboth_n4_final$frechet_deg < original_frechet_h_cf_degree,
                                "better", "worse")
#subset worse
cf_hboth_n4_final_worse <- subset(cf_hboth_n4_final, cat == "worse")
cf_hboth_n4_final_worse$sample_core <- as.factor(as.character(cf_hboth_n4_final_worse$sample_core))
cf_hboth_n4_final_worse$sample_rare <- as.factor(as.character(cf_hboth_n4_final_worse$sample_rare))
n_nodes_4_worse_core <- table(cf_hboth_n4_final_worse$sample_core)
n_nodes_4_worse_core <- data.frame(n_nodes_4_worse_core)
n_nodes_4_worse_core$which <- "worse"
n_nodes_4_worse_core$n <- 4
n_nodes_4_worse_core$Perc <- (n_nodes_4_worse_core$Freq / sum(n_nodes_4_worse_core$Freq)) * 100
n_nodes_4_worse_rare <- table(cf_hboth_n4_final_worse$sample_rare)
n_nodes_4_worse_rare <- data.frame(n_nodes_4_worse_rare)
n_nodes_4_worse_rare$which <- "worse"
n_nodes_4_worse_rare$n <- 4
n_nodes_4_worse_rare$Perc <- (n_nodes_4_worse_rare$Freq / sum(n_nodes_4_worse_rare$Freq)) * 100
# subset better
cf_hboth_n4_final_better <- subset(cf_hboth_n4_final, cat == "better")
cf_hboth_n4_final_better$sample_core <- as.factor(as.character(cf_hboth_n4_final_better$sample_core))
cf_hboth_n4_final_better$sample_rare <- as.factor(as.character(cf_hboth_n4_final_better$sample_rare))
n_nodes_4_better_core <- table(cf_hboth_n4_final_better$sample_core)
n_nodes_4_better_core <- data.frame(n_nodes_4_better_core)
n_nodes_4_better_core$which <- "better"
n_nodes_4_better_core$n <- 4
n_nodes_4_better_core$Perc <- (n_nodes_4_better_core$Freq / sum(n_nodes_4_better_core$Freq)) * 100
n_nodes_4_better_rare <- table(cf_hboth_n4_final_better$sample_rare)
n_nodes_4_better_rare <- data.frame(n_nodes_4_better_rare)
n_nodes_4_better_rare$which <- "better"
n_nodes_4_better_rare$n <- 4
n_nodes_4_better_rare$Perc <- (n_nodes_4_better_rare$Freq / sum(n_nodes_4_better_rare$Freq)) * 100

cf_hboth_n6 <- subset(cf_hboth, n_nodes == 6)
frechet_hboth_deg_n6 <- subset(frechet_attack_similarity, N == 6)
frechet_hboth_deg_n6 <- frechet_hboth_deg_n6[rep(seq_len(nrow(frechet_hboth_deg_n6)), each = 3), ]
cf_hboth_n6_final <- data.frame(cbind(cf_hboth_n6, frechet_hboth_deg_n6$cf_hboth_deg))
colnames(cf_hboth_n6_final) <- c("i", "n_nodes", "sample_core", "sample_rare", "frechet_deg")
cf_hboth_n6_final$cat <- ifelse(cf_hboth_n6_final$frechet_deg < original_frechet_h_cf_degree, "better", "worse")

#subset worse
cf_hboth_n6_final_worse <- subset(cf_hboth_n6_final, cat == "worse")
cf_hboth_n6_final_worse$sample_core <- as.factor(as.character(cf_hboth_n6_final_worse$sample_core))
cf_hboth_n6_final_worse$sample_rare <- as.factor(as.character(cf_hboth_n6_final_worse$sample_rare))
n_nodes_6_worse_core <- table(cf_hboth_n6_final_worse$sample_core)
n_nodes_6_worse_core <- data.frame(n_nodes_6_worse_core)
n_nodes_6_worse_core$which <- "worse"
n_nodes_6_worse_core$n <- 6
n_nodes_6_worse_core$Perc <- (n_nodes_6_worse_core$Freq / sum(n_nodes_6_worse_core$Freq)) * 100
n_nodes_6_worse_rare <- table(cf_hboth_n6_final_worse$sample_rare)
n_nodes_6_worse_rare <- data.frame(n_nodes_6_worse_rare)
n_nodes_6_worse_rare$which <- "worse"
n_nodes_6_worse_rare$n <- 6
n_nodes_6_worse_rare$Perc <- (n_nodes_6_worse_rare$Freq / sum(n_nodes_6_worse_rare$Freq)) * 100
# subset better
cf_hboth_n6_final_better <- subset(cf_hboth_n6_final, cat == "better")
cf_hboth_n6_final_better$sample_core <- as.factor(as.character(cf_hboth_n6_final_better$sample_core))
cf_hboth_n6_final_better$sample_rare <- as.factor(as.character(cf_hboth_n6_final_better$sample_rare))
n_nodes_6_better_core <- table(cf_hboth_n6_final_better$sample_core)
n_nodes_6_better_core <- data.frame(n_nodes_6_better_core)
n_nodes_6_better_core$which <- "better"
n_nodes_6_better_core$n <- 6
n_nodes_6_better_core$Perc <- (n_nodes_6_better_core$Freq / sum(n_nodes_6_better_core$Freq)) * 100
n_nodes_6_better_rare <- table(cf_hboth_n6_final_better$sample_rare)
n_nodes_6_better_rare <- data.frame(n_nodes_6_better_rare)
n_nodes_6_better_rare$which <- "better"
n_nodes_6_better_rare$n <- 6
n_nodes_6_better_rare$Perc <- (n_nodes_6_better_rare$Freq / sum(n_nodes_6_better_rare$Freq)) * 100

cf_hboth_n8 <- subset(cf_hboth, n_nodes == 8)
frechet_hboth_deg_n8 <- subset(frechet_attack_similarity, N == 8)
frechet_hboth_deg_n8 <- frechet_hboth_deg_n8[rep(seq_len(nrow(frechet_hboth_deg_n8)), each = 4), ]
cf_hboth_n8_final <- data.frame(cbind(cf_hboth_n8, frechet_hboth_deg_n8$cf_hboth_deg))
colnames(cf_hboth_n8_final) <- c("i", "n_nodes", "sample_core", "sample_rare", "frechet_deg")
cf_hboth_n8_final$cat <- ifelse(cf_hboth_n8_final$frechet_deg < original_frechet_h_cf_degree, "better", "worse")

#subset worse
cf_hboth_n8_final_worse <- subset(cf_hboth_n8_final, cat == "worse")
cf_hboth_n8_final_worse$sample_core <- as.factor(as.character(cf_hboth_n8_final_worse$sample_core))
cf_hboth_n8_final_worse$sample_rare <- as.factor(as.character(cf_hboth_n8_final_worse$sample_rare))
n_nodes_8_worse_core <- table(cf_hboth_n8_final_worse$sample_core)
n_nodes_8_worse_core <- data.frame(n_nodes_8_worse_core)
n_nodes_8_worse_core$which <- "worse"
n_nodes_8_worse_core$n <- 8
n_nodes_8_worse_core$Perc <- (n_nodes_8_worse_core$Freq / sum(n_nodes_8_worse_core$Freq)) * 100
n_nodes_8_worse_rare <- table(cf_hboth_n8_final_worse$sample_rare)
n_nodes_8_worse_rare <- data.frame(n_nodes_8_worse_rare)
n_nodes_8_worse_rare$which <- "worse"
n_nodes_8_worse_rare$n <- 8
n_nodes_8_worse_rare$Perc <- (n_nodes_8_worse_rare$Freq / sum(n_nodes_8_worse_rare$Freq)) * 100
# subset better
cf_hboth_n8_final_better <- subset(cf_hboth_n8_final, cat == "better")
cf_hboth_n8_final_better$sample_core <- as.factor(as.character(cf_hboth_n8_final_better$sample_core))
cf_hboth_n8_final_better$sample_rare <- as.factor(as.character(cf_hboth_n8_final_better$sample_rare))
n_nodes_8_better_core <- table(cf_hboth_n8_final_better$sample_core)
n_nodes_8_better_core <- data.frame(n_nodes_8_better_core)
n_nodes_8_better_core$which <- "better"
n_nodes_8_better_core$n <- 8
n_nodes_8_better_core$Perc <- (n_nodes_8_better_core$Freq / sum(n_nodes_8_better_core$Freq)) * 100
n_nodes_8_better_rare <- table(cf_hboth_n8_final_better$sample_rare)
n_nodes_8_better_rare <- data.frame(n_nodes_8_better_rare)
n_nodes_8_better_rare$which <- "better"
n_nodes_8_better_rare$n <- 8
n_nodes_8_better_rare$Perc <- (n_nodes_8_better_rare$Freq / sum(n_nodes_8_better_rare$Freq)) * 100


cf_hboth_n10 <- subset(cf_hboth, n_nodes == 10)
frechet_hboth_deg_n10 <- subset(frechet_attack_similarity, N == 10)
frechet_hboth_deg_n10 <- frechet_hboth_deg_n10[rep(seq_len(nrow(frechet_hboth_deg_n10)), each = 5), ]
cf_hboth_n10_final <- data.frame(cbind(cf_hboth_n10, frechet_hboth_deg_n10$cf_hboth_deg))
colnames(cf_hboth_n10_final) <- c("i", "n_nodes", "sample_core", "sample_rare", "frechet_deg")
cf_hboth_n10_final$cat <- ifelse(cf_hboth_n10_final$frechet_deg < original_frechet_h_cf_degree, "better", "worse")

#subset worse
cf_hboth_n10_final_worse <- subset(cf_hboth_n10_final, cat == "worse")
cf_hboth_n10_final_worse$sample_core <- as.factor(as.character(cf_hboth_n10_final_worse$sample_core))
cf_hboth_n10_final_worse$sample_rare <- as.factor(as.character(cf_hboth_n10_final_worse$sample_rare))
n_nodes_10_worse_core <- table(cf_hboth_n10_final_worse$sample_core)
n_nodes_10_worse_core <- data.frame(n_nodes_10_worse_core)
n_nodes_10_worse_core$which <- "worse"
n_nodes_10_worse_core$n <- 10
n_nodes_10_worse_core$Perc <- (n_nodes_10_worse_core$Freq / sum(n_nodes_10_worse_core$Freq)) * 100
n_nodes_10_worse_rare <- table(cf_hboth_n10_final_worse$sample_rare)
n_nodes_10_worse_rare <- data.frame(n_nodes_10_worse_rare)
n_nodes_10_worse_rare$which <- "worse"
n_nodes_10_worse_rare$n <- 10
n_nodes_10_worse_rare$Perc <- (n_nodes_10_worse_rare$Freq / sum(n_nodes_10_worse_rare$Freq)) * 100
# subset better
cf_hboth_n10_final_better <- subset(cf_hboth_n10_final, cat == "better")
cf_hboth_n10_final_better$sample_core <- as.factor(as.character(cf_hboth_n10_final_better$sample_core))
cf_hboth_n10_final_better$sample_rare <- as.factor(as.character(cf_hboth_n10_final_better$sample_rare))
n_nodes_10_better_core <- table(cf_hboth_n10_final_better$sample_core)
n_nodes_10_better_core <- data.frame(n_nodes_10_better_core)
n_nodes_10_better_core$which <- "better"
n_nodes_10_better_core$n <- 10
n_nodes_10_better_core$Perc <- (n_nodes_10_better_core$Freq / sum(n_nodes_10_better_core$Freq)) * 100
n_nodes_10_better_rare <- table(cf_hboth_n10_final_better$sample_rare)
n_nodes_10_better_rare <- data.frame(n_nodes_10_better_rare)
n_nodes_10_better_rare$which <- "better"
n_nodes_10_better_rare$n <- 10
n_nodes_10_better_rare$Perc <- (n_nodes_10_better_rare$Freq / sum(n_nodes_10_better_rare$Freq)) * 100


cf_hboth_n12 <- subset(cf_hboth, n_nodes == 12)
frechet_hboth_deg_n12 <- subset(frechet_attack_similarity, N == 12)
frechet_hboth_deg_n12 <- frechet_hboth_deg_n12[rep(seq_len(nrow(frechet_hboth_deg_n12)), each = 6), ]
cf_hboth_n12_final <- data.frame(cbind(cf_hboth_n12, frechet_hboth_deg_n12$cf_hboth_deg))
colnames(cf_hboth_n12_final) <- c("i", "n_nodes", "sample_core", "sample_rare", "frechet_deg")
cf_hboth_n12_final$cat <- ifelse(cf_hboth_n12_final$frechet_deg < original_frechet_h_cf_degree, "better", "worse")

#subset worse
cf_hboth_n12_final_worse <- subset(cf_hboth_n12_final, cat == "worse")
cf_hboth_n12_final_worse$sample_core <- as.factor(as.character(cf_hboth_n12_final_worse$sample_core))
cf_hboth_n12_final_worse$sample_rare <- as.factor(as.character(cf_hboth_n12_final_worse$sample_rare))
n_nodes_12_worse_core <- table(cf_hboth_n12_final_worse$sample_core)
n_nodes_12_worse_core <- data.frame(n_nodes_12_worse_core)
n_nodes_12_worse_core$which <- "worse"
n_nodes_12_worse_core$n <- 12
n_nodes_12_worse_core$Perc <- (n_nodes_12_worse_core$Freq / sum(n_nodes_12_worse_core$Freq)) * 100
n_nodes_12_worse_rare <- table(cf_hboth_n12_final_worse$sample_rare)
n_nodes_12_worse_rare <- data.frame(n_nodes_12_worse_rare)
n_nodes_12_worse_rare$which <- "worse"
n_nodes_12_worse_rare$n <- 12
n_nodes_12_worse_rare$Perc <- (n_nodes_12_worse_rare$Freq / sum(n_nodes_12_worse_rare$Freq)) * 100
# subset better
cf_hboth_n12_final_better <- subset(cf_hboth_n12_final, cat == "better")
cf_hboth_n12_final_better$sample_core <- as.factor(as.character(cf_hboth_n12_final_better$sample_core))
cf_hboth_n12_final_better$sample_rare <- as.factor(as.character(cf_hboth_n12_final_better$sample_rare))
n_nodes_12_better_core <- table(cf_hboth_n12_final_better$sample_core)
n_nodes_12_better_core <- data.frame(n_nodes_12_better_core)
n_nodes_12_better_core$which <- "better"
n_nodes_12_better_core$n <- 12
n_nodes_12_better_core$Perc <- (n_nodes_12_better_core$Freq / sum(n_nodes_12_better_core$Freq)) * 100
n_nodes_12_better_rare <- table(cf_hboth_n12_final_better$sample_rare)
n_nodes_12_better_rare <- data.frame(n_nodes_12_better_rare)
n_nodes_12_better_rare$which <- "better"
n_nodes_12_better_rare$n <- 12
n_nodes_12_better_rare$Perc <- (n_nodes_12_better_rare$Freq / sum(n_nodes_12_better_rare$Freq)) * 100


cf_hboth_n14 <- subset(cf_hboth, n_nodes == 14)
frechet_hboth_deg_n14 <- subset(frechet_attack_similarity, N == 14)
frechet_hboth_deg_n14 <- frechet_hboth_deg_n14[rep(seq_len(nrow(frechet_hboth_deg_n14)), each = 7), ]
cf_hboth_n14_final <- data.frame(cbind(cf_hboth_n14, frechet_hboth_deg_n14$cf_hboth_deg))
colnames(cf_hboth_n14_final) <- c("i", "n_nodes", "sample_core", "sample_rare", "frechet_deg")
cf_hboth_n14_final$cat <- ifelse(cf_hboth_n14_final$frechet_deg < original_frechet_h_cf_degree, "better", "worse")

#subset worse
cf_hboth_n14_final_worse <- subset(cf_hboth_n14_final, cat == "worse")
cf_hboth_n14_final_worse$sample_core <- as.factor(as.character(cf_hboth_n14_final_worse$sample_core))
cf_hboth_n14_final_worse$sample_rare <- as.factor(as.character(cf_hboth_n14_final_worse$sample_rare))
n_nodes_14_worse_core <- table(cf_hboth_n14_final_worse$sample_core)
n_nodes_14_worse_core <- data.frame(n_nodes_14_worse_core)
n_nodes_14_worse_core$which <- "worse"
n_nodes_14_worse_core$n <- 14
n_nodes_14_worse_core$Perc <- (n_nodes_14_worse_core$Freq / sum(n_nodes_14_worse_core$Freq)) * 100
n_nodes_14_worse_rare <- table(cf_hboth_n14_final_worse$sample_rare)
n_nodes_14_worse_rare <- data.frame(n_nodes_14_worse_rare)
n_nodes_14_worse_rare$which <- "worse"
n_nodes_14_worse_rare$n <- 14
n_nodes_14_worse_rare$Perc <- (n_nodes_14_worse_rare$Freq / sum(n_nodes_14_worse_rare$Freq)) * 100
# subset better
cf_hboth_n14_final_better <- subset(cf_hboth_n14_final, cat == "better")
cf_hboth_n14_final_better$sample_core <- as.factor(as.character(cf_hboth_n14_final_better$sample_core))
cf_hboth_n14_final_better$sample_rare <- as.factor(as.character(cf_hboth_n14_final_better$sample_rare))
n_nodes_14_better_core <- table(cf_hboth_n14_final_better$sample_core)
n_nodes_14_better_core <- data.frame(n_nodes_14_better_core)
n_nodes_14_better_core$which <- "better"
n_nodes_14_better_core$n <- 14
n_nodes_14_better_core$Perc <- (n_nodes_14_better_core$Freq / sum(n_nodes_14_better_core$Freq)) * 100
n_nodes_14_better_rare <- table(cf_hboth_n14_final_better$sample_rare)
n_nodes_14_better_rare <- data.frame(n_nodes_14_better_rare)
n_nodes_14_better_rare$which <- "better"
n_nodes_14_better_rare$n <- 14
n_nodes_14_better_rare$Perc <- (n_nodes_14_better_rare$Freq / sum(n_nodes_14_better_rare$Freq)) * 100


cf_hboth_n16 <- subset(cf_hboth, n_nodes == 16)
frechet_hboth_deg_n16 <- subset(frechet_attack_similarity, N == 16)
frechet_hboth_deg_n16 <- frechet_hboth_deg_n16[rep(seq_len(nrow(frechet_hboth_deg_n16)), each = 8), ]
cf_hboth_n16_final <- data.frame(cbind(cf_hboth_n16, frechet_hboth_deg_n16$cf_hboth_deg))
colnames(cf_hboth_n16_final) <- c("i", "n_nodes", "sample_core", "sample_rare", "frechet_deg")
cf_hboth_n16_final$cat <- ifelse(cf_hboth_n16_final$frechet_deg < original_frechet_h_cf_degree, "better", "worse")

#subset worse
cf_hboth_n16_final_worse <- subset(cf_hboth_n16_final, cat == "worse")
cf_hboth_n16_final_worse$sample_core <- as.factor(as.character(cf_hboth_n16_final_worse$sample_core))
cf_hboth_n16_final_worse$sample_rare <- as.factor(as.character(cf_hboth_n16_final_worse$sample_rare))
n_nodes_16_worse_core <- table(cf_hboth_n16_final_worse$sample_core)
n_nodes_16_worse_core <- data.frame(n_nodes_16_worse_core)
n_nodes_16_worse_core$which <- "worse"
n_nodes_16_worse_core$n <- 16
n_nodes_16_worse_core$Perc <- (n_nodes_16_worse_core$Freq / sum(n_nodes_16_worse_core$Freq)) * 100
n_nodes_16_worse_rare <- table(cf_hboth_n16_final_worse$sample_rare)
n_nodes_16_worse_rare <- data.frame(n_nodes_16_worse_rare)
n_nodes_16_worse_rare$which <- "worse"
n_nodes_16_worse_rare$n <- 16
n_nodes_16_worse_rare$Perc <- (n_nodes_16_worse_rare$Freq / sum(n_nodes_16_worse_rare$Freq)) * 100
# subset better
cf_hboth_n16_final_better <- subset(cf_hboth_n16_final, cat == "better")
cf_hboth_n16_final_better$sample_core <- as.factor(as.character(cf_hboth_n16_final_better$sample_core))
cf_hboth_n16_final_better$sample_rare <- as.factor(as.character(cf_hboth_n16_final_better$sample_rare))
n_nodes_16_better_core <- table(cf_hboth_n16_final_better$sample_core)
n_nodes_16_better_core <- data.frame(n_nodes_16_better_core)
n_nodes_16_better_core$which <- "better"
n_nodes_16_better_core$n <- 16
n_nodes_16_better_core$Perc <- (n_nodes_16_better_core$Freq / sum(n_nodes_16_better_core$Freq)) * 100
n_nodes_16_better_rare <- table(cf_hboth_n16_final_better$sample_rare)
n_nodes_16_better_rare <- data.frame(n_nodes_16_better_rare)
n_nodes_16_better_rare$which <- "better"
n_nodes_16_better_rare$n <- 16
n_nodes_16_better_rare$Perc <- (n_nodes_16_better_rare$Freq / sum(n_nodes_16_better_rare$Freq)) * 100


merge_hboth_core <- rbind(n_nodes_4_better_core, n_nodes_6_better_core, n_nodes_8_better_core, 
                          n_nodes_10_better_core, n_nodes_12_better_core, n_nodes_14_better_core, n_nodes_16_better_core,
                          n_nodes_4_worse_core, n_nodes_6_worse_core, n_nodes_8_worse_core, 
                          n_nodes_10_worse_core, n_nodes_12_worse_core, n_nodes_14_worse_core, n_nodes_16_worse_core)

merge_hboth_rare <- rbind(n_nodes_4_better_rare, n_nodes_6_better_rare, n_nodes_8_better_rare, 
                          n_nodes_10_better_rare, n_nodes_12_better_rare, n_nodes_14_better_rare, n_nodes_16_better_rare,
                          n_nodes_4_worse_rare, n_nodes_6_worse_rare, n_nodes_8_worse_rare, 
                          n_nodes_10_worse_rare, n_nodes_12_worse_rare, n_nodes_14_worse_rare, n_nodes_16_worse_rare)

merge_hboth_core$Var1 <- str_replace(merge_hboth_core$Var1, "_", " ")
merge_hboth_rare$Var1 <- str_replace(merge_hboth_rare$Var1, "_", " ")
merge_hboth <- data.frame(rbind(merge_hboth_core, merge_hboth_rare))

hboth_frechet_plot <- ggplot(merge_hboth) + geom_tile(aes(x=n, y=Var1, fill=Perc)) + facet_wrap(~which) +
  scale_fill_gradient(low="cyan", high="red") + ylab("") + xlab("Number of nodes added") + theme_pubr(border=TRUE, legend = "right") +
  theme(legend.title = element_blank(), axis.text.y = element_text(face = "italic"))

merge_hboth_wide <- spread(merge_hboth, key="Var1", value="Freq")
merge_hboth_wide$Perc <- NULL

merge_hboth_wide_better <- subset(merge_hboth_wide, which == "better")
merge_hboth_wide_better[is.na(merge_hboth_wide_better)] <- 0
merge_hboth_wide_better$which <- NULL
merge_hboth_wide_better <- ddply(merge_hboth_wide_better,"n",numcolwise(sum))
rownames(merge_hboth_wide_better) <- merge_hboth_wide_better$n
merge_hboth_wide_better$n <- NULL
merge_hboth_wide_better <- data.frame(merge_hboth_wide_better)
hboth_better_fisher <- vegan::fisher.alpha(merge_hboth_wide_better)
hboth_better_shannon <- vegan::diversity(merge_hboth_wide_better, index = "shannon")
hboth_better_specNum <- vegan::specnumber(merge_hboth_wide_better)
hboth_better_dominance <- vegan::diversity(merge_hboth_wide_better, index = "simpson")
hboth_better_BPindex <- microbiome::dominance(t(merge_hboth_wide_better))
div_hboth_better <- data.frame(cbind(hboth_better_fisher, hboth_better_shannon, 
                                     hboth_better_specNum, hboth_better_dominance, hboth_better_BPindex$gini))
colnames(div_hboth_better) <- c("fisher", "shannon", "specNumber", "dominance", "gini")
div_hboth_better$performance <- "better"

merge_hboth_wide_worse <- subset(merge_hboth_wide, which == "worse")
merge_hboth_wide_worse[is.na(merge_hboth_wide_worse)] <- 0
merge_hboth_wide_worse$which <- NULL
merge_hboth_wide_worse <- ddply(merge_hboth_wide_worse,"n",numcolwise(sum))
rownames(merge_hboth_wide_worse) <- merge_hboth_wide_worse$n
merge_hboth_wide_worse$n <- NULL
merge_hboth_wide_worse <- data.frame(merge_hboth_wide_worse)
hboth_worse_fisher <- vegan::fisher.alpha(merge_hboth_wide_worse)
hboth_worse_shannon <- vegan::diversity(merge_hboth_wide_worse, index = "shannon")
hboth_worse_specNum <- vegan::specnumber(merge_hboth_wide_worse)
hboth_worse_dominance <- vegan::diversity(merge_hboth_wide_worse, index = "simpson")
hboth_worse_BPindex <- microbiome::dominance(t(merge_hboth_wide_worse))
div_hboth_worse <- data.frame(cbind(hboth_worse_fisher, hboth_worse_shannon, 
                                    hboth_worse_specNum, hboth_worse_dominance, hboth_worse_BPindex$gini))
colnames(div_hboth_worse) <- c("fisher", "shannon", "specNumber", "dominance", "gini")
div_hboth_worse$performance <- "worse"

div_hboth <- data.frame(rbind(div_hboth_better, div_hboth_worse))
div_hboth$performance <- as.factor(as.character(div_hboth$performance))

# index of diversity (only a few factors dominante)
wilcox.test(shannon ~ performance, data=div_hboth)
wilcoxonR(div_hboth$shannon, g=div_hboth$performance, ci=TRUE)
# community dominance index 
wilcox.test(gini ~ performance, data=div_hboth)
wilcoxonR(div_hboth$gini, g=div_hboth$performance, ci=TRUE)

hboth_gini <- ggplot(div_hboth, aes(x=performance, y=gini)) + geom_boxplot() + 
  geom_point(size=3, alpha=0.3) + ylab("Community dominance index, gini \n") + 
  stat_compare_means( aes(label = ..p.signif..),label.x.npc = "middle", label.y = 0.75, size=6) +
  ylim(0.0, 0.8) + xlab("Performance") + theme_pubr(border=TRUE, base_size = 10)

hboth_shannon <- ggplot(div_hboth, aes(x=performance, y=shannon)) + geom_boxplot() +
  geom_point(size=3, alpha=0.3) + ylab("Shannon diversity index \n") +
  stat_compare_means(aes(label = ..p.signif..), label.x.npc = "middle", label.y = 3.9, size = 6) +
  ylim(2.5, 4.0) + xlab("Performance") + theme_pubr(border=TRUE, base_size = 10)

hboth_shannon_gini_plot <- ggarrange(hboth_shannon, hboth_gini, labels = c("a", "b"))

summary(frechet_attack_similarity$cf_hrare_rand)
summary(frechet_attack_similarity$CF_hcore_rand)
summary(frechet_attack_similarity$cf_hboth_rand)
summary(network_kernel_shortestPath$CF...Healthy..rare.)
summary(network_kernel_shortestPath$CF...Healthy..core.)
summary(network_kernel_shortestPath$CF...Healthy..both.)
summary(network_kernel_WL$CF...Healthy..rare.)
summary(network_kernel_WL$CF...Healthy..core.)
summary(network_kernel_WL$CF...Healthy..both.)

cf_hboth_persPreb <- data.frame(cbind(cf_hboth, frechet_attack_similarity$cf_hboth_deg))
colnames(cf_hboth_persPreb) <- c("i", "n_nodes", "sample_core", "sample_rare", "frechet_deg")
cf_hboth_persPreb$cat <- ifelse(cf_hboth_persPreb$frechet_deg < original_frechet_h_cf_degree, "better", "worse")
cf_hboth_persPreb_better <- subset(cf_hboth_persPreb, cat == "better")
cf_hboth_persPreb_better$cat_coreCF <- ifelse(cf_hboth_persPreb_better$sample_core %in% rownames(background_core_cf), 1, 0)

cf_hboth_persPreb_worse <- subset(cf_hboth_persPreb, cat == "worse")
cf_hboth_persPreb_worse$cat_coreCF <- ifelse(cf_hboth_persPreb_worse$sample_core %in% rownames(background_core_cf), 1, 0)
preb_better <- table(cf_hboth_persPreb_better$cat_coreCF)
in_better <- preb_better[2]
notIn_better <- preb_better[1]
in_better_perc <- in_better / (in_better + notIn_better)
in_better_perc <- round(in_better_perc,2)
preb_worse <- table(cf_hboth_persPreb_worse$cat_coreCF)
in_worse <- preb_worse[2]
notIn_worse <- preb_worse[1]
in_worse_perc <- in_worse / (in_worse + notIn_worse)
in_worse_perc <- round(in_worse_perc,2)


# Make figures ####
ggsave("image_output/Fig_4.tif", first_network_plot,
       device="tiff", width=11, height=6.9, dpi=600)

ggsave("image_output/Fig_5.tif", modulated_network_similarities, 
       device="tiff", dpi = 600, width = 8.59, height=7.98)

ggsave("image_output/Fig_6.tif", hboth_frechet_plot, dpi=600,
       device="tiff", width=11.6, height=7.81)

# Supplementary figures

ggsave("image_output/Supplementary_Fig_8.tif", stats_rare_core, dpi=300,
       device="tiff", width=10.4, height=6.9)

ggsave("image_output/Supplementary_Fig_9.tif", hboth_shannon_gini_plot, dpi=300,
       device="tiff", width=6, height=3.5)

# Make tables
#write.csv2(frechet_attack_similarity, "Frechet_attack_similarity_modulatedCF.csv", row.names = FALSE)
#write.csv2(network_kernel_WL, "kernel_WL_similarity_modulatedCF.csv", row.names = FALSE)
#write.csv2(network_kernel_shortestPath, "kernel_sp_similarity_modulatedCF.csv", row.names = FALSE)