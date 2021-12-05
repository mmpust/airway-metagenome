# title: "Step 2, Simulation runs"
# author: "Marie-Madlen Pust"
# last update: "05 December 2021"

############################################################################################################
# clean global environment
rm(list=ls())

# set working directory
setwd('C:/R')

############################################################################################################
# define global functions
# load packages and install if necessary
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, 'Package'])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)}

# packages
packages <- c('readr','viridis','dplyr','stringr','tidyr','ggplot2','factoextra','gmodels','ggpubr','plyr','purrr','Hmisc',
              'reshape','vegan','magrittr','scales','grid','reshape2','rcompanion','SimilarityMeasures','viridis','knitr', 
              'ggrepel','forcats','pheatmap', 'cowplot', 'popgraph', 'qgraph', 'sna', 'tidygraph','expm','igraph','NetSwan',
              'CINNA', 'graphkernels', 'data.table', 'microbiome', 'RVAideMemoire', 'bipartite', 'corrplot', 'ROCR')

# generate network statistics
all_indices <- function(g){
  res <- matrix(0,vcount(g),4)
  res[,1] <- igraph::degree(g)
  res[,2] <- igraph::betweenness(g)
  res[,3] <- igraph::closeness(g)
  res[,4] <- igraph::hub_score(g)$vector
  apply(res,2,function(x) round(x,8))}

# correlation parameters
weight_val_pos = 0.2 
weight_val_neg = -0.2 
sig_level = 0.01

# network parameters
directed_network=FALSE 

# network evaluation, null models
net.metric.zscore = function(obsval, nullval) {(obsval - mean(nullval))/sd(nullval)} 
net.metric.pvalue = function(x){2*pnorm(-abs(x))}

nullmodel_method = "shuffle.web" 
nullmodel_networks_n = 100 

############################################################################################################
# load packages
ipak(packages)

############################################################################################################
# import meta data table of patients
md <- read_delim('input_files/meta_data/spatial_metadata_2020_12_2.csv', ';', escape_double = FALSE, trim_ws = TRUE)
# convert to data frame object
md <- data.frame(md)
# make sample IDs row names
rownames(md) <- md$Sample 
# remove sample IDs as column
md$Sample <- NULL
# add column with BMI information
md$BMI <- round(md$Weight / ((md$Height/100) * (md$Height/100)),1)


# import pangenome data
# bcphc-normalised count table
ds_pangenome_bcphc <- read_delim('input_files/taxonomic_data/pangenome.bcphc.merged.csv', ';', escape_double = FALSE, trim_ws = TRUE, skip_empty_rows = TRUE)
# convert to data frame object
ds_pangenome_bcphc <- data.frame(ds_pangenome_bcphc)
# remove duplicate rows by getting the count sum per species
ds_pangenome_bcphc <- ddply(ds_pangenome_bcphc, 'Species', numcolwise(sum))
# set species as row names
rownames(ds_pangenome_bcphc) <- ds_pangenome_bcphc$Species
# remove non-numeric species column
ds_pangenome_bcphc$Species <- NULL

# import pangenome data
# RLE-normalised count table
ds_pangenome_rle <- read_delim('input_files/taxonomic_data/pangenome.rle.merged.csv', ';', escape_double = FALSE, trim_ws = TRUE, skip_empty_rows = TRUE)
# convert to data frame object
ds_pangenome_rle <- data.frame(ds_pangenome_rle)
# remove duplicate rows by getting the count sum per species
ds_pangenome_rle <- ddply(ds_pangenome_rle, 'Species', numcolwise(sum))
# set species as row names
rownames(ds_pangenome_rle) <- ds_pangenome_rle$Species
# remove non-numeric species column
ds_pangenome_rle$Species <- NULL

# import pangenome data
# vst-normalised count table
ds_pangenome_vst <- read_delim('input_files/taxonomic_data/pangenome.vst.merged.csv', ';', escape_double = FALSE, trim_ws = TRUE, skip_empty_rows = TRUE)
# convert to data frame object
ds_pangenome_vst <- data.frame(ds_pangenome_vst)
# remove duplicate rows by getting the count sum per species
ds_pangenome_vst <- ddply(ds_pangenome_vst, 'Species', numcolwise(sum))
# set species as row names
rownames(ds_pangenome_vst) <- ds_pangenome_vst$Species
# remove non-numeric species column
ds_pangenome_vst$Species <- NULL


# import one-strain per species data
# bcphc-normalised count table
ds_osps_bcphc <- read_delim('input_files/taxonomic_data/osps.bcphc.merged.csv', ';', escape_double = FALSE, trim_ws = TRUE, skip_empty_rows = TRUE)
# convert to data frame object
ds_osps_bcphc <- data.frame(ds_osps_bcphc)
# remove duplicate rows by getting the count sum per species
ds_osps_bcphc <- ddply(ds_osps_bcphc, 'Species', numcolwise(sum))
# set species as row names
rownames(ds_osps_bcphc) <- ds_osps_bcphc$Species
# remove non-numeric species column
ds_osps_bcphc$Species <- NULL

# import one-strain per species data
# RLE normalisation
ds_osps_rle <- read_delim('input_files/taxonomic_data/osps.rle.merged.csv', ';', escape_double = FALSE, trim_ws = TRUE, skip_empty_rows = TRUE)
# convert to data frame object
ds_osps_rle <- data.frame(ds_osps_rle)
# remove duplicate rows by getting the count sum per species
ds_osps_rle <- ddply(ds_osps_rle, 'Species', numcolwise(sum))
# set species as row names
rownames(ds_osps_rle) <- ds_osps_rle$Species
# remove non-numeric species column
ds_osps_rle$Species <- NULL

# import one-strain per species data
# vst normalisation
ds_osps_vst <- read_delim('input_files/taxonomic_data/osps.vst.merged.csv', ';', escape_double = FALSE, trim_ws = TRUE, skip_empty_rows = TRUE)
# convert to data frame object
ds_osps_vst <- data.frame(ds_osps_vst)
# remove duplicate rows by getting the count sum per species
ds_osps_vst <- ddply(ds_osps_vst, 'Species', numcolwise(sum))
# set species as row names
rownames(ds_osps_vst) <- ds_osps_vst$Species
# remove non-numeric species column
ds_osps_vst$Species <- NULL

############################################################################################################
# load background species table, BCPHC-normalised
background_species <- read_delim("background_species.csv", ";", escape_double = FALSE, trim_ws = TRUE)

# convert to dataframe object
background_species <- data.frame(background_species)
# subset background species table, extract healthy, core species, BCPHC-normalised, 25th abundance percentile
background_core_h_bcphc <- subset(background_species, State == "Healthy" & Species_type == "background_core" & Normalisation == "BCPHC" & Threshold == "25th percentile")
# remove duplicate entries
background_core_h_bcphc <- background_core_h_bcphc[!duplicated(background_core_h_bcphc$Species),]

# subset background species table, extract healthy, rare species, BCPHC-normalised, 25th abundance percentile
background_rare_h_bcphc <- subset(background_species, State == "Healthy" & Species_type == "background_rare" & Normalisation == "BCPHC" & Threshold == "25th percentile")
# remove duplicate entries
background_rare_h_bcphc <- background_rare_h_bcphc[!duplicated(background_rare_h_bcphc$Species),]

# subset background species table, extract CF, core species, BCPHC-normalised, 25th abundance percentile
background_core_cf_bcphc <- subset(background_species, State == "CF" & Species_type == "background_core" & Normalisation == "BCPHC" & Threshold == "25th percentile")
# remove duplicate entries
background_core_cf_bcphc <- background_core_cf_bcphc[!duplicated(background_core_cf_bcphc$Species),]

# subset background species table, extract CF, rare species, BCPHC-normalised, 25th abundance percentile
background_rare_cf_bcphc <- subset(background_species, State == "CF" & Species_type == "background_rare" & Normalisation == "BCPHC" & Threshold == "25th percentile")
# remove duplicate entries
background_rare_cf_bcphc <- background_rare_cf_bcphc[!duplicated(background_rare_cf_bcphc$Species),]


# load background species table, RLE-normalised
# subset background species table, extract healthy, core species, RLE-normalised, 25th abundance percentile
background_core_h_rle <- subset(background_species, State == "Healthy" & Species_type == "background_core" & Normalisation == "RLE" & Threshold == "25th percentile")
# remove duplicate entries
background_core_h_rle <- background_core_h_rle[!duplicated(background_core_h_rle$Species),]

# subset background species table, extract healthy, rare species, RLE-normalised, 25th abundance percentile
background_rare_h_rle <- subset(background_species, State == "Healthy" & Species_type == "background_rare" & Normalisation == "RLE" & Threshold == "25th percentile")
# remove duplicate entries
background_rare_h_rle <- background_rare_h_rle[!duplicated(background_rare_h_rle$Species),]

# subset background species table, extract CF, core species, RLE-normalised, 25th abundance percentile
background_core_cf_rle <- subset(background_species, State == "CF" & Species_type == "background_core" & Normalisation == "RLE" & Threshold == "25th percentile")
# remove duplicate entries
background_core_cf_rle <- background_core_cf_rle[!duplicated(background_core_cf_rle$Species),]

# subset background species table, extract CF, rare species, RLE-normalised, 25th abundance percentile
background_rare_cf_rle <- subset(background_species, State == "CF" & Species_type == "background_rare" & Normalisation == "RLE" & Threshold == "25th percentile")
# remove duplicate entries
background_rare_cf_rle <- background_rare_cf_rle[!duplicated(background_rare_cf_rle$Species),]


#load background species table, VST-normalised
# subset background species table, extract healthy, core species, VST-normalised, 25th abundance percentile
background_core_h_vst <- subset(background_species, State == "Healthy" & Species_type == "background_core" & Normalisation == "VST" & Threshold == "25th percentile")
# remove duplicate entries
background_core_h_vst <- background_core_h_vst[!duplicated(background_core_h_vst$Species),]

# subset background species table, extract healthy, rare species, VST-normalised, 25th abundance percentile
background_rare_h_vst <- subset(background_species, State == "Healthy" & Species_type == "background_rare" & Normalisation == "VST" & Threshold == "25th percentile")
# remove duplicate entries
background_rare_h_vst <- background_rare_h_vst[!duplicated(background_rare_h_vst$Species),]

# subset background species table, extract CF, core species, VST-normalised, 25th abundance percentile
background_core_cf_vst <- subset(background_species, State == "CF" & Species_type == "background_core" & Normalisation == "VST" & Threshold == "25th percentile")
# remove duplicate entries
background_core_cf_vst <- background_core_cf_vst[!duplicated(background_core_cf_vst$Species),]

# subset background species table, extract CF, rare species, VST-normalised, 25th abundance percentile
background_rare_cf_vst <- subset(background_species, State == "CF" & Species_type == "background_rare" & Normalisation == "VST" & Threshold == "25th percentile")
# remove duplicate entries
background_rare_cf_vst <- background_rare_cf_vst[!duplicated(background_rare_cf_vst$Species),]

############################################################################################################
# Correlation between rare and core species
# generate correlation matrix
# healthy, extract background rare species (as defined by pangenome and one-strain per species)
df_h_rare_bcphc <- subset(ds_pangenome_bcphc, rownames(ds_pangenome_bcphc) %in% background_core_h_bcphc$Species | rownames(ds_pangenome_bcphc) %in% background_rare_h_bcphc$Species)
# transpose dataframe
df_h_rare_bcphc <- data.frame(t(df_h_rare_bcphc))
# add column with disease state (CF vs healthy)
df_h_rare_bcphc$state <- md$State
# subset and remove CF data
df_h_rare_bcphc <- subset(df_h_rare_bcphc, state == "Healthy")
# remove non-numeric column
df_h_rare_bcphc$state <- NULL
# obtain Spearman's rank correlation of species count matrix
df_h_rare_bcphc <- rcorr(as.matrix(df_h_rare_bcphc), type = 'spearman')

# create node and edge lists, healthy
# extract and store p-values of correlation analysis
df_h_rare_bcphc_Pval <- df_h_rare_bcphc$P
# extract and store the correlation coefficient of correlation analysis
df_h_rare_bcphc_COR <- df_h_rare_bcphc$r

# generate edge list
# take data in wide format and stack columns into a single column of data for p-values
df_h_rare_bcphc_COR_edges <- reshape2::melt(df_h_rare_bcphc_Pval)
# take data in wide format and stack columns into a single column of data for correlation coefficient
df_h_rare_bcphc_coredges <- reshape2::melt(df_h_rare_bcphc_COR)
# merge data frames
df_h_rare_bcphc_COR_edges$COR <- df_h_rare_bcphc_coredges$value
# round p-value to 5 decimal places
df_h_rare_bcphc_COR_edges$value <- round(df_h_rare_bcphc_COR_edges$value, 5)
# store row names in own column
df_h_rare_bcphc_COR_edges$Label <- row.names(df_h_rare_bcphc_coredges)
# re-name columns
colnames(df_h_rare_bcphc_COR_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index row names
rownames(df_h_rare_bcphc_COR_edges) <- NULL
# remove NAs
df_h_rare_bcphc_COR_edges <- df_h_rare_bcphc_COR_edges[complete.cases(df_h_rare_bcphc_COR_edges), ]
# extract significant correlations as defined by global variables (see above)
df_h_rare_bcphc_COR_edges_short <- subset(df_h_rare_bcphc_COR_edges, pValue < sig_level & Weight >= weight_val_pos | Weight <= weight_val_neg)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
df_h_rare_bcphc_COR_edges_short$cortype <- ifelse(df_h_rare_bcphc_COR_edges_short$Weight > 0, "pos", "neg")
# clean column entries
df_h_rare_bcphc_COR_edges_short$Target <- str_replace(df_h_rare_bcphc_COR_edges_short$Target, "\\.", " ")
df_h_rare_bcphc_COR_edges_short$Source <- str_replace(df_h_rare_bcphc_COR_edges_short$Source, "\\.", " ")
# extract two columns from edge list to build node list
df_h_rare_bcphc_nodes <- select(df_h_rare_bcphc_COR_edges_short, c("Source", "cortype"))
# remove duplicates
df_h_rare_bcphc_nodes = df_h_rare_bcphc_nodes[!duplicated(df_h_rare_bcphc_nodes$Source),]
# convert to data frame object
df_h_rare_bcphc_nodes <- data.frame(df_h_rare_bcphc_nodes)
# re-index the data frame
rownames(df_h_rare_bcphc_nodes) <- NULL
# define column names of node list
colnames(df_h_rare_bcphc_nodes) <- c("Id", "CorType")
# clean entries so that they match edge list entries
df_h_rare_bcphc_nodes$Id <- str_replace(df_h_rare_bcphc_nodes$Id, "\\.", " ")
# if the species in node list are found in vector storing background rare species, add "rare", otherwise "core"
df_h_rare_bcphc_nodes$species_type <- ifelse(df_h_rare_bcphc_nodes$Id %in% background_rare_h_bcphc$Species, "rare", "core")
# add genus information
df_h_rare_bcphc_nodes$Genus <- df_h_rare_bcphc_nodes$Id
df_h_rare_bcphc_nodes$Genus <- sapply(strsplit(as.character(df_h_rare_bcphc_nodes$Genus)," "), `[`, 1)

# CF, background rare
# CF, extract background rare species (as defined by pangenome and one-strain per species)
df_cf_rare_bcphc <- subset(ds_pangenome_bcphc, rownames(ds_pangenome_bcphc) %in% background_core_cf_bcphc$Species | rownames(ds_pangenome_bcphc) %in% background_rare_cf_bcphc$Species)
# transpose dataframe
df_cf_rare_bcphc <- data.frame(t(df_cf_rare_bcphc))
# add column with disease state (CF vs healthy)
df_cf_rare_bcphc$state <- md$State
# subset and extract CF samples
df_cf_rare_bcphc <- subset(df_cf_rare_bcphc, state == "CF")
# remove non-numeric column
df_cf_rare_bcphc$state <- NULL
# perform Spearman's rank correlation analysis
df_cf_rare_bcphc <- rcorr(as.matrix(df_cf_rare_bcphc), type = 'spearman')

# create node and edge lists
# extract and store p-values of correlation analysis
df_cf_rare_bcphc_Pval <- df_cf_rare_bcphc$P
# extract and store correlation coefficient of analysis
df_cf_rare_bcphc_COR <- df_cf_rare_bcphc$r

# generate edge list
# take data in wide format and stack columns into a single column of data for p-values
df_cf_rare_bcphc_COR_edges <- reshape2::melt(df_cf_rare_bcphc_Pval)
# take data in wide format and stack columns into a single column of data for correlation coefficient
df_cf_rare_bcphc_coredges <- reshape2::melt(df_cf_rare_bcphc_COR)
# merge tables
df_cf_rare_bcphc_COR_edges$COR <- df_cf_rare_bcphc_coredges$value
# round p-value to five decimal places
df_cf_rare_bcphc_COR_edges$value <- round(df_cf_rare_bcphc_COR_edges$value, 5)
# store row names in individual column
df_cf_rare_bcphc_COR_edges$Label <- row.names(df_cf_rare_bcphc_coredges)
# re-name columns
colnames(df_cf_rare_bcphc_COR_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# remove rownames
rownames(df_cf_rare_bcphc_COR_edges) <- NULL
# remove NAs
df_cf_rare_bcphc_COR_edges <- df_cf_rare_bcphc_COR_edges[complete.cases(df_cf_rare_bcphc_COR_edges), ]
# extract only significant correlations as defined by global parameters (see above)
df_cf_rare_bcphc_COR_edges_short <- subset(df_cf_rare_bcphc_COR_edges, pValue < sig_level & Weight >= weight_val_pos | Weight <= weight_val_neg)
# if correlation is greater than 0, add "positive correlation", otherwise "negative correlation"
df_cf_rare_bcphc_COR_edges_short$cortype <- ifelse(df_cf_rare_bcphc_COR_edges_short$Weight > 0, "pos", "neg")
# clean column names
df_cf_rare_bcphc_COR_edges_short$Target <- str_replace(df_cf_rare_bcphc_COR_edges_short$Target, "\\.", " ")
df_cf_rare_bcphc_COR_edges_short$Source <- str_replace(df_cf_rare_bcphc_COR_edges_short$Source, "\\.", " ")
# take two columns to make node list
df_cf_rare_bcphc_nodes <- select(df_cf_rare_bcphc_COR_edges_short, c("Source", "cortype"))
# remove duplicate rows
df_cf_rare_bcphc_nodes = df_cf_rare_bcphc_nodes[!duplicated(df_cf_rare_bcphc_nodes$Source),]
# create dataframe object
df_cf_rare_bcphc_nodes <- data.frame(df_cf_rare_bcphc_nodes)
# re-index rownames
rownames(df_cf_rare_bcphc_nodes) <- NULL
# rename column names
colnames(df_cf_rare_bcphc_nodes) <- c("Id", "CorType")
# clean up ID name so that it matches ID names of edge list
df_cf_rare_bcphc_nodes$Id <- str_replace(df_cf_rare_bcphc_nodes$Id, "\\.", " ")
# if species in background rare list, add "rare", otherwise add "core"
df_cf_rare_bcphc_nodes$species_type <- ifelse(df_cf_rare_bcphc_nodes$Id %in% background_rare_cf_bcphc$Species,  "rare", "core")
# add genus information by splitting Species name 
df_cf_rare_bcphc_nodes$Genus <- df_cf_rare_bcphc_nodes$Id
df_cf_rare_bcphc_nodes$Genus <- sapply(strsplit(as.character(df_cf_rare_bcphc_nodes$Genus)," "), `[`, 1)


############################################################################################################
# Generate networks ####
# healthy, background core and rare species
# convert species type (core or rare) from class character to class factor
df_h_rare_bcphc_nodes$species_type <- as.factor(as.character(df_h_rare_bcphc_nodes$species_type))
# make graph object from node and edge list
h_rare_net_bcphc <- graph_from_data_frame(d=df_h_rare_bcphc_COR_edges_short, vertices=df_h_rare_bcphc_nodes, directed = directed_network)
# simplify the graph by removing multiple edges and loops
h_rare_net_bcphc <- simplify(h_rare_net_bcphc, remove.multiple = T, remove.loops = T)
# add colour code (rare species = blue, core species = green)
V(h_rare_net_bcphc)$color <- ifelse(V(h_rare_net_bcphc)$species_type == "rare", "lightsteelblue2", "springgreen4")
# choose graph layout (Fruchterman Reingold)
h.net_bcphc.plot <- layout.fruchterman.reingold(h_rare_net_bcphc)
# convert to data frame
h.net_bcphc.plot.df <- as.data.frame(h.net_bcphc.plot)
# add ID information
h.net_bcphc.plot.df$Id <- df_h_rare_bcphc_nodes$Id
# add species information
h.net_bcphc.plot.df$species_type <- df_h_rare_bcphc_nodes$species_type
# convert net to dataframe (important for plotting later on)
h.net_bcphc.plot.edges <- get.data.frame(h_rare_net_bcphc)
# set minimum x values
h.net_bcphc.plot.edges$from.x <- h.net_bcphc.plot.df$V1[match(h.net_bcphc.plot.edges$from, h.net_bcphc.plot.df$Id)]  
# set minimum y values
h.net_bcphc.plot.edges$from.y <- h.net_bcphc.plot.df$V2[match(h.net_bcphc.plot.edges$from, h.net_bcphc.plot.df$Id)]
# set maximum x values
h.net_bcphc.plot.edges$to.x <- h.net_bcphc.plot.df$V1[match(h.net_bcphc.plot.edges$to, h.net_bcphc.plot.df$Id)]  
# set maximum y values
h.net_bcphc.plot.edges$to.y <- h.net_bcphc.plot.df$V2[match(h.net_bcphc.plot.edges$to, h.net_bcphc.plot.df$Id)]
# link species to type information (core or rare species)
h.net_bcphc.plot.edges$species_type <- ifelse(h.net_bcphc.plot.edges$from %in% background_rare_h_bcphc$Species, "rare", 
                                              ifelse(h.net_bcphc.plot.edges$from %in% background_core_h_bcphc$Species, "core", "undefined"))


# cf, background core and rare 
# store edge list with new name
cf_rare_edges_bcphc <- df_cf_rare_bcphc_COR_edges_short
# store node list with new name
cf_rare_nodes_bcphc <- df_cf_rare_bcphc_nodes
# convert species type from class character to class numeric
cf_rare_nodes_bcphc$species_type <- as.factor(as.character(cf_rare_nodes_bcphc$species_type))
# make graph object from data frame
cf_rare_net_bcphc <- graph_from_data_frame(d=cf_rare_edges_bcphc, vertices=cf_rare_nodes_bcphc, directed = directed_network)
# simplify graph by removing multiple edges and loops
cf_rare_net_bcphc <- simplify(cf_rare_net_bcphc, remove.multiple = T, remove.loops = T)
# insert colour code (rare species = blue, core species = orange)
V(cf_rare_net_bcphc)$color <- ifelse(V(cf_rare_net_bcphc)$species_type == "rare", "lightsteelblue2", "springgreen4")
# introduce graph layout (Fruchterman Reingold)
cf.net_bcphc.plot <- layout.fruchterman.reingold(cf_rare_net_bcphc)
# convert to data frame
cf.net_bcphc.plot.df <- as.data.frame(cf.net_bcphc.plot)
# add species ID
cf.net_bcphc.plot.df$Id <- cf_rare_nodes_bcphc$Id
# add species type (core or rare)
cf.net_bcphc.plot.df$species_type <- cf_rare_nodes_bcphc$species_type
# add node size (fixed)
cf.net_bcphc.plot.df$node_size <- V(cf_rare_net_bcphc)$size

# generate data frame for plotting
cf.net_bcphc.plot.edges <- get.data.frame(cf_rare_net_bcphc)
# set minimum x values
cf.net_bcphc.plot.edges$from.x <- cf.net_bcphc.plot.df$V1[match(cf.net_bcphc.plot.edges$from, cf.net_bcphc.plot.df$Id)]  
# set minimum y values
cf.net_bcphc.plot.edges$from.y <- cf.net_bcphc.plot.df$V2[match(cf.net_bcphc.plot.edges$from, cf.net_bcphc.plot.df$Id)]
# set maximum x values
cf.net_bcphc.plot.edges$to.x <- cf.net_bcphc.plot.df$V1[match(cf.net_bcphc.plot.edges$to, cf.net_bcphc.plot.df$Id)]  
# set maximum y values
cf.net_bcphc.plot.edges$to.y <- cf.net_bcphc.plot.df$V2[match(cf.net_bcphc.plot.edges$to, cf.net_bcphc.plot.df$Id)]
# add species type information (core or rare)
cf.net_bcphc.plot.edges$species_type <- ifelse(cf.net_bcphc.plot.edges$from %in% background_rare_cf_bcphc$Species, "rare", 
                                         ifelse(cf.net_bcphc.plot.edges$from %in% background_core_cf_bcphc$Species, "core", "undefined"))

# generate network plot (healthy, BCPHC-normalised)
h_net_bcphc_plot <- ggplot() + 
  geom_segment(data = h.net_bcphc.plot.edges, aes(x = from.x, xend = to.x, y = from.y, yend = to.y, colour=species_type), alpha=0.6) +
  geom_point(data = h.net_bcphc.plot.df, aes(x=V1,y=V2,colour=species_type),size=2) +  
  geom_label_repel(data=h.net_bcphc.plot.df, aes(x=V1,y=V2,label=Id, colour=species_type), fontface="italic", size=2) +
  theme_pubr(border=FALSE, legend = "none")  +
  scale_colour_manual(values=c("core"="darkorange", "rare"="darkblue")) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

# generate network plot (CF, BCPHC-normalised)
cf_net_bcphc_plot <- ggplot() + 
  geom_segment(data = cf.net_bcphc.plot.edges, aes(x = from.x, xend = to.x, y = from.y, yend = to.y, colour=species_type), alpha=0.6) +
  geom_point(data = cf.net_bcphc.plot.df, aes(x=V1,y=V2,colour=species_type),size=2) +  
  geom_label_repel(data=cf.net_bcphc.plot.df, aes(x=V1,y=V2,label=Id, colour=species_type), fontface="italic", size=2) +
  theme_pubr(border=FALSE, legend = "none")  +
  scale_colour_manual(values=c("core"="darkorange", "rare"="darkblue")) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

############################################################################################################
# Correlation between rare and core species, RLE normalised
# generate correlation matrix
# healthy, extract background rare species (as defined by pangenome and one-strain per species)
df_h_rare_rle <- subset(ds_pangenome_rle, rownames(ds_pangenome_rle) %in% background_core_h_rle$Species | rownames(ds_pangenome_rle) %in% background_rare_h_rle$Species)
# transpose data frame
df_h_rare_rle <- data.frame(t(df_h_rare_rle))
# add information based on disease state (CF vs healthy)
df_h_rare_rle$state <- md$State
# subset dataframe and extract healthy samples
df_h_rare_rle <- subset(df_h_rare_rle, state == "Healthy")
# remove non-numeric column
df_h_rare_rle$state <- NULL
# undergo Spearman's rank correlation analysis
df_h_rare_rle <- rcorr(as.matrix(df_h_rare_rle), type = 'spearman')

# create node and edge lists, healthy
# extract and store p-value of correlation analysis
df_h_rare_rle_Pval <- df_h_rare_rle$P
# extract and store correlation coefficient
df_h_rare_rle_COR <- df_h_rare_rle$r

# generate edge list
# take data in wide format and stack columns into a single column of data for p-values
df_h_rare_rle_COR_edges <- reshape2::melt(df_h_rare_rle_Pval)
# take data in wide format and stack columns into a single column of data for correlation coefficient
df_h_rare_rle_coredges <- reshape2::melt(df_h_rare_rle_COR)
# merge data tables
df_h_rare_rle_COR_edges$COR <- df_h_rare_rle_coredges$value
# round p-value to five decimal places
df_h_rare_rle_COR_edges$value <- round(df_h_rare_rle_COR_edges$value, 5)
# store row names in own column
df_h_rare_rle_COR_edges$Label <- row.names(df_h_rare_rle_coredges)
# re-name columns
colnames(df_h_rare_rle_COR_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index row names
rownames(df_h_rare_rle_COR_edges) <- NULL
# remove NAs
df_h_rare_rle_COR_edges <- df_h_rare_rle_COR_edges[complete.cases(df_h_rare_rle_COR_edges), ]
# extract significant correlation as defined in global variables
df_h_rare_rle_COR_edges_short <- subset(df_h_rare_rle_COR_edges, pValue < 0.05 & Weight >= weight_val_pos | Weight <= weight_val_neg)
# if correlation coefficient is larger than 0, add "positive correlation" otherwise "negative correlation"
df_h_rare_rle_COR_edges_short$cortype <- ifelse(df_h_rare_rle_COR_edges_short$Weight > 0, "pos", "neg")
# clean names
df_h_rare_rle_COR_edges_short$Target <- str_replace(df_h_rare_rle_COR_edges_short$Target, "\\.", " ")
df_h_rare_rle_COR_edges_short$Source <- str_replace(df_h_rare_rle_COR_edges_short$Source, "\\.", " ")
# take the two columns to generate node list
df_h_rare_rle_nodes <- select(df_h_rare_rle_COR_edges_short, c("Source", "cortype"))
# remove duplicate entries
df_h_rare_rle_nodes = df_h_rare_rle_nodes[!duplicated(df_h_rare_rle_nodes$Source),]
# make data frame
df_h_rare_rle_nodes <- data.frame(df_h_rare_rle_nodes)
# re-index row names
rownames(df_h_rare_rle_nodes) <- NULL
# re-name columns
colnames(df_h_rare_rle_nodes) <- c("Id", "CorType")
# clean up IDs
df_h_rare_rle_nodes$Id <- str_replace(df_h_rare_rle_nodes$Id, "\\.", " ")
# if Species is in "rare species" list, add "rare", otherwise "core"
df_h_rare_rle_nodes$species_type <- ifelse(df_h_rare_rle_nodes$Id %in% background_rare_h_rle$Species, "rare", "core")
# get genus information
df_h_rare_rle_nodes$Genus <- df_h_rare_rle_nodes$Id
df_h_rare_rle_nodes$Genus <- sapply(strsplit(as.character(df_h_rare_rle_nodes$Genus)," "), `[`, 1)

# CF, background rare
# extract background rare species (as defined by pangenome and one-strain per species)
df_cf_rare_rle <- subset(ds_pangenome_rle, rownames(ds_pangenome_rle) %in% background_core_cf_rle$Species | rownames(ds_pangenome_rle) %in% background_rare_cf_rle$Species)
# transpose data frame
df_cf_rare_rle <- data.frame(t(df_cf_rare_rle))
# add disease state information of child (healthy, CF)
df_cf_rare_rle$state <- md$State
# subset and get CF children
df_cf_rare_rle <- subset(df_cf_rare_rle, state == "CF")
# remove non-numeric column
df_cf_rare_rle$state <- NULL
# undergo Spearman's rank correlation analysis
df_cf_rare_rle <- rcorr(as.matrix(df_cf_rare_rle), type = 'spearman')

# create node and edge lists
# extract and store p-values of correlation analysis
df_cf_rare_rle_Pval <- df_cf_rare_rle$P
# extract and store correlation coefficients of analysis
df_cf_rare_rle_COR <- df_cf_rare_rle$r

# generate edge list
# take data in wide format and stack columns into a single column of data for p-values
df_cf_rare_rle_COR_edges <- reshape2::melt(df_cf_rare_rle_Pval)
# take data in wide format and stack columns into a single column of data for correlation coefficients
df_cf_rare_rle_coredges <- reshape2::melt(df_cf_rare_rle_COR)
# merge data frames
df_cf_rare_rle_COR_edges$COR <- df_cf_rare_rle_coredges$value
# round p-values to five decimal places
df_cf_rare_rle_COR_edges$value <- round(df_cf_rare_rle_COR_edges$value, 5)
# store row names in own column
df_cf_rare_rle_COR_edges$Label <- row.names(df_cf_rare_rle_coredges)
# rename columns
colnames(df_cf_rare_rle_COR_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(df_cf_rare_rle_COR_edges) <- NULL
# remove NAs
df_cf_rare_rle_COR_edges <- df_cf_rare_rle_COR_edges[complete.cases(df_cf_rare_rle_COR_edges), ]
# extract significant correlations as defined by global variables (see above) 
df_cf_rare_rle_COR_edges_short <- subset(df_cf_rare_rle_COR_edges, pValue < sig_level & Weight >= weight_val_pos | Weight <= weight_val_neg)
# if correlation coefficient is larger than 0, add "positive correlation" otherwise "negative correlation"
df_cf_rare_rle_COR_edges_short$cortype <- ifelse(df_cf_rare_rle_COR_edges_short$Weight > 0, "pos", "neg")
# clean column entries
df_cf_rare_rle_COR_edges_short$Target <- str_replace(df_cf_rare_rle_COR_edges_short$Target, "\\.", " ")
df_cf_rare_rle_COR_edges_short$Source <- str_replace(df_cf_rare_rle_COR_edges_short$Source, "\\.", " ")
# select columns to generate node lists
df_cf_rare_rle_nodes <- select(df_cf_rare_rle_COR_edges_short, c("Source", "cortype"))
# remove duplicate rows
df_cf_rare_rle_nodes = df_cf_rare_rle_nodes[!duplicated(df_cf_rare_rle_nodes$Source),]
# convert to data frame object
df_cf_rare_rle_nodes <- data.frame(df_cf_rare_rle_nodes)
# re-index rows
rownames(df_cf_rare_rle_nodes) <- NULL
# rename columns
colnames(df_cf_rare_rle_nodes) <- c("Id", "CorType")
# clean species names
df_cf_rare_rle_nodes$Id <- str_replace(df_cf_rare_rle_nodes$Id, "\\.", " ")
# if species id in background species list, add "rare", otherwise "core"
df_cf_rare_rle_nodes$species_type <- ifelse(df_cf_rare_rle_nodes$Id %in% background_rare_cf_rle$Species, "rare", "core")
# obtain Genus information
df_cf_rare_rle_nodes$Genus <- df_cf_rare_rle_nodes$Id
df_cf_rare_rle_nodes$Genus <- sapply(strsplit(as.character(df_cf_rare_rle_nodes$Genus)," "), `[`, 1)

############################################################################################################
# Make networks, RLE-normalised
# healthy, background core and rare species
# convert column from character to factor object
df_h_rare_rle_nodes$species_type <- as.factor(as.character(df_h_rare_rle_nodes$species_type))

# make graph from data frame
h_rare_net_rle <- graph_from_data_frame(d=df_h_rare_rle_COR_edges_short, vertices=df_h_rare_rle_nodes, directed = directed_network)
# simplify graph by removing multiple edges and loops
h_rare_net_rle <- simplify(h_rare_net_rle, remove.multiple = T, remove.loops = T)
# add colour code (rare = blue, core = green)
V(h_rare_net_rle)$color <- ifelse(V(h_rare_net_rle)$species_type == "rare", "lightsteelblue2", "springgreen4")
# define layout algorithm (Fruchterman Reingold algorithm)
h.net_rle.plot <- layout.fruchterman.reingold(h_rare_net_rle)
# make data frame
h.net_rle.plot.df <- as.data.frame(h.net_rle.plot)
# add species Id
h.net_rle.plot.df$Id <- df_h_rare_rle_nodes$Id
# add species type
h.net_rle.plot.df$species_type <- df_h_rare_rle_nodes$species_type

# make dataframe for network plotting, healthy
h.net_rle.plot.edges <- get.data.frame(h_rare_net_rle)
# set minimum x value
h.net_rle.plot.edges$from.x <- h.net_rle.plot.df$V1[match(h.net_rle.plot.edges$from, h.net_rle.plot.df$Id)]  
# set minimum y values
h.net_rle.plot.edges$from.y <- h.net_rle.plot.df$V2[match(h.net_rle.plot.edges$from, h.net_rle.plot.df$Id)]
# set maximum x values
h.net_rle.plot.edges$to.x <- h.net_rle.plot.df$V1[match(h.net_rle.plot.edges$to, h.net_rle.plot.df$Id)]  
# set maximum y values
h.net_rle.plot.edges$to.y <- h.net_rle.plot.df$V2[match(h.net_rle.plot.edges$to, h.net_rle.plot.df$Id)]
# define species type, if in list rare species background, add "rare", otherwise "core"
h.net_rle.plot.edges$species_type <- ifelse(h.net_rle.plot.edges$from %in% background_rare_h_rle$Species, "rare", 
                                              ifelse(h.net_rle.plot.edges$from %in% background_core_h_rle$Species, "core", "undefined"))

# cf, background core and rare 
# rename edge list, CF, RLE-normalised
cf_rare_edges_rle <- df_cf_rare_rle_COR_edges_short
# rename node list, CF, RLE-normalised
cf_rare_nodes_rle <- df_cf_rare_rle_nodes
# convert species type from character to factor object
cf_rare_nodes_rle$species_type <- as.factor(as.character(cf_rare_nodes_rle$species_type))
# generate graph from data frame
cf_rare_net_rle <- graph_from_data_frame(d=cf_rare_edges_rle, vertices=cf_rare_nodes_rle, directed = directed_network)
# simplify graph by removing multiple edges and loops
cf_rare_net_rle <- simplify(cf_rare_net_rle, remove.multiple = T, remove.loops = T)
# insert colour code (rare species = blue, core species = green)
V(cf_rare_net_rle)$color <- ifelse(V(cf_rare_net_rle)$species_type == "rare", "lightsteelblue2", "springgreen4")
# set graph layout algorithm (Fruchterman Reingold)
cf.net_rle.plot <- layout.fruchterman.reingold(cf_rare_net_rle)
# convert to data frame
cf.net_rle.plot.df <- as.data.frame(cf.net_rle.plot)
# add species ID
cf.net_rle.plot.df$Id <- cf_rare_nodes_rle$Id
# add species type
cf.net_rle.plot.df$species_type <- cf_rare_nodes_rle$species_type
# add column with node size (was defined globally)
cf.net_rle.plot.df$node_size <- V(cf_rare_net_rle)$size

# make data frame for network plotting, CF, RLE-normalised
cf.net_rle.plot.edges <- get.data.frame(cf_rare_net_rle)
# set minimum x values
cf.net_rle.plot.edges$from.x <- cf.net_rle.plot.df$V1[match(cf.net_rle.plot.edges$from, cf.net_rle.plot.df$Id)]  
# set minimum y values
cf.net_rle.plot.edges$from.y <- cf.net_rle.plot.df$V2[match(cf.net_rle.plot.edges$from, cf.net_rle.plot.df$Id)]
# set maximum x values
cf.net_rle.plot.edges$to.x <- cf.net_rle.plot.df$V1[match(cf.net_rle.plot.edges$to, cf.net_rle.plot.df$Id)]  
# set maximum x values
cf.net_rle.plot.edges$to.y <- cf.net_rle.plot.df$V2[match(cf.net_rle.plot.edges$to, cf.net_rle.plot.df$Id)]
# if species in list rare background species, add "rare", otherwise add "core"
cf.net_rle.plot.edges$species_type <- ifelse(cf.net_rle.plot.edges$from %in% background_rare_cf_rle$Species, "rare", 
                                               ifelse(cf.net_rle.plot.edges$from %in% background_core_cf_rle$Species, "core", "undefined"))

# generate network plot (healthy, RLE-normalised)
h_net_rle_plot <- ggplot() + 
  geom_segment(data = h.net_rle.plot.edges, aes(x = from.x, xend = to.x, y = from.y, yend = to.y, colour=species_type), alpha=0.6) +
  geom_point(data = h.net_rle.plot.df, aes(x=V1,y=V2,colour=species_type),size=2) +  
  geom_label_repel(data=h.net_rle.plot.df, aes(x=V1,y=V2,label=Id, colour=species_type), fontface="italic", size=2) +
  theme_pubr(border=FALSE, legend = "none")  +
  scale_colour_manual(values=c("core"="darkorange", "rare"="darkblue")) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

# generate network plot (CF, RLE-normalised)
cf_net_rle_plot <- ggplot() + 
  geom_segment(data = cf.net_rle.plot.edges, aes(x = from.x, xend = to.x, y = from.y, yend = to.y, colour=species_type), alpha=0.6) +
  geom_point(data = cf.net_rle.plot.df, aes(x=V1,y=V2,colour=species_type),size=2) +  
  geom_label_repel(data=cf.net_rle.plot.df, aes(x=V1,y=V2,label=Id, colour=species_type), fontface="italic", size=2) +
  theme_pubr(border=FALSE, legend = "none")  +
  scale_colour_manual(values=c("core"="darkorange", "rare"="darkblue")) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

############################################################################################################

# Correlation between rare and core species, VST- normalised
# generate correlation matrix
# healthy, extract background rare species (as defined by pangenome and one-strain per species)
df_h_rare_vst <- subset(ds_pangenome_vst, rownames(ds_pangenome_vst) %in% background_core_h_vst$Species | rownames(ds_pangenome_vst) %in% background_rare_h_vst$Species)
# transpose data frame
df_h_rare_vst <- data.frame(t(df_h_rare_vst))
# add column with disease state (CF, healthy)
df_h_rare_vst$state <- md$State
# subset and extract healthy samples
df_h_rare_vst <- subset(df_h_rare_vst, state == "Healthy")
# remove non-numeric column
df_h_rare_vst$state <- NULL
# undergo Spearman's rank correlation analysis
df_h_rare_vst <- rcorr(as.matrix(df_h_rare_vst), type = 'spearman')

# create node and edge lists, healthy, VST-normalised
# extract and store p-values of correlation analysis
df_h_rare_vst_Pval <- df_h_rare_vst$P
# extract and store correlation coefficients of analysis
df_h_rare_vst_COR <- df_h_rare_vst$r

# generate edge list
# take data in wide format and stack columns into a single column of data for p-values
df_h_rare_vst_COR_edges <- reshape2::melt(df_h_rare_vst_Pval)
# take data in wide format and stack columns into a single column of data for correlation coefficients
df_h_rare_vst_coredges <- reshape2::melt(df_h_rare_vst_COR)
# merge data frames
df_h_rare_vst_COR_edges$COR <- df_h_rare_vst_coredges$value
# round p-values to five decimal places
df_h_rare_vst_COR_edges$value <- round(df_h_rare_vst_COR_edges$value, 5)
# store row names in own label
df_h_rare_vst_COR_edges$Label <- row.names(df_h_rare_vst_coredges)
# add column names
colnames(df_h_rare_vst_COR_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index row names
rownames(df_h_rare_vst_COR_edges) <- NULL
# remove NAs
df_h_rare_vst_COR_edges <- df_h_rare_vst_COR_edges[complete.cases(df_h_rare_vst_COR_edges), ]
# extract significant correlations as globally defined
df_h_rare_vst_COR_edges_short <- subset(df_h_rare_vst_COR_edges, pValue < sig_level & Weight >= weight_val_pos | Weight <= weight_val_neg)
# if correlation coefficient is larger than 0, add "positive", otherwise "negative"
df_h_rare_vst_COR_edges_short$cortype <- ifelse(df_h_rare_vst_COR_edges_short$Weight > 0, "pos", "neg")
# clean names of column items
df_h_rare_vst_COR_edges_short$Target <- str_replace(df_h_rare_vst_COR_edges_short$Target, "\\.", " ")
df_h_rare_vst_COR_edges_short$Source <- str_replace(df_h_rare_vst_COR_edges_short$Source, "\\.", " ")
# select both columns to generate node list
df_h_rare_vst_nodes <- select(df_h_rare_vst_COR_edges_short, c("Source", "cortype"))
# remove duplicates
df_h_rare_vst_nodes = df_h_rare_vst_nodes[!duplicated(df_h_rare_vst_nodes$Source),]
# make data frame
df_h_rare_vst_nodes <- data.frame(df_h_rare_vst_nodes)
# re-index rows
rownames(df_h_rare_vst_nodes) <- NULL
# rename columns
colnames(df_h_rare_vst_nodes) <- c("Id", "CorType")
# clean IDs names so that they match edge list IDs
df_h_rare_vst_nodes$Id <- str_replace(df_h_rare_vst_nodes$Id, "\\.", " ")
# if id in rare species background list, add "rare" otherwise "core"
df_h_rare_vst_nodes$species_type <- ifelse(df_h_rare_vst_nodes$Id %in% background_rare_h_vst$Species, "rare", "core")
# add genus information
df_h_rare_vst_nodes$Genus <- df_h_rare_vst_nodes$Id
df_h_rare_vst_nodes$Genus <- sapply(strsplit(as.character(df_h_rare_vst_nodes$Genus)," "), `[`, 1)


# # CF, background rare
# CF, extract background rare species (as defined by pangenome and one-strain per species)
df_cf_rare_vst <- subset(ds_pangenome_vst, rownames(ds_pangenome_vst) %in% background_core_cf_vst$Species | rownames(ds_pangenome_vst) %in% background_rare_cf_vst$Species)
# transpose data frame
df_cf_rare_vst <- data.frame(t(df_cf_rare_vst))
# add disease state (CF, healthy)
df_cf_rare_vst$state <- md$State
# subset and extract CF samples
df_cf_rare_vst <- subset(df_cf_rare_vst, state == "CF")
# remove non-numeric columns
df_cf_rare_vst$state <- NULL
# perform Spearman's rank correlation analysis
df_cf_rare_vst <- rcorr(as.matrix(df_cf_rare_vst), type = 'spearman')

# create node and edge lists
# extract and store p-values of correlation analysis
df_cf_rare_vst_Pval <- df_cf_rare_vst$P
# extract and store correlation coefficients of analysis
df_cf_rare_vst_COR <- df_cf_rare_vst$r

# generate edge list
# take data in wide format and stack columns into a single column of data for p-values
df_cf_rare_vst_COR_edges <- reshape2::melt(df_cf_rare_vst_Pval)
# take data in wide format and stack columns into a single column of data for correlation coefficients
df_cf_rare_vst_coredges <- reshape2::melt(df_cf_rare_vst_COR)
# merge tables
df_cf_rare_vst_COR_edges$COR <- df_cf_rare_vst_coredges$value
# round p-values to five decimal places
df_cf_rare_vst_COR_edges$value <- round(df_cf_rare_vst_COR_edges$value, 5)
# store row names in own column
df_cf_rare_vst_COR_edges$Label <- row.names(df_cf_rare_vst_coredges)
# rename columns
colnames(df_cf_rare_vst_COR_edges) <- c('Target', 'Source', 'pValue', 'Weight', 'Id')
# re-index rows
rownames(df_cf_rare_vst_COR_edges) <- NULL
# remove NAs
df_cf_rare_vst_COR_edges <- df_cf_rare_vst_COR_edges[complete.cases(df_cf_rare_vst_COR_edges), ]
# extract significant correlations
df_cf_rare_vst_COR_edges_short <- subset(df_cf_rare_vst_COR_edges, pValue < sig_level & Weight >= weight_val_pos | Weight <= weight_val_neg)
# if correlation coefficient is larger than 0, add "positive correlation", otherwise "negative correlation"
df_cf_rare_vst_COR_edges_short$cortype <- ifelse(df_cf_rare_vst_COR_edges_short$Weight > 0, "pos", "neg")
# clean names of column entries
df_cf_rare_vst_COR_edges_short$Target <- str_replace(df_cf_rare_vst_COR_edges_short$Target, "\\.", " ")
df_cf_rare_vst_COR_edges_short$Source <- str_replace(df_cf_rare_vst_COR_edges_short$Source, "\\.", " ")
# select both columns to create node list
df_cf_rare_vst_nodes <- select(df_cf_rare_vst_COR_edges_short, c("Source", "cortype"))
# remove duplicate entries
df_cf_rare_vst_nodes = df_cf_rare_vst_nodes[!duplicated(df_cf_rare_vst_nodes$Source),]
# make data frame
df_cf_rare_vst_nodes <- data.frame(df_cf_rare_vst_nodes)
# re-index rows
rownames(df_cf_rare_vst_nodes) <- NULL
# rename columns
colnames(df_cf_rare_vst_nodes) <- c("Id", "CorType")
# clean up id, so that they match id entries of edge list
df_cf_rare_vst_nodes$Id <- str_replace(df_cf_rare_vst_nodes$Id, "\\.", " ")
# if id in background rare species list, add "rare" otherwise add "core"
df_cf_rare_vst_nodes$species_type <- ifelse(df_cf_rare_vst_nodes$Id %in% background_rare_cf_vst$Species, "rare", "core")
# add genus information
df_cf_rare_vst_nodes$Genus <- df_cf_rare_vst_nodes$Id
df_cf_rare_vst_nodes$Genus <- sapply(strsplit(as.character(df_cf_rare_vst_nodes$Genus)," "), `[`, 1)

############################################################################################################
# Generate networks 
# healthy, background core and rare species ####
# convert species type (core, rare) from character to factor
df_h_rare_vst_nodes$species_type <- as.factor(as.character(df_h_rare_vst_nodes$species_type))
# generate graph from data frame
h_rare_net_vst <- graph_from_data_frame(d=df_h_rare_vst_COR_edges_short, vertices=df_h_rare_vst_nodes, directed = directed_network)
# simplify graph by removing multiple edges and loops
h_rare_net_vst <- simplify(h_rare_net_vst, remove.multiple = T, remove.loops = T)
# add colour code with rare = blue and core = green
V(h_rare_net_vst)$color <- ifelse(V(h_rare_net_vst)$species_type == "rare", "lightsteelblue2", "springgreen4")
# apply Fruchterman eingold graph 
h.net_vst.plot <- layout.fruchterman.reingold(h_rare_net_vst)
# convert to data frame
h.net_vst.plot.df <- as.data.frame(h.net_vst.plot)
# add column with species id
h.net_vst.plot.df$Id <- df_h_rare_vst_nodes$Id
# add column with species type (core and rare)
h.net_vst.plot.df$species_type <- df_h_rare_vst_nodes$species_type

# make data frame for network plotting
h.net_vst.plot.edges <- get.data.frame(h_rare_net_vst)
# determine minimum x
h.net_vst.plot.edges$from.x <- h.net_vst.plot.df$V1[match(h.net_vst.plot.edges$from, h.net_vst.plot.df$Id)]  
# determine minimum y
h.net_vst.plot.edges$from.y <- h.net_vst.plot.df$V2[match(h.net_vst.plot.edges$from, h.net_vst.plot.df$Id)]
# determine maximum x
h.net_vst.plot.edges$to.x <- h.net_vst.plot.df$V1[match(h.net_vst.plot.edges$to, h.net_vst.plot.df$Id)]  
# determine maximum y
h.net_vst.plot.edges$to.y <- h.net_vst.plot.df$V2[match(h.net_vst.plot.edges$to, h.net_vst.plot.df$Id)]
# if species in list of rare background species, add "rare" otherwise add "core"
h.net_vst.plot.edges$species_type <- ifelse(h.net_vst.plot.edges$from %in% background_rare_h_vst$Species, "rare", 
                                            ifelse(h.net_vst.plot.edges$from %in% background_core_h_vst$Species, "core", "undefined"))

# cf, background core and rare
#  rename edge list
cf_rare_edges_vst <- df_cf_rare_vst_COR_edges_short
# rename node list
cf_rare_nodes_vst <- df_cf_rare_vst_nodes
# convert species type (core, rare) from class character to class factor
cf_rare_nodes_vst$species_type <- as.factor(as.character(cf_rare_nodes_vst$species_type))
# convert data frame to graph
cf_rare_net_vst <- graph_from_data_frame(d=cf_rare_edges_vst, vertices=cf_rare_nodes_vst, directed = directed_network)
# simplify graph by removing multiple edges and loops
cf_rare_net_vst <- simplify(cf_rare_net_vst, remove.multiple = T, remove.loops = T)
# set colour code (rare species = blue, core species = green)
V(cf_rare_net_vst)$color <- ifelse(V(cf_rare_net_vst)$species_type == "rare", "lightsteelblue2", "springgreen4")
# apply Fruchterman Reingold layout algorithm
cf.net_vst.plot <- layout.fruchterman.reingold(cf_rare_net_vst)
# conver to data frame
cf.net_vst.plot.df <- as.data.frame(cf.net_vst.plot)
# add species ID
cf.net_vst.plot.df$Id <- cf_rare_nodes_vst$Id
# add species type (core, rare)
cf.net_vst.plot.df$species_type <- cf_rare_nodes_vst$species_type
# set node size (is globally defined)
cf.net_vst.plot.df$node_size <- V(cf_rare_net_vst)$size

# convert to dataframe for network plotting
cf.net_vst.plot.edges <- get.data.frame(cf_rare_net_vst)
# determine minimum x values
cf.net_vst.plot.edges$from.x <- cf.net_vst.plot.df$V1[match(cf.net_vst.plot.edges$from, cf.net_vst.plot.df$Id)]  
# determine minimum y values
cf.net_vst.plot.edges$from.y <- cf.net_vst.plot.df$V2[match(cf.net_vst.plot.edges$from, cf.net_vst.plot.df$Id)]
# determine maximum x values
cf.net_vst.plot.edges$to.x <- cf.net_vst.plot.df$V1[match(cf.net_vst.plot.edges$to, cf.net_vst.plot.df$Id)] 
# determine maximum y values
cf.net_vst.plot.edges$to.y <- cf.net_vst.plot.df$V2[match(cf.net_vst.plot.edges$to, cf.net_vst.plot.df$Id)]
# add species type, if id in background rare species list, add "rare", otherwise add "core"
cf.net_vst.plot.edges$species_type <- ifelse(cf.net_vst.plot.edges$from %in% background_rare_cf_vst$Species, "rare", 
                                             ifelse(cf.net_vst.plot.edges$from %in% background_core_cf_vst$Species, "core", "undefined"))

# generate network plot, healthy, VST-normalised
h_net_vst_plot <- ggplot() + geom_segment(data = h.net_vst.plot.edges, aes(x = from.x, xend = to.x, y = from.y, yend = to.y, colour=species_type), alpha=0.6) +
  geom_point(data = h.net_vst.plot.df, aes(x=V1,y=V2,colour=species_type),size=2) +  
  geom_label_repel(data=h.net_vst.plot.df, aes(x=V1,y=V2,label=Id, colour=species_type), fontface="italic", size=2) +
  theme_pubr(border=FALSE, legend = "none")  +
  scale_colour_manual(values=c("core"="darkorange", "rare"="darkblue")) +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

# generate network plot, CF, VST-normalised
cf_net_vst_plot <- ggplot() + 
  geom_segment(data = cf.net_vst.plot.edges, aes(x = from.x, xend = to.x, y = from.y, yend = to.y, colour=species_type), alpha=0.6) +
  geom_point(data = cf.net_vst.plot.df, aes(x=V1,y=V2,colour=species_type),size=2) +  
  geom_label_repel(data=cf.net_vst.plot.df, aes(x=V1,y=V2,label=Id, colour=species_type), fontface="italic", size=2) +
  theme_pubr(border=FALSE, legend = "none")  +
  scale_colour_manual(values=c("core"="darkorange", "rare"="darkblue")) +
  theme(axis.text.x = element_blank(),  axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

############################################################################################################
# robustness and vulnerability analysis (targeted attack)
# healthy, background core and rare species
# BCPHC-normalised data
set.seed(123)
# assesses networks vulnerability due either to random breakdowns or to intentional attacks with 10 iterations
h_rare_net_bcphc_attack <-swan_combinatory(h_rare_net_bcphc,10)
# convert to data frame
h_rare_net_bcphc_attack <- data.frame(h_rare_net_bcphc_attack)
# add column names
colnames(h_rare_net_bcphc_attack) <- c("fraction_removed", "betweeness", "degree", "cascading", "random")
# remove duplicate rows by taking median of duplicate values
h_rare_net_bcphc_attack <- ddply(h_rare_net_bcphc_attack,"fraction_removed",numcolwise(median))
# insert a first row to have starting point at 0
h_rare_net_bcphc_attack[1, ] <- c(0.0, 0.0, 0.0, 0.0, 0.0)

# cf, background core and rare species
set.seed(123)
# assesses networks vulnerability due either to random breakdowns or to intentional attacks with 10 iterations
cf_rare_net_bcphc_attack <-swan_combinatory(cf_rare_net_bcphc,10)
# convert to data frame
cf_rare_net_bcphc_attack <- data.frame(cf_rare_net_bcphc_attack)
# rename columns
colnames(cf_rare_net_bcphc_attack) <- c("fraction_removed", "betweeness", "degree", "cascading", "random")
# remove duplicate rows by taking median of duplicate values
cf_rare_net_bcphc_attack <- ddply(cf_rare_net_bcphc_attack,"fraction_removed",numcolwise(median))
# insert a first row to have starting point at 0
cf_rare_net_bcphc_attack <- rbind(h_rare_net_bcphc_attack[1,], cf_rare_net_bcphc_attack)

# similarity measures
# choose parameters that will be extracted for analysis
select_items_bcphc_1 <- c("fraction_removed", "degree")
select_items_bcphc_2 <- c("fraction_removed", "random")
# subset fraction_removed and degree columns, CF, BCPHC-normalised
cf_rare_net_bcphc_matrix_degree <- select(cf_rare_net_bcphc_attack, select_items_bcphc_1)
# convert data frame to matrix
cf_rare_net_bcphc_matrix_degree <- as.matrix(cf_rare_net_bcphc_matrix_degree)
# subset fraction_removed and degree columns, Healthy, BCPHC-normalised
h_rare_net_bcphc_matrix_degree <- select(h_rare_net_bcphc_attack, select_items_bcphc_1)
# convert data frame to matrix
h_rare_net_bcphc_matrix_degree <- as.matrix(h_rare_net_bcphc_matrix_degree)
# calculate Frechet distance between healthy and CF attack curves, bcphc-normalised data
original_frechet_h_cf_degree_bcphc <- Frechet(h_rare_net_bcphc_matrix_degree, cf_rare_net_bcphc_matrix_degree, testLeash = -1)
# convert degree to data frame, CF, BCPHC-normalised
cf_rare_net_bcphc_matrix_degree <- data.frame(cf_rare_net_bcphc_matrix_degree)
# add meta data
cf_rare_net_bcphc_matrix_degree$normalisation <- "BCPHC"
cf_rare_net_bcphc_matrix_degree$state <- "CF"
# convert degree to data frame, Healthy, BCPHC-normalised
h_rare_net_bcphc_matrix_degree <- data.frame(h_rare_net_bcphc_matrix_degree)
# add meta data
h_rare_net_bcphc_matrix_degree$normalisation <- "BCPHC"
h_rare_net_bcphc_matrix_degree$state <- "Healthy"

# healthy, background core and rare species
# RLE-normalised data
set.seed(123)
# assesses networks vulnerability due either to random breakdowns or to intentional attacks with 10 iterations
h_rare_net_rle_attack <-swan_combinatory(h_rare_net_rle,10)
# convert to data frame
h_rare_net_rle_attack <- data.frame(h_rare_net_rle_attack)
# rename columns
colnames(h_rare_net_rle_attack) <- c("fraction_removed", "betweeness", "degree", "cascading", "random")
# remove duplicate rows by taking median of duplicate values
h_rare_net_rle_attack <- ddply(h_rare_net_rle_attack,"fraction_removed",numcolwise(median))
# insert a first row to have starting point at 0
h_rare_net_rle_attack[1, ] <- c(0.0, 0.0, 0.0, 0.0, 0.0)

# cf, background core and rare species
set.seed(123)
# assesses networks vulnerability due either to random breakdowns or to intentional attacks with 10 iterations
cf_rare_net_rle_attack <-swan_combinatory(cf_rare_net_rle,10)
cf_rare_net_rle_attack <- data.frame(cf_rare_net_rle_attack)
# rename columns
colnames(cf_rare_net_rle_attack) <- c("fraction_removed", "betweeness", "degree", "cascading", "random")
# remove duplicate rows by taking median of duplicate values
cf_rare_net_rle_attack <- ddply(cf_rare_net_rle_attack,"fraction_removed",numcolwise(median))
# insert a first row to have starting point at 0
cf_rare_net_rle_attack <- rbind(h_rare_net_rle_attack[1,], cf_rare_net_rle_attack)

# similarity measures
# choose parameters that will be extracted for analysis
select_items_rle_1 <- c("fraction_removed", "degree")
select_items_rle_2 <- c("fraction_removed", "random")
# subset fraction_removed and degree columns, CF, RLE-normalised
cf_rare_net_rle_matrix_degree <- select(cf_rare_net_rle_attack, select_items_rle_1)
# convert to data matrix
cf_rare_net_rle_matrix_degree <- as.matrix(cf_rare_net_rle_matrix_degree)
# subset fraction_removed and degree columns, Healthy, RLE-normalised
h_rare_net_rle_matrix_degree <- select(h_rare_net_rle_attack, select_items_rle_1)
# convert to data matrix
h_rare_net_rle_matrix_degree <- as.matrix(h_rare_net_rle_matrix_degree)
# calculate Frechet distance between healthy and CF attack curve
original_frechet_h_cf_degree_rle <- Frechet(h_rare_net_rle_matrix_degree, cf_rare_net_rle_matrix_degree, testLeash = -1)
# convert degree to data frame, CF, RLE-normalised
cf_rare_net_rle_matrix_degree <- data.frame(cf_rare_net_rle_matrix_degree)
# add meta data
cf_rare_net_rle_matrix_degree$normalisation <- "RLE"
cf_rare_net_rle_matrix_degree$state <- "CF"
# convert degree to data frame, Healthy, BCPHC-normalised
h_rare_net_rle_matrix_degree <- data.frame(h_rare_net_rle_matrix_degree)
# add meta data
h_rare_net_rle_matrix_degree$normalisation <- "RLE"
h_rare_net_rle_matrix_degree$state <- "Healthy"


# healthy, background core and rare species
# VST-normalised data
set.seed(123)
# assesses networks vulnerability due either to random breakdowns or to intentional attacks with 10 iterations
h_rare_net_vst_attack <-swan_combinatory(h_rare_net_vst,10)
# generate data frame
h_rare_net_vst_attack <- data.frame(h_rare_net_vst_attack)
# rename columns
colnames(h_rare_net_vst_attack) <- c("fraction_removed", "betweeness", "degree", "cascading", "random")
# remove duplicate rows by taking median of duplicate values
h_rare_net_vst_attack <- ddply(h_rare_net_vst_attack,"fraction_removed", numcolwise(median))
# insert new first row to start at value 0
h_rare_net_vst_attack[1, ] <- c(0.0, 0.0, 0.0, 0.0, 0.0)

# cf, background core and rare species
set.seed(123)
# assesses networks vulnerability due either to random breakdowns or to intentional attacks with 10 iterations
cf_rare_net_vst_attack <-swan_combinatory(cf_rare_net_vst,10)
# convert to data frame
cf_rare_net_vst_attack <- data.frame(cf_rare_net_vst_attack)
# rename columns
colnames(cf_rare_net_vst_attack) <- c("fraction_removed", "betweeness", "degree", "cascading", "random")
# remove duplicate rows by taking median of duplicate values
cf_rare_net_vst_attack <- ddply(cf_rare_net_vst_attack,"fraction_removed",numcolwise(median))
# insert a first row to have starting point at 0
cf_rare_net_vst_attack <- rbind(h_rare_net_vst_attack[1,], cf_rare_net_vst_attack)

# similarity measures
# choose parameters that will be extracted for analysis
select_items_vst_1 <- c("fraction_removed", "degree")
select_items_vst_2 <- c("fraction_removed", "random")
# subset fraction_removed and degree columns, CF, VST-normalised
cf_rare_net_vst_matrix_degree <- select(cf_rare_net_vst_attack, select_items_vst_1)
# convert to matrix
cf_rare_net_vst_matrix_degree <- as.matrix(cf_rare_net_vst_matrix_degree)
# subset fraction_removed and degree columns, Healthy, VST-normalised
h_rare_net_vst_matrix_degree <- select(h_rare_net_vst_attack, select_items_vst_1)

# convert to matrix
h_rare_net_vst_matrix_degree <- as.matrix(h_rare_net_vst_matrix_degree)
#  calculate Frechet distance between CF and healthy network, VST-normalised data
original_frechet_h_cf_degree_vst <- Frechet(h_rare_net_vst_matrix_degree, cf_rare_net_vst_matrix_degree, testLeash = -1)

# convert degree to data frame, CF, VST-normalised
cf_rare_net_vst_matrix_degree <- data.frame(cf_rare_net_vst_matrix_degree)
# add meta data
cf_rare_net_vst_matrix_degree$normalisation <- "VST"
cf_rare_net_vst_matrix_degree$state <- "CF"
# convert degree to data frame, Healthy, VST-normalised
h_rare_net_vst_matrix_degree <- data.frame(h_rare_net_vst_matrix_degree)
# add meta data
h_rare_net_vst_matrix_degree$normalisation <- "VST"
h_rare_net_vst_matrix_degree$state <- "Healthy"

# merge all attack matrices, BCPHC, RLE and VST-normalised data
all_matrix_degree <- data.frame(rbind(cf_rare_net_bcphc_matrix_degree, cf_rare_net_rle_matrix_degree,
                                      h_rare_net_bcphc_matrix_degree, h_rare_net_rle_matrix_degree, h_rare_net_vst_matrix_degree))

############################################################################################################

# Continue robustness and vulnerability analysis for random attacks
# subset fraction_removed and random columns, CF, BCPHC-normalised
cf_rare_net_bcphc_matrix_random <- select(cf_rare_net_bcphc_attack, select_items_bcphc_2)
# convert to matrix
cf_rare_net_bcphc_matrix_random <- as.matrix(cf_rare_net_bcphc_matrix_random)
# subset fraction_removed and random columns, Healthy, BCPHC-normalised
h_rare_net_bcphc_matrix_random <- select(h_rare_net_bcphc_attack, select_items_bcphc_2)
# convert to matrix
h_rare_net_bcphc_matrix_random <- as.matrix(h_rare_net_bcphc_matrix_random)
# calculate Frechet distance between healthy and CF random attack curves
original_frechet_h_cf_random_bcphc <- Frechet(h_rare_net_bcphc_matrix_random, cf_rare_net_bcphc_matrix_random, testLeash = -1)
# convert back to data frame, cf, BCPHC-normalised
cf_rare_net_bcphc_matrix_random <- data.frame(cf_rare_net_bcphc_matrix_random)
# add meta data
cf_rare_net_bcphc_matrix_random$normalisation <- "BCPHC"
cf_rare_net_bcphc_matrix_random$state <- "CF"
# convert back to data frame, healthy, BCPHC-normalised
h_rare_net_bcphc_matrix_random <- data.frame(h_rare_net_bcphc_matrix_random)
# add meta data
h_rare_net_bcphc_matrix_random$normalisation <- "BCPHC"
h_rare_net_bcphc_matrix_random$state <- "Healthy"


# subset fraction_removed and random columns, CF, RLE-normalised
cf_rare_net_rle_matrix_random <- select(cf_rare_net_rle_attack, select_items_rle_2)
# convert to matrix
cf_rare_net_rle_matrix_random <- as.matrix(cf_rare_net_rle_matrix_random)
# subset fraction_removed and random columns, Healthy, RLE-normalised
h_rare_net_rle_matrix_random <- select(h_rare_net_rle_attack, select_items_rle_2)
# convert to matrix
h_rare_net_rle_matrix_random <- as.matrix(h_rare_net_rle_matrix_random) 
# calculate Frechet distance between healthy and CF random attack curves
original_frechet_h_cf_random_rle <- Frechet(h_rare_net_rle_matrix_random, cf_rare_net_rle_matrix_random, testLeash = -1)
# convert back to data frame, CF, RLE-normalised
cf_rare_net_rle_matrix_random <- data.frame(cf_rare_net_rle_matrix_random)
# add meta data
cf_rare_net_rle_matrix_random$normalisation <- "RLE"
cf_rare_net_rle_matrix_random$state <- "CF"
# convert back to data frame, healthy, RLE-normalised
h_rare_net_rle_matrix_random <- data.frame(h_rare_net_rle_matrix_random)
# add meta data
h_rare_net_rle_matrix_random$normalisation <- "RLE"
h_rare_net_rle_matrix_random$state <- "Healthy"


# subset fraction_removed and random columns, CF, VST-normalised
cf_rare_net_vst_matrix_random <- select(cf_rare_net_vst_attack, select_items_vst_2)
# convert to matrix
cf_rare_net_vst_matrix_random <- as.matrix(cf_rare_net_vst_matrix_random)
# subset fraction_removed and random columns, Healthy, VST-normalised
h_rare_net_vst_matrix_random <- select(h_rare_net_vst_attack, select_items_vst_2)
# convert to matrix
h_rare_net_vst_matrix_random <- as.matrix(h_rare_net_vst_matrix_random)
# calculate Frechet distance between healthy and CF random attack curves
original_frechet_h_cf_random_vst <- Frechet(h_rare_net_vst_matrix_random, cf_rare_net_vst_matrix_random, testLeash = -1)
# convert back to data frame, CF, VST-normalised
cf_rare_net_vst_matrix_random <- data.frame(cf_rare_net_vst_matrix_random)
# add meta data
cf_rare_net_vst_matrix_random$normalisation <- "VST"
cf_rare_net_vst_matrix_random$state <- "CF"
# convert back to data frame, Healthy, VST-normalised
h_rare_net_vst_matrix_random <- data.frame(h_rare_net_vst_matrix_random)
# add meta data
h_rare_net_vst_matrix_random$normalisation <- "VST"
h_rare_net_vst_matrix_random$state <- "Healthy"

# merge all random attack outputs (BCPHC, RLE and VST-normalised)
all_matrix_random <- data.frame(rbind(cf_rare_net_bcphc_matrix_random, cf_rare_net_rle_matrix_random,
                                      h_rare_net_bcphc_matrix_random, h_rare_net_rle_matrix_random, h_rare_net_vst_matrix_random))

############################################################################################################

# generate targeted plot 
targeted_attack <- ggplot(all_matrix_degree) +
  geom_vline(aes(xintercept=0.25), colour="grey", linetype="dashed") +
  geom_hline(aes(yintercept=0.44), colour="grey", linetype="dashed") +
  geom_hline(aes(yintercept=0.835), colour="grey", linetype="dashed") +
  geom_point(aes(x=fraction_removed, y=degree, colour=state, shape=normalisation), size=1.2) +
  geom_line(aes(x=fraction_removed, y=degree, colour=state, shape=normalisation), size=0.3) +
  theme_bw() + 
  scale_colour_manual(values=c("darkred", "black")) + xlab("p (targeted)") + ylab("Network efficiency") +
  theme(panel.grid = element_blank(), legend.title = element_blank(), axis.title = element_text(size=11), axis.text = element_text(size=11), 
        legend.text = element_text(size=10, family = "Arial"))

# generate random plot
random_attack <- ggplot(all_matrix_random) +
  geom_vline(aes(xintercept=0.25), colour="grey", linetype="dashed") +
  geom_hline(aes(yintercept=0.44), colour="grey", linetype="dashed") +
  geom_hline(aes(yintercept=0.54), colour="grey", linetype="dashed") +
  geom_point(aes(x=fraction_removed, y=random, colour=state, shape=normalisation), size=1.2) +
  geom_line(aes(x=fraction_removed, y=random, colour=state, shape=normalisation), size=0.3) +
  theme_bw() + 
  scale_colour_manual(values=c("darkred", "black")) + xlab("p (random)") + ylab(" ") +
  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.text = element_text(size=10, family = "Arial"), axis.title = element_text(size=11), 
        axis.text = element_text(size=11))

# merge cf and healthy plot, RLE-normalised
ab1_bcphc <- ggarrange(h_net_rle_plot, cf_net_rle_plot, nrow = 1, ncol = 2, labels = c("A", "B"), label.x = 0.01, label.y = 0.98, font.label = list(size = 12, color = "black")) 

# merge targeted and random attack plots
attack_plots <- ggarrange(targeted_attack, random_attack, labels=c("C", "D"), common.legend = TRUE, font.label = list(size = 12, color = "black"), legend = "right")

# merge cf and healthy network plots, BCPHC-normalised + attack curves of all normalisation steps
first_net_bcphc_plot <- ggarrange(ab1_bcphc, attack_plots, nrow=2, heights = c(1,1))

# merge cf and healthy network plots, VST and RLE-normalised
net_rle_vst <- ggarrange(cf_net_vst_plot, h_net_vst_plot, cf_net_bcphc_plot, h_net_bcphc_plot, heights = c(0.6,0.6), labels=c("A", "B", "C", "D"), font.label = list(size = 12, color ="black"))

############################################################################################################
# Network modulation 
# Set base seed (for reproducibility purposes)
set.seed(111)
# Generate 100 random integers between 1 and 500
random_seeds <- sample(1:500, 100, replace = FALSE) 
# number of nodes to be inserted later on
add_nodes <- c(1:20) 

# define empty variables globally for the loop assignment later on, BCPHC-normalised
frechet_attack_similarity_bcphc = NULL
network_statistics_bcphc = NULL
netwerk_kernel_randomWalk_bcphc = NULL
network_kernel_shortestPath_bcphc = NULL
network_kernel_WL_bcphc = NULL
selected_species_bcphc = NULL
network_nullmodel_bcphc = NULL

# define empty variables globally for the loop assignment later on, RLE-normalised
frechet_attack_similarity_rle = NULL
network_statistics_rle = NULL
netwerk_kernel_randomWalk_rle = NULL
network_kernel_shortestPath_rle = NULL
network_kernel_WL_rle = NULL
selected_species_rle = NULL
network_nullmodel_rle = NULL

# define empty variables globally for the loop assignment later on, VST-normalised
frechet_attack_similarity_vst = NULL
network_statistics_vst = NULL
netwerk_kernel_randomWalk_vst = NULL
network_kernel_shortestPath_vst = NULL
network_kernel_WL_vst = NULL
selected_species_vst = NULL
network_nullmodel_vst = NULL
sample_with_replacement = TRUE 

# loop over list of random seeds
for (i in random_seeds){
  # repeat with an increasing number of species
  for (n_nodes in add_nodes){
    # set seed
    set.seed(i)
    print(n_nodes)
    print(i)
    print(random_seeds)
    # start with BCPHC-normalised data
    # randomly extract species (with replacement) from healthy network, BCPHC
    sample_h_rare_0_bcphc <- sample(V(h_rare_net_bcphc)$name, n_nodes, replace = sample_with_replacement)
    # generate a new graph structure from scratch with the extracted species
    h_net_subset_for_cf_rare_bcphc <- induced_subgraph(h_rare_net_bcphc, 
                                                       which(V(h_rare_net_bcphc)$name %in% sample_h_rare_0_bcphc),
                                                       impl = "create_from_scratch")
    # make new node list
    h_net_subset_for_cf_rare_bcphc_df_nodes <- as_data_frame(h_net_subset_for_cf_rare_bcphc, what ="vertices")
    # make new edge list
    h_net_subset_for_cf_rare_bcphc_df_edges <- as_data_frame(h_net_subset_for_cf_rare_bcphc, what ="edges")
    # convert node information from original CF graph into data frame
    cf_rare_net_bcphc_df_nodes <- as_data_frame(cf_rare_net_bcphc, what ="vertices")
    # convert edge information from original CF graph into data frame
    cf_rare_net_bcphc_df_edges <- as_data_frame(cf_rare_net_bcphc, what ="edges")
    # merge sub-sampled healthy network nodes with CF network nodes
    new_cf_rare_nodes_bcphc <- data.frame(rbind(h_net_subset_for_cf_rare_bcphc_df_nodes, cf_rare_net_bcphc_df_nodes))
    # keep only unique rows from a data frame
    new_cf_rare_nodes_bcphc <- distinct(new_cf_rare_nodes_bcphc)
    # merged sub-sampled healthy network edges with CF network edges
    new_cf_rare_edges_bcphc <- data.frame(rbind(h_net_subset_for_cf_rare_bcphc_df_edges, cf_rare_net_bcphc_df_edges))
    # keep only unique rows from a data frame
    new_cf_rare_edges_bcphc <- distinct(new_cf_rare_edges_bcphc)
    
    # make a new graph from dataframe of modulated CF structure
    new_cf_rare_net_bcphc <- graph_from_data_frame(new_cf_rare_edges_bcphc, directed=directed_network, vertices=new_cf_rare_nodes_bcphc)
    
    # null model
    new_net_bcphc_asdj <- as_adjacency_matrix(new_cf_rare_net_bcphc, type="both")
    new_net_bcphc_matrix <- as.matrix(new_net_bcphc_asdj)
    # Calculate network metric extinction slope
    originalMod_net_bcphc_extinction <- networklevel(new_net_bcphc_matrix, index = "extinction slope", extinctmethod='degree')
    originalMod_net_bcphc_extinction <- originalMod_net_bcphc_extinction[2]
    # make null model
    nullmodelMod_bcphc <- nullmodel(new_net_bcphc_matrix, N=nullmodel_networks_n, method=nullmodel_method)
    nullmodelMod_bcphc_extinction <- unlist(sapply(nullmodelMod_bcphc, networklevel, index = 'extinction slope', extinctmethod='degree'))  
    nullmodelMod_bcphc_extinction <- nullmodelMod_bcphc_extinction[2,]
    # get zscore and p value between original and random network structures
    bcphc_mod_zscore <- net.metric.zscore(originalMod_net_bcphc_extinction, nullmodelMod_bcphc_extinction)
    bcphc_mod_pvalue <- net.metric.pvalue(bcphc_mod_zscore)
    
    # generate final data frame
    nullmodelMod_stats_extinction_bcphc <- data.frame(cbind(bcphc_mod_zscore, bcphc_mod_pvalue))
    nullmodelMod_stats_extinction_bcphc$normlisation <- "BCPHC"
    nullmodelMod_stats_extinction_bcphc$seed <- i
    nullmodelMod_stats_extinction_bcphc$n_nodes <- n_nodes
    nullmodelMod_stats_extinction_bcphc$index <- "Robustness"
    colnames(nullmodelMod_stats_extinction_bcphc) <- c("zscore", "p-value", "normalisation", "seed", "n_nodes", "index")
    rownames(nullmodelMod_stats_extinction_bcphc) <- NULL
    
    # assesses networks vulnerability due to targeted attacks with 10 iterations for assessing random error
    new_cf_net_attack_rare_bcphc <- swan_combinatory(new_cf_rare_net_bcphc, 10)
    # convert to data frame
    new_cf_net_attack_rare_bcphc <- data.frame(new_cf_net_attack_rare_bcphc)
    # rename columns
    colnames(new_cf_net_attack_rare_bcphc) <- c("fraction_removed", "betweeness", "degree", "cascading", "random")
    # remove duplicate rows and take median value of duplicates
    new_cf_net_attack_rare_bcphc <- ddply(new_cf_net_attack_rare_bcphc,"fraction_removed",numcolwise(median))
    # insert first row with only zero values to have zero as starting point
    new_cf_net_attack_rare_bcphc <- rbind(h_rare_net_bcphc_attack[1,], new_cf_net_attack_rare_bcphc)
    # make dataframe with sub-selection options
    select_items_bcphc <- c("fraction_removed", "degree")
    # select items from attack dataframe
    cf_modulated_net_matrix_degree_bcphc <- select(new_cf_net_attack_rare_bcphc, select_items_bcphc)
    # convert to matrix
    cf_modulated_net_matrix_degree_bcphc <- as.matrix(cf_modulated_net_matrix_degree_bcphc)
    # save selected species
    selected_species_bcphc <- rbind(selected_species_bcphc, data.frame(i, n_nodes, h_net_subset_for_cf_rare_bcphc_df_nodes))
    
    # kernel-based comparison, part 1
    # create list of graph structures
    my_list1_bcphc <- list(h_rare_net_bcphc, cf_rare_net_bcphc, new_cf_rare_net_bcphc)
    # calculate random walk kernel
    my_K_bcphc_1 <- CalculateKStepRandomWalkKernel(my_list1_bcphc, rep(1,2))
    # convert to data frame
    my_K_bcphc_df <- data.frame(my_K_bcphc_1)
    # normalise to healthy structure
    my_K_bcphc_df_1 <- my_K_bcphc_df / my_K_bcphc_df[1,1]
    # rename columns
    colnames(my_K_bcphc_df_1) <- c("Healthy", "CF", "modulated CF")
    # rename rows
    rownames(my_K_bcphc_df_1) <- c("Healthy", "CF", "modulated CF")
    # keep first row
    my_K_bcphc_df_1 <- my_K_bcphc_df_1[1,]
    
    
    # kernel-based comparison, part 2
    # calculate shortest pathway kernel
    my_K2_bcphc <- CalculateShortestPathKernel(my_list1_bcphc)
    # convert to data frame
    my_K_bcphc_df2 <- data.frame(my_K2_bcphc)
    # normalise to healthy structure
    my_K_bcphc_df2 <- my_K_bcphc_df2 / my_K_bcphc_df2[1,1]
    # rename columns
    colnames(my_K_bcphc_df2) <- c("Healthy", "CF", "modulated CF")
    # rename rows
    rownames(my_K_bcphc_df2) <- c("Healthy", "CF", "modulated CF")
    # keep first row
    my_K_bcphc_df2 <- my_K_bcphc_df2[1,]
    
    
    # kernel-based comparison, part 3
    # calculate Weisfeiler-Lehmann kernel
    my_K3_bcphc <- CalculateWLKernel(my_list1_bcphc, 5)
    # convert to data frame
    my_K_bcphc_df2 <- data.frame(my_K3_bcphc)
    # normalise to healthy structure
    my_K_bcphc_df3 <- my_K_bcphc_df2 / 	my_K_bcphc_df2[1,1]
    # rename columns
    colnames(my_K_bcphc_df3) <- c("Healthy", "CF", "modulated CF")
    # rename rows
    rownames(my_K_bcphc_df3) <- c("Healthy", "CF", "modulated CF")
    # keep first row
    my_K_bcphc_df3 <- my_K_bcphc_df3[1,]
    
    # re-set row names
    nullmodelMod_bcphc_names <- lapply(nullmodelMod_bcphc, function(x) {rownames(x) <- rownames(new_net_bcphc_matrix); x })
    # re-set column names
    nullmodelMod_bcphc_names <- lapply(nullmodelMod_bcphc_names, function(x) {colnames(x) <- colnames(new_net_bcphc_matrix); x })
    # make list of graphs from list of adjacency matrices
    nullmodelMod_bcphc_graphs <- lapply(nullmodelMod_bcphc_names, graph_from_adjacency_matrix)
    # calculate Weisfeiler-Lehmann kernel
    nullmodelMod_bcphc_spk <- CalculateWLKernel(nullmodelMod_bcphc_graphs,5)
    # normalise to healthy structure
    nullmodelMod_bcphc_spk <- nullmodelMod_bcphc_spk / my_K_bcphc_df2[1,1]
    
    # get zscore and p value between original and random network structures
    bcphc_mod_zscore_spk <- net.metric.zscore(my_K_bcphc_df3[1,3], nullmodelMod_bcphc_spk)
    bcphc_mod_pvalue_spk <- net.metric.pvalue(bcphc_mod_zscore_spk)
    # make dataframe
    nullmodelMod_stats_spk_bcphc <- data.frame(cbind(bcphc_mod_zscore_spk, bcphc_mod_pvalue_spk))
    # add meta data 
    nullmodelMod_stats_spk_bcphc$normlisation <- "BCPHC"
    nullmodelMod_stats_spk_bcphc$seed <- i
    nullmodelMod_stats_spk_bcphc$n_nodes <- n_nodes
    nullmodelMod_stats_spk_bcphc$index <- "Weisfeiler_Lehmann"
    # re-name columns
    colnames(nullmodelMod_stats_spk_bcphc) <- c("zscore", "p-value", "normalisation", "seed", "n_nodes", "index")
    # remove row names
    rownames(nullmodelMod_stats_spk_bcphc) <- NULL
        
    # export kernel output and bind with empty global variable (defined above)
    netwerk_kernel_randomWalk_bcphc <- rbind(netwerk_kernel_randomWalk_bcphc, data.frame(i, n_nodes, my_K_bcphc_df))
    network_kernel_shortestPath_bcphc <- rbind(network_kernel_shortestPath_bcphc, data.frame(i, n_nodes, my_K_bcphc_df2))
    network_kernel_WL_bcphc <- rbind(network_kernel_WL_bcphc, data.frame(i, n_nodes, my_K_bcphc_df3))
    network_nullmodel_bcphc <- rbind(network_nullmodel_bcphc, data.frame(rbind(nullmodelMod_stats_extinction_bcphc, nullmodelMod_stats_spk_bcphc)))
    
    # Frechet distance between network attack curves
    # convert healthy attack curve to matrix
    h_rare_net_bcphc_matrix_degree_2 <- as.matrix(h_rare_net_bcphc_matrix_degree[,c(1,2)])
    # convert modulated attack curve to matrix
    cf_rare_net_bcphc_matrix_degree_2 <- as.matrix(cf_modulated_net_matrix_degree_bcphc[,c(1,2)])
    # calculate Frechet distance between both matrices
    frechet_deg_rare_bcphc <- Frechet(h_rare_net_bcphc_matrix_degree_2,cf_rare_net_bcphc_matrix_degree_2)
    # export and bind with empty global variable (defined above)
    frechet_attack_similarity_bcphc = rbind(frechet_attack_similarity_bcphc, data.frame(i, n_nodes, frechet_deg_rare_bcphc))
    
    # Obtain network statistics
    # degree
    degree_cf_hrare_bcphc <- sum(degree(new_cf_rare_net_bcphc))
    # number of nodes
    nNodes_cf_hrare_bcphc <- gorder(new_cf_rare_net_bcphc) 
    # number of edges
    nEdges_cf_hrare_bcphc <- gsize(new_cf_rare_net_bcphc) 
    # network diameter
    diameter_cf_hrare_bcphc <- diameter(new_cf_rare_net_bcphc) 
    # edge density
    density_cf_hrare_bcphc <- edge_density(new_cf_rare_net_bcphc) 
    # transitivity
    trans_cf_hrare_bcphc <- transitivity(new_cf_rare_net_bcphc) 
    # eigen centrality
    eigen_cf_hrare_bcphc <- centr_eigen(new_cf_rare_net_bcphc)$centralization 
    # merge network statistics information
    network_statistics_bcphc <- rbind(network_statistics_bcphc, data.frame(i, n_nodes, degree_cf_hrare_bcphc, nNodes_cf_hrare_bcphc, nEdges_cf_hrare_bcphc, 
                                                                           diameter_cf_hrare_bcphc,  density_cf_hrare_bcphc, 
                                                                           trans_cf_hrare_bcphc, eigen_cf_hrare_bcphc))
    
    
    # continue with RLE-normalised data
    ## randomly extract species (with replacement) from healthy network, RLE
    sample_h_rare_0_rle <- sample(V(h_rare_net_rle)$name, n_nodes, replace = sample_with_replacement)
    # generate a new graph structure from scratch with the extracted species
    h_net_subset_for_cf_rare_rle <- induced_subgraph(h_rare_net_rle, 
                                                     which(V(h_rare_net_rle)$name %in% sample_h_rare_0_rle),
                                                     impl = "create_from_scratch")
    # make new node list
    h_net_subset_for_cf_rare_rle_df_nodes <- as_data_frame(h_net_subset_for_cf_rare_rle, what ="vertices")
    # make new edge list
    h_net_subset_for_cf_rare_rle_df_edges <- as_data_frame(h_net_subset_for_cf_rare_rle, what ="edges")
    # convert node information from original CF graph into data frame
    cf_rare_net_rle_df_nodes <- as_data_frame(cf_rare_net_rle, what ="vertices")
    # convert edge information from original CF graph into data frame
    cf_rare_net_rle_df_edges <- as_data_frame(cf_rare_net_rle, what ="edges")
    # merge sub-sampled healthy network nodes with CF network nodes
    new_cf_rare_nodes_rle <- data.frame(rbind(h_net_subset_for_cf_rare_rle_df_nodes, cf_rare_net_rle_df_nodes))
    # keep only unique rows from a data frame
    new_cf_rare_nodes_rle <- distinct(new_cf_rare_nodes_rle)
    # merged sub-sampled healthy network edges with CF network edges
    new_cf_rare_edges_rle <- data.frame(rbind(h_net_subset_for_cf_rare_rle_df_edges, cf_rare_net_rle_df_edges))
    # keep only unique rows from a data frame
    new_cf_rare_edges_rle <- distinct(new_cf_rare_edges_rle)
    
    # make a new graph from dataframe of modulated CF structure
    new_cf_rare_net_rle <- graph_from_data_frame(new_cf_rare_edges_rle, directed=directed_network, vertices=new_cf_rare_nodes_rle)
    
    # null model
    new_net_rle_asdj <- as_adjacency_matrix(new_cf_rare_net_rle, type="both")
    new_net_rle_matrix <- as.matrix(new_net_rle_asdj)
    # Calculate network metric extinction slope
    originalMod_net_rle_extinction <- networklevel(new_net_rle_matrix, index = "extinction slope", extinctmethod='degree')
    originalMod_net_rle_extinction <- originalMod_net_rle_extinction[2]
    # make null model
    nullmodelMod_rle <- nullmodel(new_net_rle_matrix, N=nullmodel_networks_n, method=nullmodel_method)
    nullmodelMod_rle_extinction <- unlist(sapply(nullmodelMod_rle, networklevel, index = 'extinction slope', extinctmethod='degree'))
    nullmodelMod_rle_extinction <- nullmodelMod_rle_extinction[2,]
    # get zscore and p value between original and random network structures
    rle_mod_zscore <- net.metric.zscore(originalMod_net_rle_extinction, nullmodelMod_rle_extinction)
    rle_mod_pvalue <- net.metric.pvalue(rle_mod_zscore)
    # generate final data frame
    nullmodelMod_stats_extinction_rle <- data.frame(cbind(rle_mod_zscore, rle_mod_pvalue))
    nullmodelMod_stats_extinction_rle$normlisation <- "RLE"
    nullmodelMod_stats_extinction_rle$seed <- i
    nullmodelMod_stats_extinction_rle$n_nodes <- n_nodes
    nullmodelMod_stats_extinction_rle$index <- "Robustness"
    colnames(nullmodelMod_stats_extinction_rle) <- c("zscore", "p-value", "normalisation", "seed", "n_nodes", "index")
    rownames(nullmodelMod_stats_extinction_rle) <- NULL
    
    
    # assesses networks vulnerability due to targeted attacks with 10 iterations for assessing random error
    new_cf_net_attack_rare_rle <- swan_combinatory(new_cf_rare_net_rle, 10)
    # convert to data frame
    new_cf_net_attack_rare_rle <- data.frame(new_cf_net_attack_rare_rle)
    # rename column names
    colnames(new_cf_net_attack_rare_rle) <- c("fraction_removed", "betweeness", "degree", "cascading", "random")
    # remove duplicate rows and take median value of duplicates
    new_cf_net_attack_rare_rle <- ddply(new_cf_net_attack_rare_rle,"fraction_removed",numcolwise(median))
    # insert first row with only zero values to have zero as starting point
    new_cf_net_attack_rare_rle <- rbind(h_rare_net_rle_attack[1,], new_cf_net_attack_rare_rle)
    # make dataframe with sub-selection options
    select_items_rle <- c("fraction_removed", "degree")
    # select items from attack dataframe
    cf_modulated_net_matrix_degree_rle <- select(new_cf_net_attack_rare_rle, select_items_rle)
    # convert to matrix
    cf_modulated_net_matrix_degree_rle <- as.matrix(cf_modulated_net_matrix_degree_rle)
    # save selected species
    selected_species_rle <- rbind(selected_species_rle, data.frame(i, n_nodes, h_net_subset_for_cf_rare_rle_df_nodes))
    
    # kernel-based comparison, part 1
    # create list of graph structures
    my_list1_rle <- list(h_rare_net_rle, cf_rare_net_rle, new_cf_rare_net_rle)
    # calculate random walk kernel
    my_K_rle_1 <- CalculateKStepRandomWalkKernel(my_list1_rle, rep(1,2))
    # convert to data frame
    my_K_rle_df <- data.frame(my_K_rle_1)
    # normalise to healthy structure
    my_K_rle_df_1 <- my_K_rle_df / my_K_rle_df[1,1]
    # rename columns
    colnames(my_K_rle_df_1) <- c("Healthy", "CF", "modulated CF")
    # rename rows
    rownames(my_K_rle_df_1) <- c("Healthy", "CF", "modulated CF")
    # keep first row
    my_K_rle_df_1 <- my_K_rle_df_1[1,]
    
    
    # kernel-based comparison, part 2
    # calculate shortest pathway kernel
    my_K2_rle <- CalculateShortestPathKernel(my_list1_rle)
    # convert to data frame
    my_K_rle_df2 <- data.frame(my_K2_rle)
    # normalise to healthy structure
    my_K_rle_df2 <- my_K_rle_df2 / my_K_rle_df2[1,1]
    # rename columns
    colnames(my_K_rle_df2) <- c("Healthy", "CF", "modulated CF")
    # rename rows
    rownames(my_K_rle_df2) <- c("Healthy", "CF", "modulated CF")
    # keep first row
    my_K_rle_df2 <- my_K_rle_df2[1,]
    
    
    # kernel-based comparison, part 3
    # calculate Weisfeiler-Lehmann kernel
    my_K3_rle <- CalculateWLKernel(my_list1_rle, 5)
    # convert to data frame
    my_K_rle_df2 <- data.frame(my_K3_rle)
    # normalise to healthy structure
    my_K_rle_df3 <- my_K_rle_df2 / 	my_K_rle_df2[1,1]
    # rename columns
    colnames(my_K_rle_df3) <- c("Healthy", "CF", "modulated CF")
    # rename rows
    rownames(my_K_rle_df3) <- c("Healthy", "CF", "modulated CF")
    # keep first row
    my_K_rle_df3 <- my_K_rle_df3[1,]
    
    # re-set row names
    nullmodelMod_rle_names <- lapply(nullmodelMod_rle, function(x) {rownames(x) <- rownames(new_net_rle_matrix); x })
    # re-set column names
    nullmodelMod_rle_names <- lapply(nullmodelMod_rle_names, function(x) {colnames(x) <- colnames(new_net_rle_matrix); x })
    # make list of graphs from list of adjacency matrices
    nullmodelMod_rle_graphs <- lapply(nullmodelMod_rle_names, graph_from_adjacency_matrix)
    # calculate Weisfeiler-Lehmann kernel
    nullmodelMod_rle_spk <- CalculateWLKernel(nullmodelMod_rle_graphs,5)
    # normalise to healthy structure
    nullmodelMod_rle_spk <- nullmodelMod_rle_spk / my_K_rle_df2[1,1]
    
    # get zscore and p value between original and random network structures
    rle_mod_zscore_spk <- net.metric.zscore(my_K_rle_df3[1,3], nullmodelMod_rle_spk)
    rle_mod_pvalue_spk <- net.metric.pvalue(rle_mod_zscore_spk)
    # make dataframe
    nullmodelMod_stats_spk_rle <- data.frame(cbind(rle_mod_zscore_spk, rle_mod_pvalue_spk))
    # add meta data 
    nullmodelMod_stats_spk_rle$normlisation <- "RLE"
    nullmodelMod_stats_spk_rle$seed <- i
    nullmodelMod_stats_spk_rle$n_nodes <- n_nodes
    nullmodelMod_stats_spk_rle$index <- "Weisfeiler_Lehmann"
    # re-name columns
    colnames(nullmodelMod_stats_spk_rle) <- c("zscore", "p-value", "normalisation", "seed", "n_nodes", "index")
    # remove row names
    rownames(nullmodelMod_stats_spk_rle) <- NULL
    
    # export kernel output and bind with empty global variable (defined above)
    netwerk_kernel_randomWalk_rle <- rbind(netwerk_kernel_randomWalk_rle, data.frame(i, n_nodes, my_K_rle_df))
    network_kernel_shortestPath_rle <- rbind(network_kernel_shortestPath_rle, data.frame(i, n_nodes, my_K_rle_df2))
    network_kernel_WL_rle <- rbind(network_kernel_WL_rle, data.frame(i, n_nodes, my_K_rle_df3))
    network_nullmodel_rle <- rbind(network_nullmodel_rle, data.frame(rbind(nullmodelMod_stats_extinction_rle, nullmodelMod_stats_spk_rle)))
    
    # Frechet distance between network attack curves
    # convert healthy attack curve to matrix
    h_rare_net_rle_matrix_degree_2 <- as.matrix(h_rare_net_rle_matrix_degree[,c(1,2)])
    # convert modulated attack curve to matrix
    cf_rare_net_rle_matrix_degree_2 <- as.matrix(cf_modulated_net_matrix_degree_rle[,c(1,2)])
    # calculate Frechet distance between both matrices
    frechet_deg_rare_rle <- Frechet(h_rare_net_rle_matrix_degree_2,cf_rare_net_rle_matrix_degree_2)
    # export and bind with empty global variable (defined above)
    frechet_attack_similarity_rle = rbind(frechet_attack_similarity_rle, data.frame(i, n_nodes, frechet_deg_rare_rle))
    
    # Obtain network statistics
    # degree
    degree_cf_hrare_rle <- sum(degree(new_cf_rare_net_rle))
    # number of nodes
    nNodes_cf_hrare_rle <- gorder(new_cf_rare_net_rle) 
    # number of edges
    nEdges_cf_hrare_rle <- gsize(new_cf_rare_net_rle) 
    # network diameter
    diameter_cf_hrare_rle <- diameter(new_cf_rare_net_rle) 
    # edge density
    density_cf_hrare_rle <- edge_density(new_cf_rare_net_rle) 
    # transitivity
    trans_cf_hrare_rle <- transitivity(new_cf_rare_net_rle) 
    # eigen centrality
    eigen_cf_hrare_rle <- centr_eigen(new_cf_rare_net_rle)$centralization 
    # merge network statistics information
    network_statistics_rle <- rbind(network_statistics_rle, data.frame(i, n_nodes, degree_cf_hrare_rle,  nNodes_cf_hrare_rle,
                                                                       nEdges_cf_hrare_rle, diameter_cf_hrare_rle, 
                                                                       density_cf_hrare_rle, trans_cf_hrare_rle, eigen_cf_hrare_rle))
    
    
    # continue with VST-normalised data
    # randomly extract species (with replacement) from healthy network, BCPHC
    sample_h_rare_0_vst <- sample(V(h_rare_net_vst)$name, n_nodes, replace = sample_with_replacement)
    # generate a new graph structure from scratch with the extracted species
    h_net_subset_for_cf_rare_vst <- induced_subgraph(h_rare_net_vst, 
                                                     which(V(h_rare_net_vst)$name %in% sample_h_rare_0_vst),
                                                     impl = "create_from_scratch")
    # make new node list
    h_net_subset_for_cf_rare_vst_df_nodes <- as_data_frame(h_net_subset_for_cf_rare_vst, what ="vertices")
    # make new edge list
    h_net_subset_for_cf_rare_vst_df_edges <- as_data_frame(h_net_subset_for_cf_rare_vst, what ="edges")
    # convert node information from original CF graph into data frame
    cf_rare_net_vst_df_nodes <- as_data_frame(cf_rare_net_vst, what ="vertices")
    # convert edge information from original CF graph into data frame
    cf_rare_net_vst_df_edges <- as_data_frame(cf_rare_net_vst, what ="edges")
    # merge sub-sampled healthy network nodes with CF network nodes
    new_cf_rare_nodes_vst <- data.frame(rbind(h_net_subset_for_cf_rare_vst_df_nodes, cf_rare_net_vst_df_nodes))
    # keep only unique rows from a data frame
    new_cf_rare_nodes_vst <- distinct(new_cf_rare_nodes_vst)
    # merged sub-sampled healthy network edges with CF network edges
    new_cf_rare_edges_vst <- data.frame(rbind(h_net_subset_for_cf_rare_vst_df_edges, cf_rare_net_vst_df_edges))
    # keep only unique rows from a data frame
    new_cf_rare_edges_vst <- distinct(new_cf_rare_edges_vst)
    
    # make a new graph from dataframe of modulated CF structure
    new_cf_rare_net_vst <- graph_from_data_frame(new_cf_rare_edges_vst, directed=directed_network, vertices=new_cf_rare_nodes_vst)
    
    # null model
    new_net_vst_asdj <- as_adjacency_matrix(new_cf_rare_net_vst, type="both")
    new_net_vst_matrix <- as.matrix(new_net_vst_asdj)
    # Calculate network metric extinction slope
    originalMod_net_vst_extinction <- networklevel(new_net_vst_matrix, index = "extinction slope", extinctmethod='degree')
    originalMod_net_vst_extinction <- originalMod_net_vst_extinction[2]
    # make null model
    nullmodelMod_vst <- nullmodel(new_net_vst_matrix, N=nullmodel_networks_n, method=nullmodel_method)
    nullmodelMod_vst_extinction <- unlist(sapply(nullmodelMod_vst, networklevel, index = 'extinction slope', extinctmethod='degree')) 
    nullmodelMod_vst_extinction <- nullmodelMod_vst_extinction[2,]
    # get zscore and p value between original and random network structures
    vst_mod_zscore <- net.metric.zscore(originalMod_net_vst_extinction, nullmodelMod_vst_extinction)
    vst_mod_pvalue <- net.metric.pvalue(vst_mod_zscore)
    # generate final data frame
    nullmodelMod_stats_extinction_vst <- data.frame(cbind(vst_mod_zscore, vst_mod_pvalue))
    nullmodelMod_stats_extinction_vst$normlisation <- "VST"
    nullmodelMod_stats_extinction_vst$seed <- i
    nullmodelMod_stats_extinction_vst$n_nodes <- n_nodes
    nullmodelMod_stats_extinction_vst$index <- "Robustness"
    colnames(nullmodelMod_stats_extinction_vst) <- c("zscore", "p-value", "normalisation", "seed", "n_nodes", "index")
    rownames(nullmodelMod_stats_extinction_vst) <- NULL
    
    
    # assesses networks vulnerability due to targeted attacks with 10 iterations for assessing random error
    new_cf_net_attack_rare_vst <- swan_combinatory(new_cf_rare_net_vst,10)
    # convert to data frame
    new_cf_net_attack_rare_vst <- data.frame(new_cf_net_attack_rare_vst)
    # rename columns
    colnames(new_cf_net_attack_rare_vst) <- c("fraction_removed", "betweeness", "degree", "cascading", "random")
    # remove duplicate rows and take median value of duplicates
    new_cf_net_attack_rare_vst <- ddply(new_cf_net_attack_rare_vst,"fraction_removed",numcolwise(median))
    # insert first row with only zero values to have zero as starting point
    new_cf_net_attack_rare_vst <- rbind(h_rare_net_vst_attack[1,], new_cf_net_attack_rare_vst)
    # make dataframe with sub-selection options
    select_items_vst <- c("fraction_removed", "degree")
    # select items from attack dataframe
    cf_modulated_net_matrix_degree_vst <- select(new_cf_net_attack_rare_vst, select_items_vst)
    # convert to matrix
    cf_modulated_net_matrix_degree_vst <- as.matrix(cf_modulated_net_matrix_degree_vst)
    
    # save selected species
    selected_species_vst <- rbind(selected_species_vst, data.frame(i, n_nodes, h_net_subset_for_cf_rare_vst_df_nodes))
    
    # kernel-based comparison, part 1
    # create list of graph structures
    my_list1_vst <- list(h_rare_net_vst, cf_rare_net_vst, new_cf_rare_net_vst)
    # calculate random walk kernel
    my_K_vst_1 <- CalculateKStepRandomWalkKernel(my_list1_vst, rep(1,2))
    # convert to data frame
    my_K_vst_df <- data.frame(my_K_vst_1)
    # normalise to healthy structure
    my_K_vst_df_1 <- my_K_vst_df / my_K_vst_df[1,1]
    # rename columns
    colnames(my_K_vst_df_1) <- c("Healthy", "CF", "modulated CF")
    # rename rows
    rownames(my_K_vst_df_1) <- c("Healthy", "CF", "modulated CF")
    # keep first row
    my_K_vst_df_1 <- my_K_vst_df_1[1,]
    
    
    # kernel-based comparison, part 2
    # calculate shortest path kernel
    my_K2_vst <- CalculateShortestPathKernel(my_list1_vst)
    # convert to data frame
    my_K_vst_df2 <- data.frame(my_K2_vst)
    # normalise to healthy structure
    my_K_vst_df2 <- my_K_vst_df2 / my_K_vst_df2[1,1]
    # rename columns
    colnames(my_K_vst_df2) <- c("Healthy", "CF", "modulated CF")
    # rename rows
    rownames(my_K_vst_df2) <- c("Healthy", "CF", "modulated CF")
    # keep first row
    my_K_vst_df2 <- my_K_vst_df2[1,]
    
    
    # kernel-based comparison, part 3
    # calculate Weisfeiler-Lehmann kernel
    my_K3_vst <- CalculateWLKernel(my_list1_vst, 5)
    # convert to data frame
    my_K_vst_df2 <- data.frame(my_K3_vst)
    # normalise to healthy structure
    my_K_vst_df3 <- my_K_vst_df2 / 	my_K_vst_df2[1,1]
    # rename columns
    colnames(my_K_vst_df3) <- c("Healthy", "CF", "modulated CF")
    # rename rows
    rownames(my_K_vst_df3) <- c("Healthy", "CF", "modulated CF")
    # keep first row
    my_K_vst_df3 <- my_K_vst_df3[1,]
    
    # re-set row names
    nullmodelMod_vst_names <- lapply(nullmodelMod_vst, function(x) {rownames(x) <- rownames(new_net_vst_matrix); x })
    # re-set column names
    nullmodelMod_vst_names <- lapply(nullmodelMod_vst_names, function(x) {colnames(x) <- colnames(new_net_vst_matrix); x })
    # make list of graphs from list of adjacency matrices
    nullmodelMod_vst_graphs <- lapply(nullmodelMod_vst_names, graph_from_adjacency_matrix)
    # calculate Weisfeiler-Lehmann kernel
    nullmodelMod_vst_spk <- CalculateWLKernel(nullmodelMod_vst_graphs,5)
    # normalise to healthy structure
    nullmodelMod_vst_spk <- nullmodelMod_vst_spk / my_K_vst_df2[1,1]
    
    # get zscore and p value between original and random network structures
    vst_mod_zscore_spk <- net.metric.zscore(my_K_vst_df3[1,3], nullmodelMod_vst_spk)
    vst_mod_pvalue_spk <- net.metric.pvalue(vst_mod_zscore_spk)
    # make dataframe
    nullmodelMod_stats_spk_vst <- data.frame(cbind(vst_mod_zscore_spk, vst_mod_pvalue_spk))
    # add meta data 
    nullmodelMod_stats_spk_vst$normlisation <- "VST"
    nullmodelMod_stats_spk_vst$seed <- i
    nullmodelMod_stats_spk_vst$n_nodes <- n_nodes
    nullmodelMod_stats_spk_vst$index <- "Weisfeiler_Lehmann"
    # re-name columns
    colnames(nullmodelMod_stats_spk_vst) <- c("zscore", "p-value", "normalisation", "seed", "n_nodes", "index")
    # remove row names
    rownames(nullmodelMod_stats_spk_vst) <- NULL
    
    # export kernel output and bind with empty global variable (defined above)
    netwerk_kernel_randomWalk_vst <- rbind(netwerk_kernel_randomWalk_vst, data.frame(i, n_nodes, my_K_vst_df))
    network_kernel_shortestPath_vst <- rbind(network_kernel_shortestPath_vst, data.frame(i, n_nodes, my_K_vst_df2))
    network_kernel_WL_vst <- rbind(network_kernel_WL_vst, data.frame(i, n_nodes, my_K_vst_df3))
    network_nullmodel_vst <- rbind(network_nullmodel_vst, data.frame(rbind(nullmodelMod_stats_extinction_vst, nullmodelMod_stats_spk_vst)))
    
    # Frechet distance between network attack curves
    # convert healthy attack curve to matrix
    h_rare_net_vst_matrix_degree_2 <- as.matrix(h_rare_net_vst_matrix_degree[,c(1:2)])
    # convert modulated attack curve to matrix
    cf_rare_net_vst_matrix_degree_2 <- as.matrix(cf_modulated_net_matrix_degree_vst[,c(1:2)])
    # calculate Frechet distance between both matrices
    frechet_deg_rare_vst <- Frechet(h_rare_net_vst_matrix_degree_2,cf_rare_net_vst_matrix_degree_2)
    # export and bind with empty global variable (defined above)
    frechet_attack_similarity_vst = rbind(frechet_attack_similarity_vst, data.frame(i, n_nodes, frechet_deg_rare_vst))
    
    # Obtain network statistics
    # degree
    degree_cf_hrare_vst <- sum(degree(new_cf_rare_net_vst))
    # number of nodes
    nNodes_cf_hrare_vst <- gorder(new_cf_rare_net_vst) 
    # number of edges
    nEdges_cf_hrare_vst <- gsize(new_cf_rare_net_vst) 
    # network diameter
    diameter_cf_hrare_vst <- diameter(new_cf_rare_net_vst) 
    # edge density
    density_cf_hrare_vst <- edge_density(new_cf_rare_net_vst) 
    # transitivity
    trans_cf_hrare_vst <- transitivity(new_cf_rare_net_vst) 
    # eigen centrality
    eigen_cf_hrare_vst <- centr_eigen(new_cf_rare_net_vst)$centralization 
    # merge network statistics information and export to global variable
    network_statistics_vst <- rbind(network_statistics_vst, data.frame(i, n_nodes, degree_cf_hrare_vst, nNodes_cf_hrare_vst, nEdges_cf_hrare_vst, 
                                                                       diameter_cf_hrare_vst, density_cf_hrare_vst, trans_cf_hrare_vst, eigen_cf_hrare_vst))    
    
  }
}


############################################################################################################
# Evaluation of attack curve simulations and Frechet distances
# BCPHC-normalised data
# calculate median Frechet distance per number of nodes that was inserted into CF network
frechet_attack_similarity_bcphc_median <- ddply(frechet_attack_similarity_bcphc,"n_nodes",numcolwise(median))
# rename columns
colnames(frechet_attack_similarity_bcphc_median) <- c("n_nodes_med", "i_med","cf_hrare_deg_frechet_med")
# add a column with Frechet distance between original CF network and original healthy network
frechet_attack_similarity_bcphc_median$cf_vs_healthy_degree <- original_frechet_h_cf_degree_bcphc
# calculate minimum Frechet distance per number of nodes that was inserted into CF network
frechet_attack_similarity_bcphc_minimum <- ddply(frechet_attack_similarity_bcphc,"n_nodes",numcolwise(min))
# rename columns
colnames(frechet_attack_similarity_bcphc_minimum) <- c("n_nodes_min", "i_min","cf_hrare_deg_frechet_min")
# calculate maximum Frechet distance per number of nods that was inserted into CF network
frechet_attack_similarity_bcphc_maximum <- ddply(frechet_attack_similarity_bcphc,"n_nodes",numcolwise(max))
# rename columns
colnames(frechet_attack_similarity_bcphc_maximum) <- c("n_nodes_max", "i_max","cf_hrare_deg_frechet_max")
# bind information into one data frame
frechet_attack_similartiy_degree <- data.frame(cbind(frechet_attack_similarity_bcphc_median, 
                                                     frechet_attack_similarity_bcphc_minimum, 
                                                     frechet_attack_similarity_bcphc_maximum))


# RLE-normalised data
# calculate median Frechet distance per number of nodes that was inserted into CF network
frechet_attack_similarity_rle_median <- ddply(frechet_attack_similarity_rle,"n_nodes",numcolwise(median))
# rename columns
colnames(frechet_attack_similarity_rle_median) <- c("n_nodes_med", "i_med","cf_hrare_deg_frechet_med")
# add a column with Frechet distance between original CF network and original healthy network
frechet_attack_similarity_rle_median$cf_vs_healthy_degree <- original_frechet_h_cf_degree_rle
# calculate minimum Frechet distance per number of nodes that was inserted into CF network
frechet_attack_similarity_rle_minimum <- ddply(frechet_attack_similarity_rle,"n_nodes",numcolwise(min))
# rename columns
colnames(frechet_attack_similarity_rle_minimum) <- c("n_nodes_min", "i_min","cf_hrare_deg_frechet_min")
# calculate maximum Frechet distance per number of nodes that was inserted into CF network
frechet_attack_similarity_rle_maximum <- ddply(frechet_attack_similarity_rle,"n_nodes",numcolwise(max))
# rename columns
colnames(frechet_attack_similarity_rle_maximum) <- c("n_nodes_max", "i_max","cf_hrare_deg_frechet_max")
# bind information into one data frame
frechet_attack_similartiy_degree_rle <- data.frame(cbind(frechet_attack_similarity_rle_median, 
                                                         frechet_attack_similarity_rle_minimum, 
                                                         frechet_attack_similarity_rle_maximum))


# VST-normalised data
# calculate median Frechet distance per number of nodes that was inserted into CF network
frechet_attack_similarity_vst_median <- ddply(frechet_attack_similarity_vst,"n_nodes",numcolwise(median))
# rename columns
colnames(frechet_attack_similarity_vst_median) <- c("n_nodes_med", "i_med","cf_hrare_deg_frechet_med")
# add a column with Frechet distance between original CF network and original healthy network
frechet_attack_similarity_vst_median$cf_vs_healthy_degree <- original_frechet_h_cf_degree_vst
# calculate minimum Frechet distance per number of nodes that was inserted into CF network
frechet_attack_similarity_vst_minimum <- ddply(frechet_attack_similarity_vst,"n_nodes",numcolwise(min))
# rename columns
colnames(frechet_attack_similarity_vst_minimum) <- c("n_nodes_min", "i_min","cf_hrare_deg_frechet_min")
# calculate maximum Frechet distance per number of nodes that was inserted into CF network
frechet_attack_similarity_vst_maximum <- ddply(frechet_attack_similarity_vst,"n_nodes",numcolwise(max))
# rename columns
colnames(frechet_attack_similarity_vst_maximum) <- c("n_nodes_max", "i_max","cf_hrare_deg_frechet_max")
# bind information into one data frame
frechet_attack_similartiy_degree_vst <- data.frame(cbind(frechet_attack_similarity_vst_median,
                                                         frechet_attack_similarity_vst_minimum, 
                                                         frechet_attack_similarity_vst_maximum))


# rename column names for all normalisation strategies so that they are matching
colnames(frechet_attack_similarity_vst) <- c("i", "n_nodes", "frechet")
colnames(frechet_attack_similarity_rle) <- c("i", "n_nodes", "frechet")
colnames(frechet_attack_similarity_bcphc) <- c("i", "n_nodes", "frechet")
# add meta data
frechet_attack_similarity_vst$normalisation <- "VST"
frechet_attack_similarity_bcphc$normalisation <- "BCPHC"
frechet_attack_similarity_rle$normalisation <- "RLE"
# bind all data frames into one
frechet_attack_similarity <- data.frame(rbind(frechet_attack_similarity_vst, 
                                              frechet_attack_similarity_rle, 
                                              frechet_attack_similarity_bcphc))
# add meta data
frechet_attack_similartiy_degree_vst$normalisation <- "VST"
frechet_attack_similartiy_degree$normalisation <- "BCPHC"
frechet_attack_similartiy_degree_rle$normalisation <- "RLE"
# bind data frames into one
frechet_attack_df <- data.frame(rbind(frechet_attack_similartiy_degree_vst, 
                                      frechet_attack_similartiy_degree, 
                                      frechet_attack_similartiy_degree_rle))

############################################################################################################
# Kernel-based evaluation of simulation runs
# Shortest pathway kernel, BCPHC
# calculate the median shortest path kernel per number of nodes that was inserted into CF network
kernel_shortestPahway_median_bcphc <- ddply(network_kernel_shortestPath_bcphc,"n_nodes",numcolwise(median))
# rename columns
colnames(kernel_shortestPahway_median_bcphc) <- c("n_nodes_med", "i_med","Healthy_med", "CF_med", "CF_modulated_med")
# calculate the minimum shortest path kernel per number of nodes that was inserted into CF network
kernel_shortestPahway_min_bcphc <- ddply(network_kernel_shortestPath_bcphc,"n_nodes",numcolwise(min))
# rename columns
colnames(kernel_shortestPahway_min_bcphc) <- c("n_nodes_min", "i_med","Healthy_min", "CF_min", "CF_modulated_min")
# calculate the maximum shortest path kernel per number of nodes that was inserted into CF network
kernel_shortestPahway_max_bcphc <- ddply(network_kernel_shortestPath_bcphc,"n_nodes",numcolwise(max))
# rename columns
colnames(kernel_shortestPahway_max_bcphc) <- c("n_nodes_max", "i_med","Healthy_max", "CF_max", "CF_modulated_max")
# merge median, minimum and maximum values into one data frame
kernel_shortestPahway_df_bcphc <- data.frame(cbind(kernel_shortestPahway_median_bcphc, kernel_shortestPahway_min_bcphc,kernel_shortestPahway_max_bcphc))
# get median value relative to healthy median value
kernel_shortestPahway_df_bcphc$CF_modulated_med_2 <- kernel_shortestPahway_df_bcphc$CF_modulated_med / kernel_shortestPahway_df_bcphc$Healthy_med
# get maximum value relative to healthy maximum value
kernel_shortestPahway_df_bcphc$CF_modulated_max_2 <- kernel_shortestPahway_df_bcphc$CF_modulated_max / kernel_shortestPahway_df_bcphc$Healthy_max
# get minimum value relative to healthy minimum value
kernel_shortestPahway_df_bcphc$CF_modulated_min_2 <- kernel_shortestPahway_df_bcphc$CF_modulated_min / kernel_shortestPahway_df_bcphc$Healthy_min
# get original CF value relative to healthy value
kernel_shortestPahway_df_bcphc$CF_med_2 <- kernel_shortestPahway_df_bcphc$CF_med / kernel_shortestPahway_df_bcphc$Healthy_med


# Continue shortest path kernel evaluations with RLE-normalised data
# calculate the median shortest path kernel per number of nodes that was inserted into CF network
kernel_shortestPahway_median_rle <- ddply(network_kernel_shortestPath_rle,"n_nodes",numcolwise(median))
# rename columns
colnames(kernel_shortestPahway_median_rle) <- c("n_nodes_med", "i_med","Healthy_med", "CF_med", "CF_modulated_med")
# calculate the minimum shortest path kernel per number of nodes that was inserted into CF network
kernel_shortestPahway_min_rle <- ddply(network_kernel_shortestPath_rle,"n_nodes",numcolwise(min))
# rename columns
colnames(kernel_shortestPahway_min_rle) <- c("n_nodes_min", "i_med","Healthy_min", "CF_min", "CF_modulated_min")
# calculate the maximum shortest path kernel per number of nodes that was inserted into CF network
kernel_shortestPahway_max_rle <- ddply(network_kernel_shortestPath_rle,"n_nodes",numcolwise(max))
# rename columns
colnames(kernel_shortestPahway_max_rle) <- c("n_nodes_max", "i_med","Healthy_max", "CF_max", "CF_modulated_max")
# merge median, minimum and maximum values into one data frame
kernel_shortestPahway_df_rle <- data.frame(cbind(kernel_shortestPahway_median_rle, kernel_shortestPahway_min_rle,kernel_shortestPahway_max_rle))
# get median value relative to healthy median value
kernel_shortestPahway_df_rle$CF_modulated_med_2 <- kernel_shortestPahway_df_rle$CF_modulated_med / kernel_shortestPahway_df_rle$Healthy_med
# get maximum value relative to healthy maximum value
kernel_shortestPahway_df_rle$CF_modulated_max_2 <- kernel_shortestPahway_df_rle$CF_modulated_max / kernel_shortestPahway_df_rle$Healthy_max
# get minimum value relative to healthy minimum value
kernel_shortestPahway_df_rle$CF_modulated_min_2 <- kernel_shortestPahway_df_rle$CF_modulated_min / kernel_shortestPahway_df_rle$Healthy_min
# get original CF value relative to healthy value
kernel_shortestPahway_df_rle$CF_med_2 <- kernel_shortestPahway_df_rle$CF_med / kernel_shortestPahway_df_rle$Healthy_med


# Continue shortest path kernel evaluations with VST-normalised data
# calculate the median shortest path kernel per number of nodes that was inserted into CF network
kernel_shortestPahway_median_vst <- ddply(network_kernel_shortestPath_vst,"n_nodes",numcolwise(median))
# rename columns
colnames(kernel_shortestPahway_median_vst) <- c("n_nodes_med", "i_med","Healthy_med", "CF_med", "CF_modulated_med")
# calculate the minimum shortest path kernel per number of nodes that was inserted into CF network
kernel_shortestPahway_min_vst <- ddply(network_kernel_shortestPath_vst,"n_nodes",numcolwise(min))
# rename columns
colnames(kernel_shortestPahway_min_vst) <- c("n_nodes_min", "i_med","Healthy_min", "CF_min", "CF_modulated_min")
# calculate the maximum shortest path kernel per number of nodes that was inserted into CF network
kernel_shortestPahway_max_vst <- ddply(network_kernel_shortestPath_vst,"n_nodes",numcolwise(max))
# rename columns
colnames(kernel_shortestPahway_max_vst) <- c("n_nodes_max", "i_med","Healthy_max", "CF_max", "CF_modulated_max")
# merge median, minimum and maximum values into one data frame
kernel_shortestPahway_df_vst <- data.frame(cbind(kernel_shortestPahway_median_vst, kernel_shortestPahway_min_vst,kernel_shortestPahway_max_vst))
# get median value relative to healthy median value
kernel_shortestPahway_df_vst$CF_modulated_med_2 <- kernel_shortestPahway_df_vst$CF_modulated_med / kernel_shortestPahway_df_vst$Healthy_med
# get maximum value relative to healthy maximum value
kernel_shortestPahway_df_vst$CF_modulated_max_2 <- kernel_shortestPahway_df_vst$CF_modulated_max / kernel_shortestPahway_df_vst$Healthy_max
# get minimum value relative to healthy minimum value
kernel_shortestPahway_df_vst$CF_modulated_min_2 <- kernel_shortestPahway_df_vst$CF_modulated_min / kernel_shortestPahway_df_vst$Healthy_min
# get original CF value relative to healthy value
kernel_shortestPahway_df_vst$CF_med_2 <- kernel_shortestPahway_df_vst$CF_med / kernel_shortestPahway_df_vst$Healthy_med


# Weisfeiler-Lehman Graph Kernels, BCPHC
# calculate the median shortest path kernel per number of nodes that was inserted into CF network
kernel_WL_median_bcphc <- ddply(network_kernel_WL_bcphc,"n_nodes",numcolwise(median))
# rename columns
colnames(kernel_WL_median_bcphc) <- c("n_nodes_med", "i_med","Healthy_med", "CF_med", "CF_modulated_med")
# calculate the minimum shortest path kernel per number of nodes that was inserted into CF network
kernel_WL_min_bcphc <- ddply(network_kernel_WL_bcphc,"n_nodes",numcolwise(min))
# rename columns
colnames(kernel_WL_min_bcphc) <- c("n_nodes_min", "i_med","Healthy_min", "CF_min", "CF_modulated_min")
# calculate the maximum shortest path kernel per number of nodes that was inserted into CF network
kernel_WL_max_bcphc <- ddply(network_kernel_WL_bcphc,"n_nodes",numcolwise(max))
# rename columns
colnames(kernel_WL_max_bcphc) <- c("n_nodes_max", "i_med","Healthy_max", "CF_max", "CF_modulated_max")
# merge median, minimum and maximum values into one data frame
kernel_WL_df_bcphc <- data.frame(cbind(kernel_WL_median_bcphc, kernel_WL_min_bcphc, kernel_WL_max_bcphc))
# get median value relative to healthy median value
kernel_WL_df_bcphc$CF_modulated_med_2 <- kernel_WL_df_bcphc$CF_modulated_med / kernel_WL_df_bcphc$Healthy_med
# get maximum value relative to healthy maximum value
kernel_WL_df_bcphc$CF_modulated_max_2 <- kernel_WL_df_bcphc$CF_modulated_max / kernel_WL_df_bcphc$Healthy_max
# get minimum value relative to healthy minimum value
kernel_WL_df_bcphc$CF_modulated_min_2 <- kernel_WL_df_bcphc$CF_modulated_min / kernel_WL_df_bcphc$Healthy_min
# get original CF value relative to healthy value
kernel_WL_df_bcphc$CF_med_2 <- kernel_WL_df_bcphc$CF_med / kernel_WL_df_bcphc$Healthy_med


# Continue Weisfeiler-Lehman Graph Kernel evaluations with RLE-normalised data
# calculate the median shortest path kernel per number of nodes that was inserted into CF network
kernel_WL_median_rle <- ddply(network_kernel_WL_rle,"n_nodes",numcolwise(median))
# rename columns
colnames(kernel_WL_median_rle) <- c("n_nodes_med", "i_med","Healthy_med", "CF_med", "CF_modulated_med")
# calculate the minimum shortest path kernel per number of nodes that was inserted into CF network
kernel_WL_min_rle <- ddply(network_kernel_WL_rle,"n_nodes",numcolwise(min))
# rename columns
colnames(kernel_WL_min_rle) <- c("n_nodes_min", "i_med","Healthy_min", "CF_min", "CF_modulated_min")
# calculate the maximum shortest path kernel per number of nodes that was inserted into CF network
kernel_WL_max_rle <- ddply(network_kernel_WL_rle,"n_nodes",numcolwise(max))#
# rename columns
colnames(kernel_WL_max_rle) <- c("n_nodes_max", "i_med","Healthy_max", "CF_max", "CF_modulated_max")
# merge median, minimum and maximum values into one data frame
kernel_WL_df_rle <- data.frame(cbind(kernel_WL_median_rle, kernel_WL_min_rle, kernel_WL_max_rle))
# get median value relative to healthy median value
kernel_WL_df_rle$CF_modulated_med_2 <- kernel_WL_df_rle$CF_modulated_med / kernel_WL_df_rle$Healthy_med
# get maximum value relative to healthy maximum value
kernel_WL_df_rle$CF_modulated_max_2 <- kernel_WL_df_rle$CF_modulated_max / kernel_WL_df_rle$Healthy_max
# get minimum value relative to healthy minimum value
kernel_WL_df_rle$CF_modulated_min_2 <- kernel_WL_df_rle$CF_modulated_min / kernel_WL_df_rle$Healthy_min
# get original CF value relative to healthy value
kernel_WL_df_rle$CF_med_2 <- kernel_WL_df_rle$CF_med / kernel_WL_df_rle$Healthy_med


# Continue Weisfeiler-Lehman Graph Kernel evaluations with VST-normalised data
# calculate the median shortest path kernel per number of nodes that was inserted into CF network
kernel_WL_median_vst <- ddply(network_kernel_WL_vst,"n_nodes",numcolwise(median))
# rename columns
colnames(kernel_WL_median_vst) <- c("n_nodes_med", "i_med","Healthy_med", "CF_med", "CF_modulated_med")
# calculate the minimum shortest path kernel per number of nodes that was inserted into CF network
kernel_WL_min_vst <- ddply(network_kernel_WL_vst,"n_nodes",numcolwise(min))
# rename columns
colnames(kernel_WL_min_vst) <- c("n_nodes_min", "i_med","Healthy_min", "CF_min", "CF_modulated_min")
# calculate the maximum shortest path kernel per number of nodes that was inserted into CF network
kernel_WL_max_vst <- ddply(network_kernel_WL_vst,"n_nodes",numcolwise(max))
# rename columns
colnames(kernel_WL_max_vst) <- c("n_nodes_max", "i_med","Healthy_max", "CF_max", "CF_modulated_max")
# merge median, minimum and maximum values into one data frame
kernel_WL_df_vst <- data.frame(cbind(kernel_WL_median_vst, kernel_WL_min_vst, kernel_WL_max_vst))
# get median value relative to healthy median value
kernel_WL_df_vst$CF_modulated_med_2 <- kernel_WL_df_vst$CF_modulated_med / kernel_WL_df_vst$Healthy_med
# get maximum value relative to healthy maximum value
kernel_WL_df_vst$CF_modulated_max_2 <- kernel_WL_df_vst$CF_modulated_max / kernel_WL_df_vst$Healthy_max
# get minimum value relative to healthy minimum value
kernel_WL_df_vst$CF_modulated_min_2 <- kernel_WL_df_vst$CF_modulated_min / kernel_WL_df_vst$Healthy_min
# get original CF value relative to healthy value
kernel_WL_df_vst$CF_med_2 <- kernel_WL_df_vst$CF_med / kernel_WL_df_vst$Healthy_med

# add meta data 
network_kernel_WL_vst$normalisation <- "VST"
network_kernel_WL_bcphc$normalisation <- "BCPHC"
network_kernel_WL_rle$normalisation <- "RLE"
# merge VST, BCPHC, RLE-normalised data frames
network_kernel_WL <- data.frame(rbind(network_kernel_WL_vst, network_kernel_WL_rle, network_kernel_WL_bcphc))
# add meta data
kernel_WL_df_vst$normalisation <- "VST"
kernel_WL_df_rle$normalisation <- "RLE"
kernel_WL_df_bcphc$normalisation <- "BCPHC"
# merge VST, BCPHC, RLE-normalised data frames
kernel_WL_df <- data.frame(rbind(kernel_WL_df_vst, kernel_WL_df_rle, kernel_WL_df_bcphc))

network_kernel_shortestPath_vst$normalisation <- "VST"
network_kernel_shortestPath_bcphc$normalisation <- "BCPHC"
network_kernel_shortestPath_rle$normalisation <- "RLE"
network_kernel_shortestPath <- data.frame(rbind(network_kernel_shortestPath_vst, network_kernel_shortestPath_rle, network_kernel_shortestPath_bcphc))

kernel_shortestPahway_df_vst$normalisation <- "VST"
kernel_shortestPahway_df_rle$normalisation <- "RLE"
kernel_shortestPahway_df_bcphc$normalisation <- "BCPHC"
kernel_shortestPathway_df <- data.frame(rbind(kernel_shortestPahway_df_vst, kernel_shortestPahway_df_rle, kernel_shortestPahway_df_bcphc))

############################################################################################################
# generate simulation plots
# Weisfeiler-Lehman plot
kernel_WL_plot <-ggplot() + 
  geom_jitter(data=network_kernel_WL, aes(x=n_nodes, y=modulated.CF, shape=normalisation), colour="grey", width = 0.3, alpha=0.5, size=0.8) +
  geom_pointrange(data=kernel_WL_df, aes(x=n_nodes_med, y=CF_modulated_med_2, ymin=CF_modulated_min_2, ymax=CF_modulated_max_2, color=normalisation)) +
  geom_line(data=kernel_WL_df, aes(x=n_nodes_med, y=CF_modulated_med_2, color=normalisation), size=1, alpha=0.5) +
  geom_line(data=kernel_WL_df, aes(x=n_nodes_med, y=CF_med_2, color=normalisation), size=1, alpha=1, linetype="dotted") +
  xlab("Number of nodes") + 
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8), limits=c(0.0, 0.9)) + 
  scale_colour_manual(values=c("black", "red", "dodgerblue4")) + theme_bw(base_size=10) +
  theme(legend.title = element_blank(), panel.grid = element_blank(), axis.text = element_text(size=11), 
        axis.title = element_text(size=11), legend.text = element_text(size=11)) +
  ylab("1-dim Weisfeiler-Lehman kernel") 

# Frechet distance plot
frechet_plot <- ggplot() + 
  geom_jitter(data=frechet_attack_similarity, aes(x=n_nodes, y=frechet_2, shape=normalisation), color="grey", width = 0.3, alpha=0.5, size=0.8) +
  geom_pointrange(data=frechet_attack_df, aes(x=n_nodes_med, y=cf_hrare_deg_frechet_med_2, ymin=cf_hrare_deg_frechet_min_2, ymax=cf_hrare_deg_frechet_max_2, color=normalisation)) +
  geom_line(data=frechet_attack_df, aes(x=n_nodes_med, y=cf_hrare_deg_frechet_med_2, color=normalisation), size=1, alpha=0.5) +
  geom_line(data=frechet_attack_df, aes(x=n_nodes_med, y=round(cf_vs_healthy_degree,2), color=normalisation), size=1, alpha=1, linetype="dotted") +
  xlab("Number of nodes") + 
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8), limits=c(0.0, 0.5)) + 
  scale_colour_manual(values=c("black", "red", "dodgerblue4")) + 
  theme_bw() +
  theme(legend.title = element_blank(), panel.grid = element_blank(), legend.position = c(0.8,0.8), axis.text = element_text(size=11), 
        axis.title = element_text(size=11), legend.text = element_text(size=11)) +
  ylab("Network vulnerability") 

############################################################################################################
# Obtain correlation statistics
# Weisfeiler-lehmann graph kernel of modulated CF network and number of nodes
# VST-normalised data
cor.test(network_kernel_WL_vst$modulated.CF, network_kernel_WL_vst$n_nodes, method = "spearman") 
spearman.ci(network_kernel_WL_vst$modulated.CF, network_kernel_WL_vst$n_nodes, nrep=1000)
# BCPHC-normalised data
cor.test(network_kernel_WL_bcphc$modulated.CF, network_kernel_WL_bcphc$n_nodes, method = "spearman")
spearman.ci(network_kernel_WL_bcphc$modulated.CF, network_kernel_WL_bcphc$n_nodes, nrep=1000)
# RLE-normalised data
cor.test(network_kernel_WL_rle$modulated.CF, network_kernel_WL_rle$n_nodes, method = "spearman")
spearman.ci(network_kernel_WL_rle$modulated.CF, network_kernel_WL_rle$n_nodes, nrep=1000)


# Frechet distance of modulated CF network to healthy network and number of nodes
# VST-normalised data
cor.test(frechet_attack_similarity_vst$frechet, frechet_attack_similarity_vst$n_nodes, method = "spearman")
spearman.ci(frechet_attack_similarity_vst$frechet, frechet_attack_similarity_vst$n_nodes, nrep=1000)
# BCPHC-normalised data
cor.test(frechet_attack_similarity_bcphc$frechet, frechet_attack_similarity_bcphc$n_nodes, method = "spearman")
spearman.ci(frechet_attack_similarity_bcphc$frechet, frechet_attack_similarity_bcphc$n_nodes, nrep=1000)
# RLE-normalised data
cor.test(frechet_attack_similarity_rle$frechet, frechet_attack_similarity_rle$n_nodes, method = "spearman")
spearman.ci(frechet_attack_similarity_rle$frechet, frechet_attack_similarity_rle$n_nodes, nrep=1000)

############################################################################################################
# Evaluation of simulation outcome based on species number, species diversity and species dominance
# BCPHC-normalised data, join information of Frechet distance outputs and selected species
frechet_attack_similarity_bcphc <- frechet_attack_similarity_bcphc[,-c(4:5)]
info_bcphc <- selected_species_bcphc %>% right_join(frechet_attack_similarity_bcphc, by=c("i","n_nodes"))
# add meta data
info_bcphc$normalisation <- "BCPHC"

# RLE-normalised data, join information of Frechet distance outputs and selected species
frechet_attack_similarity_rle <- frechet_attack_similarity_rle[,-c(4:5)]
info_rle <- selected_species_rle %>% right_join(frechet_attack_similarity_rle, by=c("i","n_nodes"))
# add meta data
info_rle$normalisation <- "RLE"

# VST-normalised data, join information of Frechet distance outputs and selected species
frechet_attack_similarity_vst <- frechet_attack_similarity_vst[,-c(4:5)]
info_vst <- selected_species_vst %>% right_join(frechet_attack_similarity_vst, by=c("i","n_nodes"))
# add meta data
info_vst$normalisation <- "VST"

info_all <- data.frame(rbind(info_bcphc, info_rle, info_vst))

# add outcome column, if Frechet distance between modulated CF network and original healthy network is smaller than distance between
# original CF network and original healthy network, add "improved" otherwise add "worse"
info_all$Outcome <- 
  with(info_all, 
       ifelse(normalisation == "BCPHC" & frechet > original_frechet_h_cf_degree_bcphc, "Destabilisation", 
              ifelse(normalisation == "BCPHC" & frechet < original_frechet_h_cf_degree_bcphc, "Stabilisation", 
                     ifelse(normalisation == "RLE" & frechet > original_frechet_h_cf_degree_rle, "Destabilisation", 
                            ifelse(normalisation == "RLE" & frechet < original_frechet_h_cf_degree_rle, "Stabilisation",
                                   ifelse(normalisation == "VST" & frechet > original_frechet_h_cf_degree_vst, "Destabilisation",
                                          ifelse(normalisation == "VST" & frechet < original_frechet_h_cf_degree_vst, "Stabilisation", "Destabilisation")))))))

network_nullmodels <- data.frame(rbind(network_nullmodel_bcphc, network_nullmodel_rle, network_nullmodel_vst))
colnames(network_nullmodels) <- c("zscore", "p.value", "normalisation", "i", "n_nodes", "index" )

# generate data table with robustness information
network_nullmodels_robustness <- subset(network_nullmodels, index == "Robustness")
# combine data frames
network_nullmodels_robustness_1 <- network_nullmodels_robustness %>% right_join(info_all, by=c("i","n_nodes", "normalisation"))
# remove NA entries
network_nullmodels_robustness_2 <- na.omit(network_nullmodels_robustness_1)

# add columns based on p-value
network_nullmodels_robustness_2$p.value.cat <- ifelse(network_nullmodels_robustness_2$p.value >= 0.05, "p-value >= 0.05", "p-value < 0.05")
network_nullmodels_robustness_2$p.value.cat.normalisation <- paste(network_nullmodels_robustness_2$p.value.cat,";",network_nullmodels_robustness_2$normalisation)
network_nullmodels_robustness_2$p.value.cat.normalisation <- str_replace(network_nullmodels_robustness_2$p.value.cat.normalisation, " ;", ";")
network_nullmodels_robustness_2$p.value.size <- ifelse(network_nullmodels_robustness_2$p.value > 0.05, "not significant", "significant")

 # plot null model results based on robustness score
nullmodel_robustness <-
  ggplot(network_nullmodels_robustness_2) +
  geom_hline(yintercept = -1.96, linetype="dashed") +
  geom_hline(yintercept = 1.96, linetype="dashed") + 
  ylim(-30, 20) +
  geom_jitter(aes(x=n_nodes, y=zscore, colour=p.value.cat.normalisation, shape=p.value.cat.normalisation, size=p.value.size), alpha=0.4,width=1, height = 1) +
  facet_wrap(~Outcome, nrow=1) + theme_bw() + 
  xlab("Number of nodes") + ylab("Z-score (Null model)") +
  scale_colour_manual(values=c("black", "red", "dodgerblue4", "black", "red", "dodgerblue4"),
                      name = "Normalisation",
                      labels = c("p-value < 0.05; BCPHC", "p-value < 0.05; RLE", "p-value < 0.05; VST",
                                 "p-value > 0.05; BCPHC", "p-value > 0.05; RLE", "p-value > 0.05; VST")) +
  scale_shape_manual(values=c(19,19,19,17,17,17),
                     name = "Normalisation",
                     labels = c("p-value < 0.05; BCPHC", "p-value < 0.05; RLE", "p-value < 0.05; VST",
                                "p-value > 0.05; BCPHC", "p-value > 0.05; RLE", "p-value > 0.05; VST")) +
  scale_size_manual(values=c(1, 0.001), guide = "none") +
  scale_x_continuous(breaks = c(1,5,10,15,20), labels = c(1,5,10,15,20), limits = c(0,22)) + 
  guides(colour = guide_legend(override.as = list(alpha=1, size=3))) + 
  theme(panel.grid = element_blank(), legend.position = "none", strip.background = element_rect(fill="white"), strip.text = element_text(size=11), axis.text = element_text(size=11), axis.title = element_text(size=11), legend.text = element_text(size=11), legend.title = element_blank()) 

 
# generate data table with Weisfeiler Lehmann information
network_nullmodels_wlk <- subset(network_nullmodels, index == "Weisfeiler_Lehmann")
# combine data frames
new_figure_LW <- data.frame(rbind(network_kernel_WL_bcphc, network_kernel_WL_rle, network_kernel_WL_vst))
# remove row-names
rownames(new_figure_LW) <- NULL
# make character entries numeric
new_figure_LW$CF <- as.numeric(new_figure_LW$CF)
new_figure_LW$modulated.CF <- as.numeric(new_figure_LW$modulated.CF)
new_figure_LW$outcome <- ifelse(new_figure_LW$CF < new_figure_LW$modulated.CF, "Improvement", "Disimprovement")
# join dataframe
network_nullmodels_wlk_2 <- network_nullmodels_wlk %>% right_join(new_figure_LW, by=c("i","n_nodes", "normalisation"))
# remove NA values
network_nullmodels_wlk_2 <- na.omit(network_nullmodels_wlk_2)
# add columns based on p-value
network_nullmodels_wlk_2$p.value.cat <- ifelse(network_nullmodels_wlk_2$p.value >= 0.05, "p-value >= 0.05", "p-value < 0.05")
network_nullmodels_wlk_2$p.value.cat.normalisation <- paste(network_nullmodels_wlk_2$p.value.cat,network_nullmodels_wlk_2$normalisation)
network_nullmodels_wlk_2$p.value.size <- ifelse(network_nullmodels_wlk_2$p.value > 0.05, 0.35, 0.1)

# plot null model results based on topology score
nullmodel_wlk <- 
  ggplot(network_nullmodels_wlk_2) +
  geom_hline(yintercept = -1.96, linetype="dashed") +
  geom_hline(yintercept = 1.96, linetype="dashed") + 
  ylim(-10, 80) +
  geom_jitter(aes(x=n_nodes, y=zscore, colour=p.value.cat.normalisation, shape=p.value.cat.normalisation), size=0.01, alpha=0.5,width=1, height=0.5) + # size = 0.3
  facet_wrap(~outcome, nrow=1) + theme_bw() + 
  xlab("Number of nodes") + ylab("Z-score (Null model)") +
  scale_colour_manual(values=c("black", "red", "dodgerblue4", "black", "red", "dodgerblue4"),
                      name = "Normalisation",
                      labels = c("p-value < 0.05; BCPHC", "p-value < 0.05; RLE", "p-value < 0.05; VST",
                                 "p-value > 0.05; BCPHC", "p-value > 0.05; RLE", "p-value > 0.05; VST")) +
  scale_shape_manual(values=c(19,19,19,17,17,17),
                     name = "Normalisation",
                     labels = c("p-value < 0.05; BCPHC", "p-value < 0.05; RLE", "p-value < 0.05; VST",
                                "p-value > 0.05; BCPHC", "p-value > 0.05; RLE", "p-value > 0.05; VST")) +
  scale_x_continuous(breaks = c(1,5,10,15,20), labels = c(1,5,10,15,20), limits = c(0,22)) + scale_size(guide = "none") +
  theme(panel.grid = element_blank(), legend.position = "none", strip.background = element_rect(fill="white"), strip.text = element_text(size=11), axis.text = element_text(size=11), axis.title = element_text(size=11), legend.text = element_text(size=11), legend.title = element_blank()) 

# combine plots
rob_lwk_plot <- ggarrange(frechet_plot, nullmodel_robustness, labels = c("A", "B"), nrow=1, widths = c(1,1))
nullmodel_merged <- ggarrange(kernel_WL_plot, nullmodel_wlk, labels = c("C", "D"), legend = "none", nrow=1, widths = c(1,1))
rob_lwk_nullmodel <- ggarrange(rob_lwk_plot, nullmodel_merged, nrow=2)

#########################################################################################
# remove networks, which were not significantly different from the random structures
network_nullmodels_robustness_3 <- subset(network_nullmodels_robustness_2, p.value < 0.05 & zscore < 0)
# obtain background rare species for all normalisation strategies
background_rare <- c(background_rare_h_rle$Species, background_rare_h_bcphc$Species,background_rare_h_vst$Species)
# remove duplicates
background_rare <- background_rare[!duplicated(background_rare)]
# obtain background core species for all normalisation strategies
background_core <- c(background_core_h_rle$Species, background_core_h_bcphc$Species,background_core_h_vst$Species)
# remove duplicates
background_core <- background_core[!duplicated(background_core)]
# select variables 
network_table <- select(network_nullmodels_robustness_3, name, frechet, i, n_nodes, Outcome)
# generate new column with merged information
network_table$seed_node <- paste(network_table$i,"_",network_table$n_nodes,"_",network_table$Outcome)
# make data frame
network_table_2 <- data.frame(table(network_table$name, network_table$Outcome, network_table$seed_node))
# store data frame in new variable
network_table_3 <- network_table_2
# define species as "core" or "rare 
network_table_3$species_type <- ifelse(network_table_3$Var1 %in% background_rare, "Rare", ifelse(network_table_3$Var1 %in% background_core, "Core", "Core"))
# convert long data frame to format wide
network_table_2_long <- spread(network_table_3, key=Var3, value=Freq)
# store in new variable
network_table_3_long <- network_table_2_long
# make data frame with genus information
network_table_3_genus <- network_table_3_long
network_table_3_genus$Var1 <- as.character(network_table_3_genus$Var1)
network_table_3_genus2 <- sapply(strsplit(network_table_3_genus$Var1," "), `[`, 1)
network_table_3_genus$Genus <- network_table_3_genus2
# remove columns
network_table_3_genus$Var1 <- NULL
network_table_3_genus$Var2 <- NULL
network_table_3_genus$species_type <- NULL
# remove duplicate rows and sum values
network_table_3_genus <- ddply(network_table_3_genus, "Genus", numcolwise(sum))
rownames(network_table_3_genus) <- network_table_3_genus$Genus
network_table_3_genus$Genus <- NULL
# transpose data frame
network_table_3_genus_t <- data.frame(t(network_table_3_genus))
# make data frame with species type (core vs rare) information
network_table_2_core_rare <- network_table_2_long
network_table_2_core_rare$Var1 <- NULL
network_table_2_core_rare$Var2 <- NULL
network_table_2_core_rare <- ddply(network_table_2_core_rare, "species_type", numcolwise(sum))
rownames(network_table_2_core_rare) <- network_table_2_core_rare$species_type
network_table_2_core_rare$species_type <- NULL
# transpose data frame
network_table_2_core_rare_t <- data.frame(t(network_table_2_core_rare))

# remove columns
network_table_3_long$Var2 <- NULL
network_table_3_long$species_type <- NULL
network_table_3_long <- ddply(network_table_3_long, "Var1", numcolwise(sum))
rownames(network_table_3_long) <- network_table_3_long$Var1
network_table_3_long$Var1 <- NULL
network_table_3_long_t <- data.frame(t(network_table_3_long))
network_table_3_long_t <- network_table_3_long_t[order(rownames(network_table_3_long_t)),]
network_table_3_genus_t <- network_table_3_genus_t[order(rownames(network_table_3_genus_t)),]

# merge data frames
network_table_3_t <- data.frame(cbind(network_table_3_long_t, network_table_3_genus_t))
# built species diversity index based on species information
network_table_3_t$Shannon <- vegan::diversity(network_table_3_t, index="shannon")
network_table_3_t$Simpson <- vegan::diversity(network_table_3_t, index="simpson")

# add column based on stabilization/destabilization effect
network_table_3_t$Outcome <- ifelse(grepl("Destabilisation", rownames(network_table_3_long_t)),"Stabilisation", "Destabilisation")

# log-transform core and rare species counts and add pseudo count of 0.01
network_table_3_t$core_species <- log2(network_table_2_core_rare_t$Core+0.01)
network_table_3_t$rare_species <- log2(network_table_2_core_rare_t$Rare+0.01)
network_table_3_t <- network_table_3_t[,-c(1:29)]
network_table_3_t$Outcome <- factor(network_table_3_t$Outcome, levels = c("Stabilisation", "Destabilisation"))

############################################################################################################
# binomial regression analysis
# split dataframe into train and test dataset (70:30)
train <- network_table_3_t[1:1972,]
test <- network_table_3_t[1973:nrow(network_table_3_t),]

# run regression model based on train dataset
model <- glm(Outcome ~.,family=binomial(link='logit'),data=train)

# extract model summary
model_df <- summary(model)
model_df_2 <- data.frame(model_df$coefficients)
model_df_2$z.value <- ifelse(model_df_2$Pr...z.. > 0.05, 0, model_df_2$z.value)
model_df_2$variables <- rownames(model_df_2)
model_df_2$Category <- ifelse(model_df_2$z.value < 0, "Stabilisation", "Destabilisation")
model_df_3 <- model_df_2[-1,]
 
# rename column entries
model_df_3$variables <- str_replace(model_df_3$variables, "core_species", "Core species")
model_df_3$variables <- str_replace(model_df_3$variables, "rare_species", "Rare species")
model_df_3$variables <- str_replace(model_df_3$variables, "Shannon", "Shannon diversity")
model_df_3$variables <- str_replace(model_df_3$variables, "Simpson", "Simpson diversity")

# evaluate variance
anova(model, test="Chisq")
fitted.results <- predict(model,newdata=test,type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)

# model validation with test dataset
p <- predict(model, newdata=test, type="response")
pr <- prediction(p, test$Outcome, label.ordering = c("Stabilisation","Destabilisation"))
# true positive/false positive rate evaluation
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
# area under the curve
auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]] # 0.9358931

# make new data frame and rename/mutate columns
model_df_4 <- model_df_3 %>%
  mutate(
    variables = factor(variables, levels = variables[order(z.value, decreasing = TRUE)]),
    label_y = ifelse(z.value < 0, 0.2, -0.2),
    label_hjust = ifelse(z.value < 0, 0, 1))

# generate z-score plot
my_plot <- ggplot(model_df_4, aes(x = variables, y = z.value, fill = Category)) +
  geom_bar(stat = "identity", col = "black") +
  geom_text(aes(y = label_y, label = variables, hjust = label_hjust)) +
  coord_flip() +
  scale_fill_manual(values = c(Stabilisation = "darkgreen", Destabilisation = "darkred")) +
  theme_bw() + ylab("Z-score") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(size=12), axis.title.x = element_text(size=12),
        panel.grid = element_blank(), legend.text = element_text(size=12), legend.title = element_blank(), legend.position = "none") +
  scale_y_continuous(limits = c(-12, 12))

# generate area under the curve validation
auc_plot <- ggplot(data=NULL) +
  geom_line(aes(x=prf@x.values[[1]], y=prf@y.values[[1]])) +
  theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(size=12), 
                     axis.title = element_text(size=12)) + 
  geom_label(aes(x=0.7,0.1),label=paste("AUC =",round(auc,2))) +
  xlab("False positive rate") + ylab("True positive rate")

# merge both plots
plot_auc <- ggarrange(my_plot, auc_plot, widths = c(1,0.8))

############################################################################################################
# save figures
ggsave(filename = "output_figures/Figure_04.pdf",first_net_bcphc_plot, width=20, height=20, device=cairo_pdf, units="cm")
ggsave(filename="output_figures/Figure_05.pdf", rob_lwk_nullmodel, width = 25, height= 20, units="cm")
ggsave(filename = "output_figures/Figure_06.pdf", plot_auc, width = 20, height = 10, units = "cm")
ggsave(filename="output_figures/Supplementary_Figure_S2.pdf",net_rle_vst, width=20, height=20, units="cm")
