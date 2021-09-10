# title: "Step 3, Simulation runs"
# author: "Marie-Madlen Pust"
# last update: "10 September 2021"

############################################################################################################
# clean global environment
rm(list=ls())

# set working directory
setwd('C:/Simulations_early_infant_microbiome/R')

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
              'CINNA', 'graphkernels', 'data.table', 'microbiome', 'RVAideMemoire')

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

############################################################################################################
# load packages
ipak(packages)
############################################################################################################

# import meta data table of patients
md <- read_delim('input_files/meta_data/spatial_metadata_2020_12_2.csv', ';', escape_double = FALSE, trim_ws = TRUE)
# convert to data frame object
md <- data.frame(md)
# make sample IDs rownames
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
# Generate networks 
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
  geom_point(data = h.net_bcphc.plot.df, aes(x=V1,y=V2,colour=species_type),size=0.9) +  
  theme_pubr(border=FALSE, legend = "none")  +
  scale_colour_manual(values=c("springgreen4", "darkblue")) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

# generate network plot (CF, BCPHC-normalised)
cf_net_bcphc_plot <- ggplot() + 
  geom_segment(data = cf.net_bcphc.plot.edges, aes(x = from.x, xend = to.x, y = from.y, yend = to.y, colour=species_type), alpha=0.6) +
  geom_point(data = cf.net_bcphc.plot.df, aes(x=V1,y=V2,colour=species_type),size=0.9) +  
  theme_pubr(border=FALSE, legend = "none")  +
  scale_colour_manual(values=c("springgreen4", "darkblue")) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

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
  geom_point(data = h.net_rle.plot.df, aes(x=V1,y=V2,colour=species_type),size=0.9) +  
  theme_pubr(border=FALSE, legend = "none")  +
  scale_colour_manual(values=c("springgreen4", "darkblue")) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

# generate network plot (CF, RLE-normalised)
cf_net_rle_plot <- ggplot() + 
  geom_segment(data = cf.net_rle.plot.edges, aes(x = from.x, xend = to.x, y = from.y, yend = to.y, colour=species_type), alpha=0.6) +
  geom_point(data = cf.net_rle.plot.df, aes(x=V1,y=V2,colour=species_type),size=0.9) +  
  theme_pubr(border=FALSE, legend = "none")  +
  scale_colour_manual(values=c("springgreen4", "darkblue")) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

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
# apply fruchterman reingold graph 
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
  geom_point(data = h.net_vst.plot.df, aes(x=V1,y=V2,colour=species_type),size=0.9) +  
  theme_pubr(border=FALSE, legend = "none")  +
  scale_colour_manual(values=c("springgreen4", "darkblue")) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

# generate network plot, CF, VST-normalised
cf_net_vst_plot <- ggplot() + 
  geom_segment(data = cf.net_vst.plot.edges, aes(x = from.x, xend = to.x, y = from.y, yend = to.y, colour=species_type), alpha=0.6) +
  geom_point(data = cf.net_vst.plot.df, aes(x=V1,y=V2,colour=species_type),size=0.9) +  
  theme_pubr(border=FALSE, legend = "none")  +
  scale_colour_manual(values=c("springgreen4", "darkblue")) +
  theme(axis.text.x = element_blank(),  
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

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
  theme_bw(base_size=10) + 
  scale_colour_manual(values=c("darkred", "black")) + xlab("p (targeted)") + ylab("Network efficiency") +
  theme(panel.grid = element_blank(), 
        legend.title = element_blank(),
        axis.title = element_text(size=9))

# generate random plot
random_attack <- ggplot(all_matrix_random) +
  geom_vline(aes(xintercept=0.25), colour="grey", linetype="dashed") +
  geom_hline(aes(yintercept=0.44), colour="grey", linetype="dashed") +
  geom_hline(aes(yintercept=0.54), colour="grey", linetype="dashed") +
  geom_point(aes(x=fraction_removed, y=random, colour=state, shape=normalisation), size=1.2) +
  geom_line(aes(x=fraction_removed, y=random, colour=state, shape=normalisation), size=0.3) +
  theme_bw(base_size=10) + 
  scale_colour_manual(values=c("darkred", "black")) + xlab("p (random)") + ylab(" ") +
  theme(panel.grid = element_blank(), 
        legend.title = element_blank(),
        axis.title = element_text(size=9))

# merge cf and healthy plot, BCPHC-normalised
ab1_bcphc <- ggarrange(h_net_bcphc_plot, cf_net_bcphc_plot, nrow = 1, ncol = 2, labels = c("A", "B"), 
                       label.x = 0.029, label.y = 0.98, font.label = list(size = 10, color = "black")) 

# merge targeted and random attack plots
attack_plots <- ggarrange(targeted_attack, random_attack, labels=c("C", "D"), common.legend = TRUE, 
                          font.label = list(size = 10, color = "black"))

# merge cf and healthy network plots, BCPHC-normalised + attack curves of all normalisation steps
first_net_bcphc_plot <- ggarrange(ab1_bcphc, attack_plots, nrow=2, heights = c(0.6,1))

# merge cf and healthy network plots, VST and RLE-normalised
net_rle_vst <- ggarrange(cf_net_vst_plot, h_net_vst_plot, cf_net_rle_plot, h_net_rle_plot, heights = c(0.6,0.6))


############################################################################################################
# Network modulation 
# Set base seed (for reproducibility purposes)
set.seed(111)
# Generate random integers between 1 and 200
random_seeds <- sample(1:500, 200, replace = FALSE) 
# number of nodes to be inserted later on
add_nodes <- c(1:20)

# define empty variables globally for the loop assignment later on, BCPHC-normalised
frechet_attack_similarity_bcphc = NULL
network_statistics_bcphc = NULL
netwerk_kernel_randomWalk_bcphc = NULL
network_kernel_shortestPath_bcphc = NULL
network_kernel_WL_bcphc = NULL
selected_species_bcphc = NULL

# define empty variables globally for the loop assignment later on, RLE-normalised
frechet_attack_similarity_rle = NULL
network_statistics_rle = NULL
netwerk_kernel_randomWalk_rle = NULL
network_kernel_shortestPath_rle = NULL
network_kernel_WL_rle = NULL
selected_species_rle = NULL

# define empty variables globally for the loop assignment later on, VST-normalised
frechet_attack_similarity_vst = NULL
network_statistics_vst = NULL
netwerk_kernel_randomWalk_vst = NULL
network_kernel_shortestPath_vst = NULL
network_kernel_WL_vst = NULL
selected_species_vst = NULL
sample_with_replacement = TRUE 

# loop over list of random seeds
for (i in random_seeds){
  # repeat with an increasing number of species
  for (n_nodes in add_nodes){
    # set seed
    set.seed(i)
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
    new_cf_rare_net_bcphc <- graph_from_data_frame(new_cf_rare_edges_bcphc, 
                                                   directed=directed_network, vertices=new_cf_rare_nodes_bcphc)
    # assesses networks vulnerability due to targeted attacks with 10 iterations for assessing random error
    new_cf_net_attack_rare_bcphc <- swan_combinatory(new_cf_rare_net_bcphc,10)
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
    my_K_bcphc_df3 <- data.frame(my_K3_bcphc)
    # normalise to healthy structure
    my_K_bcphc_df3 <- my_K_bcphc_df3 / 	my_K_bcphc_df3[1,1]
    # rename columns
    colnames(my_K_bcphc_df3) <- c("Healthy", "CF", "modulated CF")
    # rename rows
    rownames(my_K_bcphc_df3) <- c("Healthy", "CF", "modulated CF")
    # keep first row
    my_K_bcphc_df3 <- my_K_bcphc_df3[1,]
    
    
    # export kernel output and bind with empty global variable (defined above)
    netwerk_kernel_randomWalk_bcphc <- rbind(netwerk_kernel_randomWalk_bcphc, data.frame(i, n_nodes, my_K_bcphc_df))
    network_kernel_shortestPath_bcphc <- rbind(network_kernel_shortestPath_bcphc, data.frame(i, n_nodes, my_K_bcphc_df2))
    network_kernel_WL_bcphc <- rbind(network_kernel_WL_bcphc, data.frame(i, n_nodes, my_K_bcphc_df3))
    
    
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
                                                                           diameter_cf_hrare_bcphc,  density_cf_hrare_bcphc, trans_cf_hrare_bcphc, eigen_cf_hrare_bcphc))
    
    
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
    new_cf_rare_net_rle <- graph_from_data_frame(new_cf_rare_edges_rle, 
                                                 directed=directed_network, vertices=new_cf_rare_nodes_rle)
    # assesses networks vulnerability due to targeted attacks with 10 iterations for assessing random error
    new_cf_net_attack_rare_rle <- swan_combinatory(new_cf_rare_net_rle,10)
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
    my_K_rle_df3 <- data.frame(my_K3_rle)
    # normalise to healthy structure
    my_K_rle_df3 <- my_K_rle_df3 / 	my_K_rle_df3[1,1]
    # rename columns
    colnames(my_K_rle_df3) <- c("Healthy", "CF", "modulated CF")
    # rename rows
    rownames(my_K_rle_df3) <- c("Healthy", "CF", "modulated CF")
    # keep first row
    my_K_rle_df3 <- my_K_rle_df3[1,]
    
    # export kernel output and bind with empty global variable (defined above)
    netwerk_kernel_randomWalk_rle <- rbind(netwerk_kernel_randomWalk_rle, data.frame(i, n_nodes, my_K_rle_df))
    network_kernel_shortestPath_rle <- rbind(network_kernel_shortestPath_rle, data.frame(i, n_nodes, my_K_rle_df2))
    network_kernel_WL_rle <- rbind(network_kernel_WL_rle, data.frame(i, n_nodes, my_K_rle_df3))
    
    
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
    my_K_vst_df3 <- data.frame(my_K3_vst)
    # normalise to healthy structure
    my_K_vst_df3 <- my_K_vst_df3 / 	my_K_vst_df3[1,1]
    # rename columns
    colnames(my_K_vst_df3) <- c("Healthy", "CF", "modulated CF")
    # rename rows
    rownames(my_K_vst_df3) <- c("Healthy", "CF", "modulated CF")
    # keep first row
    my_K_vst_df3 <- my_K_vst_df3[1,]
    
    # export kernel output and bind with empty global variable (defined above)
    netwerk_kernel_randomWalk_vst <- rbind(netwerk_kernel_randomWalk_vst, data.frame(i, n_nodes, my_K_vst_df))
    network_kernel_shortestPath_vst <- rbind(network_kernel_shortestPath_vst, data.frame(i, n_nodes, my_K_vst_df2))
    network_kernel_WL_vst <- rbind(network_kernel_WL_vst, data.frame(i, n_nodes, my_K_vst_df3))
    
    
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
frechet_attack_similartiy_degree <- data.frame(cbind(frechet_attack_similarity_bcphc_median, frechet_attack_similarity_bcphc_minimum, frechet_attack_similarity_bcphc_maximum))


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
frechet_attack_similartiy_degree_rle <- data.frame(cbind(frechet_attack_similarity_rle_median, frechet_attack_similarity_rle_minimum, frechet_attack_similarity_rle_maximum))


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
frechet_attack_similartiy_degree_vst <- data.frame(cbind(frechet_attack_similarity_vst_median,frechet_attack_similarity_vst_minimum, frechet_attack_similarity_vst_maximum))
# rename column names for all normalisation strategies so that they are matching
colnames(frechet_attack_similarity_vst) <- c("i", "n_nodes", "frechet")
colnames(frechet_attack_similarity_rle) <- c("i", "n_nodes", "frechet")
colnames(frechet_attack_similarity_bcphc) <- c("i", "n_nodes", "frechet")
# add meta data
frechet_attack_similarity_vst$normalisation <- "VST"
frechet_attack_similarity_bcphc$normalisation <- "BCPHC"
frechet_attack_similarity_rle$normalisation <- "RLE"
# bind all data frames into one
frechet_attack_similarity <- data.frame(rbind(frechet_attack_similarity_vst, frechet_attack_similarity_rle, frechet_attack_similarity_bcphc))
# add meta data
frechet_attack_similartiy_degree_vst$normalisation <- "VST"
frechet_attack_similartiy_degree$normalisation <- "BCPHC"
frechet_attack_similartiy_degree_rle$normalisation <- "RLE"
# bind data frames into one
frechet_attack_df <- data.frame(rbind(frechet_attack_similartiy_degree_vst, frechet_attack_similartiy_degree, frechet_attack_similartiy_degree_rle))

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
  geom_jitter(data=network_kernel_WL,
              aes(x=n_nodes, y=modulated.CF, shape=normalisation), 
              colour="grey", width = 0.3, alpha=0.5, size=0.8) +
  geom_pointrange(data=kernel_WL_df,
                  aes(x=n_nodes_med, y=CF_modulated_med_2, 
                      ymin=CF_modulated_min_2, 
                      ymax=CF_modulated_max_2, 
                      color=normalisation)) +
  geom_line(data=kernel_WL_df,
            aes(x=n_nodes_med, y=CF_modulated_med_2, color=normalisation), 
            size=1, alpha=0.5) +
  geom_line(data=kernel_WL_df,
            aes(x=n_nodes_med, y=CF_med_2, color=normalisation), 
            size=1, alpha=1, linetype="dotted") +
  xlab("Number of nodes") + 
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8), limits=c(0.0, 0.9)) + 
  scale_colour_manual(values=c("black", "darkgreen", "blue")) + theme_bw(base_size=10) +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size=9)) +
  ylab("Weisfeiler-Lehman kernel") 

# Shortest path kernel
kernel_SP_plot <-ggplot() + 
  geom_jitter(data=network_kernel_shortestPath,
              aes(x=n_nodes, y=modulated.CF, shape=normalisation), 
              colour="grey", width = 0.3, alpha=0.5, size=0.8) +
  geom_pointrange(data=kernel_shortestPathway_df,
                  aes(x=n_nodes_med, y=CF_modulated_med_2, 
                      ymin=CF_modulated_min_2, 
                      ymax=CF_modulated_max_2, 
                      color=normalisation)) +
  geom_line(data=kernel_shortestPathway_df,
            aes(x=n_nodes_med, y=CF_modulated_med_2, color=normalisation), 
            size=1, alpha=0.5) +
  geom_line(data=kernel_shortestPathway_df,
            aes(x=n_nodes_med, y=CF_med_2, color=normalisation), 
            size=1, alpha=1, linetype="dotted") +
  xlab("Number of nodes") + 
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8), limits=c(0.0, 0.9)) + 
  scale_colour_manual(values=c("black", "darkgreen", "blue")) + theme_bw(base_size=10) +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size=9)) +
  ylab("Shortest pathway kernel") 

# Frechet distance plot
frechet_plot <- ggplot() + 
  geom_jitter(data=frechet_attack_similarity, aes(x=n_nodes, y=frechet, shape=normalisation), 
              color="grey", width = 0.3, alpha=0.5, size=0.8) +
  geom_pointrange(data=frechet_attack_df,
                  aes(x=n_nodes_med, y=cf_hrare_deg_frechet_med, 
                      ymin=cf_hrare_deg_frechet_min,
                      ymax=cf_hrare_deg_frechet_max,
                      color=normalisation)) +
  geom_line(data=frechet_attack_df,
            aes(x=n_nodes_med, y=cf_hrare_deg_frechet_med, color=normalisation), 
            size=1, alpha=0.5) +
  geom_line(data=frechet_attack_df,
            aes(x=n_nodes_med, y=round(cf_vs_healthy_degree,2), color=normalisation), 
            size=1, alpha=1, linetype="dotted") +
  xlab("Number of nodes") + 
  scale_y_continuous(breaks = c(0.0, 0.2, 0.4, 0.6, 0.8), limits=c(0.0, 0.9)) + 
  scale_colour_manual(values=c("black", "darkgreen", "blue")) + theme_bw(base_size=10) +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size=9)) +
  ylab("Network vulnerability") 

# combine simulation plots
simulation_plots <- ggarrange(kernel_WL_plot, kernel_SP_plot, frechet_plot, common.legend = TRUE, 
                              nrow=1, labels = c("A", "B", "C"), font.label = list(size = 9, color = "black"))

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

# Shortest path graph kernel of modulated CF network and number of nodes
# BCPHC-normalised data
cor.test(network_kernel_shortestPath_bcphc$modulated.CF, network_kernel_shortestPath_bcphc$n_nodes, method = "spearman")
spearman.ci(network_kernel_shortestPath_bcphc$modulated.CF, network_kernel_shortestPath_bcphc$n_nodes, nrep=1000)
# VST-normalised data
cor.test(network_kernel_shortestPath_vst$modulated.CF, network_kernel_shortestPath_vst$n_nodes, method = "spearman")
spearman.ci(network_kernel_shortestPath_vst$modulated.CF, network_kernel_shortestPath_vst$n_nodes, nrep=1000)
# RLE-normalised data
cor.test(network_kernel_shortestPath_rle$modulated.CF, network_kernel_shortestPath_rle$n_nodes, method = "spearman")
spearman.ci(network_kernel_shortestPath_rle$modulated.CF, network_kernel_shortestPath_rle$n_nodes, nrep=1000)

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

# join information of all three normalisation strategies (VST, RLE and BCPHC)
info_all <- data.frame(rbind(info_vst, info_bcphc, info_rle))

# add outcome column, if Frechet distance between modulated CF network and original healthy network is smaller than distance between
# original CF network and original healthy network, add "improved" otherwise add "worse"
info_all$Outcome <- 
  with(info_all, 
       ifelse(normalisation == "BCPHC" & frechet > original_frechet_h_cf_degree_bcphc, "Worse", 
              ifelse(normalisation == "BCPHC" & frechet < original_frechet_h_cf_degree_bcphc, "Improved", 
                     ifelse(normalisation == "RLE" & frechet > original_frechet_h_cf_degree_rle, "Worse", 
                            ifelse(normalisation == "RLE" & frechet < original_frechet_h_cf_degree_rle, "Improved",
                                   ifelse(normalisation == "VST" & frechet > original_frechet_h_cf_degree_vst, "Worse",
                                          ifelse(normalisation == "VST" & frechet < original_frechet_h_cf_degree_vst, "Improved", "unknown")))))))
# count number of worse and improved occurrences per node, species and normalisation strategy
table_outcome_name <- as.data.frame(table(info_all$Outcome, info_all$name, info_all$n_nodes, info_all$normalisation))
# subset and extract all runs that improved the CF outcome
table_outcome_name_improved <- subset(table_outcome_name, Var1 == "Improved")
# convert frequency to percentage
table_outcome_name_improved$Per <- (table_outcome_name_improved$Freq / sum(table_outcome_name_improved$Freq)) * 100
# subset and extract all runs that destabilized the CF outcome
table_outcome_name_worse <- subset(table_outcome_name, Var1 == "Worse")
# convert frequency to percentage
table_outcome_name_worse$Per <- (table_outcome_name_worse$Freq / sum(table_outcome_name_worse$Freq)) * 100

# bind improved and worse output tables
table_outcome <- data.frame(rbind(table_outcome_name_worse, table_outcome_name_improved))
# convert count table from character to numeric object
table_outcome$Var3 <- as.numeric(as.character(table_outcome$Var3))
# scale the percentage column
table_outcome$Per_scale <- rescale(table_outcome$Per, to=c(-2,2))

# obtain background rare species for all normalisation strategies
background_rare <- c(background_rare_h_rle$Species, background_rare_h_bcphc$Species,background_rare_h_vst$Species)
# remove duplicates
background_rare <- background_rare[!duplicated(background_rare)]
# obtain background core species for all normalisation strategies
background_core <- c(background_core_h_rle$Species, background_core_h_bcphc$Species,background_core_h_vst$Species)
# remove duplicates
background_core <- background_core[!duplicated(background_core)]
# add core/rare information to count table
table_outcome$species_type <- ifelse(table_outcome$Var2 %in% background_rare, "Rare", ifelse(table_outcome$Var2 %in% background_core, "Core", "Core"))
# extract runs that improved modulated CF networks
table_outcome_improved <- subset(table_outcome, Var1 == "Improved")
# convert dataframe from long to wide format
merge_hboth_wide <- spread(table_outcome, key="Var2", value="Freq")
# remove percentage and scale columns
merge_hboth_wide$Per <- NULL
merge_hboth_wide$Per_scale <- NULL

# extract runs that improved the network
merge_hboth_wide_better <- subset(merge_hboth_wide, Var1 == "Improved")
# set NAs to 0
merge_hboth_wide_better[is.na(merge_hboth_wide_better)] <- 0
# remove non-numeric species type (core, rare) column
merge_hboth_wide_better$species_type <- NULL

# subset and extract BCPHC-normalised data
merge_hboth_wide_better_bcphc <- subset(merge_hboth_wide_better, Var4 == "BCPHC")
# remove duplicate rows and sum those entries
merge_hboth_wide_better_bcphc <- ddply(merge_hboth_wide_better_bcphc,"Var3",numcolwise(sum))
# rename rows
rownames(merge_hboth_wide_better_bcphc) <- merge_hboth_wide_better_bcphc$Var3
# remove column
merge_hboth_wide_better_bcphc$Var3 <- NULL
# generate data frame object
merge_hboth_wide_better_bcphc <- data.frame(merge_hboth_wide_better_bcphc)
# get alpha diversity indices
hboth_better_fisher_bcphc <- vegan::fisher.alpha(merge_hboth_wide_better_bcphc)
# Shannon diversity index
hboth_better_shannon_bcphc <- vegan::diversity(merge_hboth_wide_better_bcphc, index = "shannon")
# Species number
hboth_better_specNum_bcphc <- vegan::specnumber(merge_hboth_wide_better_bcphc)
# Simpson diversity index
hboth_better_dominance_bcphc <- vegan::diversity(merge_hboth_wide_better_bcphc, index = "simpson")
# Dominance indices
hboth_better_BPindex_bcphc <- microbiome::dominance(t(merge_hboth_wide_better_bcphc))
# store all indices in one data frame
div_hboth_better_bcphc <- data.frame(cbind(hboth_better_fisher_bcphc, hboth_better_shannon_bcphc, hboth_better_specNum_bcphc, hboth_better_dominance_bcphc, 
                                           hboth_better_BPindex_bcphc$gini))
# rename columns
colnames(div_hboth_better_bcphc) <- c("fisher", "shannon", "specNumber", "dominance", "gini")
# add meta data
div_hboth_better_bcphc$performance <- "Improved"
div_hboth_better_bcphc$normalisation <- "BCPHC"


# subset and extract RLE-normalised data
merge_hboth_wide_better_rle <- subset(merge_hboth_wide_better, Var4 == "RLE")
# remove duplicate rows and sum those entries
merge_hboth_wide_better_rle <- ddply(merge_hboth_wide_better_rle,"Var3",numcolwise(sum))
# rename rows
rownames(merge_hboth_wide_better_rle) <- merge_hboth_wide_better_rle$Var3
# remove column
merge_hboth_wide_better_rle$Var3 <- NULL
# generate data frame object
merge_hboth_wide_better_rle <- data.frame(merge_hboth_wide_better_rle)
# get alpha diversity indices
hboth_better_fisher_rle <- vegan::fisher.alpha(merge_hboth_wide_better_rle)
# Shannon diversity index
hboth_better_shannon_rle <- vegan::diversity(merge_hboth_wide_better_rle, index = "shannon")
# Species number
hboth_better_specNum_rle <- vegan::specnumber(merge_hboth_wide_better_rle)
# Simpson diversity index
hboth_better_dominance_rle <- vegan::diversity(merge_hboth_wide_better_rle, index = "simpson")
# Dominance indices
hboth_better_BPindex_rle <- microbiome::dominance(t(merge_hboth_wide_better_rle))
# store all indices in one data frame
div_hboth_better_rle <- data.frame(cbind(hboth_better_fisher_rle, hboth_better_shannon_rle, hboth_better_specNum_rle, hboth_better_dominance_rle, 
                                         hboth_better_BPindex_rle$gini))
# rename columns
colnames(div_hboth_better_rle) <- c("fisher", "shannon", "specNumber", "dominance", "gini")
# add meta data
div_hboth_better_rle$performance <- "Improved"
div_hboth_better_rle$normalisation <- "RLE"

# subset and extract VST-normalised data
merge_hboth_wide_better_vst <- subset(merge_hboth_wide_better, Var4 == "VST")
# remove duplicate rows and sum those entries
merge_hboth_wide_better_vst <- ddply(merge_hboth_wide_better_vst,"Var3",numcolwise(sum))
# rename rows
rownames(merge_hboth_wide_better_vst) <- merge_hboth_wide_better_vst$Var3
# remove column
merge_hboth_wide_better_vst$Var3 <- NULL
# generate data frame object
merge_hboth_wide_better_vst <- data.frame(merge_hboth_wide_better_vst)
# get alpha diversity indices
hboth_better_fisher_vst <- vegan::fisher.alpha(merge_hboth_wide_better_vst)
# Shannon diversity index
hboth_better_shannon_vst <- vegan::diversity(merge_hboth_wide_better_vst, index = "shannon")
# Species number
hboth_better_specNum_vst <- vegan::specnumber(merge_hboth_wide_better_vst)
# Simpson diversity index
hboth_better_dominance_vst <- vegan::diversity(merge_hboth_wide_better_vst, index = "simpson")
# Dominance indices
hboth_better_BPindex_vst <- microbiome::dominance(t(merge_hboth_wide_better_vst))
# store all indices in one data frame
div_hboth_better_vst <- data.frame(cbind(hboth_better_fisher_vst, hboth_better_shannon_vst, 
                                         hboth_better_specNum_vst, hboth_better_dominance_vst, 
                                         hboth_better_BPindex_vst$gini))
# rename columns
colnames(div_hboth_better_vst) <- c("fisher", "shannon", "specNumber", "dominance", "gini")
# add metadata 
div_hboth_better_vst$performance <- "Improved"
div_hboth_better_vst$normalisation <- "VST"

# combine all normalisation outputs in one data frame
div_hboth_better <- data.frame(rbind(div_hboth_better_bcphc,div_hboth_better_rle,div_hboth_better_vst))

# extract runs that destabilised the network
merge_hboth_wide_worse <- subset(merge_hboth_wide, Var1 == "Worse")

# subset and extract BCPHC-normalised data
merge_hboth_wide_worse_bcphc <- subset(merge_hboth_wide_worse, Var4 == "BCPHC")
# set NAs to 0
merge_hboth_wide_worse_bcphc[is.na(merge_hboth_wide_worse_bcphc)] <- 0
# remove non-numeric species type (core, rare) column
merge_hboth_wide_worse_bcphc$species_type <- NULL
# remove duplicate rows and sum those entries
merge_hboth_wide_worse_bcphc <- ddply(merge_hboth_wide_worse_bcphc,"Var3",numcolwise(sum))
# rename rows
rownames(merge_hboth_wide_worse_bcphc) <- merge_hboth_wide_worse_bcphc$Var3
# remove columns
merge_hboth_wide_worse_bcphc$Var3 <- NULL
merge_hboth_wide_worse_bcphc$Var4 <- NULL
# generate data frame object
merge_hboth_wide_worse_bcphc <- data.frame(merge_hboth_wide_worse_bcphc)
# get alpha diversity indices
hboth_worse_fisher_bcphc <- vegan::fisher.alpha(merge_hboth_wide_worse_bcphc)
# Shannon diversity index
hboth_worse_shannon_bcphc <- vegan::diversity(merge_hboth_wide_worse_bcphc, index = "shannon")
# Species number
hboth_worse_specNum_bcphc <- vegan::specnumber(merge_hboth_wide_worse_bcphc)
# Simpson diversity index
hboth_worse_dominance_bcphc <- vegan::diversity(merge_hboth_wide_worse_bcphc, index = "simpson")
# Dominance indices
hboth_worse_BPindex_bcphc <- microbiome::dominance(t(merge_hboth_wide_worse_bcphc))
# store all indices in one data frame
div_hboth_worse_bcphc <- data.frame(cbind(hboth_worse_fisher_bcphc, hboth_worse_shannon_bcphc, hboth_worse_specNum_bcphc, hboth_worse_dominance_bcphc, 
                                          hboth_worse_BPindex_bcphc$gini))
# rename columns
colnames(div_hboth_worse_bcphc) <- c("fisher", "shannon", "specNumber", "dominance", "gini")
# add metadata
div_hboth_worse_bcphc$performance <- "Worse"
div_hboth_worse_bcphc$normalisation <- "BCPHC"

# subset and extract RLE-normalised data
merge_hboth_wide_worse_rle <- subset(merge_hboth_wide_worse, Var4 == "RLE")
# set NAs to zero
merge_hboth_wide_worse_rle[is.na(merge_hboth_wide_worse_rle)] <- 0
# remove column
merge_hboth_wide_worse_rle$species_type <- NULL
# remove duplicate rows and sum those entries
merge_hboth_wide_worse_rle <- ddply(merge_hboth_wide_worse_rle,"Var3",numcolwise(sum))
# rename rows
rownames(merge_hboth_wide_worse_rle) <- merge_hboth_wide_worse_rle$Var3
# remove columns
merge_hboth_wide_worse_rle$Var3 <- NULL
merge_hboth_wide_worse_rle$Var4 <- NULL
# generate data frame object
merge_hboth_wide_worse_rle <- data.frame(merge_hboth_wide_worse_rle)
# get alpha diversity indices
hboth_worse_fisher_rle <- vegan::fisher.alpha(merge_hboth_wide_worse_rle)
# Shannon diversity index
hboth_worse_shannon_rle <- vegan::diversity(merge_hboth_wide_worse_rle, index = "shannon")
# Species number
hboth_worse_specNum_rle <- vegan::specnumber(merge_hboth_wide_worse_rle)
# Simpson diversity
hboth_worse_dominance_rle <- vegan::diversity(merge_hboth_wide_worse_rle, index = "simpson")
# Dominance index
hboth_worse_BPindex_rle <- microbiome::dominance(t(merge_hboth_wide_worse_rle))
# store all indices in one data frame
div_hboth_worse_rle <- data.frame(cbind(hboth_worse_fisher_rle, hboth_worse_shannon_rle, hboth_worse_specNum_rle, hboth_worse_dominance_rle, 
                                        hboth_worse_BPindex_rle$gini))
# add column names
colnames(div_hboth_worse_rle) <- c("fisher", "shannon", "specNumber", "dominance", "gini")
# add meta data
div_hboth_worse_rle$performance <- "Worse"
div_hboth_worse_rle$normalisation <- "RLE"

# subset and extract VST-normalised data
merge_hboth_wide_worse_vst <- subset(merge_hboth_wide_worse, Var4 == "VST")
# set NAs to zero
merge_hboth_wide_worse_vst[is.na(merge_hboth_wide_worse_vst)] <- 0
# remove column
merge_hboth_wide_worse_vst$species_type <- NULL
# remove duplicate rows and sum those entries
merge_hboth_wide_worse_vst <- ddply(merge_hboth_wide_worse_vst,"Var3",numcolwise(sum))
# rename rows
rownames(merge_hboth_wide_worse_vst) <- merge_hboth_wide_worse_vst$Var3
# remove columns
merge_hboth_wide_worse_vst$Var3 <- NULL
merge_hboth_wide_worse_vst$Var4 <- NULL
# generate data frame object
merge_hboth_wide_worse_vst <- data.frame(merge_hboth_wide_worse_vst)
# get alpha diversity indices
hboth_worse_fisher_vst <- vegan::fisher.alpha(merge_hboth_wide_worse_vst)
# Shannon diversity index
hboth_worse_shannon_vst <- vegan::diversity(merge_hboth_wide_worse_vst, index = "shannon")
# Species number
hboth_worse_specNum_vst <- vegan::specnumber(merge_hboth_wide_worse_vst)
# Simpson diversity index
hboth_worse_dominance_vst <- vegan::diversity(merge_hboth_wide_worse_vst, index = "simpson")
# Dominance index
hboth_worse_BPindex_vst <- microbiome::dominance(t(merge_hboth_wide_worse_vst))
# store all indices in one data frame
div_hboth_worse_vst <- data.frame(cbind(hboth_worse_fisher_vst, hboth_worse_shannon_vst, 
                                        hboth_worse_specNum_vst, hboth_worse_dominance_vst, 
                                        hboth_worse_BPindex_vst$gini))
# rename columns
colnames(div_hboth_worse_vst) <- c("fisher", "shannon", "specNumber", "dominance", "gini")
# add meta data
div_hboth_worse_vst$performance <- "Worse"
div_hboth_worse_vst$normalisation <- "VST"

# combine all normalisation outputs in one data frame
div_hboth_worse <- data.frame(rbind(div_hboth_worse_bcphc,div_hboth_worse_rle, div_hboth_worse_vst))
# combine improved and worse output into one data frame
div_hboth <- data.frame(rbind(div_hboth_better, div_hboth_worse))
# convert performance ("improved", "worse) from class character to class factor
div_hboth$performance <- as.factor(as.character(div_hboth$performance))

############################################################################################################
# Statistics
# Is the Shannon diversity index significantly different between good and bad performance of modulated CF network?
wilcox.test(shannon ~ performance, data=div_hboth)
# Calculate effect size with confidence intervals
wilcoxonR(div_hboth$shannon, g=div_hboth$performance, ci=TRUE)

# Is the Community dominance index significantly different between good and bad performance of modulated CF network?
wilcox.test(gini ~ performance, data=div_hboth)
# Calculate effect size with confidence intervals
wilcoxonR(div_hboth$gini, g=div_hboth$performance, ci=TRUE)

# Frequency plot
# Are core species more important than rare speces in improving the CF network?
wilcox.test(Per ~ species_type, data=table_outcome_improved)
wilcoxonR(table_outcome_improved$Per, g=table_outcome_improved$species_type, ci=TRUE)

############################################################################################################
# generate plots for simulation runs
# plot heatmap
heatmap_plot <- ggplot(table_outcome) +
  geom_tile(aes(x=Var3, y=Var2, fill=Per_scale)) +
  scale_fill_gradientn(colours=c("white", "peachpuff", "firebrick1", "firebrick2", "black")) +
  facet_wrap(~Var1 ~Var4, nrow=1) + xlab("Number of nodes") +  ylab("") +
  scale_x_continuous(breaks=c(3,6,9), limits = c(1.5,10.5)) +
  theme_bw(base_size=10) + theme(panel.grid = element_blank(), panel.background = element_blank()) +
  theme(legend.title = element_blank(), axis.text.y = element_text(face = "italic"),
        strip.background = element_rect(fill = "white"))

freq_plot <-  ggplot(data=table_outcome_improved, aes(x=species_type, y=scale(Per))) +
  geom_boxplot(width=0.4, colour="darkred", outlier.alpha = 0) +
  geom_jitter(width=0.05, colour="black", alpha=0.4) +
  stat_summary(fun=median, geom="point", color="darkred") +
  stat_compare_means(label = "p.signif", label.x.npc = "centre", label.y = 1.4, size=3) +
  ylab("Scaled frequency") + theme_bw(base_size=10) + xlab("Species type") + coord_flip() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=8))

hboth_gini <- ggplot(div_hboth, aes(x=performance, y=gini)) +
  geom_boxplot(width=0.4, colour="darkred", outlier.size = 0, outlier.alpha = 0) +
  geom_jitter(width=0.05, colour="black", alpha=0.4) +
  stat_summary(fun=median, geom="point", color="darkred") +
  stat_compare_means(label = "p.signif", label.x = 1.4, label.y = 0.77 ,size=5) +
  ylab("Community dominance index") + theme_bw(base_size=10) + xlab("Performance") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=8)) + ylim(0.1, 0.8)

hboth_shannon <- ggplot(div_hboth, aes(x=performance, y=shannon)) + 
  geom_boxplot(width=0.4, colour="darkred", outlier.size = 0, outlier.alpha = 0) +
  geom_jitter(width=0.05, colour="black", alpha=0.4) +
  stat_summary(fun=median, geom="point", color="darkred") +
  stat_compare_means(label = "p.signif", label.x = 1.4, label.y = 3.90, size=5) +
  ylab("Shannon diversity index") + theme_bw(base_size=10) + xlab("Performance") +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=8)) + ylim(1, 4)

hboth_shannon_gini_plot <- ggarrange(hboth_shannon, hboth_gini, labels = c("A", "B"), font.label = list(size = 10, color = "black"))
freq_plot_2 <- ggarrange(freq_plot, labels = c("C"), font.label = list(size = 10, color = "black"))
all_stats_plot <- ggarrange(hboth_shannon_gini_plot, freq_plot_2, nrow=2, heights=c(1,0.6))

############################################################################################################
# save networks and attack simulations
# save BCPHC-normalised networks
tiff(filename="output_figures/Figure_04.tif", res=600, units="in", width=6, height = 6)
first_net_bcphc_plot
dev.off()

# save VST and RLE-normalised networks
tiff(filename="output_figures/other_nets.tif", res=600, units="in", width=6, height = 6)
net_rle_vst
dev.off()

# save simulation plots
tiff(filename="output_figures/Figure_05.tif", res=300, units="in", width=9, height = 4)
simulation_plots
dev.off()

# save simulation statistic plots
tiff(filename="output_figures/Figure_07.tif", res=600, units="in", width=4.5, height = 5)
all_stats_plot
dev.off()

# save heatmap
tiff(filename="output_figures/Figure_06.tif", res=600, units="in", width=9, height=6)
heatmap_plot
dev.off()
############################################################################################################
        
