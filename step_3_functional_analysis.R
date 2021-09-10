# title: "Functional analysis"
# author: "Marie-Madlen Pust"
# date: "09 September 2021"


############################################################################################################
# clean global environment
rm(list=ls())

# set working directory
setwd('C:/Users/marie/Desktop/Simulations_early_infant_microbiome/R')

############################################################################################################
# define global functions
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, 'Package'])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)}

# load packages
packages <- c('readr', 'plyr', 'purrr', 'dplyr', 'ggplot2', 'vegan', 'ggpubr', 'factoextra', 'stringr',
              'pheatmap', 'RColorBrewer', 'ggdendro', 'tidyr', 'viridis', 'rstatix')

set_diff = 10 
set_statistics = mean 
score_value = 50
############################################################################################################

# load required R packages
ipak(packages)

############################################################################################################
# import data sets
# import table with background species
bg_species <- read_delim("input_files/taxonomic_data/background_species.csv", 
                                 delim = ";", escape_double = FALSE, trim_ws = TRUE)
bg_species <- data.frame(bg_species)
# extract species from the 25th abundance percentile definition
bg_species_25 <- subset(bg_species, Threshold == "25th percentile")
# extract species from the RLE normalisation data set
bg_species_25_rle <- subset(bg_species_25, Normalisation == "RLE")
# divide data frame into healthy and CF species
bg_species_25_rle_h <- subset(bg_species_25_rle, State == "Healthy")
bg_species_25_rle_cf <- subset(bg_species_25_rle, State == "CF")
# extract species from the VST normalisation data set
bg_species_25_vst <- subset(bg_species_25, Normalisation == "VST")
# divide data frame into healthy and CF species
bg_species_25_vst_h <- subset(bg_species_25_vst, State == "Healthy")
bg_species_25_vst_cf <- subset(bg_species_25_vst, State == "CF")
# extract species from the BCPHC normalisation data set
bg_species_25_bcphc <- subset(bg_species_25, Normalisation == "BCPHC")
# divide data frame into healthy and CF species
bg_species_25_bcphc_h <- subset(bg_species_25_bcphc, State == "Healthy")
bg_species_25_bcphc_cf <- subset(bg_species_25_bcphc, State == "CF")

# import the FM scores obtained from the PAO1 reference sequence
Pao1 <- read_delim("input_files/functional_data/PA01_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
Pao1 <- data.frame(Pao1)
Pao1$pathway_key_reactions <- NULL
Pao1$pathway_completeness <- NULL
Pao1$pathway_candidate_reaction <- NULL
# rename columns
colnames(Pao1) <- c("id", "pathway", "Pao1")

# import the FM scores obtained from the Staphylococcus aureus reference sequence
SA <- read_delim("input_files/functional_data/SA_pathways_reactions.csv", ";", 
                   escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                           pathway_completeness = col_number(), 
                                                           pathway_completeness_perc = col_number(), 
                                                           pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
SA <- data.frame(SA)
SA$pathway_key_reactions <- NULL
SA$pathway_completeness <- NULL
SA$pathway_candidate_reaction <- NULL
# rename columns
colnames(SA) <- c("id", "pathway", "SA")

# import the FM scores obtained from the Actinomyces israelii reference sequence
ActIsr <- read_delim("input_files/functional_data/ActIsr_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
ActIsr <- data.frame(ActIsr)
ActIsr$pathway_key_reactions <- NULL
ActIsr$pathway_completeness <- NULL
ActIsr$pathway_candidate_reaction <- NULL
# rename columns
colnames(ActIsr) <- c("id", "pathway", "ActIsr")


# import the FM scores obtained from the Capnocytophaga endodontalis reference sequence
CapEnd <- read_delim("input_files/functional_data/CapEnd_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
CapEnd <- data.frame(CapEnd)
CapEnd$pathway_key_reactions <- NULL
CapEnd$pathway_completeness <- NULL
CapEnd$pathway_candidate_reaction <- NULL
# rename columns
colnames(CapEnd) <- c("id", "pathway", "CapEnd")

# import the FM scores obtained from the Capnocytophaga sputigena reference sequence
CapSpu <- read_delim("input_files/functional_data/CapSpu_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
CapSpu <- data.frame(CapSpu)
CapSpu$pathway_key_reactions <- NULL
CapSpu$pathway_completeness <- NULL
CapSpu$pathway_candidate_reaction <- NULL
# rename columns
colnames(CapSpu) <- c("id", "pathway", "CapSpu")

# import the FM scores obtained from the Fusobacterium nucleatum reference sequence
FusNuc <- read_delim("input_files/functional_data/FusNuc_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
FusNuc <- data.frame(FusNuc)
FusNuc$pathway_key_reactions <- NULL
FusNuc$pathway_completeness <- NULL
FusNuc$pathway_candidate_reaction <- NULL
# rename columns
colnames(FusNuc) <- c("id", "pathway", "FusNuc")

# import the FM scores obtained from the Gemella haemolysans reference sequence
GemHae <- read_delim("input_files/functional_data/GemHae_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
GemHae <- data.frame(GemHae)
GemHae$pathway_key_reactions <- NULL
GemHae$pathway_completeness <- NULL
GemHae$pathway_candidate_reaction <- NULL
# rename columns
colnames(GemHae) <- c("id", "pathway", "GemHae")

# import the FM scores obtained from the Gemella morbillorum reference sequence
GemMor <- read_delim("input_files/functional_data/GemMor_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
GemMor <- data.frame(GemMor)
GemMor$pathway_key_reactions <- NULL
GemMor$pathway_completeness <- NULL
GemMor$pathway_candidate_reaction <- NULL
# rename columns
colnames(GemMor) <- c("id", "pathway", "GemMor")

# import the FM scores obtained from the Gemella sanguinis reference sequence
GemSan <- read_delim("input_files/functional_data/GemSan_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
GemSan <- data.frame(GemSan)
GemSan$pathway_key_reactions <- NULL
GemSan$pathway_completeness <- NULL
GemSan$pathway_candidate_reaction <- NULL
# rename columns
colnames(GemSan) <- c("id", "pathway", "GemSan")


# import the FM scores obtained from the Haemophilus parainflunzae reference sequence
HaePar <- read_delim("input_files/functional_data/HaePar_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
HaePar <- data.frame(HaePar)
HaePar$pathway_key_reactions <- NULL
HaePar$pathway_completeness <- NULL
HaePar$pathway_candidate_reaction <- NULL
# rename columns
colnames(HaePar) <- c("id", "pathway", "HaePar")

# import the FM scores obtained from the Leptotrichia buccalis reference sequence
LepBuc <- read_delim("input_files/functional_data/LepBuc_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
LepBuc <- data.frame(LepBuc)
LepBuc$pathway_key_reactions <- NULL
LepBuc$pathway_completeness <- NULL
LepBuc$pathway_candidate_reaction <- NULL
# rename columns
colnames(LepBuc) <- c("id", "pathway", "LepBuc")

# import the FM scores obtained from the Leptotrichia hofstadii reference sequence
LepHof <- read_delim("input_files/functional_data/LepHof_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
LepHof <- data.frame(LepHof)
LepHof$pathway_key_reactions <- NULL
LepHof$pathway_completeness <- NULL
LepHof$pathway_candidate_reaction <- NULL
# rename columns
colnames(LepHof) <- c("id", "pathway", "LepHof")

# import the FM scores obtained from the Leptotrichia hongkongensis reference sequence
LepHon <- read_delim("input_files/functional_data/LepHon_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
LepHon <- data.frame(LepHon)
LepHon$pathway_key_reactions <- NULL
LepHon$pathway_completeness <- NULL
LepHon$pathway_candidate_reaction <- NULL
# rename columns
colnames(LepHon) <- c("id", "pathway", "LepHon")


# import the FM scores obtained from the Neisseria cinerea reference sequence
NeiCin <- read_delim("input_files/functional_data/NeiCin_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
NeiCin <- data.frame(NeiCin)
NeiCin$pathway_key_reactions <- NULL
NeiCin$pathway_completeness <- NULL
NeiCin$pathway_candidate_reaction <- NULL
# rename columns
colnames(NeiCin) <- c("id", "pathway", "NeiCin")

# import the FM scores obtained from the Neisseria polysaccharea reference sequence
NeiPol <- read_delim("input_files/functional_data/NeiPol_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
NeiPol <- data.frame(NeiPol)
NeiPol$pathway_key_reactions <- NULL
NeiPol$pathway_completeness <- NULL
NeiPol$pathway_candidate_reaction <- NULL
# rename columns
colnames(NeiPol) <- c("id", "pathway", "NeiPol")

# import the FM scores obtained from the Neisseria sicca reference sequence
NeiSic <- read_delim("input_files/functional_data/NeiSic_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
NeiSic <- data.frame(NeiSic)
NeiSic$pathway_key_reactions <- NULL
NeiSic$pathway_completeness <- NULL
NeiSic$pathway_candidate_reaction <- NULL
# rename columns
colnames(NeiSic) <- c("id", "pathway", "NeiSic")

# import the FM scores obtained from the Neisseria subflava reference sequence
NeiSub <- read_delim("input_files/functional_data/NeiSub_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
NeiSub <- data.frame(NeiSub)
NeiSub$pathway_key_reactions <- NULL
NeiSub$pathway_completeness <- NULL
NeiSub$pathway_candidate_reaction <- NULL
# rename columns
colnames(NeiSub) <- c("id", "pathway", "NeiSub")

# import the FM scores obtained from the Porphyromonas asaccharolytica reference sequence
PorAsa <- read_delim("input_files/functional_data/PorAsa_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
PorAsa <- data.frame(PorAsa)
PorAsa$pathway_key_reactions <- NULL
PorAsa$pathway_completeness <- NULL
PorAsa$pathway_candidate_reaction <- NULL
# rename columns
colnames(PorAsa) <- c("id", "pathway", "PorAsa")

# import the FM scores obtained from the Prevotella jejuni reference sequence
PreJej_II <- read_delim("input_files/functional_data/PreJej_II_pathways_reactions.csv", ";", 
                        escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                                pathway_completeness = col_number(), 
                                                                pathway_completeness_perc = col_number(), 
                                                                pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
PreJej_II <- data.frame(PreJej_II)
PreJej_II$pathway_key_reactions <- NULL
PreJej_II$pathway_completeness <- NULL
PreJej_II$pathway_candidate_reaction <- NULL
# rename columns
colnames(PreJej_II) <- c("id", "pathway", "PreJej_II")

# import the FM scores obtained from the Prevotella jejuni reference sequence
PreJej <- read_delim("input_files/functional_data/PreJej_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
PreJej <- data.frame(PreJej)
PreJej$pathway_key_reactions <- NULL
PreJej$pathway_completeness <- NULL
PreJej$pathway_candidate_reaction <- NULL
# rename columns
colnames(PreJej) <- c("id", "pathway", "PreJej")
# Prevotella jejuni has two chromosomes, so here we merge the information into one column
PreJej_II$PreJej_III <- PreJej_II$PreJej_II + PreJej$PreJej

# import the FM scores obtained from the Prevotella melaninogenica reference sequence
PreMel_II <- read_delim("input_files/functional_data/PreMel_II_pathways_reactions.csv", ";", 
                        escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                                pathway_completeness = col_number(), 
                                                                pathway_completeness_perc = col_number(), 
                                                                pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
PreMel_II <- data.frame(PreMel_II)
PreMel_II$pathway_key_reactions <- NULL
PreMel_II$pathway_completeness <- NULL
PreMel_II$pathway_candidate_reaction <- NULL
# rename columns
colnames(PreMel_II) <- c("id", "pathway", "PreMel_II")

# import the FM scores obtained from the Prevotella melaninogenica reference sequence
PreMel <- read_delim("input_files/functional_data/PreMel_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
PreMel <- data.frame(PreMel)
PreMel$pathway_key_reactions <- NULL
PreMel$pathway_completeness <- NULL
PreMel$pathway_candidate_reaction <- NULL
# rename columns
colnames(PreMel) <- c("id", "pathway", "PreMel")
# Prevotella melaninogenica has two chromosomes, so here we merge the information into one column
PreMel_II$PreMel_III <- PreMel_II$PreMel_II + PreMel$PreMel

# import the FM scores obtained from the Prevotella oris reference sequence
PreOri <- read_delim("input_files/functional_data/PreOri_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
PreOri <- data.frame(PreOri)
PreOri$pathway_key_reactions <- NULL
PreOri$pathway_completeness <- NULL
PreOri$pathway_candidate_reaction <- NULL
# rename columns
colnames(PreOri) <- c("id", "pathway", "PreOri")

# import the FM scores obtained from the Prevotella ruminicola reference sequence
PreRum <- read_delim("input_files/functional_data/PreRum_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
PreRum <- data.frame(PreRum)
PreRum$pathway_key_reactions <- NULL
PreRum$pathway_completeness <- NULL
PreRum$pathway_candidate_reaction <- NULL
# rename columns
colnames(PreRum) <- c("id", "pathway", "PreRum")

# import the FM scores obtained from the Pseudoleptotrichia goodfellowii reference sequence
PseGoo <- read_delim("input_files/functional_data/PseGoo_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
PseGoo <- data.frame(PseGoo)
PseGoo$pathway_key_reactions <- NULL
PseGoo$pathway_completeness <- NULL
PseGoo$pathway_candidate_reaction <- NULL
# rename columns
colnames(PseGoo) <- c("id", "pathway", "PseGoo")

# import the FM scores obtained from the Rothia mucilaginosa reference sequence
RotMuc <- read_delim("input_files/functional_data/RotMuc_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
RotMuc <- data.frame(RotMuc)
RotMuc$pathway_key_reactions <- NULL
RotMuc$pathway_completeness <- NULL
RotMuc$pathway_candidate_reaction <- NULL
# rename columns
colnames(RotMuc) <- c("id", "pathway", "RotMuc")

# import the FM scores obtained from the Schaalia odntolytica reference sequence
SchOdo <- read_delim("input_files/functional_data/SchOdo_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
SchOdo <- data.frame(SchOdo)
SchOdo$pathway_key_reactions <- NULL
SchOdo$pathway_completeness <- NULL
SchOdo$pathway_candidate_reaction <- NULL
# rename columns
colnames(SchOdo) <- c("id", "pathway", "SchOdo")

# import the FM scores obtained from the Streptococcus australis reference sequence
StrAus <- read_delim("input_files/functional_data/StrAus_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
StrAus <- data.frame(StrAus)
StrAus$pathway_key_reactions <- NULL
StrAus$pathway_completeness <- NULL
StrAus$pathway_candidate_reaction <- NULL
# rename columns
colnames(StrAus) <- c("id", "pathway", "StrAus")


# import the FM scores obtained from the Streptococcus gordonii reference sequence
StrGor <- read_delim("input_files/functional_data/StrGor_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
StrGor <- data.frame(StrGor)
StrGor$pathway_key_reactions <- NULL
StrGor$pathway_completeness <- NULL
StrGor$pathway_candidate_reaction <- NULL
# rename columns
colnames(StrGor) <- c("id", "pathway", "StrGor")

# import the FM scores obtained from the Streptococcus koorensis reference sequence
StrKor <- read_delim("input_files/functional_data/StrKor_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
StrKor <- data.frame(StrKor)
StrKor$pathway_key_reactions <- NULL
StrKor$pathway_completeness <- NULL
StrKor$pathway_candidate_reaction <- NULL
# rename columns
colnames(StrKor) <- c("id", "pathway", "StrKor")

# import the FM scores obtained from the Streptococcus mitis reference sequence
StrMit <- read_delim("input_files/functional_data/StrMit_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
StrMit <- data.frame(StrMit)
StrMit$pathway_key_reactions <- NULL
StrMit$pathway_completeness <- NULL
StrMit$pathway_candidate_reaction <- NULL
# rename columns
colnames(StrMit) <- c("id", "pathway", "StrMit")


# import the FM scores obtained from the Streptococcus parasanguinis reference sequence
StrPar <- read_delim("input_files/functional_data/StrPar_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
StrPar <- data.frame(StrPar)
StrPar$pathway_key_reactions <- NULL
StrPar$pathway_completeness <- NULL
StrPar$pathway_candidate_reaction <- NULL
# rename columns
colnames(StrPar) <- c("id", "pathway", "StrPar")

# import the FM scores obtained from the Streptococcus pseudopneumoniae reference sequence
StrPsePneu <- read_delim("input_files/functional_data/StrPsePneu_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
StrPsePneu <- data.frame(StrPsePneu)
StrPsePneu$pathway_key_reactions <- NULL
StrPsePneu$pathway_completeness <- NULL
StrPsePneu$pathway_candidate_reaction <- NULL
# rename columns
colnames(StrPsePneu) <- c("id", "pathway", "StrPsePneu")

# import the FM scores obtained from the Veillonella atypica reference sequence
VeiAty <- read_delim("input_files/functional_data/VeiAty_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
VeiAty <- data.frame(VeiAty)
VeiAty$pathway_key_reactions <- NULL
VeiAty$pathway_completeness <- NULL
VeiAty$pathway_candidate_reaction <- NULL
# rename columns
colnames(VeiAty) <- c("id", "pathway", "VeiAty")

# import the FM scores obtained from the Veillonella dispar reference sequence
VeiDis <- read_delim("input_files/functional_data/VeiDis_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
VeiDis <- data.frame(VeiDis)
VeiDis$pathway_key_reactions <- NULL
VeiDis$pathway_completeness <- NULL
VeiDis$pathway_candidate_reaction <- NULL
# rename columns
colnames(VeiDis) <- c("id", "pathway", "VeiDis")

# import the FM scores obtained from the Veillonella rodentium reference sequence
VeiRod <- read_delim("input_files/functional_data/VeiRod_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
# clean the data frame
VeiRod <- data.frame(VeiRod)
VeiRod$pathway_key_reactions <- NULL
VeiRod$pathway_completeness <- NULL
VeiRod$pathway_candidate_reaction <- NULL
# rename columns
colnames(VeiRod) <- c("id", "pathway", "VeiRod")

# merge all tables into one
merge_species_0 <- data.frame(cbind(ActIsr, CapEnd$CapEnd, CapSpu$CapSpu, FusNuc$FusNuc, 
                                  GemHae$GemHae, GemMor$GemMor, GemSan$GemSan, HaePar$HaePar,
                                  LepBuc$LepBuc, LepHof$LepHof, LepHon$LepHon,
                                  NeiCin$NeiCin, NeiPol$NeiPol, NeiSic$NeiSic, NeiSub$NeiSub,
                                  PorAsa$PorAsa, PreJej_II$PreJej_III, PreMel_II$PreMel_III, PreOri$PreOri, 
                                  PreRum$PreRum, PseGoo$PseGoo, RotMuc$RotMuc, SchOdo$SchOdo, 
                                  StrAus$StrAus, StrGor$StrGor, StrKor$StrKor, StrMit$StrMit, 
                                  StrPar$StrPar, StrPsePneu$StrPsePneu,
                                  VeiAty$VeiAty, VeiDis$VeiDis, VeiRod$VeiRod))
merge_species_2 <- merge_species_0
merge_species_2$pathway <- NULL
# merge all pathways with same id and keep the maximum bit-score value that was obtained per species.
merge_species_2 <- ddply(merge_species_2,"id",numcolwise(max))
rownames(merge_species_2) <- merge_species_2$id
merge_species_2$id <- NULL
# transpose data frame
merge_species_3 <- data.frame(t(merge_species_2))

# rename row names with appropriate species labels
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "ActIsr", "Actinomyces israelii")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "CapEnd.CapEnd", "Capnocytophaga endodontalis")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "CapSpu.CapSpu", "Capnocytophaga sputigena")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "FusNuc.FusNuc", "Fusobacterium nucleatum")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "GemHae.GemHae", "Gemella haemolysans")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "GemMor.GemMor", "Gemella morbillorum")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "GemSan.GemSan", "Gemella sanguinis")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "HaePar.HaePar", "Haemophilus parainfluenzae")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "LepBuc.LepBuc", "Leptotrichia buccalis")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "LepHof.LepHof", "Leptotrichia hofstadii")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "LepHon.LepHon", "Leptotrichia hongkongensis")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "NeiCin.NeiCin", "Neisseria cinerea")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "NeiPol.NeiPol", "Neisseria polysaccharea")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "NeiSic.NeiSic", "Neisseria sicca")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "NeiSub.NeiSub", "Neisseria subflava")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "PorAsa.PorAsa", "Porphyromonas asaccharolytica")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "PreJej_II.PreJej_III", "Prevotella jejuni") 
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "PreMel_II.PreMel_III", "Prevotella melaninogenica")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "PreOri.PreOri", "Prevotella oris")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "PreRum.PreRum", "Prevotella ruminicola")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "PseGoo.PseGoo", "Pseudoleptotrichia goodfellowii")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "RotMuc.RotMuc", "Rothia mucilaginosa")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "SchOdo.SchOdo", "Schaalia odontolytica")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "StrAus.StrAus", "Streptococcus australis")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "StrGor.StrGor", "Streptococcus gordonii")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "StrKor.StrKor", "Streptococcus koreensis")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "StrMit.StrMit", "Streptococcus mitis")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "StrPar.StrPar", "Streptococcus parasanguinis")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "StrPsePneu.StrPsePneu", "Streptococcus pseudopneumoniae")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "VeiAty.VeiAty", "Veillonella atypica")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "VeiDis.VeiDis", "Veillonella dispar")
rownames(merge_species_3) <- str_replace(rownames(merge_species_3), "VeiRod.VeiRod", "Veillonella rodentium")

############################################################################################################
# identify FM score similarity between species with Hierarchical clustering
d_species = as.dendrogram(hclust(dist(merge_species_3, method = "canberra"), method="ward.D2"))
d_strepto_2 = dendro_data(d_species)
labs_strepto <- label(d_strepto_2)

# plot the dendrogram
strep_hier <-
  ggplot(segment(d_strepto_2)) +
  geom_segment(aes(x=x,y=y, xend=xend, yend=yend)) +
  geom_text(data=labs_strepto, aes(label=label, x=x, y=-500), size=2.5, fontface='italic') +
  ylim(-850, 2000) + ylab("                                                       Canberra distance") +
  scale_colour_manual(values=c("darkblue", "orange3")) +
  theme_bw(base_size=8) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank()) +
  coord_flip()

############################################################################################################
# find the FM score difference between healthy core and rare species biosphere
# rename the original data frame
bg_species_bcphc_0 <- merge_species_3
bg_species_rle_0 <- merge_species_3
bg_species_vst_0 <- merge_species_3

# subset the original data frame with FM scores based on the background species
# detected in healthy children, 25th abundance percentile, bcphc-normalised count data
bg_species_bcphc_h_0 <- bg_species_bcphc_0[rownames(bg_species_bcphc_0) %in% bg_species_25_bcphc_h$Species,]
bg_species_25_bcphc_h <- bg_species_25_bcphc_h[order(bg_species_25_bcphc_h$Species),]
bg_species_bcphc_h_0 <- bg_species_bcphc_h_0[order(rownames(bg_species_bcphc_h_0)),]
bg_species_bcphc_h_0$Species_type <- bg_species_25_bcphc_h$Species_type
number_species_bcphc_h <- table(bg_species_bcphc_h_0$Species_type)
# get the mean value based on core and rare species biosphere
bg_species_bcphc_h_1 <- ddply(bg_species_bcphc_h_0,"Species_type",numcolwise(set_statistics))
rownames(bg_species_bcphc_h_1) <- bg_species_bcphc_h_1$Species_type
bg_species_bcphc_h_1$Species_type <- NULL
# remove all rows that sum to zero
bg_species_bcphc_h_2 <- bg_species_bcphc_h_1[, colSums(bg_species_bcphc_h_1) > 0]
# transpose the data frame
bg_species_bcphc_h_3 <- data.frame(t(bg_species_bcphc_h_2))
# find those rows that are different between core and rare species biosphere
bg_species_bcphc_h_4 <- subset(bg_species_bcphc_h_3, background_core != background_rare)
# add a column that quantifies the difference
bg_species_bcphc_h_4$diff <- abs(bg_species_bcphc_h_4$background_core - bg_species_bcphc_h_4$background_rare)
# keep only the rows with more than 10% difference between core and rare species biosphere
bg_species_bcphc_h_5 <- subset(bg_species_bcphc_h_4, diff > set_diff) 
bg_species_bcphc_h_5$diff <- NULL
# clean the row names
rownames(bg_species_bcphc_h_5) <- str_replace(rownames(bg_species_bcphc_h_5), "^X", "")
rownames(bg_species_bcphc_h_5) <- str_replace_all(rownames(bg_species_bcphc_h_5), "\\.", "-")
colnames(bg_species_bcphc_h_5) <- c("Core species", "Rare species")
bg_species_bcphc_h_5$Pathway <- rownames(bg_species_bcphc_h_5)
rownames(bg_species_bcphc_h_5) <- NULL
# convert the data frame into the long format
bg_species_bcphc_h_6 <- gather(bg_species_bcphc_h_5, key="Type", value="FMs", -"Pathway")
# add corresponding meta data to the table
bg_species_bcphc_h_6$State <- "Healthy"
bg_species_bcphc_h_6$Normalisation <- "BCPHC"



# subset the original data frame with FM scores based on the background species
# detected in healthy children, 25th abundance percentile, RLE-normalised count data
bg_species_rle_h_0 <- bg_species_rle_0[rownames(bg_species_rle_0) %in% bg_species_25_rle_h$Species,]
bg_species_25_rle_h <- bg_species_25_rle_h[order(bg_species_25_rle_h$Species),]
bg_species_rle_h_0 <- bg_species_rle_h_0[order(rownames(bg_species_rle_h_0)),]
bg_species_rle_h_0$Species_type <- bg_species_25_rle_h$Species_type
number_species_rle_h <- table(bg_species_rle_h_0$Species_type)
# get the mean value based on core and rare species biosphere
bg_species_rle_h_1 <- ddply(bg_species_rle_h_0,"Species_type",numcolwise(set_statistics))
rownames(bg_species_rle_h_1) <- bg_species_rle_h_1$Species_type
bg_species_rle_h_1$Species_type <- NULL
# remove all rows that sum to zero
bg_species_rle_h_2 <- bg_species_rle_h_1[, colSums(bg_species_rle_h_1) > 0]
# transpose the data frame
bg_species_rle_h_3 <- data.frame(t(bg_species_rle_h_2))
# find those rows that are different between core and rare species biosphere
bg_species_rle_h_4 <- subset(bg_species_rle_h_3, background_core != background_rare)
# add a column that quantifies the difference
bg_species_rle_h_4$diff <- abs(bg_species_rle_h_4$background_core - bg_species_rle_h_4$background_rare)
# keep only the rows with more than 10% difference between core and rare species biosphere
bg_species_rle_h_5 <- subset(bg_species_rle_h_4, diff > set_diff) 
bg_species_rle_h_5$diff <- NULL
# clean the row names
rownames(bg_species_rle_h_5) <- str_replace(rownames(bg_species_rle_h_5), "^X", "")
rownames(bg_species_rle_h_5) <- str_replace_all(rownames(bg_species_rle_h_5), "\\.", "-")
colnames(bg_species_rle_h_5) <- c("Core species", "Rare species")
bg_species_rle_h_5$Pathway <- rownames(bg_species_rle_h_5)
rownames(bg_species_rle_h_5) <- NULL
# convert the data frame into the long format
bg_species_rle_h_6 <- gather(bg_species_rle_h_5, key="Type", value="FMs", -"Pathway")
# add corresponding meta data to the table
bg_species_rle_h_6$State <- "Healthy"
bg_species_rle_h_6$Normalisation <- "RLE"


# subset the original data frame with FM scores based on the background species
# detected in healthy children, 25th abundance percentile, VST-normalised count data
bg_species_vst_h_0 <- bg_species_vst_0[rownames(bg_species_vst_0) %in% bg_species_25_vst_h$Species,]
bg_species_25_vst_h <- bg_species_25_vst_h[order(bg_species_25_vst_h$Species),]
bg_species_vst_h_0 <- bg_species_vst_h_0[order(rownames(bg_species_vst_h_0)),]
bg_species_vst_h_0$Species_type <- bg_species_25_vst_h$Species_type
number_species_vst_h <- table(bg_species_vst_h_0$Species_type)
# get the mean value based on core and rare species biosphere
bg_species_vst_h_1 <- ddply(bg_species_vst_h_0,"Species_type",numcolwise(set_statistics))
rownames(bg_species_vst_h_1) <- bg_species_vst_h_1$Species_type
bg_species_vst_h_1$Species_type <- NULL
# remove all rows that sum to zero
bg_species_vst_h_2 <- bg_species_vst_h_1[, colSums(bg_species_vst_h_1) > 0]
# transpose the data frame
bg_species_vst_h_3 <- data.frame(t(bg_species_vst_h_2))
# find those rows that are different between core and rare species biosphere
bg_species_vst_h_4 <- subset(bg_species_vst_h_3, background_core != background_rare)
# add a column that quantifies the difference
bg_species_vst_h_4$diff <- abs(bg_species_vst_h_4$background_core - bg_species_vst_h_4$background_rare)
# keep only the rows with more than 10% difference between core and rare species biosphere
bg_species_vst_h_5 <- subset(bg_species_vst_h_4, diff > set_diff) 
bg_species_vst_h_5$diff <- NULL
# clean the row names
rownames(bg_species_vst_h_5) <- str_replace(rownames(bg_species_vst_h_5), "^X", "")
rownames(bg_species_vst_h_5) <- str_replace_all(rownames(bg_species_vst_h_5), "\\.", "-")
colnames(bg_species_vst_h_5) <- c("Core species", "Rare species")
bg_species_vst_h_5$Pathway <- rownames(bg_species_vst_h_5)
rownames(bg_species_vst_h_5) <- NULL
# convert the data frame into the long format
bg_species_vst_h_6 <- gather(bg_species_vst_h_5, key="Type", value="FMs", -"Pathway")
# add corresponding meta data to the table
bg_species_vst_h_6$State <- "Healthy"
bg_species_vst_h_6$Normalisation <- "VST"

# find the matching elements between normalisation methods
Pathway_h <- Reduce(intersect, list(bg_species_vst_h_6$Pathway, bg_species_rle_h_6$Pathway, bg_species_bcphc_h_6$Pathway))

# merge all data frames
bg_species_h_6 <- data.frame(rbind(bg_species_bcphc_h_6, bg_species_rle_h_6, bg_species_vst_h_6))
# subset the dataframe to display only items that match between normlisation strategies
bg_species_h_6_short <- bg_species_h_6[bg_species_h_6$Pathway %in% Pathway_h,]

# exclude 'shadow pathways' that are identical but have different metacyc IDs
exclude_pwy <- c("PWY-7969", "PWY-7968", "PWY-7967", "PWY-7301", "PWY-7966", "PWY3O-450")
bg_species_h_6_short <- bg_species_h_6_short[!bg_species_h_6_short$Pathway %in% exclude_pwy,]

# add an extra column with increasing numbers
bg_species_h_6_short$Pathway_as_number <- rep(seq.int(1,(nrow(bg_species_h_6_short) / 6)),6)

# plot the heatmap of FM score differences between core and rare species biosphere
FM_plot_H <- ggplot(bg_species_h_6_short) +
  geom_tile(aes(x=Type, y=Pathway_as_number, fill=FMs), 
            width=0.9, height=1, size=0.1, colour="black") +
  theme_bw(base_size=10) +
  theme(axis.text.y = element_text(size = 5), 
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(angle=45, size=10, hjust=1),
        strip.background =element_rect(fill="white"),
        strip.text = element_text(size=10),
        legend.title = element_blank(), legend.position = "top",
        panel.grid = element_blank()) +
  ylab("Metabolic pathways") + xlab(" ") + 
  facet_wrap(~State ~Normalisation , nrow=1, scales = "free_y") +
  scale_y_continuous(breaks = c(seq.int(from=5,to=105, by=5)), limits = c(1,110)) +
  scale_fill_gradientn(colours = c("white", "powderblue", "steelblue4", 
                                   "lightgreen", "forestgreen", "yellow", "orange"))

# print out a data frame with the metabolic pathway information per row.
write.csv(FM_plot_H$data, "Metabolic_pathway.csv")

############################################################################################################
# Find differences between healthy core and rare species FM scores and PAO1 FM score
# Therefore, the CF FM scores per core and rare species biosphere have to be generated

# subset the original data frame with FM scores based on the background species
# detected in CF children, 25th abundance percentile, bcphc-normalised count data
bg_species_bcphc_cf_0 <- bg_species_bcphc_0[rownames(bg_species_bcphc_0) %in% bg_species_25_bcphc_cf$Species,]
bg_species_25_bcphc_cf <- bg_species_25_bcphc_cf[order(bg_species_25_bcphc_cf$Species),]
bg_species_bcphc_cf_0 <- bg_species_bcphc_cf_0[order(rownames(bg_species_bcphc_cf_0)),]
bg_species_bcphc_cf_0$Species_type <- bg_species_25_bcphc_cf$Species_type
number_species_bcphc_cf <- table(bg_species_bcphc_cf_0$Species_type)
# get the mean value based on core and rare species biosphere
bg_species_bcphc_cf_1 <- ddply(bg_species_bcphc_cf_0,"Species_type",numcolwise(set_statistics))
rownames(bg_species_bcphc_cf_1) <- bg_species_bcphc_cf_1$Species_type
bg_species_bcphc_cf_1$Species_type <- NULL
# remove all rows that sum to zero
bg_species_bcphc_cf_2 <- bg_species_bcphc_cf_1[, colSums(bg_species_bcphc_cf_1) > 0]
# transpose the data frame
bg_species_bcphc_cf_3 <- data.frame(t(bg_species_bcphc_cf_2))
# clean up row names
rownames(bg_species_bcphc_cf_3) <- str_replace(rownames(bg_species_bcphc_cf_3), "^X", "")
rownames(bg_species_bcphc_cf_3) <- str_replace_all(rownames(bg_species_bcphc_cf_3), "\\.", "-")
colnames(bg_species_bcphc_cf_3) <- c("Core species", "Rare species")
bg_species_bcphc_cf_3$Pathway <- rownames(bg_species_bcphc_cf_3)
rownames(bg_species_bcphc_cf_3) <- NULL
# convert data frame to long format
bg_species_bcphc_cf_6 <- gather(bg_species_bcphc_cf_3, key="Type", value="FMs", -"Pathway")
# add meta data information
bg_species_bcphc_cf_6$State <- "CF"
bg_species_bcphc_cf_6$Normalisation <- "BCPHC"



# subset the original data frame with FM scores based on the background species
# detected in CF children, 25th abundance percentile, rle-normalised count data
bg_species_rle_cf_0 <- bg_species_rle_0[rownames(bg_species_rle_0) %in% bg_species_25_rle_cf$Species,]
bg_species_25_rle_cf <- bg_species_25_rle_cf[order(bg_species_25_rle_cf$Species),]
bg_species_rle_cf_0 <- bg_species_rle_cf_0[order(rownames(bg_species_rle_cf_0)),]
bg_species_rle_cf_0$Species_type <- bg_species_25_rle_cf$Species_type
number_species_rle_cf <- table(bg_species_rle_cf_0$Species_type)
# get the mean value based on core and rare species biosphere
bg_species_rle_cf_1 <- ddply(bg_species_rle_cf_0,"Species_type",numcolwise(set_statistics))
rownames(bg_species_rle_cf_1) <- bg_species_rle_cf_1$Species_type
bg_species_rle_cf_1$Species_type <- NULL
# remove all rows that sum to zero
bg_species_rle_cf_2 <- bg_species_rle_cf_1[, colSums(bg_species_rle_cf_1) > 0]
# transpose the data frame
bg_species_rle_cf_3 <- data.frame(t(bg_species_rle_cf_2))
# find those rows that are different between core and rare species biosphere
rownames(bg_species_rle_cf_3) <- str_replace(rownames(bg_species_rle_cf_3), "^X", "")
rownames(bg_species_rle_cf_3) <- str_replace_all(rownames(bg_species_rle_cf_3), "\\.", "-")
colnames(bg_species_rle_cf_3) <- c("Core species", "Rare species")
bg_species_rle_cf_3$Pathway <- rownames(bg_species_rle_cf_3)
rownames(bg_species_rle_cf_3) <- NULL
# convert data frame to long format
bg_species_rle_cf_6 <- gather(bg_species_rle_cf_3, key="Type", value="FMs", -"Pathway")
# add meta data information
bg_species_rle_cf_6$State <- "CF"
bg_species_rle_cf_6$Normalisation <- "RLE"



# subset the original data frame with FM scores based on the background species
# detected in CF children, 25th abundance percentile, vst-normalised count data
bg_species_vst_cf_0 <- bg_species_vst_0[rownames(bg_species_vst_0) %in% bg_species_25_vst_cf$Species,]
bg_species_25_vst_cf <- bg_species_25_vst_cf[order(bg_species_25_vst_cf$Species),]
bg_species_vst_cf_0 <- bg_species_vst_cf_0[order(rownames(bg_species_vst_cf_0)),]
bg_species_vst_cf_0$Species_type <- bg_species_25_vst_cf$Species_type
number_species_vst_cf <- table(bg_species_vst_cf_0$Species_type)
# get the mean value based on core and rare species biosphere
bg_species_vst_cf_1 <- ddply(bg_species_vst_cf_0,"Species_type",numcolwise(set_statistics))
rownames(bg_species_vst_cf_1) <- bg_species_vst_cf_1$Species_type
bg_species_vst_cf_1$Species_type <- NULL
# remove all rows that sum to zero
bg_species_vst_cf_2 <- bg_species_vst_cf_1[, colSums(bg_species_vst_cf_1) > 0]
# transpose the data frame
bg_species_vst_cf_3 <- data.frame(t(bg_species_vst_cf_2))
# find those rows that are different between core and rare species biosphere
rownames(bg_species_vst_cf_3) <- str_replace(rownames(bg_species_vst_cf_3), "^X", "")
rownames(bg_species_vst_cf_3) <- str_replace_all(rownames(bg_species_vst_cf_3), "\\.", "-")
# in this case, we only add one column name, because there are no rare species 
colnames(bg_species_vst_cf_3) <- c("Core species")
bg_species_vst_cf_3$Pathway <- rownames(bg_species_vst_cf_3)
rownames(bg_species_vst_cf_3) <- NULL
# convert data frame to long format
bg_species_vst_cf_6 <- gather(bg_species_vst_cf_3, key="Type", value="FMs", -"Pathway")
# add meta data information
bg_species_vst_cf_6$State <- "CF"
bg_species_vst_cf_6$Normalisation <- "VST"

# find matching pathways between CF normalisations
Pathway_cf <- Reduce(intersect, list(bg_species_vst_cf_6$Pathway, bg_species_rle_cf_6$Pathway, bg_species_bcphc_cf_6$Pathway))
# merge the normalised CF data frames
bg_species_cf_6 <- data.frame(rbind(bg_species_bcphc_cf_6, bg_species_rle_cf_6, bg_species_vst_cf_6))
# extract the overlapping pathways
bg_species_cf_6_short <- bg_species_cf_6[bg_species_cf_6$Pathway %in% Pathway_h,]

# Prepare PAO1 dataframe
Pao1_df <- data.frame(cbind(Pao1$id))
# Add meta data
Pao1_df$Type <- "PAO1"
Pao1_df$FM <- Pao1$Pao1
Pao1_df$State <- "Both"
Pao1_df$Normalisation <- "All"

# Compare CF, BCPHC_normalised data with PAO1 data
Pao1_df_bcphc_cf <- Pao1_df
colnames(Pao1_df_bcphc_cf) <- colnames(bg_species_bcphc_cf_6)
Pao_bg_species_bcphc_cf_6 <- data.frame(rbind(bg_species_bcphc_cf_6, Pao1_df_bcphc_cf))
# obtained the square-rooted FM score per pathway by adding a pseudo count of 1
Pao_bg_species_bcphc_cf_6$FMs <- Pao_bg_species_bcphc_cf_6$FMs + 1
Pao_bg_species_bcphc_cf_6$FMs_sqrt <- sqrt(Pao_bg_species_bcphc_cf_6$FMs)
# subset the data frame into core and rare species biosphere
Pao_bg_species_bcphc_cf_6_core <- subset(Pao_bg_species_bcphc_cf_6, Type != "Rare species")
Pao_bg_species_bcphc_cf_6_rare <- subset(Pao_bg_species_bcphc_cf_6, Type != "Core species")
# approach a statistical comparison between core species and PAO1 and store output in data frame
Pao_bcphc_cf_core <- rcompanion::wilcoxonR(x=Pao_bg_species_bcphc_cf_6_core$FMs_sqrt, 
                                           g=Pao_bg_species_bcphc_cf_6_core$Type, ci=TRUE)
Pao_bcphc_cf_core$Name <- "Core-BCPHC"
Pao_bcphc_cf_core$State <- "CF"
# approach a statistical comparison between rare species and PAO1 and store output in data frame
Pao_bcphc_cf_rare <- rcompanion::wilcoxonR(x=Pao_bg_species_bcphc_cf_6_rare$FMs_sqrt, 
                                           g=Pao_bg_species_bcphc_cf_6_rare$Type, ci=TRUE)
Pao_bcphc_cf_rare$Name <- "Rare-BCPHC"
Pao_bcphc_cf_rare$State <- "CF"


# Compare Healthy, BCPHC_normalised data with PAO1 data
Pao1_df_bcphc_h <- Pao1_df
colnames(Pao1_df_bcphc_h) <- colnames(bg_species_bcphc_h_6)
Pao_bg_species_bcphc_h_6 <- data.frame(rbind(bg_species_bcphc_h_6, Pao1_df_bcphc_h))
# obtained the square-rooted FM score per pathway by adding a pseudo count of 1
Pao_bg_species_bcphc_h_6$FMs <- Pao_bg_species_bcphc_h_6$FMs + 1
Pao_bg_species_bcphc_h_6$FMs_sqrt <- sqrt(Pao_bg_species_bcphc_h_6$FMs)
# subset the data frame into core and rare species biosphere
Pao_bg_species_bcphc_h_6_core <- subset(Pao_bg_species_bcphc_h_6, Type != "Rare species")
Pao_bg_species_bcphc_h_6_rare <- subset(Pao_bg_species_bcphc_h_6, Type != "Core species")
# approach a statistical comparison between core species and PAO1 and store output in data frame
Pao_bcphc_h_core <- rcompanion::wilcoxonR(x=Pao_bg_species_bcphc_h_6_core$FMs_sqrt, 
                                          g=Pao_bg_species_bcphc_h_6_core$Type, ci=TRUE)
Pao_bcphc_h_core$Name <- "Core-BCPHC"
Pao_bcphc_h_core$State <- "Healthy"
# approach a statistical comparison between rare species and PAO1 and store output in data frame
Pao_bcphc_h_rare <- rcompanion::wilcoxonR(x=Pao_bg_species_bcphc_h_6_rare$FMs_sqrt, 
                                          g=Pao_bg_species_bcphc_h_6_rare$Type, ci=TRUE)
Pao_bcphc_h_rare$Name <- "Rare-BCPHC"
Pao_bcphc_h_rare$State <- "Healthy"


# Compare CF, RLE_normalised data with PAO1 data
Pao1_df_rle_cf <- Pao1_df
colnames(Pao1_df_rle_cf) <- colnames(bg_species_rle_cf_6)
Pao_bg_species_rle_cf_6 <- data.frame(rbind(bg_species_rle_cf_6, Pao1_df_rle_cf))
# obtained the square-rooted FM score per pathway by adding a pseudo count of 1
Pao_bg_species_rle_cf_6$FMs <- Pao_bg_species_rle_cf_6$FMs + 1
Pao_bg_species_rle_cf_6$FMs_sqrt <- sqrt(Pao_bg_species_rle_cf_6$FMs)
# subset the data frame into core and rare species biosphere
Pao_bg_species_rle_cf_6_core <- subset(Pao_bg_species_rle_cf_6, Type != "Rare species")
Pao_bg_species_rle_cf_6_rare <- subset(Pao_bg_species_rle_cf_6, Type != "Core species")
# approach a statistical comparison between core species and PAO1 and store output in data frame
Pao_rle_cf_core <- rcompanion::wilcoxonR(x=Pao_bg_species_rle_cf_6_core$FMs_sqrt, 
                                           g=Pao_bg_species_rle_cf_6_core$Type, ci=TRUE)
Pao_rle_cf_core$Name <- "Core-RLE"
Pao_rle_cf_core$State <- "CF"
# approach a statistical comparison between rare species and PAO1 and store output in data frame
Pao_rle_cf_rare <- rcompanion::wilcoxonR(x=Pao_bg_species_rle_cf_6_rare$FMs_sqrt, 
                                           g=Pao_bg_species_rle_cf_6_rare$Type, ci=TRUE)
Pao_rle_cf_rare$Name <- "Rare-RLE"
Pao_rle_cf_rare$State <- "CF"


# Compare healthy, RLE_normalised data with PAO1 data
Pao1_df_rle_h <- Pao1_df
colnames(Pao1_df_rle_h) <- colnames(bg_species_rle_h_6)
Pao_bg_species_rle_h_6 <- data.frame(rbind(bg_species_rle_h_6, Pao1_df_rle_h))
# obtained the square-rooted FM score per pathway by adding a pseudo count of 1
Pao_bg_species_rle_h_6$FMs <- Pao_bg_species_rle_h_6$FMs + 1
Pao_bg_species_rle_h_6$FMs_sqrt <- sqrt(Pao_bg_species_rle_h_6$FMs)
# subset the data frame into core and rare species biosphere
Pao_bg_species_rle_h_6_core <- subset(Pao_bg_species_rle_h_6, Type != "Rare species")
Pao_bg_species_rle_h_6_rare <- subset(Pao_bg_species_rle_h_6, Type != "Core species")
# approach a statistical comparison between core species and PAO1 and store output in data frame
Pao_rle_h_core <- rcompanion::wilcoxonR(x=Pao_bg_species_rle_h_6_core$FMs_sqrt, 
                                          g=Pao_bg_species_rle_h_6_core$Type, ci=TRUE)
Pao_rle_h_core$Name <- "Core-RLE"
Pao_rle_h_core$State <- "Healthy"
# approach a statistical comparison between rare species and PAO1 and store output in data frame
Pao_rle_h_rare <- rcompanion::wilcoxonR(x=Pao_bg_species_rle_h_6_rare$FMs_sqrt, 
                                          g=Pao_bg_species_rle_h_6_rare$Type, ci=TRUE)
Pao_rle_h_rare$Name <- "Rare-RLE"
Pao_rle_h_rare$State <- "Healthy"


# Compare CF, VST_normalised data with PAO1 data
Pao1_df_vst_cf <- Pao1_df
colnames(Pao1_df_vst_cf) <- colnames(bg_species_vst_cf_6)
Pao_bg_species_vst_cf_6 <- data.frame(rbind(bg_species_vst_cf_6, Pao1_df_vst_cf))
# obtained the square-rooted FM score per pathway by adding a pseudo count of 1
Pao_bg_species_vst_cf_6$FMs <- Pao_bg_species_vst_cf_6$FMs + 1
Pao_bg_species_vst_cf_6$FMs_sqrt <- sqrt(Pao_bg_species_vst_cf_6$FMs)
# subset the data frame into core and rare species biosphere
Pao_bg_species_vst_cf_6_core <- subset(Pao_bg_species_vst_cf_6, Type != "Rare species")
Pao_bg_species_vst_cf_6_rare <- subset(Pao_bg_species_vst_cf_6, Type != "Core species")
# approach a statistical comparison between core species and PAO1 and store output in data frame
Pao_vst_cf_core <- rcompanion::wilcoxonR(x=Pao_bg_species_vst_cf_6_core$FMs_sqrt, 
                                         g=Pao_bg_species_vst_cf_6_core$Type, ci=TRUE)
Pao_vst_cf_core$Name <- "Core-VST"
Pao_vst_cf_core$State <- "CF"
#  a statistical comparison is not approached, because there are not
# enough observations for rare species in VST-normalised data
#Pao_vst_cf_rare <- rcompanion::wilcoxonR(x=Pao_bg_species_vst_cf_6_rare$FMs_sqrt, 
   #                                      g=Pao_bg_species_vst_cf_6_rare$Type, ci=TRUE)
#Pao_vst_cf_rare$Name <- "Rare-VST"
#Pao_vst_cf_rare$State <- "CF"


# Compare healthy, VST_normalised data with PAO1 data
Pao1_df_vst_h <- Pao1_df
colnames(Pao1_df_vst_h) <- colnames(bg_species_vst_h_6)
Pao_bg_species_vst_h_6 <- data.frame(rbind(bg_species_vst_h_6, Pao1_df_vst_h))
# obtained the square-rooted FM score per pathway by adding a pseudo count of 1
Pao_bg_species_vst_h_6$FMs <- Pao_bg_species_vst_h_6$FMs + 1
Pao_bg_species_vst_h_6$FMs_sqrt <- sqrt(Pao_bg_species_vst_h_6$FMs)
# subset the data frame into core and rare species biosphere
Pao_bg_species_vst_h_6_core <- subset(Pao_bg_species_vst_h_6, Type != "Rare species")
Pao_bg_species_vst_h_6_rare <- subset(Pao_bg_species_vst_h_6, Type != "Core species")
# approach a statistical comparison between core species and PAO1 and store output in data frame
Pao_vst_h_core <- rcompanion::wilcoxonR(x=Pao_bg_species_vst_h_6_core$FMs_sqrt, 
                                        g=Pao_bg_species_vst_h_6_core$Type, ci=TRUE)
Pao_vst_h_core$Name <- "Core-VST"
Pao_vst_h_core$State <- "Healthy"
# approach a statistical comparison between rare species and PAO1 and store output in data frame
Pao_vst_h_rare <- rcompanion::wilcoxonR(x=Pao_bg_species_vst_h_6_rare$FMs_sqrt, 
                                        g=Pao_bg_species_vst_h_6_rare$Type, ci=TRUE)
Pao_vst_h_rare$Name <- "Rare-VST"
Pao_vst_h_rare$State <- "Healthy"


# merge all data frames into one dataframe
merge_all <- data.frame(rbind(Pao_bcphc_cf_core, Pao_bcphc_h_core, Pao_rle_cf_core, Pao_rle_h_core,
                   Pao_vst_cf_core, Pao_vst_h_core, Pao_bcphc_cf_rare, Pao_bcphc_h_rare,
                   Pao_rle_cf_rare, Pao_rle_h_rare,Pao_vst_h_rare))

# plot the graph
PAO1_stats_plot <-
  ggplot(merge_all) +
  geom_linerange(aes(y=Name, xmin=upper.ci, xmax=lower.ci), size=0.3) +
  geom_point(aes(y=Name, x=r), size=0.5) + facet_wrap(~State) + 
  geom_vline(aes(xintercept=0), colour="red", size=0.3) +
  theme_bw(base_size = 5) + theme(panel.grid = element_blank(),
                                   strip.background = element_rect(fill="white")) + 
  xlab("Effect size r") + ylab(" ") +
  scale_x_continuous(breaks = c(-0.1, 0, 0.1), limits = c(-0.15, 0.15)) 


############################################################################################################
# Adhesion pattern analysis
# Import data sets
# Import Actinomyces israelii adhesin results
ActIsr_adhesion <- read_delim("input_files/functional_data/ActIsr_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
ActIsr_adhesion <- data.frame(ActIsr_adhesion)
# rename columns
colnames(ActIsr_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove non-significant matches
ActIsr_adhesion <- subset(ActIsr_adhesion, Evalue < 0.01)
# extract all significant scores
ActIsr_adhesion <- subset(ActIsr_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
ActIsr_adhesion <- ddply(ActIsr_adhesion,"Adhesion",numcolwise(max))
# add species name to column
ActIsr_adhesion$Species <- "Actinomyces israelii"


# Import Capnocytophaga endodontalis adhesin results
CapEnd_adhesion <- read_delim("input_files/functional_data/CapEnd_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
CapEnd_adhesion <- data.frame(CapEnd_adhesion)
# rename columns
colnames(CapEnd_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove non-significant matches
CapEnd_adhesion <- subset(CapEnd_adhesion, Evalue < 0.01)
# extract all significant scores
CapEnd_adhesion <- subset(CapEnd_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
CapEnd_adhesion <- ddply(CapEnd_adhesion,"Adhesion",numcolwise(max))
# add species name to column
CapEnd_adhesion$Species <- "Capnocytophaga endodontalis"


# Import Capnocytophaga sputigena adhesin results
CapSpu_adhesion <- read_delim("input_files/functional_data/CapSpu_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
CapSpu_adhesion <- data.frame(CapSpu_adhesion)
# rename columns
colnames(CapSpu_adhesion) <- c("Adhesion", "Score", "Evalue")
# keep significant matches
CapSpu_adhesion <- subset(CapSpu_adhesion, Evalue < 0.01)
# extract all significant scores
CapSpu_adhesion <- subset(CapSpu_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
CapSpu_adhesion <- ddply(CapSpu_adhesion,"Adhesion",numcolwise(max))
# add species name to column
CapSpu_adhesion$Species <- "Capnocytophaga sputigena"


# Import Fusobacterium nucleatum adhesin results
FusNuc_adhesion <- read_delim("input_files/functional_data/FusNuc_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
FusNuc_adhesion <- data.frame(FusNuc_adhesion)
# rename columns
colnames(FusNuc_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
FusNuc_adhesion$Evalue <- str_replace(FusNuc_adhesion$Evalue,";", "")
FusNuc_adhesion$Evalue <- as.numeric(FusNuc_adhesion$Evalue)
# keep significant matches
FusNuc_adhesion <- subset(FusNuc_adhesion, Evalue < 0.01)
# extract all significant scores
FusNuc_adhesion <- subset(FusNuc_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
FusNuc_adhesion <- ddply(FusNuc_adhesion,"Adhesion",numcolwise(max))
# add species name to column
FusNuc_adhesion$Species <- "Fusobacterium nucleatum"


# Import Gemella haemolysans adhesin results
GemHae_adhesion <- read_delim("input_files/functional_data/GemHae_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
GemHae_adhesion <- data.frame(GemHae_adhesion)
# rename columns
colnames(GemHae_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
GemHae_adhesion$Evalue <- str_replace(GemHae_adhesion$Evalue,";", "")
GemHae_adhesion$Evalue <- as.numeric(GemHae_adhesion$Evalue)
# keep significant matches
GemHae_adhesion <- subset(GemHae_adhesion, Evalue < 0.01)
# extract all significant scores
GemHae_adhesion <- subset(GemHae_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
GemHae_adhesion <- ddply(GemHae_adhesion,"Adhesion",numcolwise(max))
# add species name to column
GemHae_adhesion$Species <- "Gemella haemolysans"


# Import Gemella morbillorum adhesin results
GemMor_adhesion <- read_delim("input_files/functional_data/GemMor_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
# rename columns
GemMor_adhesion <- data.frame(GemMor_adhesion)
colnames(GemMor_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
GemMor_adhesion$Evalue <- str_replace(GemMor_adhesion$Evalue,";", "")
GemMor_adhesion$Evalue <- as.numeric(GemMor_adhesion$Evalue)
# keep significant matches
GemMor_adhesion <- subset(GemMor_adhesion, Evalue < 0.01)
# extract all significant scores
GemMor_adhesion <- subset(GemMor_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
GemMor_adhesion <- ddply(GemMor_adhesion,"Adhesion",numcolwise(max))
# add species name to column
GemMor_adhesion$Species <- "Gemella morbillorum"


# Import Gemella sanguinis adhesin results
GemSan_adhesion <- read_delim("input_files/functional_data/GemSan_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
# rename columns
GemSan_adhesion <- data.frame(GemSan_adhesion)
colnames(GemSan_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
GemSan_adhesion$Evalue <- str_replace(GemSan_adhesion$Evalue,";", "")
GemSan_adhesion$Evalue <- as.numeric(GemSan_adhesion$Evalue)
# keep significant matches
GemSan_adhesion <- subset(GemSan_adhesion, Evalue < 0.01)
# extract all significant scores
GemSan_adhesion <- subset(GemSan_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
GemSan_adhesion <- ddply(GemSan_adhesion,"Adhesion",numcolwise(max))
# add species name to column
GemSan_adhesion$Species <- "Gemella sanguinis"


# Import Haemophilus parainfluenzae adhesin results
HaePar_adhesion <- read_delim("input_files/functional_data/HaePar_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
HaePar_adhesion <- data.frame(HaePar_adhesion)
# rename columns
colnames(HaePar_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
HaePar_adhesion$Evalue <- str_replace(HaePar_adhesion$Evalue,";", "")
HaePar_adhesion$Evalue <- as.numeric(HaePar_adhesion$Evalue)
# keep significant matches
HaePar_adhesion <- subset(HaePar_adhesion, Evalue < 0.01)
# extract all significant scores
HaePar_adhesion <- subset(HaePar_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
HaePar_adhesion <- ddply(HaePar_adhesion,"Adhesion",numcolwise(max))
# add species name to column
HaePar_adhesion$Species <- "Haemophilus parainfluenzae"


# Import Leptotrichia hofstadii adhesin results
LepHof_adhesion <- read_delim("input_files/functional_data/LepHof_adhesion.csv", 
                              delim = ";", escape_double = FALSE, col_names = FALSE, 
                              trim_ws = TRUE)
LepHof_adhesion <- data.frame(LepHof_adhesion)
# rename columns
colnames(LepHof_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
LepHof_adhesion$Evalue <- str_replace(LepHof_adhesion$Evalue,";", "")
LepHof_adhesion$Evalue <- as.numeric(LepHof_adhesion$Evalue)
# keep significant matches
LepHof_adhesion <- subset(LepHof_adhesion, Evalue < 0.01)
# extract all significant scores
LepHof_adhesion <- subset(LepHof_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
LepHof_adhesion <- ddply(LepHof_adhesion,"Adhesion",numcolwise(max))
# add species name to column
LepHof_adhesion$Species <- "Leptotrichia hofstadii"


# Import Leptotrichia hongkongensis adhesin results
LepHon_adhesion <- read_delim("input_files/functional_data/LepHon_adhesion.csv", 
                              delim = ";", escape_double = FALSE, col_names = FALSE, 
                              trim_ws = TRUE)
LepHon_adhesion <- data.frame(LepHon_adhesion)
# rename columns
colnames(LepHon_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
LepHon_adhesion$Evalue <- str_replace(LepHon_adhesion$Evalue,";", "")
LepHon_adhesion$Evalue <- as.numeric(LepHon_adhesion$Evalue)
# keep significant matches
LepHon_adhesion <- subset(LepHon_adhesion, Evalue < 0.01)
LepHon_adhesion <- subset(LepHon_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
LepHon_adhesion <- ddply(LepHon_adhesion,"Adhesion",numcolwise(max))
# add species name to column
LepHon_adhesion$Species <- "Leptotrichia hongkongensis"


# Import Neisseria cinerea adhesin results
NeiCin_adhesion <- read_delim("input_files/functional_data/NeiCin_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
NeiCin_adhesion <- data.frame(NeiCin_adhesion)
# rename columns
colnames(NeiCin_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
NeiCin_adhesion$Evalue <- str_replace(NeiCin_adhesion$Evalue,";", "")
NeiCin_adhesion$Evalue <- as.numeric(NeiCin_adhesion$Evalue)
# keep significant matches
NeiCin_adhesion <- subset(NeiCin_adhesion, Evalue < 0.01)
NeiCin_adhesion <- subset(NeiCin_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
NeiCin_adhesion <- ddply(NeiCin_adhesion,"Adhesion",numcolwise(max))
# add species name to column
NeiCin_adhesion$Species <- "Neisseria cinerea"


# Import Neisseria polysaccharea adhesin results
NeiPol_adhesion <- read_delim("input_files/functional_data/NeiPol_adhesion.csv", 
                    delim = ";", escape_double = FALSE, col_names = FALSE, 
                    trim_ws = TRUE)
NeiPol_adhesion <- data.frame(NeiPol_adhesion)
# rename columns
colnames(NeiPol_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
NeiPol_adhesion$Evalue <- str_replace(NeiPol_adhesion$Evalue,";", "")
NeiPol_adhesion$Evalue <- as.numeric(NeiPol_adhesion$Evalue)
# keep significant matches
NeiPol_adhesion <- subset(NeiPol_adhesion, Evalue < 0.01)
NeiPol_adhesion <- subset(NeiPol_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
NeiPol_adhesion <- ddply(NeiPol_adhesion,"Adhesion",numcolwise(max))
# add species name to column
NeiPol_adhesion$Species <- "Neisseria polysaccharea"


# Import Neisseria sicca adhesin results
NeiSic_adhesion <- read_delim("input_files/functional_data/NeiSic_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
NeiSic_adhesion <- data.frame(NeiSic_adhesion)
# rename columns
colnames(NeiSic_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
NeiSic_adhesion$Evalue <- str_replace(NeiSic_adhesion$Evalue,";", "")
NeiSic_adhesion$Evalue <- as.numeric(NeiSic_adhesion$Evalue)
# keep significant matches
NeiSic_adhesion <- subset(NeiSic_adhesion, Evalue < 0.01)
NeiSic_adhesion <- subset(NeiSic_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
NeiSic_adhesion <- ddply(NeiSic_adhesion,"Adhesion",numcolwise(max))
# add species name to column
NeiSic_adhesion$Species <- "Neisseria sicca"


# Import Neisseria subflava adhesin results
NeiSub_adhesion <- read_delim("input_files/functional_data/NeiSub_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
NeiSub_adhesion <- data.frame(NeiSub_adhesion)
# rename columns
colnames(NeiSub_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
NeiSub_adhesion$Evalue <- str_replace(NeiSub_adhesion$Evalue,";", "")
NeiSub_adhesion$Evalue <- as.numeric(NeiSub_adhesion$Evalue)
# keep significant matches
NeiSub_adhesion <- subset(NeiSub_adhesion, Evalue < 0.01)
NeiSub_adhesion <- subset(NeiSub_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
NeiSub_adhesion <- ddply(NeiSub_adhesion,"Adhesion",numcolwise(max))
# add species name to column
NeiSub_adhesion$Species <- "Neisseria subflava"


# Import Porphyromonas asaccharolytica adhesin results
PorAsa_adhesion <- read_delim("input_files/functional_data/PorAsa_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
PorAsa_adhesion <- data.frame(PorAsa_adhesion)
# rename columns
colnames(PorAsa_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
PorAsa_adhesion$Evalue <- str_replace(PorAsa_adhesion$Evalue,";", "")
PorAsa_adhesion$Evalue <- as.numeric(PorAsa_adhesion$Evalue)
# keep significant matches
PorAsa_adhesion <- subset(PorAsa_adhesion, Evalue < 0.01)
PorAsa_adhesion <- subset(PorAsa_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
PorAsa_adhesion <- ddply(PorAsa_adhesion,"Adhesion",numcolwise(max))
# add species name to column
PorAsa_adhesion$Species <- "Porphyromonas asaccharolytica"


# Import Prevotella jejuni adhesin results
PreJej_adhesion <- read_delim("input_files/functional_data/PreJej_adhesion.csv", 
                        delim = ";", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE)
PreJej_adhesion <- data.frame(PreJej_adhesion)
# rename columns
colnames(PreJej_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
PreJej_adhesion$Evalue <- str_replace(PreJej_adhesion$Evalue,";", "")
PreJej_adhesion$Evalue <- as.numeric(PreJej_adhesion$Evalue)
# keep significant matches
PreJej_adhesion <- subset(PreJej_adhesion, Evalue < 0.01)
PreJej_adhesion <- subset(PreJej_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
PreJej_adhesion <- ddply(PreJej_adhesion,"Adhesion",numcolwise(max))
# add species name to column
PreJej_adhesion$Species <- "Prevotella jejuni"


# Import Prevotella melaninogenica adhesin results
PreMel_adhesion <- read_delim("input_files/functional_data/PreMel_adhesion.csv", 
                        delim = ";", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE)
PreMel_adhesion <- data.frame(PreMel_adhesion)
# rename columns
colnames(PreMel_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
PreMel_adhesion$Evalue <- str_replace(PreMel_adhesion$Evalue,";", "")
PreMel_adhesion$Evalue <- as.numeric(PreMel_adhesion$Evalue)
# keep significant matches
PreMel_adhesion <- subset(PreMel_adhesion, Evalue < 0.01)
PreMel_adhesion <- subset(PreMel_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
PreMel_adhesion <- ddply(PreMel_adhesion,"Adhesion",numcolwise(max))
# add species name to column
PreMel_adhesion$Species <- "Prevotella melaninogenica"


# Import Prevotella oris adhesin results
PreOri_adhesion <- read_delim("input_files/functional_data/PreOri_adhesion.csv", 
                              delim = ";", escape_double = FALSE, col_names = FALSE, 
                              trim_ws = TRUE)
PreOri_adhesion <- data.frame(PreOri_adhesion)
# rename columns
colnames(PreOri_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
PreOri_adhesion$Evalue <- str_replace(PreOri_adhesion$Evalue,";", "")
PreOri_adhesion$Evalue <- as.numeric(PreOri_adhesion$Evalue)
# keep significant matches
PreOri_adhesion <- subset(PreOri_adhesion, Evalue < 0.01)
PreOri_adhesion <- subset(PreOri_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
PreOri_adhesion <- ddply(PreOri_adhesion,"Adhesion",numcolwise(max))
# add species name to column
PreOri_adhesion$Species <- "Prevotella oris"


# Import Prevotella ruminicola adhesin results
PreRum_adhesion <- read_delim("input_files/functional_data/PreRum_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
PreRum_adhesion <- data.frame(PreRum_adhesion)
# rename columns
colnames(PreRum_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
PreRum_adhesion$Evalue <- str_replace(PreRum_adhesion$Evalue,";", "")
PreRum_adhesion$Evalue <- as.numeric(PreRum_adhesion$Evalue)
# keep significant matches
PreRum_adhesion <- subset(PreRum_adhesion, Evalue < 0.01)
PreRum_adhesion <- subset(PreRum_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
PreRum_adhesion <- ddply(PreRum_adhesion,"Adhesion",numcolwise(max))
# add species name to column
PreRum_adhesion$Species <- "Prevotella ruminicola"


# Import Pseudoleptotrichia goodfellowii adhesin results
PseGoo_adhesion <- read_delim("input_files/functional_data/PseGoo_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
PseGoo_adhesion <- data.frame(PseGoo_adhesion)
# rename columns
colnames(PseGoo_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
PseGoo_adhesion$Evalue <- str_replace(PseGoo_adhesion$Evalue,";", "")
PseGoo_adhesion$Evalue <- as.numeric(PseGoo_adhesion$Evalue)
# keep significant matches
PseGoo_adhesion <- subset(PseGoo_adhesion, Evalue < 0.01)
PseGoo_adhesion <- subset(PseGoo_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
PseGoo_adhesion <- ddply(PseGoo_adhesion,"Adhesion",numcolwise(max))
# add species name to column
PseGoo_adhesion$Species <- "Pseudoleptotrichia goodfellowii"


# Import Rothia mucilaginosa adhesin results
RotMuc_adhesion <- read_delim("input_files/functional_data/RotMuc_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
RotMuc_adhesion <- data.frame(RotMuc_adhesion)
# rename columns
colnames(RotMuc_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
RotMuc_adhesion$Evalue <- str_replace(RotMuc_adhesion$Evalue,";", "")
RotMuc_adhesion$Evalue <- as.numeric(RotMuc_adhesion$Evalue)
# keep significant matches
RotMuc_adhesion <- subset(RotMuc_adhesion, Evalue < 0.01)
RotMuc_adhesion <- subset(RotMuc_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
RotMuc_adhesion <- ddply(RotMuc_adhesion,"Adhesion",numcolwise(max))
# add species name to column
RotMuc_adhesion$Species <- "Rothia mucilaginosa"


# Import Schaalia odontolytica adhesin results
SchOdo_adhesion <- read_delim("input_files/functional_data/SchOdo_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
SchOdo_adhesion <- data.frame(SchOdo_adhesion)
# rename columns
colnames(SchOdo_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
SchOdo_adhesion$Evalue <- str_replace(SchOdo_adhesion$Evalue,";", "")
SchOdo_adhesion$Evalue <- as.numeric(SchOdo_adhesion$Evalue)
# keep significant matches
SchOdo_adhesion <- subset(SchOdo_adhesion, Evalue < 0.01)
SchOdo_adhesion <- subset(SchOdo_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
SchOdo_adhesion <- ddply(SchOdo_adhesion,"Adhesion",numcolwise(max))
# add species name to column
SchOdo_adhesion$Species <- "Schaalia odontolytica"


# Import Streptococcus australis adhesin results
StrAus_adhesion <- read_delim("input_files/functional_data/StrAus_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
StrAus_adhesion <- data.frame(StrAus_adhesion)
# rename columns
colnames(StrAus_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
StrAus_adhesion$Evalue <- str_replace(StrAus_adhesion$Evalue,";", "")
StrAus_adhesion$Evalue <- as.numeric(StrAus_adhesion$Evalue)
# keep significant matches
StrAus_adhesion <- subset(StrAus_adhesion, Evalue < 0.01)
StrAus_adhesion <- subset(StrAus_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
StrAus_adhesion <- ddply(StrAus_adhesion,"Adhesion",numcolwise(max))
# add species name to column
StrAus_adhesion$Species <- "Streptococcus australis"


# Import Streptococcus gordonii adhesin results
StrGor_adhesion <- read_delim("input_files/functional_data/StrGor_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
StrGor_adhesion <- data.frame(StrGor_adhesion)
# rename columns
colnames(StrGor_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
StrGor_adhesion$Evalue <- str_replace(StrGor_adhesion$Evalue,";", "")
StrGor_adhesion$Evalue <- as.numeric(StrGor_adhesion$Evalue)
# keep significant matches
StrGor_adhesion <- subset(StrGor_adhesion, Evalue < 0.01)
StrGor_adhesion <- subset(StrGor_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
StrGor_adhesion <- ddply(StrGor_adhesion,"Adhesion",numcolwise(max))
# add species name to column
StrGor_adhesion$Species <- "Streptococcus gordonii"


# Import Streptococcus koreensis adhesin results
StrKor_adhesion <- read_delim("input_files/functional_data/StrKor_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
StrKor_adhesion <- data.frame(StrKor_adhesion)
# rename columns
colnames(StrKor_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
StrKor_adhesion$Evalue <- str_replace(StrKor_adhesion$Evalue,";", "")
StrKor_adhesion$Evalue <- as.numeric(StrKor_adhesion$Evalue)
# keep significant matches
StrKor_adhesion <- subset(StrKor_adhesion, Evalue < 0.01)
StrKor_adhesion <- subset(StrKor_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
StrKor_adhesion <- ddply(StrKor_adhesion,"Adhesion",numcolwise(max))
# add species name to column
StrKor_adhesion$Species <- "Streptococcus koreensis"


# Import Streptococcus mitis adhesin results
StrMit_adhesion <- read_delim("input_files/functional_data/StrMit_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
StrMit_adhesion <- data.frame(StrMit_adhesion)
# rename columns
colnames(StrMit_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
StrMit_adhesion$Evalue <- str_replace(StrMit_adhesion$Evalue,";", "")
StrMit_adhesion$Evalue <- as.numeric(StrMit_adhesion$Evalue)
# keep significant matches
StrMit_adhesion <- subset(StrMit_adhesion, Evalue < 0.01)
StrMit_adhesion <- subset(StrMit_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
StrMit_adhesion <- ddply(StrMit_adhesion,"Adhesion",numcolwise(max))
# add species name to column
StrMit_adhesion$Species <- "Streptococcus mitis"


# Import Streptococcus parasanguinis adhesin results
StrPar_adhesion <- read_delim("input_files/functional_data/StrPar_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
StrPar_adhesion <- data.frame(StrPar_adhesion)
# rename columns
colnames(StrPar_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
StrPar_adhesion$Evalue <- str_replace(StrPar_adhesion$Evalue,";", "")
StrPar_adhesion$Evalue <- as.numeric(StrPar_adhesion$Evalue)
# keep significant matches
StrPar_adhesion <- subset(StrPar_adhesion, Evalue < 0.01)
StrPar_adhesion <- subset(StrPar_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
StrPar_adhesion <- ddply(StrPar_adhesion,"Adhesion",numcolwise(max))
# add species name to column
StrPar_adhesion$Species <- "Streptococcus parasanguinis"


# Import Streptococcus pseudopneumoniae adhesin results
StrPsePneu_adhesion <- read_delim("input_files/functional_data/StrPsePneu_adhesion.csv", 
                         delim = ";", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)
StrPsePneu_adhesion <- data.frame(StrPsePneu_adhesion)
# rename columns
colnames(StrPsePneu_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
StrPsePneu_adhesion$Evalue <- str_replace(StrPsePneu_adhesion$Evalue,";", "")
StrPsePneu_adhesion$Evalue <- as.numeric(StrPsePneu_adhesion$Evalue)
# keep significant matches
StrPsePneu_adhesion <- subset(StrPsePneu_adhesion, Evalue < 0.01)
StrPsePneu_adhesion <- subset(StrPsePneu_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
StrPsePneu_adhesion <- ddply(StrPsePneu_adhesion,"Adhesion",numcolwise(max))
# add species name to column
StrPsePneu_adhesion$Species <- "Streptococcus pseudopneumoniae"


# Import Veillonella atypica adhesin results
VeiAty_adhesion <- read_delim("input_files/functional_data/VeiAty_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
VeiAty_adhesion <- data.frame(VeiAty_adhesion)
# rename columns
colnames(VeiAty_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
VeiAty_adhesion$Evalue <- str_replace(VeiAty_adhesion$Evalue,";", "")
VeiAty_adhesion$Evalue <- as.numeric(VeiAty_adhesion$Evalue)
# keep significant matches
VeiAty_adhesion <- subset(VeiAty_adhesion, Evalue < 0.01)
VeiAty_adhesion <- subset(VeiAty_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
VeiAty_adhesion <- ddply(VeiAty_adhesion,"Adhesion",numcolwise(max))
# add species name to column
VeiAty_adhesion$Species <- "Veillonella atypica"


# Import Veillonella dispar adhesin results
VeiDis_adhesion <- read_delim("input_files/functional_data/VeiDis_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
VeiDis_adhesion <- data.frame(VeiDis_adhesion)
# rename columns
colnames(VeiDis_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
VeiDis_adhesion$Evalue <- str_replace(VeiDis_adhesion$Evalue,";", "")
VeiDis_adhesion$Evalue <- as.numeric(VeiDis_adhesion$Evalue)
# keep significant matches
VeiDis_adhesion <- subset(VeiDis_adhesion, Evalue < 0.01)
VeiDis_adhesion <- subset(VeiDis_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
VeiDis_adhesion <- ddply(VeiDis_adhesion,"Adhesion",numcolwise(max))
# add species name to column
VeiDis_adhesion$Species <- "Veillonella dispar"


# Import Veillonella rodentium adhesin results
VeiRod_adhesion <- read_delim("input_files/functional_data/VeiRod_adhesion.csv", 
                     delim = ";", escape_double = FALSE, col_names = FALSE, 
                     trim_ws = TRUE)
VeiRod_adhesion <- data.frame(VeiRod_adhesion)
# rename columns
colnames(VeiRod_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
VeiRod_adhesion$Evalue <- str_replace(VeiRod_adhesion$Evalue,";", "")
VeiRod_adhesion$Evalue <- as.numeric(VeiRod_adhesion$Evalue)
# keep significant matches
VeiRod_adhesion <- subset(VeiRod_adhesion, Evalue < 0.01)
VeiRod_adhesion <- subset(VeiRod_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
VeiRod_adhesion <- ddply(VeiRod_adhesion,"Adhesion",numcolwise(max))
# add species name to column
VeiRod_adhesion$Species <- "Veillonella rodentium"


# Import Pseudomonas aeruginosa adhesin results
Pao1_adhesion <- read_delim("input_files/functional_data/PAO1_adhesion.csv", 
                              delim = ";", escape_double = FALSE, col_names = FALSE, 
                              trim_ws = TRUE)
Pao1_adhesion <- data.frame(Pao1_adhesion)
# rename columns
colnames(Pao1_adhesion) <- c("Adhesion", "Score", "Evalue")
# remove error in data frame
Pao1_adhesion$Evalue <- str_replace(Pao1_adhesion$Evalue,";", "")
Pao1_adhesion$Evalue <- as.numeric(Pao1_adhesion$Evalue)
# keep significant matches
Pao1_adhesion <- subset(Pao1_adhesion, Evalue < 0.01)
Pao1_adhesion <- subset(Pao1_adhesion, Score > score_value)
# remove duplicate rows and keep the maximum score that was obtained for an adhesin
Pao1_adhesion <- ddply(Pao1_adhesion,"Adhesion",numcolwise(max))
# add species name to column
Pao1_adhesion$Species <- "Pseudomonas aeruginosa"
Pao1_adhesion$State <- "Pathogen"
Pao1_adhesion$Species_type <- "Pathogen"
Pao1_adhesion$Normalisation <- "not relevant"


# merge all data frames
adhesion_all <- data.frame(rbind(ActIsr_adhesion, CapEnd_adhesion, CapSpu_adhesion, FusNuc_adhesion,
                                 GemHae_adhesion, GemMor_adhesion, GemSan_adhesion, HaePar_adhesion, 
                                 LepHof_adhesion, LepHon_adhesion, NeiCin_adhesion, NeiPol_adhesion,
                                 NeiSic_adhesion, NeiSub_adhesion, PorAsa_adhesion, PreJej_adhesion,
                                 PreMel_adhesion, PreOri_adhesion, PreRum_adhesion, PseGoo_adhesion,
                                 RotMuc_adhesion, SchOdo_adhesion, StrAus_adhesion,
                                 StrGor_adhesion, StrKor_adhesion, StrMit_adhesion, StrPar_adhesion, 
                                 StrPsePneu_adhesion, VeiAty_adhesion, VeiDis_adhesion, VeiRod_adhesion))

############################################################################################################
# Subset data frame

# BCPHC, healthy
# Subset core and rare species
bg_species_25_bcphc_h_rare <- subset(bg_species_25_bcphc_h, Species_type == "background_rare")
bg_species_25_bcphc_h_core <- subset(bg_species_25_bcphc_h, Species_type == "background_core")
# extract adhesin pattern of rare species
adhesion_bcphc_25_h_rare <- adhesion_all[adhesion_all$Species %in% bg_species_25_bcphc_h_rare$Species,]
# add meta data
adhesion_bcphc_25_h_rare$State <- "Healthy"
adhesion_bcphc_25_h_rare$Species_type <- "Rare"
adhesion_bcphc_25_h_rare$Normalisation <- "BCPHC"
# extract adhesin pattern of core species
adhesion_bcphc_25_h_core <- adhesion_all[adhesion_all$Species %in% bg_species_25_bcphc_h_core$Species,]
adhesion_bcphc_25_h_core$State <- "Healthy"
adhesion_bcphc_25_h_core$Species_type <- "Core"
adhesion_bcphc_25_h_core$Normalisation <- "BCPHC"


# BCPHC, CF
# Subset core and rare species
bg_species_25_bcphc_cf_rare <- subset(bg_species_25_bcphc_cf, Species_type == "background_rare")
bg_species_25_bcphc_cf_core <- subset(bg_species_25_bcphc_cf, Species_type == "background_core")
# extract adhesin pattern of rare species
adhesion_bcphc_25_cf_rare <- adhesion_all[adhesion_all$Species %in% bg_species_25_bcphc_cf_rare$Species,]
# add meta data
adhesion_bcphc_25_cf_rare$State <- "CF"
adhesion_bcphc_25_cf_rare$Species_type <- "Rare"
adhesion_bcphc_25_cf_rare$Normalisation <- "BCPHC"
# extract adhesin pattern of core species
adhesion_bcphc_25_cf_core <- adhesion_all[adhesion_all$Species %in% bg_species_25_bcphc_cf_core$Species,]
adhesion_bcphc_25_cf_core$State <- "CF"
adhesion_bcphc_25_cf_core$Species_type <- "Core"
adhesion_bcphc_25_cf_core$Normalisation <- "BCPHC"


# RLE, healthy
# Subset core and rare species
bg_species_25_rle_h_rare <- subset(bg_species_25_rle_h, Species_type == "background_rare")
bg_species_25_rle_h_core <- subset(bg_species_25_rle_h, Species_type == "background_core")
# extract adhesin pattern of rare species
adhesion_rle_25_h_rare <- adhesion_all[adhesion_all$Species %in% bg_species_25_rle_h_rare$Species,]
# add meta data
adhesion_rle_25_h_rare$State <- "Healthy"
adhesion_rle_25_h_rare$Species_type <- "Rare"
adhesion_rle_25_h_rare$Normalisation <- "RLE"
# extract adhesin pattern of core species
adhesion_rle_25_h_core <- adhesion_all[adhesion_all$Species %in% bg_species_25_rle_h_core$Species,]
adhesion_rle_25_h_core$State <- "Healthy"
adhesion_rle_25_h_core$Species_type <- "Core"
adhesion_rle_25_h_core$Normalisation <- "RLE"

# rle, CF
# Subset core and rare species
bg_species_25_rle_cf_rare <- subset(bg_species_25_rle_cf, Species_type == "background_rare")
bg_species_25_rle_cf_core <- subset(bg_species_25_rle_cf, Species_type == "background_core")
# extract adhesin pattern of rare species
adhesion_rle_25_cf_rare <- adhesion_all[adhesion_all$Species %in% bg_species_25_rle_cf_rare$Species,]
# add meta data
adhesion_rle_25_cf_rare$State <- "CF"
adhesion_rle_25_cf_rare$Species_type <- "Rare"
adhesion_rle_25_cf_rare$Normalisation <- "RLE"
# extract adhesin pattern of core species
adhesion_rle_25_cf_core <- adhesion_all[adhesion_all$Species %in% bg_species_25_rle_cf_core$Species,]
adhesion_rle_25_cf_core$State <- "CF"
adhesion_rle_25_cf_core$Species_type <- "Core"
adhesion_rle_25_cf_core$Normalisation <- "RLE"


# VST, healthy
# Subset core and rare species
bg_species_25_vst_h_rare <- subset(bg_species_25_vst_h, Species_type == "background_rare")
bg_species_25_vst_h_core <- subset(bg_species_25_vst_h, Species_type == "background_core")
# extract adhesin pattern of rare species
adhesion_vst_25_h_rare <- adhesion_all[adhesion_all$Species %in% bg_species_25_vst_h_rare$Species,]
# add meta data
adhesion_vst_25_h_rare$State <- "Healthy"
adhesion_vst_25_h_rare$Species_type <- "Rare"
adhesion_vst_25_h_rare$Normalisation <- "VST"
# extract adhesin pattern of core species
adhesion_vst_25_h_core <- adhesion_all[adhesion_all$Species %in% bg_species_25_vst_h_core$Species,]
adhesion_vst_25_h_core$State <- "Healthy"
adhesion_vst_25_h_core$Species_type <- "Core"
adhesion_vst_25_h_core$Normalisation <- "VST"


# VST, CF
# Subset core species (there are no rare species for CF, VST-normalised data)
bg_species_25_vst_cf_core <- subset(bg_species_25_vst_cf, Species_type == "background_core")
# extract adhesin pattern of rare species
adhesion_vst_25_cf_core <- adhesion_all[adhesion_all$Species %in% bg_species_25_vst_cf_core$Species,]
# add meta data
adhesion_vst_25_cf_core$State <- "CF"
adhesion_vst_25_cf_core$Species_type <- "Core"
adhesion_vst_25_cf_core$Normalisation <- "VST"


# merge all adhesion dataframes
adhesion_25 <- data.frame(rbind(adhesion_bcphc_25_h_core, adhesion_bcphc_25_cf_core,
                                adhesion_bcphc_25_h_rare, adhesion_bcphc_25_cf_rare,
                                adhesion_rle_25_h_core, adhesion_rle_25_cf_core,
                                adhesion_rle_25_h_rare, adhesion_rle_25_cf_rare,
                                adhesion_vst_25_h_core, adhesion_vst_25_h_rare,
                                adhesion_vst_25_cf_core))

# clean up the adhesion name
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "1", "")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "\\/2", "")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, " 2", "")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "\\/3", "")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, " 3", "")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, " 1", "")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "\\/4", "")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "MafA ", "MafA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "CdiA", "Contact-dependent growth inhibition CdiA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Accessory Sec system prot", "Accessory Sec system")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Adhesin ", "")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "BmaC autotrans", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Adhesin/invasin TibA autot", "Adhesin/invasin TibA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Adhesion and penetration p", "Adhesion/penetration protein")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "ATP-dependent zinc metallo", "Zinc metalloprotease FtsH")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Autotransporter adhesin Nh", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Autotransporter adhesin Sa", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Autotransporter adhesin Eh", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Autotransporter adhesi", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Trimeric autotransport", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Immunoglobulin-binding pro", "Immunoglobulin-binding protein")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "HTH-type transcriptional r", "TcaR")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Toxin ", "")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Bone sialoprotein-binding p", "Bone sialoprotein-binding protein")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Deoxyribonuclease CdiA", "Contact-dependent growth inhibition CdiA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Cyclic pyranopterin monoph", "Pyranopterin monophosphate")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Extracellular matrix-bindin", "Extracellular matrix binding protein Ebh")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "High-affinity zinc uptake", "High-affinity zinc uptake ZnuA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Histidine protein kinase S", "Histidine protein kinase SaeS")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "IgA-specific serine endopep", "IgA-specific autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Manganese import ATP-bindi", "ABC transport system")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Manganese import system pe", "ABC transport system")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Metal ABC transporter subs", "ABC transport system")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Outer membrane protei", "Outer membrane protein")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Outer membrane porin F", "Outer membrane porin")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Outer membrane usher prote", "Outer membrane protein")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Plasmin and fibronectin-bi", "Fibronectin-binding protein B, FnBPB")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Pneumococcal serine-ri", "Pneumococcal serine-rich repeat protein, PsrP")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Poly-beta-,6-N-acetyl-D-g", "Poly-beta-1,6-N-acetyl-D-glucosamine")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Probable cell-surface antig", "Cell-surface antigen")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Probable manganese-depende", "Inorganic pyrophosphatase, ppaC")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Probable site-specific rec", "Site-specific recombinase")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Protein translocase subun", "Protein translocase SecA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Putat", "ABC transport system")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Putative metal ABC transpor", "ABC transport system")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Putative zinc metalloprote", "Putative zinc metalloprotease")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Ribosomal RNA small subuni", "Translational regulator CsrA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Serine-aspartate repeat-co", "Bone sialoprotein-binding protein")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Serine-rich adhesin for pl", "Serine-rich adhesin for platelets, sraP")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Serine protease HtrA-like", "Serine protease HtrA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Serine protease SepA autot", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Subtilase cytotoxin subuni", "Subtilase cytotoxin, subA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Temperature-sensitive hemag", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Trimeric autotransporterer", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Trimeric autotransporterer a", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Trimeric autotransporterern Ba", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "tRNA nuclease CdiA Contact-dependent growth inhibition CdiA", "Contact-dependent growth inhibition CdiA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "tRNA nuclease CdiA-2 Contact-dependent growth inhibition CdiA", "Contact-dependent growth inhibition CdiA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Contact-dependent growth inhibition CdiA-2", "Contact-dependent growth inhibition CdiA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Cell division protein SepF", "Cell division protein SepF")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Two-component sensor PprA", "Two-component sensor PprA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Chaperone protein DnaK", "Chaperone protein DnaK")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Cell surface antigen I/II", "Cell surface antigen")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Response regulator SaeR", "Response regulator protein")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Hemagglutinin tsh autotransporter", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Two-component response reg", "Two-component response")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "UDP-N-acetylglucosamine--p", "UDP-N-acetylglucosamine")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Uncharacterized outer memb", "Outer membrane protein")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "ABC transport systemive metal ABC transpor", "ABC transport system")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "ABC transport systemive zinc metalloprote", "ABC transport system")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Uncharacterized protein YcgV", "Uncharacterized protein YcgV")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Trimeric autotransportern Ba", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Trimeric autotransporter a", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Sensor protein VraS", "Sensor protein VraS")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Ata autotransporter", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Manganese ABC transporter", "ABC transport system")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Outer membrane proteinn Ics", "Outer membrane protein")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Uro-adherence factor A", "Uro-adherence factor A")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Cell surface antigen", "Cell-surface antigen")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Probable autotransporter", "Trimeric autotransporter")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "YadA", "Collagen-binding outer membrane protein, YadA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "6S rRNA endonuclease Contact-dependent growth inhibition CdiA", 
                                        "Contact-dependent growth inhibition CdiA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "Deoxyribonuclease Contact-dependent growth inhibition CdiA", 
                                        "Contact-dependent growth inhibition CdiA")
adhesion_25$Adhesion <- str_replace_all(adhesion_25$Adhesion, "tRNA nuclease Contact-dependent growth inhibition CdiA", 
                                        "Contact-dependent growth inhibition CdiA")

# generate data frame with adhesion count data per species and normalisation type
adhesion_25_df <- as.data.frame(table(adhesion_25$Adhesion, adhesion_25$State, adhesion_25$Species_type, adhesion_25$Normalisation))
# add column with merged name
adhesion_25_df$Merge <- paste(adhesion_25_df$Var3,"-",adhesion_25_df$Var4)
# add pseudo count of 1 to frequency value
adhesion_25_df$Freq_2 <- adhesion_25_df$Freq + 1
# Log transformation of frequency value
adhesion_25_df$Freq_3 <- log(adhesion_25_df$Freq_2)
# set negative value to 0
adhesion_25_df$Freq_3 <- ifelse(adhesion_25_df$Freq_3 < 0, 0, adhesion_25_df$Freq_3)

# generate heatmap
adhesion_heatmap <-
  ggplot(adhesion_25_df, aes(x=Merge, y=Var1, fill=Freq_3)) +
  geom_tile(colour="white", size=0.1, width=0.9) + 
  xlab(" ") + ylab(" ") +
  scale_fill_gradientn(colours=c("white", "beige", "yellow", "darkslategray4", "darkslategrey")) +
  facet_wrap(~Var2, scales = "free_x") + theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle=45, size=12, hjust=1),
        axis.text.y = element_text(size=10, vjust=0.5, hjust=0.5),
        strip.text = element_text(size=12),
        legend.position = "bottom",
        panel.grid = element_blank(),
        strip.background = element_rect(fill="white"),
        legend.title = element_blank()) 

############################################################################################################
# Output, figures
tiff(filename = "output_figures/adhesion_heatmap.tif", width=30, height = 32, 
     res=300, bg="white", units="cm", pointsize = 14)
adhesion_heatmap
dev.off()

tiff(filename="output_figures/fm_plot.tif", width=12, height=22, 
     units="cm", res=600, bg="white")
FM_plot_H
dev.off()

tiff(filename="output_figures/dendrogram.tif", res=300, 
     units="in", width=6, height=5)
strep_hier
dev.off()

tiff(filename="output_figures/pao1_plots.tif", width=3.4, height=1.4, 
     units="in", res=300, bg="white")
PAO1_stats_plot
dev.off()

############################################################################################################