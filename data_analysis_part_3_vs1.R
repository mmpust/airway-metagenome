# title: "Functional analysis of the human airway rare species"
# author: "Marie-Madlen Pust"
# date: "13 04 2021"

# clean global environment
rm(list=ls())

# set working directory
setwd("C:/Users/marie/Desktop/Phd studies3/e_rare_biosphere/R")

# load packages
library('readr')
library('plyr')
library('purrr')
library('dplyr')
library('ggplot2')
library('vegan')
library('ggpubr')
library('factoextra')
library('stringr')
library('pheatmap')
library('RColorBrewer')
library('ggdendro')

ActIsr <- read_delim("data_input/ActIsr_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
ActIsr <- data.frame(ActIsr)
ActIsr$pathway_key_reactions <- NULL
ActIsr$pathway_completeness <- NULL
ActIsr$pathway_candidate_reaction <- NULL
colnames(ActIsr) <- c("id", "pathway", "ActIsr")


CamCon <- read_delim("data_input/CamCon_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
CamCon <- data.frame(CamCon)
CamCon$pathway_key_reactions <- NULL
CamCon$pathway_completeness <- NULL
CamCon$pathway_candidate_reaction <- NULL
colnames(CamCon) <- c("id", "pathway", "CamCon")


CapEnd <- read_delim("data_input/CapEnd_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
CapEnd <- data.frame(CapEnd)
CapEnd$pathway_key_reactions <- NULL
CapEnd$pathway_completeness <- NULL
CapEnd$pathway_candidate_reaction <- NULL
colnames(CapEnd) <- c("id", "pathway", "CapEnd")


CapSpu <- read_delim("data_input/CapSpu_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
CapSpu <- data.frame(CapSpu)
CapSpu$pathway_key_reactions <- NULL
CapSpu$pathway_completeness <- NULL
CapSpu$pathway_candidate_reaction <- NULL
colnames(CapSpu) <- c("id", "pathway", "CapSpu")


FusNuc <- read_delim("data_input/FusNuc_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
FusNuc <- data.frame(FusNuc)
FusNuc$pathway_key_reactions <- NULL
FusNuc$pathway_completeness <- NULL
FusNuc$pathway_candidate_reaction <- NULL
colnames(FusNuc) <- c("id", "pathway", "FusNuc")


GemHae <- read_delim("data_input/GemHae_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
GemHae <- data.frame(GemHae)
GemHae$pathway_key_reactions <- NULL
GemHae$pathway_completeness <- NULL
GemHae$pathway_candidate_reaction <- NULL
colnames(GemHae) <- c("id", "pathway", "GemHae")


GemMor <- read_delim("data_input/GemMor_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
GemMor <- data.frame(GemMor)
GemMor$pathway_key_reactions <- NULL
GemMor$pathway_completeness <- NULL
GemMor$pathway_candidate_reaction <- NULL
colnames(GemMor) <- c("id", "pathway", "GemMor")


GemSan <- read_delim("data_input/GemSan_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
GemSan <- data.frame(GemSan)
GemSan$pathway_key_reactions <- NULL
GemSan$pathway_completeness <- NULL
GemSan$pathway_candidate_reaction <- NULL
colnames(GemSan) <- c("id", "pathway", "GemSan")


HaeHae <- read_delim("data_input/HaeHae_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
HaeHae <- data.frame(HaeHae)
HaeHae$pathway_key_reactions <- NULL
HaeHae$pathway_completeness <- NULL
HaeHae$pathway_candidate_reaction <- NULL
colnames(HaeHae) <- c("id", "pathway", "HaeHae")

HaePar <- read_delim("data_input/HaePar_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
HaePar <- data.frame(HaePar)
HaePar$pathway_key_reactions <- NULL
HaePar$pathway_completeness <- NULL
HaePar$pathway_candidate_reaction <- NULL
colnames(HaePar) <- c("id", "pathway", "HaePar")


LepBuc <- read_delim("data_input/LepBuc_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
LepBuc <- data.frame(LepBuc)
LepBuc$pathway_key_reactions <- NULL
LepBuc$pathway_completeness <- NULL
LepBuc$pathway_candidate_reaction <- NULL
colnames(LepBuc) <- c("id", "pathway", "LepBuc")


LepHof <- read_delim("data_input/LepHof_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
LepHof <- data.frame(LepHof)
LepHof$pathway_key_reactions <- NULL
LepHof$pathway_completeness <- NULL
LepHof$pathway_candidate_reaction <- NULL
colnames(LepHof) <- c("id", "pathway", "LepHof")


LepHon <- read_delim("data_input/LepHon_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
LepHon <- data.frame(LepHon)
LepHon$pathway_key_reactions <- NULL
LepHon$pathway_completeness <- NULL
LepHon$pathway_candidate_reaction <- NULL
colnames(LepHon) <- c("id", "pathway", "LepHon")


LepWad <- read_delim("data_input/LepWad_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
LepWad <- data.frame(LepWad)
LepWad$pathway_key_reactions <- NULL
LepWad$pathway_completeness <- NULL
LepWad$pathway_candidate_reaction <- NULL
colnames(LepWad) <- c("id", "pathway", "LepWad")


NeiMuc <- read_delim("data_input/NeiMuc_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
NeiMuc <- data.frame(NeiMuc)
NeiMuc$pathway_key_reactions <- NULL
NeiMuc$pathway_completeness <- NULL
NeiMuc$pathway_candidate_reaction <- NULL
colnames(NeiMuc) <- c("id", "pathway", "NeiMuc")


NeiPol <- read_delim("data_input/NeiPol_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
NeiPol <- data.frame(NeiPol)
NeiPol$pathway_key_reactions <- NULL
NeiPol$pathway_completeness <- NULL
NeiPol$pathway_candidate_reaction <- NULL
colnames(NeiPol) <- c("id", "pathway", "NeiPol")


NeiSic <- read_delim("data_input/NeiSic_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
NeiSic <- data.frame(NeiSic)
NeiSic$pathway_key_reactions <- NULL
NeiSic$pathway_completeness <- NULL
NeiSic$pathway_candidate_reaction <- NULL
colnames(NeiSic) <- c("id", "pathway", "NeiSic")


NeiSub <- read_delim("data_input/NeiSub_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
NeiSub <- data.frame(NeiSub)
NeiSub$pathway_key_reactions <- NULL
NeiSub$pathway_completeness <- NULL
NeiSub$pathway_candidate_reaction <- NULL
colnames(NeiSub) <- c("id", "pathway", "NeiSub")


PorAsa <- read_delim("data_input/PorAsa_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
PorAsa <- data.frame(PorAsa)
PorAsa$pathway_key_reactions <- NULL
PorAsa$pathway_completeness <- NULL
PorAsa$pathway_candidate_reaction <- NULL
colnames(PorAsa) <- c("id", "pathway", "PorAsa")


PreJej_II <- read_delim("data_input/PreJej_II_pathways_reactions.csv", ";", 
                        escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                                pathway_completeness = col_number(), 
                                                                pathway_completeness_perc = col_number(), 
                                                                pathway_key_reactions = col_character()), trim_ws = TRUE)
PreJej_II <- data.frame(PreJej_II)
PreJej_II$pathway_key_reactions <- NULL
PreJej_II$pathway_completeness <- NULL
PreJej_II$pathway_candidate_reaction <- NULL
colnames(PreJej_II) <- c("id", "pathway", "PreJej_II")


PreJej <- read_delim("data_input/PreJej_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
PreJej <- data.frame(PreJej)
PreJej$pathway_key_reactions <- NULL
PreJej$pathway_completeness <- NULL
PreJej$pathway_candidate_reaction <- NULL
colnames(PreJej) <- c("id", "pathway", "PreJej")
PreJej_II$PreJej_III <- PreJej_II$PreJej_II + PreJej$PreJej


PreMel_II <- read_delim("data_input/PreMel_II_pathways_reactions.csv", ";", 
                        escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                                pathway_completeness = col_number(), 
                                                                pathway_completeness_perc = col_number(), 
                                                                pathway_key_reactions = col_character()), trim_ws = TRUE)
PreMel_II <- data.frame(PreMel_II)
PreMel_II$pathway_key_reactions <- NULL
PreMel_II$pathway_completeness <- NULL
PreMel_II$pathway_candidate_reaction <- NULL
colnames(PreMel_II) <- c("id", "pathway", "PreMel_II")


PreMel <- read_delim("data_input/PreMel_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
PreMel <- data.frame(PreMel)
PreMel$pathway_key_reactions <- NULL
PreMel$pathway_completeness <- NULL
PreMel$pathway_candidate_reaction <- NULL
colnames(PreMel) <- c("id", "pathway", "PreMel")
PreMel_II$PreMel_III <- PreMel_II$PreMel_II + PreMel$PreMel


PreOri <- read_delim("data_input/PreOri_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
PreOri <- data.frame(PreOri)
PreOri$pathway_key_reactions <- NULL
PreOri$pathway_completeness <- NULL
PreOri$pathway_candidate_reaction <- NULL
colnames(PreOri) <- c("id", "pathway", "PreOri")


PreRum <- read_delim("data_input/PreRum_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
PreRum <- data.frame(PreRum)
PreRum$pathway_key_reactions <- NULL
PreRum$pathway_completeness <- NULL
PreRum$pathway_candidate_reaction <- NULL
colnames(PreRum) <- c("id", "pathway", "PreRum")


PrevEno <- read_delim("data_input/PrevEno_pathways_reactions.csv", ";", 
                      escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                              pathway_completeness = col_number(), 
                                                              pathway_completeness_perc = col_number(), 
                                                              pathway_key_reactions = col_character()), trim_ws = TRUE)
PrevEno <- data.frame(PrevEno)
PrevEno$pathway_key_reactions <- NULL
PrevEno$pathway_completeness <- NULL
PrevEno$pathway_candidate_reaction <- NULL
colnames(PrevEno) <- c("id", "pathway", "PrevEno")


PseGoo <- read_delim("data_input/PseGoo_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
PseGoo <- data.frame(PseGoo)
PseGoo$pathway_key_reactions <- NULL
PseGoo$pathway_completeness <- NULL
PseGoo$pathway_candidate_reaction <- NULL
colnames(PseGoo) <- c("id", "pathway", "PseGoo")


RotMuc <- read_delim("data_input/RotMuc_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
RotMuc <- data.frame(RotMuc)
RotMuc$pathway_key_reactions <- NULL
RotMuc$pathway_completeness <- NULL
RotMuc$pathway_candidate_reaction <- NULL
colnames(RotMuc) <- c("id", "pathway", "RotMuc")


SchMey <- read_delim("data_input/SchMey_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
SchMey <- data.frame(SchMey)
SchMey$pathway_key_reactions <- NULL
SchMey$pathway_completeness <- NULL
SchMey$pathway_candidate_reaction <- NULL
colnames(SchMey) <- c("id", "pathway", "SchMey")


SchOdo <- read_delim("data_input/SchOdo_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
SchOdo <- data.frame(SchOdo)
SchOdo$pathway_key_reactions <- NULL
SchOdo$pathway_completeness <- NULL
SchOdo$pathway_candidate_reaction <- NULL
colnames(SchOdo) <- c("id", "pathway", "SchOdo")


StrAus <- read_delim("data_input/StrAus_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
StrAus <- data.frame(StrAus)
StrAus$pathway_key_reactions <- NULL
StrAus$pathway_completeness <- NULL
StrAus$pathway_candidate_reaction <- NULL
colnames(StrAus) <- c("id", "pathway", "StrAus")


StrEqu <- read_delim("data_input/StrEqu_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
StrEqu <- data.frame(StrEqu)
StrEqu$pathway_key_reactions <- NULL
StrEqu$pathway_completeness <- NULL
StrEqu$pathway_candidate_reaction <- NULL
colnames(StrEqu) <- c("id", "pathway", "StrEqu")


StrGor <- read_delim("data_input/StrGor_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
StrGor <- data.frame(StrGor)
StrGor$pathway_key_reactions <- NULL
StrGor$pathway_completeness <- NULL
StrGor$pathway_candidate_reaction <- NULL
colnames(StrGor) <- c("id", "pathway", "StrGor")


StrHim <- read_delim("data_input/StrHim_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
StrHim <- data.frame(StrHim)
StrHim$pathway_key_reactions <- NULL
StrHim$pathway_completeness <- NULL
StrHim$pathway_candidate_reaction <- NULL
colnames(StrHim) <- c("id", "pathway", "StrHim")


StrMit <- read_delim("data_input/StrMit_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
StrMit <- data.frame(StrMit)
StrMit$pathway_key_reactions <- NULL
StrMit$pathway_completeness <- NULL
StrMit$pathway_candidate_reaction <- NULL
colnames(StrMit) <- c("id", "pathway", "StrMit")


StrOra <- read_delim("data_input/StrOra_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
StrOra <- data.frame(StrOra)
StrOra$pathway_key_reactions <- NULL
StrOra$pathway_completeness <- NULL
StrOra$pathway_candidate_reaction <- NULL
colnames(StrOra) <- c("id", "pathway", "StrOra")


StrPar <- read_delim("data_input/StrPar_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
StrPar <- data.frame(StrPar)
StrPar$pathway_key_reactions <- NULL
StrPar$pathway_completeness <- NULL
StrPar$pathway_candidate_reaction <- NULL
colnames(StrPar) <- c("id", "pathway", "StrPar")


StrPne <- read_delim("data_input/StrPne_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
StrPne <- data.frame(StrPne)
StrPne$pathway_key_reactions <- NULL
StrPne$pathway_completeness <- NULL
StrPne$pathway_candidate_reaction <- NULL
colnames(StrPne) <- c("id", "pathway", "StrPne")


StrSal <- read_delim("data_input/StrSal_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
StrSal <- data.frame(StrSal)
StrSal$pathway_key_reactions <- NULL
StrSal$pathway_completeness <- NULL
StrSal$pathway_candidate_reaction <- NULL
colnames(StrSal) <- c("id", "pathway", "StrSal")


StrVes <- read_delim("data_input/StrVes_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
StrVes <- data.frame(StrVes)
StrVes$pathway_key_reactions <- NULL
StrVes$pathway_completeness <- NULL
StrVes$pathway_candidate_reaction <- NULL
colnames(StrVes) <- c("id", "pathway", "StrVes")


VeiAty <- read_delim("data_input/VeiAty_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
VeiAty <- data.frame(VeiAty)
VeiAty$pathway_key_reactions <- NULL
VeiAty$pathway_completeness <- NULL
VeiAty$pathway_candidate_reaction <- NULL
colnames(VeiAty) <- c("id", "pathway", "VeiAty")


VeiDis <- read_delim("data_input/VeiDis_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
VeiDis <- data.frame(VeiDis)
VeiDis$pathway_key_reactions <- NULL
VeiDis$pathway_completeness <- NULL
VeiDis$pathway_candidate_reaction <- NULL
colnames(VeiDis) <- c("id", "pathway", "VeiDis")


VeiPar <- read_delim("data_input/VeiPar_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
VeiPar <- data.frame(VeiPar)
VeiPar$pathway_key_reactions <- NULL
VeiPar$pathway_completeness <- NULL
VeiPar$pathway_candidate_reaction <- NULL
colnames(VeiPar) <- c("id", "pathway", "VeiPar")


VeiRod <- read_delim("data_input/VeiRod_pathways_reactions.csv", ";", 
                     escape_double = FALSE, col_types = cols(pathway_candidate_reaction = col_character(), 
                                                             pathway_completeness = col_number(), 
                                                             pathway_completeness_perc = col_number(), 
                                                             pathway_key_reactions = col_character()), trim_ws = TRUE)
VeiRod <- data.frame(VeiRod)
VeiRod$pathway_key_reactions <- NULL
VeiRod$pathway_completeness <- NULL
VeiRod$pathway_candidate_reaction <- NULL
colnames(VeiRod) <- c("id", "pathway", "VeiRod")

merge_species <- data.frame(cbind(ActIsr, CamCon$CamCon, CapEnd$CapEnd, CapSpu$CapSpu, FusNuc$FusNuc, 
                                  GemHae$GemHae, GemMor$GemMor, GemSan$GemSan, HaeHae$HaeHae, HaePar$HaePar,
                                  LepBuc$LepBuc, LepHof$LepHof, LepHon$LepHon ,LepWad$LepWad,
                                  NeiMuc$NeiMuc, NeiPol$NeiPol, NeiSic$NeiSic, NeiSub$NeiSub,
                                  PorAsa$PorAsa, PreJej_II$PreJej_III, PreMel_II$PreMel_III, PreOri$PreOri, 
                                  PreRum$PreRum, PrevEno$PrevEno, PseGoo$PseGoo, RotMuc$RotMuc, 
                                  SchMey$SchMey, SchOdo$SchOdo,
                                  VeiAty$VeiAty, VeiDis$VeiDis, VeiPar$VeiPar, VeiRod$VeiRod))
merge_species_2 <- merge_species
merge_species_2$pathway <- NULL
merge_species_2 <- ddply(merge_species_2,"id",numcolwise(max))
rownames(merge_species_2) <- merge_species_2$id
merge_species_2$id <- NULL
merge_species_3 <- data.frame(t(merge_species_2))

merge_species_3$species <- c("rare", "core", "rare", "rare", "rare", "core", "rare", "core", "core", "rare",
                             "rare", "rare", "rare", "core", "core", "rare", "rare", "core","rare", "core",
                             "core", "rare", "rare", "rare", "rare", "core", "core", 
                             "rare", "core", "core", "rare", "rare")
merge_species_3 <- ddply(merge_species_3,"species",numcolwise(median))
rownames(merge_species_3) <- merge_species_3$species
merge_species_3$species <- NULL
merge_species_4 <- merge_species_3[, colSums(merge_species_3) > 0]
merge_species_5 <- data.frame(t(merge_species_4))

merge_species_6 <- subset(merge_species_5, core != rare)
merge_species_6$diff <- merge_species_6$core - merge_species_6$rare
merge_species_6$diff2 <- ifelse(merge_species_6$diff < 0, (merge_species_6$diff * -1), merge_species_6$diff)
merge_species_7 <- subset(merge_species_6, diff2 > 24.0) # 24
merge_species_7$diff2 <- NULL
merge_species_7$diff <- NULL
rownames(merge_species_7) <- str_replace(rownames(merge_species_7), "^X", "")
rownames(merge_species_7) <- str_replace_all(rownames(merge_species_7), "\\.", "-")
colnames(merge_species_7) <- c("Core species", "Rare species")

plot_heat_all <- pheatmap(merge_species_7, scale = "none", cutree_rows = 5, cutree_cols = 2,
                          treeheight_col = 0, cluster_rows = TRUE, show_rownames = TRUE, 
                          show_colnames = TRUE, angle_col = 0, fontsize_row = 6)

#all_sum <- subset(merge_species, id %in% important_ids_3)
pheatmap_rows <- rownames(merge_species_7[plot_heat_all$tree_row[["order"]],])
pheatmap_rows <- str_replace(pheatmap_rows, "^X", "")
pheatmap_rows <- str_replace_all(pheatmap_rows, "\\.", "_")
pheatmap_rows <- str_replace_all(pheatmap_rows, "-", "_")
merge_species$id <- str_replace(merge_species$id, "^X", "")
merge_species$id <- str_replace_all(merge_species$id, "\\.", "_")
merge_species$id <- str_replace_all(merge_species$id, "-", "_")
merge_species_sub_1 <- subset(merge_species, id %in% pheatmap_rows)
merge_species_sub_2 <- merge_species_sub_1$id
names(merge_species_sub_2) <- merge_species_sub_1$pathway

get_names = NULL
for(items in pheatmap_rows){
  a = (names(merge_species_sub_2)[merge_species_sub_2 == items])
  get_names <- rbind(get_names, data.frame(a))}
pheatmap_rows_df <- data.frame(cbind(pheatmap_rows, get_names$a))

pheatmap_rows_df$pheatmap_rows <- str_replace_all(pheatmap_rows_df$pheatmap_rows, "_", "-")
colnames(pheatmap_rows_df) <- c("Pathway_id", "Pathway")

pheatmap_rows_df$Superclasses <-
  c("Nucleoside and Nucleotide Metabolism",              # PWY-7199
    "Fatty Acid and Lipid Biosynthesis",                 # PWY-7755
    "Amino Acid Metabolism",                             # THREONINE-DEG2-PWY
    "Amino Acid Metabolism",                             # PWY0-823
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY0-1507
    "Carbohydrate Metabolism",                           # PWY-7581
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-7380
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-7147
    "Inorganic Nutrient Metabolism",                     # PWY-6593
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-6578
    "Aromatic Compound Metabolism",                      # PWY-6077
    "Amino Acid Metabolism",                             # GLUTAMINDEG-PWY
    "Amine and Polyamine Metabolism",                    # PWY-40
    "Carbohydrate Metabolism",                           # FUCCAT-PWY
    "Carbohydrate Metabolism",                           # PWY-7178
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-8128
    "Generation of Precursor Metabolites and Energy",    # PWY3O-440
    "Carbohydrate Metabolism",                           # PWY0-1324 !!!!!
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-7970
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-7969, Cobamide biosynthesis
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-7968, Cobamide biosynthesis
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-7967
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-7966
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-7965
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-7964 
    "Nucleoside and Nucleotide Metabolism",              # PWY-7193
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-5509
    "Carboxylate Degradation",                           # PWY-5177
    "Amine and Polyamine Metabolism",                    # PWY-43
    "Fatty Acid and Lipid Degradation",                  # PWY-4261
    "Carboxylate Degradation",                           # PROPIONMET-PWY
    "Amino Acid Metabolism",                             # ARGDEG-III-PWY
    "Generation of Precursor Metabolites and Energy",    # ENTNER-DOUDOROFF-PWY
    "Detoxification",                                    # PWY0-1299 !!!!
    "Amino Acid Metabolism",                             # ASPARAGINESYN-PWY
    "Generation of Precursor Metabolites and Energy",    # PWY-7980
    "Fatty Acid and Lipid Biosynthesis",                 # PWY-4381
    "Generation of Precursor Metabolites and Energy",    # PWY-6938-1
    "Generation of Precursor Metabolites and Energy",    # PWY-8136
    "Amino Acid Metabolism",                             # PWY-5436
    "Nucleoside and Nucleotide Metabolism",              # PWY-7205
    "Cell Structure Biosynthesis",                       # PWY-8133
    "Nucleoside and Nucleotide Metabolism",              # PWY-6609
    "Carbohydrate Metabolism",                           # DTDPRHAMSYN-PWY
    "Carbohydrate Metabolism",                           # PWY-5659
    "Amine and Polyamine Metabolism",                    # PWY-3641
    "Amino Acid Metabolism",                             # SERSYN-PWY
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-2161
    "Generation of Precursor Metabolites and Energy",    # PWY-3781
    "Carbohydrate Metabolism",                           # PWY-5661
    "Nucleoside and Nucleotide Metabolism",              # SALVPURINE2-PWY
    "Nucleic Acid Processing",                           # PWY0-1587
    "Nucleoside and Nucleotide Metabolism",              # PWY-7226
    "Amino Acid Metabolism",                             # CYSTSYN-PWY
    "Amino Acid Metabolism",                             # GLUTDEG-PWY 
    "Nucleic Acid Processing",                           # PWY-6700, Queuosine Biosynthesis and Salvage
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-6614, folate biosynthesis
    "Generation of Precursor Metabolites and Energy",    # PYRUVDEHYD-PWY
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-6910
    "Generation of Precursor Metabolites and Energy",    # PWY-6785
    "Generation of Precursor Metabolites and Energy",    # PWY-6759
    "Amino Acid Metabolism",                             # ALANINE-SYN2-PWY
    "Detoxification",                                    # FORMASS-PWY
    "Nucleoside and Nucleotide Metabolism",              # PWY-8131
    "Carbohydrate Metabolism",                           # PWY-7774
    "Aromatic Compound Metabolism",                      # PWY-7773
    "Aromatic Compound Metabolism",                      # PWY-7772
    "Generation of Precursor Metabolites and Energy",    # PWY-7564, toxin biosynthesis !!
    "Aromatic Compound Metabolism",                      # PWY-7081
    "Aromatic Compound Metabolism",                      # PWY-6210
    "Aromatic Compound Metabolism",                      # PWY-6087
    "Aromatic Compound Metabolism",                      # PWY-6041
    "Chlorinated Compound Degradation",                  # 12DICHLORETHDEG-PWY
    "Carboxylate Degradation",                           # PWY-5654
    "Nucleoside and Nucleotide Metabolism",              # PWY0-1471
    "Carbohydrate Metabolism",                           # PWY-8121
    "Carbohydrate Metabolism",                           # PWY-8089
    "Detoxification",                                    # PWY-8043
    "Amino Acid Metabolism",                             # PWY-7870
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-6987
    "Inorganic Nutrient Metabolism",                     # PWY-6523
    "Cell Structure Biosynthesis",                       # PWY-6454, Vancomycin Resistance !!!
    "Nucleoside and Nucleotide Metabolism",              # PWY-6430
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-5523
    "Amino Acid Metabolism",                             # PWY-3982
    "Aromatic Compound Metabolism",                      # P343-PWY
    "Aromatic Compound Metabolism",                      # 4TOLCARBDEG-PWY
    "Amino Acid Metabolism",                             # ARGDEG-V-PWY
    "Aromatic Compound Metabolism",                      # PWY-5487
    "Aromatic Compound Metabolism",                      # PWY-7002
    "Aromatic Compound Metabolism",                      # TOLSULFDEG-PWY
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # HEME-BIOSYNTHESIS-II
    "Carbohydrate Metabolism",                           # PWY-7775
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-7356
    "Generation of Precursor Metabolites and Energy",    # RUMP-PWY
    "Generation of Precursor Metabolites and Energy",    # OXIDATIVEPENT-PWY
    "Amino Acid Metabolism",                             # TRPKYNCAT-PWY
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY0-501
    "Generation of Precursor Metabolites and Energy",    # PWY0-1576, anerobic respiration
    "Generation of Precursor Metabolites and Energy",    # PWY0-1336, anerobic respiration
    "Amino Acid Metabolism",                             # PWY-8120
    "Carbohydrate Metabolism",                           # PWY-7335
    "Carbohydrate Metabolism",                           # PWY-7334
    "Carbohydrate Metabolism",                           # PWY-7333
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-6894
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-6893
    "Cofactor, Carrier, and Vitamin Biosynthesis",       # PWY-6890
    "Aromatic Compound Metabolism",                      # PWY-6406, Salicylate
    "Fatty Acid and Lipid Biosynthesis",                 # PWY-5669
    "Amino Acid Metabolism",                             # CITRULLINE-DEG-PWY
    "Cofactor, Carrier, and Vitamin Biosynthesis")       # GLUTATHIONESYN-PWY

pheatmap_rows_df <- pheatmap_rows_df[order(pheatmap_rows_df$Pathway_id),]
merge_species_7 <- merge_species_7[order(rownames(merge_species_7)),]
core_rare_metabolic_differences <- data.frame(cbind(pheatmap_rows_df, merge_species_7))
rownames(core_rare_metabolic_differences) <- NULL

# Show Prevotella example
prev_names <- c("PreJej_II.PreJej_III", "PreMel_II.PreMel_III", "PreOri.PreOri", "PreRum.PreRum", "PrevEno.PrevEno")
merge_prev_1 <- select(merge_species_2, all_of(prev_names))
merge_prev_2 <- data.frame(t(merge_prev_1))
merge_prev_3 <- merge_prev_2[, colSums(merge_prev_2) > 0]
merge_prev_3$type <- c("core", "core", "rare", "rare", "rare")
rownames(merge_prev_3) <- c("Prevotella jejuni", "Prevotella melaninogenica", "Prevotella oris", "Prevotella ruminicola", 
                            "Prevotella enoeca")

pca_patternRec_prev <- prcomp(merge_prev_3[,1:(ncol(merge_prev_3)-1)], center = TRUE, scale. = FALSE)

pcaInd_patternRec_prev <- fviz_pca_ind(pca_patternRec_prev, geom.ind = c("text", "point"), pointsize = 2,
                                       palette = c("orange", "darkblue"), fill.ind = merge_prev_3$type, 
                                       col.ind = merge_prev_3$type, addEllipses = FALSE, ellipse.type = "confidence",
                                       ellipse.level=0.95, repel = TRUE, legend.title = "") + ggtitle("") +
  theme_pubr(border = TRUE, base_size = 12, legend = "none")

pcaVAR_patternRec_prev <- fviz_pca_var(pca_patternRec_prev, select.var = list(contrib = 50), col.var = "contrib", 
                                       geom = c("point", "text"), labelsize = 1, gradient.cols = c(
                                         "#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) +
  theme_pubr(border=TRUE, base_size = 12, legend = "none") +
  ggtitle(" ") + xlab("Dim1") + ylab("Dim2") 

# Extract variable contributions of PCA
varContribution_50_prev <- data.frame(pcaVAR_patternRec_prev$data[1:50,])
# rename columns
colnames(varContribution_50_prev) <- c("variable", "Dim1", "Dim2", "coord", "cos2", "contrib")
# remove rownames
rownames(varContribution_50_prev) <- NULL
# order table
varContribution_50_prev <- varContribution_50_prev[order(-varContribution_50_prev$contrib),]
# reset row numbering
rownames(varContribution_50_prev) <- NULL
varContribution_50_prev$variable
varContribution_50_prev$superclass <-
  c("Carbohydrate degradation",                       # BGALACT-PWY 
    "Carboxylate Degradation",                        # GLUCONSUPER.PWY
    "Polyprenyl Biosynthesis",                        # PWY.5122
    "Cofactor, Carrier, and Vitamin Biosynthesis",    # PWY.5964 
    "Inorganic Nutrient Metabolism",                  # PWY.6348 !!!!
    "Pyruvate fermentation to 2,3-butanediol",        # PWY.6392
    "Cofactor, Carrier, and Vitamin Biosynthesis",    # PWY.6466
    "Nucleoside and Nucleotide Degradation",          # PWY.6756
    "Carboxylate Degradation",                        # PWY.7247 
    "Amino Acid Degradation",                         # PWY1.2
    "Amino Acid Degradation",                         # GLUTAMINDEG.PWY
    "Amino Acid Biosynthesis",                        # PWY.5155
    "Cofactor, Carrier, and Vitamin Biosynthesis",    # PWY.8097
    "Amino Acid Biosynthesis",                        # ASPARAGINE.BIOSYNTHESIS
    "Amino Acid Degradation",                         # GLUTAMINEFUM.PWY
    "Amino Acid Biosynthesis",                        # GLUTSYN.PWY
    "Amino Acid Biosynthesis",                        # PWY.4341
    "Carbohydrate Biosynthesis",                      # PWY-4861
    "Generation of Precursor Metabolites and Energy", # PWY.5535 
    "Carbohydrate Biosynthesis",                      # PWY.63
    "Carbohydrate Degradation",                       # PWY.6717 !!!
    "Carbohydrate Biosynthesis",                      # PWY.7127D
    "Carbohydrate Biosynthesis",                      # PWY.7330
    "Carbohydrate Degradation",                       # PWY.7568
    "Carboxylate Degradation",                        # PWY.7686  
    "Protein Modification",                           # PWY.7801
    "Carboxylate Degradation",                        # PWY0.1306
    "Nucleoside and Nucleotide Biosynthesis",         # PWY.7226
    "Nucleoside and Nucleotide Biosynthesis",         # PWY.7227
    "Amino Acid Biosynthesis",                        # HISTSYN.PWY
    "Amino Acid Degradation",                         # SERDEG.PWY
    "Nucleoside and Nucleotide Degradation",          # PWY.7177
    "Carboxylate Degradation",                        # IDNCAT.PWY 
    "Secondary Metabolite Degradation",               # PWY.6507
    "Metabolic Regulator Biosynthesis",               # PPGPPMET.PWY
    "Generation of Precursor Metabolites and Energy", # PWY.8015
    "Generation of Precursor Metabolites and Energy", # PWY0.1312
    "Amino Acid Degradation",                         # HISDEG-PWY
    "Carbohydrate Degradation",                       # RHAMCAT-PWY 
    "Amino Acid Biosynthesis",                        # ASPARAGINESYN.PWY
    "Amino Acid Degradation",                         # ASPARTATE.DEG1.PWY 
    "Amino Acid Biosynthesis",                        # ASPARTATESYN.PWY
    "Amino Acid Biosynthesis",                        # CYSTSYN.PWY
    "Amino Acid Biosynthesis",                        # GLNSYN.PWY
    "Amino Acid Biosynthesis",                        # PGLYSYN.ALA.PWY
    "Detoxification",                                 # P401.PWY
    "Generation of Precursor Metabolites and Energy", # PWY.7980
    "Inorganic Nutrient Metabolism",                  # PWY.6964
    "Amino Acid Biosynthesis",                        # GLUTORN.PWY
    "Carboxylate Degradation")                        # GALACTUROCAT.PWY
varContribution_50_prev$genus <- "Prevotella"


# Show Leptotrichia example
lepto_names <- c("LepBuc.LepBuc", "LepHof.LepHof", "LepHon.LepHon", "LepWad.LepWad")
merge_lepto_1 <- select(merge_species_2, all_of(lepto_names))
merge_lepto_2 <- data.frame(t(merge_lepto_1))
merge_lepto_3 <- merge_lepto_2[, colSums(merge_lepto_2) > 0]
merge_lepto_3$type <- c("rare", "rare", "rare", "core")
rownames(merge_lepto_3) <- c("Leptotrichia buccalis", "Leptotrichia hofstadii", 
                             "Leptotrichia hongkongensis", "Leptotrichia wadei")

pca_patternRec_lepto <- prcomp(merge_lepto_3[,1:(ncol(merge_lepto_3)-1)], center = TRUE, scale. = FALSE)

pcaInd_patternRec_lepto <- fviz_pca_ind(pca_patternRec_lepto, geom.ind = c("text", "point"), pointsize = 2,
                                        palette = c("orange", "darkblue"), fill.ind = merge_lepto_3$type, 
                                        col.ind = merge_lepto_3$type, addEllipses = FALSE, ellipse.type = "confidence",
                                        ellipse.level=0.95, repel = TRUE, legend.title = "") + ggtitle("") +
  theme_pubr(border = TRUE, base_size = 12, legend = "none")

pcaVAR_patternRec_lepto <- fviz_pca_var(pca_patternRec_lepto, select.var = list(contrib = 50), col.var = "contrib", 
                                        geom = c("point", "text"), labelsize = 1, gradient.cols = c(
                                          "#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) +
  theme_pubr(border=TRUE, base_size = 12, legend = "none") +
  ggtitle(" ") + xlab("Dim1") + ylab("Dim2") 

# Extract variable contributions of PCA
varContribution_50_lepto <- data.frame(pcaVAR_patternRec_lepto$data[1:50,])
# rename columns
colnames(varContribution_50_lepto) <-
  c("variable", "Dim1", "Dim2", "coord", "cos2", "contrib")
# remove rownames
rownames(varContribution_50_lepto) <- NULL
# order table
varContribution_50_lepto <- varContribution_50_lepto[order(-varContribution_50_lepto$contrib),]
# reset row numbering
rownames(varContribution_50_lepto) <- NULL
varContribution_50_lepto$superclass <- c(
  "Amino Acid Degradation",                      # ASPARAGINE-DEG1-PWY 
  "Carbohydrate Degradation",                    # MANNCAT-PWY         
  "Cofactor, Carrier, and Vitamin Biosynthesis", # PWY-6466            
  "Amino Acid Degradation",                      # GLUTAMINEFUM-PWY    
  "Amino Acid Biosynthesis",                     # GLUTSYN-PWY         
  "Amino Acid Biosynthesis",                     # PWY-4341
  "Cofactor, Carrier, and Vitamin Biosynthesis", # PWY-8097
  "Aromatic Compound Degradation",               # PWY-6781
  "Nucleoside and Nucleotide Biosynthesis",      # PWY-7183
  "Carbohydrate Degradation",                    # PWY0-1301
  "Carbohydrate Degradation",                    # FUCCAT-PWY
  "Carbohydrate Degradation",                    # LYXMET-PWY
  "Carboxylate Degradation",                     # PWY0-301
  "Carbohydrate Degradation",                    # ARABCAT-PWY
  "Carboxylate Degradation",                     # PWY-6961
  "Carbohydrate Degradation",                    # DARABCATK12-PWY
  "Carbohydrate Degradation",                    # TREDEGLOW-PWY
  "Carbohydrate Biosynthesis",                   # GDPRHAMSYN-PWY
  "Carbohydrate Degradation",                    # PWY0-1182
  "Amino Acid Biosynthesis",                     # HISTSYN-PWY
  "Nucleoside and Nucleotide Degradation",       # PWY-7177
  "Carbohydrate Degradation",                    # DARABCAT-PWY
  "Aromatic Compound Degradation",               # PWY-6077
  "Nucleoside and Nucleotide Biosynthesis",      # SALVPURINE2-PWY
  "Protein Modification",                        # PWY-7942
  "Carbohydrate Degradation",                    # PWY-2722
  "Secondary Metabolite Biosynthesis",           # PWY-3161, hormone biosynthesis
  "Secondary Metabolite Biosynthesis",           # PWY-5025
  "Methylglyoxal Detoxification",                # PWY-5453
  "Carbohydrate Biosynthesis",                   # PWY-5661
  "Carbohydrate Biosynthesis",                   # PWY-5739
  "Carbohydrate Biosynthesis",                   # PWY-6478
  "Cofactor, Carrier, and Vitamin Biosynthesis", # PWY-6890
  "Cofactor, Carrier, and Vitamin Biosynthesis", # PWY-6893
  "Cofactor, Carrier, and Vitamin Biosynthesis", # PWY-6894
  "Inorganic Nutrient Metabolism",               # PWY-6964
  "Acrylonitrile Degradation",                   # PWY-7308
  "Cofactor, Carrier, and Vitamin Biosynthesis", # PWY-8096
  "Cofactor, Carrier, and Vitamin Biosynthesis", # PWY1G-126
  "Amino Acid Degradation",                      # GLUTAMINDEG-PWY
  "Amino Acid Biosynthesis",                     # GLYSYN-THR-PWY
  "Nucleoside and Nucleotide Biosynthesis",      # PWY-6599
  "Fatty Acid and Lipid Biosynthesis",           # PWY-6795
  "Cofactor, Carrier, and Vitamin Biosynthesis", # PWY-7380
  "Polymeric Compound Degradation",              # PWY-7794
  "Cofactor, Carrier, and Vitamin Biosynthesis", # PWY0-1507
  "Fatty Acid and Lipid Degradation",            # PWY-5533
  "Carbohydrate Degradation",                    # PWY-6527
  "Amino Acid Biosynthesis",                     # GLUTORN-PWY
  "Cofactor, Carrier, and Vitamin Biosynthesis") # PLPSAL-PWY

varContribution_50_lepto$genus <- "Leptotrichia"

# Show Veillonella example
veill_names <- c("VeiRod.VeiRod", "VeiPar.VeiPar", "VeiDis.VeiDis", "VeiAty.VeiAty")
merge_veill_1 <- select(merge_species_2, all_of(veill_names))
merge_veill_2 <- data.frame(t(merge_veill_1))
merge_veill_3 <- merge_veill_2[, colSums(merge_veill_2) > 0]
merge_veill_3$type <- c("rare", "rare", "core", "core")
rownames(merge_veill_3) <- c("Veillonella rodentium", "Veillonella parvula", 
                             "Veillonella dispar", "Veillonella atypica")

pca_patternRec_veill <- prcomp(merge_veill_3[,1:(ncol(merge_veill_3)-1)], center = TRUE, scale. = FALSE)

pcaInd_patternRec_veill <-fviz_pca_ind(pca_patternRec_veill, geom.ind = c("text", "point"), pointsize = 2,
                                       palette = c("orange", "darkblue"),
                                       fill.ind = merge_veill_3$type, col.ind = merge_veill_3$type,
                                       addEllipses = FALSE, ellipse.type = "confidence", ellipse.level=0.95,
                                       repel = TRUE, legend.title = "") + ggtitle("") + 
  theme_pubr(border = TRUE, base_size = 12, legend = "none")

pcaVAR_patternRec_veill <- fviz_pca_var(pca_patternRec_veill, select.var = list(contrib = 50), col.var = "contrib", 
                                        geom = c("point", "text"), labelsize = 1, gradient.cols = c(
                                          "#00AFBB","#E7B800", "#FC4E07"), repel = TRUE) +
  theme_pubr(border=TRUE, base_size = 12, legend = "none") +
  ggtitle(" ") + xlab("Dim1") + ylab("Dim2") 

# Extract variable contributions of PCA
varContribution_50_veill <- data.frame(pcaVAR_patternRec_veill$data[1:50,])
# rename columns
colnames(varContribution_50_veill) <- c("variable", "Dim1", "Dim2", "coord", "cos2", "contrib")
# remove rownames
rownames(varContribution_50_veill) <- NULL
# order table
varContribution_50_veill <- varContribution_50_veill[order(-varContribution_50_veill$contrib),]
# reset row numbering
rownames(varContribution_50_veill) <- NULL

varContribution_50_veill$superclass <- c(
  "Cofactor, Carrier, and Vitamin Biosynthesis",           #PWY-8097
  "Amino Acid Biosynthesis",                               #GLYSYN-ALA-PWY
  "Lactaldehyde Degradation",                              #PWY0-1315
  "Carbohydrate Biosynthesis",                             #PWY-5512
  "Carbohydrate Biosynthesis",                             #PWY-63
  "Carbohydrate Biosynthesis",                             #PWY-7344
  "Generation of Precursor Metabolites and Energy",        #PWY-6130
  "Amine and Polyamine Degradation",                       #PWY-7167
  "Reactive Oxygen Species Degradation",                   #DETOX1-PWY
  "Protein Modification",                                  #PWY-7942
  "Inorganic Nutrient Metabolism",                         #PWY-8062
  "Alcohol Degradation",                                   #ETOH-ACETYLCOA-ANA-PWY
  "Cofactor, Carrier, and Vitamin Biosynthesis",           #PWY-6890
  "Cofactor, Carrier, and Vitamin Biosynthesis",           #PWY-8096
  "Amine and Polyamine Biosynthesis",                      #ARGSPECAT-PWY
  "Secondary Metabolite Biosynthesis",                     #PWY-5394
  "Alcohol Degradation",                                   #PWY0-1280
  "Amino Acid Degradation",                                #THRDLCTCAT-PWY
  "Amino Acid Degradation",                                #THREONINE-DEG2-PWY
  "Amino Acid Degradation",                                #CITRULLINE-DEG-PWY
  "Generation of Precursor Metabolites and Energy",        #PWY-7013
  "Amino Acid Biosynthesis",                               #PWY-3982
  "Nucleoside and Nucleotide Degradation",                 #PWY-6430
  "Carbohydrate Degradation",                              #PWY-7180
  "Inorganic Nutrient Metabolism",                         #PHOSPHONOTASE-PWY
  "Carboxylate Degradation",                               #PWY-5162
  "Amino Acid Degradation",                                #PWY-5436
  "Generation of Precursor Metabolites and Energy",        #PWY-5480
  "Amine and Polyamine Degradation",                       #PWY-5698
  "Generation of Precursor Metabolites and Energy",        #PWY-6587
  "Generation of Precursor Metabolites and Energy ",       #PWY-7541
  "Carbohydrate Degradation",                              #PWY-8060
  "Cofactor, Carrier, and Vitamin Biosynthesis",           #PWY-7761
  "Cofactor, Carrier, and Vitamin Biosynthesis",           #PWY-6268
  "Amino Acid Degradation",                                #GLUTAMINDEG-PWY
  "Amino Acid Degradation",                                #GLYCLEAV-PWY
  "Generation of Precursor Metabolites and Energy",        #PWY-5046
  "Generation of Precursor Metabolites and Energy",        #PWY-5084               
  "Nucleoside and Nucleotide Biosynthesis",                #PWY-6610               
  "Carbohydrate Degradation",                              #PWY-8089
  "Cofactor, Carrier, and Vitamin Biosynthesis",           #PWY-8095
  "Nucleic Acid Processing",                               #PWY0-1554
  "Generation of Precursor Metabolites and Energy",        #PYRUVDEHYD-PWY
  "Fatty Acid and Lipid Biosynthesis",                     #PWY-5994
  "Amino Acid Biosynthesis",                               #PWY-3461
  "Amino Acid Biosynthesis",                               #PWY-3462
  "Amino Acid Biosynthesis",                               #PWY-6120
  "Cofactor, Carrier, and Vitamin Biosynthesis",           #PWY-7378
  "Cofactor, Carrier, and Vitamin Biosynthesis",           #PWY-5523
  "Cofactor, Carrier, and Vitamin Biosynthesis")           #PWY-5653
varContribution_50_veill$genus <- "Veillonella"

# merge pca dataframes
varContribution_150_merge <- data.frame(rbind(
  varContribution_50_prev, varContribution_50_veill, varContribution_50_lepto))
varContribution_150_merge$Dim1 <- round(varContribution_150_merge$Dim1,1)
varContribution_150_merge$Dim2 <- round(varContribution_150_merge$Dim2,1)
varContribution_150_merge$cos2 <- round(varContribution_150_merge$cos2,1)
varContribution_150_merge$contrib <- round(varContribution_150_merge$contrib,1)
varContribution_150_merge$coord <- NULL
varContribution_150_merge$variable <- str_replace_all(varContribution_150_merge$variable, "\\.", "-")

#dendrogram
merge_species_all <- cbind(
  StrAus, StrEqu$StrEqu,StrGor$StrGor, StrHim$StrHim, StrMit$StrMit, StrOra$StrOra,
  StrPar$StrPar, StrPne$StrPne, StrSal$StrSal, StrVes$StrVes,
  VeiAty$VeiAty, VeiDis$VeiDis, VeiPar$VeiPar, VeiRod$VeiRod, SchOdo$SchOdo, SchMey$SchMey,
  RotMuc$RotMuc, PseGoo$PseGoo, PreJej_II$PreJej_III, PreMel_II$PreMel_III,
  PrevEno$PrevEno, PreRum$PreRum, PreOri$PreOri, PorAsa$PorAsa,
  NeiSub$NeiSub, NeiSic$NeiSic, NeiPol$NeiPol, NeiMuc$NeiMuc, LepWad$LepWad, LepHof$LepHof, LepHon$LepHon,
  LepBuc$LepBuc, HaePar$HaePar, HaeHae$HaeHae, GemSan$GemSan, GemMor$GemMor, GemHae$GemHae, FusNuc$FusNuc,
  CapSpu$CapSpu, CapEnd$CapEnd, CamCon$CamCon, ActIsr$ActIsr)

#merge_species_all$pathway_id <- NULL
merge_species_all$pathway_id <- paste(rownames(merge_species_all), "_", merge_species_all$id)
rownames(merge_species_all) <- merge_species_all$pathway_id
merge_species_all$pathway_id <- NULL
merge_species_all$id <- NULL
merge_species_all$pathway <- NULL
merge_species_all <- merge_species_all[rowSums(merge_species_all) != 0,]
merge_species_all2 <- merge_species_all / 100
colnames(merge_species_all2) <- map(strsplit(colnames(merge_species_all), split = "\\$"), 1)

d_all = as.dendrogram(hclust(dist(t(merge_species_all2), method = "canberra"), method="ward.D2"))
d_all_2 = dendro_data(d_all, type="rectangle")
labs_all <- label(d_all_2)

labs_all$label <- str_replace(labs_all$label, "CapSpu", "Capnocytophaga sputigena")
labs_all$label <- str_replace(labs_all$label, "CapEnd", "Capnocytophaga endodontalis")
labs_all$label <- str_replace(labs_all$label, "PreJej_II", "Prevotella jejuni")
labs_all$label <- str_replace(labs_all$label, "PreMel_II", "Prevotella melaninogenica")
labs_all$label <- str_replace(labs_all$label, "PorAsa", "Porphyromonas asaccharolytica")
labs_all$label <- str_replace(labs_all$label, "PreRum", "Prevotella ruminicola")
labs_all$label <- str_replace(labs_all$label, "PrevEno", "Prevotella enoeca")
labs_all$label <- str_replace(labs_all$label, "PreOri", "Prevotella oris")
labs_all$label <- str_replace(labs_all$label, "HaePar", "Haemophilus parainfluenzae")
labs_all$label <- str_replace(labs_all$label, "HaeHae", "Haemophilus haemolyticus")
labs_all$label <- str_replace(labs_all$label, "ActIsr", "Actinomyces israelii")
labs_all$label <- str_replace(labs_all$label, "SchOdo", "Schaalia odontolytica")
labs_all$label <- str_replace(labs_all$label, "SchMey", "Schaalia meyeri")
labs_all$label <- str_replace(labs_all$label, "GemHae", "Gemella haemolysans")
labs_all$label <- str_replace(labs_all$label, "GemSan", "Gemella sanguinis")
labs_all$label <- str_replace(labs_all$label, "GemMor", "Gemella morbillorum")
labs_all$label <- str_replace(labs_all$label, "StrGor", "Streptococcus gordonii")
labs_all$label <- str_replace(labs_all$label, "StrMit", "Streptococcus mitis")
labs_all$label <- str_replace(labs_all$label, "StrAus", "Streptococcus australis")
labs_all$label <- str_replace(labs_all$label, "StrPar", "Streptococcus parasanguinis")
labs_all$label <- str_replace(labs_all$label, "StrHim", "Streptococcus himalayensis")
labs_all$label <- str_replace(labs_all$label, "StrEqu", "Streptococcus equinus")
labs_all$label <- str_replace(labs_all$label, "StrSal", "Streptococcus salivarius")
labs_all$label <- str_replace(labs_all$label, "StrVes", "Streptococcus vestibularis")
labs_all$label <- str_replace(labs_all$label, "StrOra", "Streptococcus oralis")
labs_all$label <- str_replace(labs_all$label, "StrPne", "Streptococcus pneumoniae")
labs_all$label <- str_replace(labs_all$label, "NeiPol", "Neisseria polysaccharea")
labs_all$label <- str_replace(labs_all$label, "NeiSic", "Neisseria sicca")
labs_all$label <- str_replace(labs_all$label, "NeiSub", "Neisseria subflava")
labs_all$label <- str_replace(labs_all$label, "NeiMuc", "Neisseria mucosa")
labs_all$label <- str_replace(labs_all$label, "RotMuc", "Rothia mucilaginosa")
labs_all$label <- str_replace(labs_all$label, "CamCon", "Campylobacter concisus")
labs_all$label <- str_replace(labs_all$label, "VeiRod", "Veillonella rodentium")
labs_all$label <- str_replace(labs_all$label, "VeiAty", "Veillonella atypica")
labs_all$label <- str_replace(labs_all$label, "VeiDis", "Veillonella dispar")
labs_all$label <- str_replace(labs_all$label, "VeiPar", "Veillonella parvula")
labs_all$label <- str_replace(labs_all$label, "FusNuc", "Fusobacterium nucleatum")
labs_all$label <- str_replace(labs_all$label, "PseGoo", "Pseudoleptotrichia goodfellowii")
labs_all$label <- str_replace(labs_all$label, "LepWad", "Leptotrichia wadei")
labs_all$label <- str_replace(labs_all$label, "LepHon", "Leptotrichia hongkongensis")
labs_all$label <- str_replace(labs_all$label, "LepHof", "Leptotrichia hofstadii")
labs_all$label <- str_replace(labs_all$label, "LepBuc", "Leptotrichia buccalis")

labs_all$group <- with(
  labs_all, 
  ifelse(label == "Streptococcus gordonii", "rare",
         ifelse(label == "Streptococcus mitis", "core",
                ifelse(label == "Streptococcus australis", "core",
                       ifelse(label == "Streptococcus parasanguinis","core",
                              ifelse(label == "Streptococcus himalayensis","rare",
                                     ifelse(label == "Streptococcus equinus", "rare",
                                            ifelse(label == "Streptococcus salivarius", "core",
                                                   ifelse(label == "Streptococcus vestibularis", "rare",
                                                          ifelse(label == "Streptococcus oralis", "core",
                                                                 ifelse(label == "Streptococcus pneumoniae", "rare",
                                                                        ifelse(label == "Actinomyces israelii", "rare",
                                                                               ifelse(label == "Capnocytophaga sputigena", "rare",
                                                                                      ifelse(label == "Capnocytophaga endodontalis", "rare",
                                                                                             ifelse(label == "Fusobacterium nucleatum", "rare",
                                                                                                    ifelse(label == "Gemella morbillorum", "rare",
                                                                                                           ifelse(label == "Haemophilus haemolyticus", "rare",
                                                                                                                  ifelse(label == "Leptotrichia buccalis", "rare",
                                                                                                                         ifelse(label == "Leptotrichia hofstadii", "rare",
                                                                                                                                ifelse(label == "Leptotrichia hongkongensis", "rare",
                                                                                                                                       ifelse(label == "Neisseria polysaccharea", "rare",
                                                                                                                                              ifelse(label == "Neisseria sicca", "rare",
                                                                                                                                                     ifelse(label == "Porphyromonas asaccharolytica", "rare",
                                                                                                                                                            ifelse(label == "Prevotella enoeca", "rare",
                                                                                                                                                                   ifelse(label == "Prevotella oris", "rare",
                                                                                                                                                                          ifelse(label == "Prevotella ruminicola", "rare",
                                                                                                                                                                                 ifelse(label == "Pseudoleptotrichia goodfellowii", "rare",
                                                                                                                                                                                        ifelse(label == "Schaalia meyeri", "rare",
                                                                                                                                                                                               ifelse(label == "Veillonella parvula", "rare",
                                                                                                                                                                                                      ifelse(label == "Veillonella rodentium", "rare",
                                                                                                                                                                                                             ifelse(label == "Campylobacter concisus", "core",
                                                                                                                                                                                                                    ifelse(label == "Gemella haemolysans", "core", 
                                                                                                                                                                                                                           ifelse(label == "Gemella sanguinis", "core",
                                                                                                                                                                                                                                  ifelse(label == "Haemophilus parainfluenzae", "core",
                                                                                                                                                                                                                                         ifelse(label == "Leptotrichia wadei", "core",
                                                                                                                                                                                                                                                ifelse(label == "Neisseria mucosa", "core",
                                                                                                                                                                                                                                                       ifelse(label == "Neisseria subflava", "core",
                                                                                                                                                                                                                                                              ifelse(label == "Prevotella jejuni", "core",
                                                                                                                                                                                                                                                                     ifelse(label == "Prevotella melaninogenica", "core",
                                                                                                                                                                                                                                                                            ifelse(label == "Rothia mucilaginosa", "core",
                                                                                                                                                                                                                                                                                   ifelse( label == "Schaalia odontolytica", "core",
                                                                                                                                                                                                                                                                                           ifelse(label == "Veillonella atypica", "core",
                                                                                                                                                                                                                                                                                                  ifelse(label == "Veillonella dispar", "core", NA
                                                                                                                                                                                                                                                                                                  )))))))))))))))))))))))))))))))))))))))))))

all_hier <- ggplot(segment(d_all_2)) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = labs_all, aes(label = label, x = x, y = -230, colour=group), size=3, fontface='italic') +
  ylim(-400, 1600) +
  scale_colour_manual(values=c("orange3", "darkblue")) +
  theme_pubr(border=TRUE, legend = "none") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab("Canberra distance") + coord_flip()

# Show Streptococcus example
merge_species_strep <- data.frame(cbind(StrAus, StrEqu$StrEqu, StrGor$StrGor, StrHim$StrHim, StrMit$StrMit, 
                                        StrOra$StrOra, StrPar$StrPar, StrPne$StrPne, StrSal$StrSal, StrVes$StrVes))
merge_species_strep_2 <- merge_species_strep
merge_species_strep_2$pathway <- NULL
merge_species_strep_2 <- ddply(merge_species_strep_2,"id",numcolwise(max))
rownames(merge_species_strep_2) <- merge_species_strep_2$id
merge_species_strep_2$id <- NULL

merge_strepto_2 <- data.frame(t(merge_species_strep_2))
merge_strepto_3 <- merge_strepto_2[, colSums(merge_strepto_2) > 0]
merge_strepto_3$type <- c("core", "rare", "rare", "rare", "core", 
                          "core", "core", "rare", "core", "rare")
merge_strepto_3$cluster <- c("A", "B", "A", "C", "C", 
                             "C", "A", "C", "B", "B")
rownames(merge_strepto_3) <- c("Streptococcus australis", 
                               "Streptococcus equinus", 
                               "Streptococcus gordonii", 
                               "Streptococcus himalayensis", 
                               "Streptococcus mitis", 
                               "Streptococcus oralis",
                               "Streptococcus parasanguinis", 
                               "Streptococcus pneumoniae", 
                               "Streptococcus salivarius", 
                               "Streptococcus vestibularis")

pca_patternRec_strepto <- prcomp(merge_strepto_3[,1:(ncol(merge_strepto_3)-2)],center = TRUE, scale. = FALSE)

pcaInd_patternRec_strepto <- fviz_pca_ind(pca_patternRec_strepto, geom.ind = c("text", "point"), 
                                          pointsize = 2, palette = c("darkgreen", "black", "lightsalmon4"),
                                          fill.ind = merge_strepto_3$cluster, col.ind = merge_strepto_3$cluster,
                                          addEllipses = FALSE, ellipse.type = "confidence", ellipse.level=0.95,
                                          repel = TRUE, legend.title = "") +
  ggtitle("") + theme_pubr(border = TRUE, base_size = 12, legend = "none")

pcaVAR_patternRec_strepto <- fviz_pca_var(pca_patternRec_strepto, select.var = list(contrib = 50), col.var = "contrib", 
                                          geom = c("point", "text"), labelsize = 1, gradient.cols = c(
                                            "#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) +
  theme_pubr(border=TRUE, base_size = 12, legend = "none") +
  ggtitle(" ") + xlab("Dim1") + ylab("Dim2") 

pca_plots <- ggarrange(pcaInd_patternRec_prev, pcaInd_patternRec_veill,
                       pcaInd_patternRec_lepto, pcaInd_patternRec_strepto,
                       labels = c("a", "b", "c", "d"))

StrAus_stress <- read.csv("data_input/StrAus_stress.tbl", sep = "\t")
StrAus_stress$Species <- "Streptococcus australis"
StrAus_stress <- subset(StrAus_stress, Value < 1e-50)
StrAus_stress$Sequences1 <- map(strsplit(StrAus_stress$Sequences, split="\\|"), 1)

StrPar_stress <- read.csv("data_input/StrPar_stress.tbl", sep="\t")
StrPar_stress$Species <- "Streptococcus parasanguinis"
StrPar_stress <- subset(StrPar_stress, Value < 1e-50)
StrPar_stress$Sequences1 <- map(strsplit(StrPar_stress$Sequences, split="\\|"), 1)

StrAus_stress_vector <- unlist(StrAus_stress$Sequences1)
StrPar_stress_vector <- unlist(StrPar_stress$Sequences1)
unique_values_StrAus <- setdiff(StrAus_stress_vector, StrPar_stress_vector)
unique_values_StrPar <- setdiff(StrPar_stress_vector, StrAus_stress_vector)
unique_values <- c(unique_values_StrAus, unique_values_StrPar)
Str_stress <- data.frame(rbind(StrAus_stress, StrPar_stress))
Str_stress_2 <- Str_stress[Str_stress$Sequences1 %in% unique_values,]

ggplot(Str_stress_2, aes(x=Species, y=Sequences, fill=Score)) + geom_tile() +
  theme_pubr(border=TRUE, legend = "right") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

StrAus_adhesion <- read.csv("data_input/StrAus_adhesion.tbl", sep = "\t")
StrAus_adhesion$Species <- "S. australis"
StrAus_adhesion <- subset(StrAus_adhesion, Value < 0.01)
StrAus_adhesion$Sequences1 <- map(strsplit(StrAus_adhesion$Sequences, split="\\|"), 1)

StrPar_adhesion <- read.csv("data_input/StrPar_adhesion.tbl", sep="\t")
StrPar_adhesion$Species <- "S. parasanguinis"
StrPar_adhesion <- subset(StrPar_adhesion, Value < 0.01)
StrPar_adhesion$Sequences1 <- map(strsplit(StrPar_adhesion$Sequences, split="\\|"), 1)

StrAus_adhesion_vector <- unlist(StrAus_adhesion$Sequences1)
StrPar_adhesion_vector <- unlist(StrPar_adhesion$Sequences1)
unique_values_StrAus_adhesion <- setdiff(StrAus_adhesion_vector, StrPar_adhesion_vector)
unique_values_StrPar_adhesion <- setdiff(StrPar_adhesion_vector, StrAus_adhesion_vector)
unique_values_adhesion <- c(unique_values_StrAus_adhesion, unique_values_StrPar_adhesion)
Str_adhesion <- data.frame(rbind(StrAus_adhesion, StrPar_adhesion))
Str_adhesion_2 <- Str_adhesion[Str_adhesion$Sequences1 %in% unique_values_adhesion,]
Str_adhesion_2$Sequences3 <- c(
  "Bone sialoprotein-binding protein", #Q6GJA6.1|BBP_STAAR|Full=Bone sialoprotein-binding p
  "Bone sialoprotein-binding protein", # Q6GJA6.1|BBP_STAAR|Full=Bone sialoprotein-binding p
  "Bone sialoprotein-binding protein", # Q6GJA6.1|BBP_STAAR|Full=Bone sialoprotein-binding p
  "Bone sialoprotein-binding protein", # Q14U76.1|BBP_STAAU|Full=Bone sialoprotein-binding p
  "Bone sialoprotein-binding protein", # Q14U76.1|BBP_STAAU|Full=Bone sialoprotein-binding p
  "Bone sialoprotein-binding protein", # Q14U76.1|BBP_STAAU|Full=Bone sialoprotein-binding p
  "Glucosyltransferase 3", # B5A7L9.1|GTF3_STRPA|Full=Glucosyltransferase 3
  "Glucosyltransferase 3", # B5A7L9.1|GTF3_STRPA|Full=Glucosyltransferase 3     
  "Glucosyltransferase 3", # B5A7L9.1|GTF3_STRPA|Full=Glucosyltransferase 3
  "Sensor protein VraS", # Q6G849.1|VRAS_STAAS|Full=Sensor protein VraS
  "Sensor protein VraS", # Q6GFH2.1|VRAS_STAAR|Full=Sensor protein VraS
  "Sensor protein VraS", # Q5HEN9.1|VRAS_STAAC|Full=Sensor protein VraS
  "Sensor protein VraS", # Q7A0H9.1|VRAS_STAAW|Full=Sensor protein VraS
  "Sensor protein VraS", # Q99SZ7.1|VRAS_STAAN|Full=Sensor protein VraS
  "Sensor protein VraS", # Q6G849.1|VRAS_STAAS|Full=Sensor protein VraS
  "Sensor protein VraS", # Q6GFH2.1|VRAS_STAAR|Full=Sensor protein VraS
  "Sensor protein VraS", # Q5HEN9.1|VRAS_STAAC|Full=Sensor protein VraS
  "Sensor protein VraS", # Q7A0H9.1|VRAS_STAAW|Full=Sensor protein VraS
  "Sensor protein VraS", # Q99SZ7.1|VRAS_STAAN|Full=Sensor protein VraS
  "Sensor protein VraS", # Q6G849.1|VRAS_STAAS|Full=Sensor protein VraS
  "Sensor protein VraS", # Q6GFH2.1|VRAS_STAAR|Full=Sensor protein VraS
  "Sensor protein VraS", # Q5HEN9.1|VRAS_STAAC|Full=Sensor protein VraS
  "Sensor protein VraS", # Q7A0H9.1|VRAS_STAAW|Full=Sensor protein VraS
  "Sensor protein VraS", # Q99SZ7.1|VRAS_STAAN|Full=Sensor protein VraS 
  "Sensor protein VraS", # Q9KWK8.1|VRAS_STAA1|Full=Sensor protein VraS
  "Sensor protein VraS", # Q7A2Q0.1|VRAS_STAAM|Full=Sensor protein VraS 
  "Sensor protein VraS", # Q9KWK8.1|VRAS_STAA1|Full=Sensor protein VraS
  "Sensor protein VraS", # Q7A2Q0.1|VRAS_STAAM|Full=Sensor protein VraS
  "Sensor protein VraS", # Q9KWK8.1|VRAS_STAA1|Full=Sensor protein VraS
  "Sensor protein VraS", # Q7A2Q0.1|VRAS_STAAM|Full=Sensor protein VraS
  "Agglutinin receptor", # P16952.2|SSP5_STRGN|Full=Agglutinin receptor|AltNa
  "Agglutinin receptor", # P16952.2|SSP5_STRGN|Full=Agglutinin receptor|AltNa
  "Agglutinin receptor", # P16952.2|SSP5_STRGN|Full=Agglutinin receptor|AltNa
  "Probable cell-surface antigen I/II", # Q9LBG3.1|PAA_STRCG|Full=Probable cell-surface antig
  "Probable cell-surface antigen I/II", # Q9LBG3.1|PAA_STRCG|Full=Probable cell-surface antig
  "Probable cell-surface antigen I/II", # Q9LBG3.1|PAA_STRCG|Full=Probable cell-surface antig
  "Protease PrtH", # P46071.1|PRTH_PORGI|Full=Protease PrtH
  "Protease PrtH", # P46071.1|PRTH_PORGI|Full=Protease PrtH
  "Protease PrtH", # P46071.1|PRTH_PORGI|Full=Protease PrtH
  "Probable cell-surface antigen I/II", # Q9KW51.2|PAS_STRIT|Full=Probable cell-surface antig
  "Probable cell-surface antigen I/II", # Q9KW51.2|PAS_STRIT|Full=Probable cell-surface antig
  "Probable cell-surface antigen I/II", # Q9KW51.2|PAS_STRIT|Full=Probable cell-surface antig
  "Internalin J", # Q8Y3L4.1|INLJ_LISMO|Full=Internalin J|Flags: Precu
  "Internalin J", # Q8Y3L4.1|INLJ_LISMO|Full=Internalin J|Flags: Precu
  "Internalin J", # Q8Y3L4.1|INLJ_LISMO|Full=Internalin J|Flags: Precu
  "Transcriptional regulator IcaR", # A6QKF4.1|ICAR_STAAE|Full=Biofilm operon icaADBC HTH
  "Transcriptional regulator IcaR", # A6QKF4.1|ICAR_STAAE|Full=Biofilm operon icaADBC HTH
  "Transcriptional regulator IcaR", # A6QKF4.1|ICAR_STAAE|Full=Biofilm operon icaADBC HTH
  "CFA/I fimbrial subunit D", # P25393.2|CFAD_ECOLX|Full=CFA/I fimbrial subunit D|  
  "CFA/I fimbrial subunit D", # P25393.2|CFAD_ECOLX|Full=CFA/I fimbrial subunit D| 
  "CFA/I fimbrial subunit D", # P25393.2|CFAD_ECOLX|Full=CFA/I fimbrial subunit D|
  "Cell surface antigen I/II", # P21979.2|SPAA_STRDO|Full=Cell surface antigen I/II;
  "Cell surface antigen I/II", # P21979.2|SPAA_STRDO|Full=Cell surface antigen I/II;
  "Cell surface antigen I/II", # P21979.2|SPAA_STRDO|Full=Cell surface antigen I/II;
  "Cell surface antigen I/II", # P23504.2|SPAP_STRMU|Full=Cell surface antigen I/II;
  "Cell surface antigen I/II", # P23504.2|SPAP_STRMU|Full=Cell surface antigen I/II;
  "Cell surface antigen I/II", # P23504.2|SPAP_STRMU|Full=Cell surface antigen I/II;
  "Serine-rich repeat protein PsrP", # A0A0H2URK1.1|PSRP_STRPN|Full=Pneumococcal serine-ri
  "Serine-rich repeat protein PsrP", # A0A0H2URK1.1|PSRP_STRPN|Full=Pneumococcal serine-ri
  "Serine-rich repeat protein PsrP", # A0A0H2URK1.1|PSRP_STRPN|Full=Pneumococcal serine-ri
  "Serine-aspartate repeat-containing protein D", # O86488.1|SDRD_STAAE|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # O86488.1|SDRD_STAAE|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # O86488.1|SDRD_STAAE|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # Q5HIB3.1|SDRD_STAAC|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # Q2FJ78.1|SDRD_STAA3|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # Q5HIB3.1|SDRD_STAAC|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # Q2FJ78.1|SDRD_STAA3|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # Q5HIB3.1|SDRD_STAAC|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # Q2FJ78.1|SDRD_STAA3|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # Q2G0L4.1|SDRD_STAA8|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # Q2G0L4.1|SDRD_STAA8|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # Q2G0L4.1|SDRD_STAA8|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # Q6GBS5.1|SDRD_STAAS|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # Q6GBS5.1|SDRD_STAAS|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # Q6GBS5.1|SDRD_STAAS|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # Q8NXX6.1|SDRD_STAAW|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # Q8NXX6.1|SDRD_STAAW|Full=Serine-aspartate repeat-co
  "Serine-aspartate repeat-containing protein D", # Q8NXX6.1|SDRD_STAAW|Full=Serine-aspartate repeat-co
  "Fibrinogen-binding protein", # O70022.1|FBE_STAEP|Full=Fibrinogen-binding protein
  "Fibrinogen-binding protein", # O70022.1|FBE_STAEP|Full=Fibrinogen-binding protein
  "Fibrinogen-binding protein") # O70022.1|FBE_STAEP|Full=Fibrinogen-binding protein)

Str_adhesion_2_SA <- subset(Str_adhesion_2, Species == "S. australis")
Str_adhesion_3_SA <- ddply(Str_adhesion_2_SA, "Sequences3", numcolwise(median)) 

Str_adhesion_2_SP <- subset(Str_adhesion_2, Species == "S. parasanguinis")
Str_adhesion_3_SP <- ddply(Str_adhesion_2_SP, "Sequences3", numcolwise(median)) 

Str_adhesion_3 <- data.frame(rbind(Str_adhesion_2_SA, Str_adhesion_2_SP))

SA_vs_SP <-
  ggplot(Str_adhesion_3, aes(x=Species, y=Sequences3, fill=Score)) + geom_tile() +
  theme_pubr(border=TRUE, legend = "right", base_size = 10) + 
  scale_fill_gradientn(colours = c("gray68", "gray58", "gray48", "gray38",
                                   "gray28", "gray18")) + xlab("") + ylab("") +
  theme(axis.text.x = element_text(face="italic"))

# Export tables
# write.csv2(core_rare_metabolic_differences, "metabolic_differences.csv", row.names = FALSE)
# write.csv2(varContribution_150_merge, "metabolic_diff_genus.csv", row.names = FALSE)

# Export plots
ggsave("image_output/Supplementary_Fig_6.tif", plot_heat_all, device = "tiff", dpi= 600,
       width = 8.8, height = 12)

ggsave("image_output/Supplementary_Fig_7.tif", pca_plots, device = "tiff", dpi = 300, 
       width = 12, height = 8.36)

ggsave("image_output/Fig_3.tif", all_hier, device = "tiff", dpi = 300,
       width = 8.46, height = 8.52)

ggsave("image_output/Supplementary_Fig_2.tif", SA_vs_SP, device = "tiff", dpi = 300,
       width = 7.51, height = 4.22)