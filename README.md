## Bacterial low-abundant taxa are key determinants of a healthy airway metagenome in the early years of human life
<br>
Marie-Madlen Pust<sup>1,2</sup>, Burkhard Tümmler<sup>1,2</sup> <br>
<sup>1</sup>Department of Paediatric Pneumology, Allergology, and Neonatology, Hannover Medical School (MHH), Germany <br>
<sup>2</sup>Biomedical Research in Endstage and Obstructive Lung Disease Hannover (BREATH), German Center for Lung Research, Hannover Medical School, Germany <br><br><br/>


**Reference databases** <br/>
The one-strain per species multi-FASTA file

```bash

# Download database
wget https://sync.academiccloud.de/index.php/s/h1it8NhwGSMaKGe/download

# Unzip the reference database 
gunzip complete_bacterialRefSeqs_201910_3.fasta.gz

# Generate an index of the reference fasta depending on the alignment tool of your choice
samtools complete_bacterialRefSeqs_201910_3.fasta
bwa index complete_bacterialRefSeqs_201910_3.fasta
```

<br/>
The pangenome multi-FASTA file 

```bash
# Download database
wget https://sync.academiccloud.de/index.php/s/vOTDJ9qDR6tvn0w/download

# Unzip the reference database 
gunzip 2020_09_reference.fa.gz

# Generate an index of the reference fasta depending on the alignment tool of your choice
samtools faidx 2020_09_reference.fa
bwa index 2020_09_reference.fa
```
<br/>


The adhesin protein sequence multi-FASTA file

```bash
# Download database
wget https://sync.academiccloud.de/index.php/s/YaL6NXMEavuZWkd/download
```
<br/>

**R files with in-text comments** <br/>
1. step_1_bootstrapping_aggregations.R
2. step_2_simulation_runs.R
3. step_3_functional_analysis.R
<br/>

**R session, information** <br/>
```bash
R version 4.1.1 (2021-08-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19042)

Matrix products: default

locale:
LC_COLLATE=German_Germany.1252
LC_CTYPE=German_Germany.1252   
LC_MONETARY=German_Germany.1252
LC_NUMERIC=C
LC_TIME=German_Germany.1252   

attached base packages:
grid, stats, graphics, grDevices, utils, datasets, methods, base   

other attached packages:
[01] Boruta_7.0.0           ggVennDiagram_1.1.4    randomForest_4.6-14    RVAideMemoire_0.9-80   microbiome_1.14.0      phyloseq_1.36.0        data.table_1.14.0      graphkernels_1.6       CINNA_1.1.54          
[10] NetSwan_0.1            igraph_1.2.6           expm_0.999-6           Matrix_1.3-4           tidygraph_1.2.0        sna_2.6                network_1.17.1         statnet.common_4.5.0   qgraph_1.6.9          
[19] popgraph_1.5.2         cowplot_1.1.1          forcats_0.5.1          ggrepel_0.9.1          knitr_1.33             SimilarityMeasures_1.4 rcompanion_2.4.1       reshape2_1.4.4         scales_1.1.1          
[28] magrittr_2.0.1         reshape_0.8.8          Hmisc_4.5-0            Formula_1.2-4          survival_3.2-11        gmodels_2.18.1         rstatix_0.7.0          viridis_0.6.1          viridisLite_0.4.0     
[37] tidyr_1.1.3            ggdendro_0.1.22        RColorBrewer_1.1-2     pheatmap_1.0.12        stringr_1.4.0          factoextra_1.0.7       ggpubr_0.4.0           vegan_2.5-7            lattice_0.20-44       
[46] permute_0.9-5          ggplot2_3.3.5          dplyr_1.0.7            purrr_0.3.4            plyr_1.8.6             readr_2.0.1           

loaded via a namespace (and not attached):
  [01] GGally_2.1.2           qdapTools_1.3.5        lavaan_0.6-9           coda_0.19-4            ragg_1.1.3             bit64_4.0.5            multcomp_1.4-17        rpart_4.1-15           RCurl_1.98-1.4        
 [10] generics_0.1.0         BiocGenerics_0.38.0    TH.data_1.0-10         proxy_0.4-26           chron_2.3-56           bit_4.0.4              tzdb_0.1.2             assertthat_0.2.1       xfun_0.25             
 [19] hms_1.1.0              fansi_0.5.0            dendextend_1.15.1      readxl_1.3.1           DBI_1.1.1              tmvnsim_1.0-2          htmlwidgets_1.5.3      stats4_4.1.1           ellipsis_0.3.2        
 [28] corrplot_0.90          backports_1.2.1        FactoMineR_2.4         pbivnorm_0.6.0         libcoin_1.0-8          vctrs_0.3.8            Biobase_2.52.0         abind_1.4-5            withr_2.4.2           
 [37] RVenn_1.1.0            checkmate_2.0.0        vroom_1.5.4            fdrtool_1.2.16         mnormt_2.0.2           cluster_2.1.2          ape_5.5                crayon_1.4.1           pkgconfig_2.0.3       
 [46] labeling_0.4.2         units_0.7-2            GenomeInfoDb_1.28.1    nlme_3.1-152           nnet_7.3-16            rlang_0.4.11           lifecycle_1.0.0        sandwich_3.0-1         cellranger_1.1.0      
 [55] matrixStats_0.60.0     lmtest_0.9-38          carData_3.0-4          Rhdf5lib_1.14.2        boot_1.3-28            zoo_1.8-9              base64enc_0.1-3        png_0.1-7              rootSolve_1.8.2.2     
 [64] bitops_1.0-7           KernSmooth_2.23-20     rhdf5filters_1.4.0     Biostrings_2.60.2      classInt_0.4-3         multcompView_0.1-8     coin_1.4-1             jpeg_0.1-9             S4Vectors_0.30.0      
 [73] ggsignif_0.6.2         leaps_3.1              lpSolve_5.6.15         gdata_2.18.0           zlibbioc_1.38.0        compiler_4.1.1         intergraph_2.0-2       cli_3.0.1              ade4_1.7-17           
 [82] XVector_0.32.0         pbapply_1.4-3          htmlTable_2.2.1        MASS_7.3-54            mgcv_1.8-36            tidyselect_1.1.1       stringi_1.7.3          textshaping_0.3.5      latticeExtra_0.6-29   
 [91] tools_4.1.1            lmom_2.8               parallel_4.1.1         rio_0.5.27             rstudioapi_0.13        foreach_1.5.1          foreign_0.8-81         gridExtra_2.3          gld_2.6.2             
[100] scatterplot3d_0.3-41   farver_2.1.0           Rtsne_0.15             digest_0.6.27          nortest_1.0-4          Rcpp_1.0.7             car_3.0-11             broom_0.7.9            sf_1.0-2              
[109] psych_2.1.6            colorspace_2.0-2       ranger_0.13.1          IRanges_2.26.0         splines_4.1.1          centiserve_1.0.0       sp_1.4-5               multtest_2.48.0        Exact_2.1             
[118] systemfonts_1.0.2      jsonlite_1.7.2         corpcor_1.6.9          glasso_1.11            flashClust_1.01-2      modeltools_0.2-23      R6_2.5.1               pillar_1.6.2           htmltools_0.5.1.1     
[127] glue_1.4.2             DT_0.18                class_7.3-19           codetools_0.2-18       mvtnorm_1.1-2          utf8_1.2.2             tibble_3.1.3           curl_4.3.2             DescTools_0.99.42     
[136] gtools_3.9.2           zip_2.2.0              openxlsx_4.2.4         sampling_2.9           biomformat_1.20.0      munsell_0.5.0          e1071_1.7-8            rhdf5_2.36.0           GenomeInfoDbData_1.2.6
[145] iterators_1.0.13       haven_2.4.3            gtable_0.3.0    
```
