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
