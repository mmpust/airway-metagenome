## Bacterial low-abundant taxa are key determinants of a healthy airway metagenome in the early years of human life
<br>
Marie-Madlen Pust<sup>1,2</sup>, Burkhard Tümmler<sup>1,2</sup> <br>
<sup>1</sup>Department of Paediatric Pneumology, Allergology, and Neonatology, Hannover Medical School (MHH), Germany <br>
<sup>2</sup>Biomedical Research in Endstage and Obstructive Lung Disease Hannover (BREATH), German Center for Lung Research, Hannover Medical School, Germany <br><br>


**Reference databases** <br/>
The one-strain per species multi-FASTA file can be obtained (see instructions below).

```bash
https://sync.academiccloud.de/index.php/s/wSpNHavgEX4CTeX/download

# Unzip the reference database 
gunzip 2020_09_reference.fa.gz

# Generate an index of the reference fasta depending on the alignment tool of your choice
samtools faidx 2020_09_reference.fa
bwa index 2020_09_reference.fa
```

**Reference databases** <br/>
The pangenome multi-FASTA file can be obtained (see instructions below).

```bash
https://sync.academiccloud.de/index.php/s/wSpNHavgEX4CTeX/download

# Unzip the reference database 
gunzip 2020_09_reference.fa.gz

# Generate an index of the reference fasta depending on the alignment tool of your choice
samtools faidx 2020_09_reference.fa
bwa index 2020_09_reference.fa
```
