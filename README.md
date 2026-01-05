**WORK IN PROGRESS**

**Introduction**

Many students and scientist who would like to start learning genomic analysis independently have computers, e.g. laptops, with limited computational power to run primary and secondary level of bioinformatic analysis. One excellent alternative is to
set up a cloud instance, e.g AWS EC2, with 64GB RAM and create and install conda environments and dependencies, e.g. fastQC, etc., to run analyses of fastq files. However, using cloud services imply financial costs. Then, for those who would like to stick to the idea of learning genomic analysis in their local computer, with limited 8 to 16 GB RAM and less than 60GB memory available, would be important to think what kind of sequencing (fastq) files to download, how to create a conda environment, how to set out a folder structure and pipelines for running the analyses.

**Aims**

I.	Create conda environment

II.	Find & download small size FASTQ files

III.	Find & download an indexed reference human genome

IV.	Create folder structure (example)

V.	Bash Scripting




I would suggest the following: start with small fastq files.
While small fastq files can be < 1GB, which contain sequences of few genes, big fastq files can be between >10 to > 100 GB (See table)




# Reference Genome

**Genome:** GRCh38 / hg38  
**Source:** Broad Institute – GATK Resource Bundle  
**Bucket:** gcp-public-data--broad-references/hg38/v0

## Files and checksums

Homo_sapiens_assembly38.fasta  
MD5: 7ff134953dcca8c8997453bbb80b6b5e

Homo_sapiens_assembly38.fasta.fai  
MD5: f76371b113734a56cde236bc0372de0a

Homo_sapiens_assembly38.dict  
MD5: 3884c62eb0e53fa92459ed9bff133ae6


Variant annotation was performed using Ensembl Variant Effect Predictor (VEP, GRCh38). Due to hardware and storage limitations, annotation was executed via the Ensembl web interface, enabling gnomAD allele frequencies, MANE/APPRIS transcripts, phenotype annotations, and pathogenicity predictors (REVEL, ClinPred). Command-line execution was tested but limited by local resource constraints.
