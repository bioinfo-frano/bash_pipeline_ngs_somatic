**WORK IN PROGRESS**

# Introduction

Many students and scientists who want to learn genomic data analysis independently often rely on personal computers (e.g., laptops) with limited computational resources. These constraints can make it challenging to run secondary and tertiary NGS analyses starting from FASTQ files through variant annotation.
One alternative is to use cloud-based resources (e.g., AWS EC2 instances with large amounts of RAM), where software dependencies can be installed via Conda. However, cloud services incur financial costs. 
Therefore, for learners who prefer to work locally — typically with 8–16 GB RAM and limited disk storage (<60 GB) — it becomes essential to carefully select small sequencing datasets and design lightweight Conda environments. Equally important is the creation of a simple and efficient analysis pipeline and a well-defined folder structure.

**Creating the computing environment for NGS - DNA analysis**

I.	Create conda environment

II.	Find & download small size FASTQ files

III.	Find & download an indexed reference human genome

IV.	Create folder structure (example)

V.	Bash Scripting


**Purpose**

This repository provides a minimal, reproducible guide for running a DNA-seq (NGS) analysis pipeline locally, from FASTQ files to variant calling and annotation, using modest computational resources.

**Disclaimer**

This guide was developed and tested on macOS running on Intel processors. Users on Apple Silicon (M1/M2/…/M5) or Linux systems should adapt the present NGS workflow accordingly.

> **IMPORTANT: Conda prerequisites**
I assume that you have already installed miniconda3 in your computer. If not, please find the documentation on how to install miniconda in the link below.

<https://docs.conda.io/projects/conda/en/stable/user-guide/install/macos.html>

<https://www.youtube.com/watch?v=OH0E7FIHyQo>

If Miniconda is already installed, you should see the (base) environment activated in your Terminal.


I. Create a specific conda environment called "DNA"

a) Create environment with dependencies

** IMPORTANT **
You don't need to deactivate from 'base' when creating a new env.

`conda create -n DNA \
  -c conda-forge -c bioconda -c defaults \
  python=3.9 \
  openjdk=17 \
  perl=5.32 \
  fastqc \
  multiqc \
  cutadapt \
  bwa \
  samtools \
  picard \
  htslib \
  gatk4 \
  bcftools \
  vcftools \
  snpeff \
  ensembl-vep \
  bedtools \
  coreutils \
  pigz \
  pbzip2 \
  pandas \
  numpy \
  matplotlib \
  seaborn \
  -y`


  b) Sanity checks - post installation

  1. Check the list of env. The new 'DNA' env should appear:

  conda env list

  2. Activate the new env.

  conda activate DNA

  3. Check the dependencies in 'DNA'

  conda list

  gatk --version
  java -version
  bwa
  samtools --version
  vep --help
  snpEff -version


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
