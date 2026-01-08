**WORK IN PROGRESS**

# Introduction

Many students and scientists who want to learn genomic data analysis independently often rely on personal computers (e.g., laptops) with limited computational resources. These constraints can make it challenging to run secondary and tertiary NGS analyses starting from FASTQ files through variant annotation.
One alternative is to use cloud-based resources (e.g., AWS EC2 instances with large amounts of RAM), where software dependencies can be installed via Conda. However, cloud services incur financial costs. 
Therefore, for learners who prefer to work locally — typically with 8–16 GB RAM and limited disk storage (<60 GB) — it becomes essential to carefully select small sequencing datasets and design lightweight Conda environments. Equally important is the creation of a simple and efficient analysis pipeline and a well-defined folder structure.

In practice, FASTQ file sizes vary widely depending on the sequencing strategy and target region. Small targeted sequencing panels may generate FASTQ files smaller than 1 GB, whereas larger experiments such as whole-exome sequencing (WES) or whole-genome sequencing (WGS) can easily produce tens to hundreds of gigabytes per sample. For users working on local machines with limited RAM and disk space, this variability has a direct impact on dataset selection and pipeline feasibility. Table 1 provides representative examples of FASTQ file sizes from cancer-related datasets deposited in the NCBI Sequence Read Archive (SRA).

**Table 1. Representative FASTQ file sizes from cancer-related SRA datasets**

**Table 1. Representative FASTQ file sizes from cancer-related SRA datasets**

| **Size Category** | **Estimated FASTQ Download Size (GB)** | **Sequencing Strategy** | **Example SRA Run (approx.)** | **Cancer Type / Comments** |
|------------------|----------------------------------------|-------------------------|-------------------------------|---------------------------|
| Small            | ~0.05 – 0.2                            | Targeted gene panel     | [SRX11805868](https://www.ncbi.nlm.nih.gov/sra/SRX11805868) (~0.2 GB) | Panel targeting ~95 cancer genes |
| Small            | ~0.05 – 0.1                            | Targeted gene panel     | [SRX18078826](https://www.ncbi.nlm.nih.gov/sra/SRX18078826) (~0.05 GB) | Adult breast cancer |
| Small            | ~0.06 – 0.07                           | Targeted gene panel     | [SRX18078666](https://www.ncbi.nlm.nih.gov/sra/SRX18078666) (~0.065 GB) | Adult breast cancer |
| Medium           | ~3 – 7                                 | Targeted panel / DNA-seq| [SRX28185140](https://www.ncbi.nlm.nih.gov/sra/SRX28185140) (~6 GB) | SCLC tumor sample |
| Medium / Large   | ~3 – 12                                | WES                     | [SRX29598354](https://www.ncbi.nlm.nih.gov/sra/SRX29598354) (~3.4 GB raw) | TNBC tumor |
| Large            | >50                                    | WGS                     | [SRX160842](https://www.ncbi.nlm.nih.gov/sra/SRX160842) (>>50 GB) | Acute Lymphoblastic Leukemia (ALL) Phase II |

For this reason, the examples and pipelines in this repository focus on small targeted sequencing datasets that can be processed comfortably on standard workstations or laptops.


## Creating the computing environment for NGS - DNA analysis

I.	Create conda environment

II.	Create folder structure (example)

III.	Find & download small size FASTQ files

IV.	Find & download an indexed reference human genome

V.	Bash shell Scripting


## Purpose

This repository provides a minimal, reproducible guide for running a DNA-seq (NGS) analysis pipeline locally, from FASTQ files to variant calling and annotation, using modest computational resources.

> **Note:**  
> This guide was developed and tested on macOS running on Intel processors. Users on Apple Silicon (M1/M2/…/M5) or Linux systems may need to adapt the workflow accordingly.


> **IMPORTANT – Conda prerequisites:**  
> I assume that you have already installed Miniconda3 on your computer. If not, please find the [documentation](https://docs.conda.io/projects/conda/en/stable/user-guide/install/macos.html) on how to install Miniconda and also check out this [YouTube](https://www.youtube.com/watch?v=OH0E7FIHyQo) video.


When Miniconda is already installed, you should see the `(base)` environment activated in your Terminal.


## I. Create a specific conda environment called "DNA"

### a) Create the environment with these dependencies.

> **IMPORTANT**
You don't need to deactivate from `(base)` when creating a new environment.

```bash
conda create -n DNA \
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
  -y
```

### b) Sanity checks - post installation

  1. Check the list of env. The new 'DNA' env should appear:

  ```bash
  conda env list
  ```

  2. Activate the new env.

```bash
  conda activate DNA
```

  3. Check the dependencies in 'DNA' and the version of each of these dependencies:

```bash
  conda list
```

  gatk --version
  
  java -version
  
  bwa
  
  samtools --version
  
  vep --help
  
  snpEff -version
  

It could happens that `snpEff -version` shows and error like this: 

```bash
Error: LinkageError
UnsupportedClassVersionError
class file version 65.0
Java runtime only recognizes up to 61.0
```
Meaning: 
snpEff 5.3 was compiled with Java 21
You are running Java 17 (`openjdk=17`)
Java 17 cannot run Java 21 bytecode

Recommendation: Downgrade snpEff (this will guarantee stability)

Then, in bash:

`conda install -n DNA -c bioconda snpeff=5.1`

Check again:

`snpEff -version`

Expected output:

`SnpEff	5.1d	2022-04-19`


## II. Create folder structure

All FASTQ files and others like the reference (e.g. human) genome and scripts should be located in specific folders. Here there's a proposed example of a folder structure:

```bash
Genomics_simple_2026/
├── reference/                 # Reference genomes and known sites
│   └── GRCh38/
│       ├── fasta/             # Reference FASTA files
│       └── known_sites/       # Known variant sites (e.g. dbSNP, Mills)
├── data/
│   └── SRA_ID/                # Sample-specific directory (e.g. SRX11805868)
│       ├── raw_fastq/         # Original FASTQ files
│       ├── qc/                # FastQC / MultiQC reports
│       ├── trimmed/           # Adapter- and quality-trimmed FASTQ files
│       ├── aligned/           # BAM files and indexes
│       ├── variants/          # VCF files
│       └── annotation/        # Annotated variants
├── scripts/                   # scripts
└── logs/                      # Log files from pipeline execution
```

Multiple samples can be processed by creating one directory per SRA accession under `data/`.

In bash:

```bash
mkdir -p Genomics_simple_2026/{reference/GRCh38/{fasta,known_sites},data/SRA_ID/{raw_fastq,qc,trimmed,aligned,variants,annotation},scripts,logs}
```

## III. Find & download small size FASTQ files of gene panels for cancer diagnostics





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


> **Note:**  
> Table values are approximate and intended for educational purposes. Dataset sizes may vary depending on sequencing depth, platform, and compression.

