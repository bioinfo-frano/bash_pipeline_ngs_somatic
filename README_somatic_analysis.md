**WORK IN PROGRESS**

#  Part II – Somatic analysis (Bash pipelines)

## Introduction

Building on [Part I – Preparation & setup](README_setup.md), this tutorial demonstrates how to analyze DNA NGS data on a **local workstation or laptop**.
The workflow is designed to run under **limited computational resources**, and has been tested on macOS (Intel-based MacBook Pro). With minor adaptations, the same approach can be applied to Linux systems.

In **Part I**, steps I–IV were completed.

### Completed steps

I.	Create conda environment

II.	Create a folder structure

III.	Find and download small-sized FASTQ datasets

IV. Download a reference human genome and indexes

V.	Bash shell scripting (covered in this part)

---

## Types of DNA-NGS analysis

When investigating genomic alterations, there are two main types of DNA-NGS analysis:

### 1. Somatic analysis 

Somatic analysis focuses on mutations acquired by a tumor. These variants are **not inherited** and are typically absent from normal tissue. 

Typical applications include:

   - Identifying actionable mutations for targeted therapies (e.g. *PIK3CA, ERBB2, ESR1, AKT1* in breast cancer)

   - Monitoring disease progression and minimal residual disease (MRD), often via liquid biopsy and circulating tumor DNA (ctDNA)

   - Clinical trials: Matching the patient to mutation-specific trials

   - Understanding tumor biology: Subtyping the tumor more precisely, identifying driver vs passenger mutations, and detecting resistance mutations
    
---

### 2. Germline analysis
 
Germline analysis aims to detect **inherited variants** present in all cells of the body. These variants may predispose individuals to cancer or other genetic diseases.

Common use cases include:

  - Assessing hereditary cancer risk: early-onset cancer, strong family history, bilateral breast cancer, triple-negative breast cancer < 60 years

  - Clinical decision-making when germline status influences treatment
    
  - Genetic counseling and family screening

**Sample types:** blood or saliva (non-tumor tissue)

---
    

### Somatic vs germline (summary)

- **Somatic NGS analysis:** understanding and treating the tumor  
- **Germline NGS analysis:** assessing inherited cancer risk


In many clinical settings, somatic (tumor) and germline (normal) DNA samples can be sequenced simultaneously, allowing both somatic and germline analyses from the same patient. This matched approach provides the most complete, accurate, and clinically actionable genomic portrait for a cancer patient and it is considered the best practice in modern genomic and precision oncology, **although** cost and access may limit its use.

## Bioinformatic DNA-NGS workflows

DNA-NGS analysis can be automated using several workflow approaches, including:
- Bash scripting
- Nextflow
- Snakemake
- Galaxy
- Others (Cromwell/WDL, Toil, Airflow, Argo Workflows, Make, Ruffus, Bpipe)


**Table 1**: Advantages and disadvantages of selected workflow approaches.

| Tool | Advantages | Disadvantages |
|------|------------|---------------|
| **Bash scripting** | • Maximum transparency and control<br>• No new syntax to learn<br>• Simple to start<br>• Universally available<br>• Ideal for learning and prototyping | • Limited scalability<br>• Manual error handling<br>• No native workflow management |
| **Snakemake** | • Python-based and readable<br>• Explicit file dependencies<br>• Excellent for complex workflows | • Scaling requires configuration<br>• Less cloud-native by default |
| **Nextflow** | • Highly portable (“write once, run anywhere”)<br>• Built-in parallelization and resume<br>• Strong community (nf-core) | • Learning curve (Groovy DSL)<br>• Debugging channels can be challenging |



In this **Part II**, bash scripting will be used to perform a **somatic DNA-NGS analysis** on the small FASTQ dataset introduced in [Part I – Preparation & setup](README_setup.md).


## Bioinformatics overview of step-by-step somatic DNA-NGS pipeline

### Preprocessing

```bash
Quality control (QC) FastQC + MultiQC
 ↓
Trimming - Adapter/UMI trimming             # FastQC + MultiQC
 ↓
QC                                          # FastQC + MultiQC
 ↓
BWA-MEM (with Read Groups) → SAM		        # Alignment step (🔴 GATK requires read groups)
 ↓
Convert SAM → BAM
 ↓
Sort BAM (samtools / Picard)                # 🔴 CRITICAL (Sorting before MarkDuplicates. Picard requires coordinate-sorted BAMs.)
 ↓
MarkDuplicates (Picard) (or UMI collapsing)
 ↓
Index BAM                                   # 🔴 REQUIRED before GATK
```

>**Note**: Base Quality Score Recalibration (BQSR) is often omitted for small targeted panels or UMI-based datasets and is therefore not included in this tutorial.

### Somatic variant analysis pipeline 

```bash
[Preprocessed BAM]
 ↓
Mutect2 (somatic variant caller)
 ↓
FilterMutectCalls
 ↓
Variant Annotation (VEP, ClinVar, COSMIC, SnpEff)
```
---

## V. Bash scripting

According to [freeCodeCamp](https://www.freecodecamp.org/news/bash-scripting-tutorial-linux-shell-script-and-command-line-for-beginners/), a bash script "is a file containing a sequence of commands that are executed by the bash program line by line. It allows you to perform a series of actions, such as navigating to a specific directory, creating a folder, and launching a process using the command line. By saving these commands in a script, you can repeat the same sequence of steps multiple times and execute them by running the script". You can watch this tutorial on [YouTube](https://www.youtube.com/watch?v=tK9Oc6AEnR4&t=614s), also from **freeCodeCamp**, which explains step-wise how to create and apply bash scripting. 

Each step in the somatic DNA-NGS analysis will have its own script. It's not necessary that you figure out how to write the scripts; they will be linked so that you can download and run them on your computer.

### Quality control (QC) 👉 [01_qc.sh](01_qc.sh)

#### FastQC & MultiQC

To get the documentation on how to interpret FastQC reports, check the documentation from the [Babraham Bioinformatics](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- FASTQ file path: ~/Genomics_cancer/data/SRA_ID/raw_fastq
- QC reports path: ~/Genomics_cancer/data/SRA_ID/qc
- REPORTS: `SRR30536566_1_fastqc.html` (R1), `SRR30536566_2_fastqc.html` (R2), `multiqc_report.html`
- Observations: 
      - High content of duplicated reads (~60% and ~66% in R1 and R2, respectively) and high GC content. 
      - Overall base quality is good across reads.
      - The first and last ~5 bp show slightly reduced quality and will be trimmed (see **Figure 1**).
      
      **Figure 1.** Per-base sequence quality plots for R1 and R2 show high-quality base calls across most read positions.
      ![Figure 1: Per-base sequence quality](images/FastQ_Per_base_sequence_quality)
  
  
      
- Regarding the **high content of duplicated reads**: The NCBI-SRA information of the dataset (RUN: SRR30536566) states: `Strategy: AMPLICON` and `Selection: PCR`. This is a **targeted amplicon sequencing**. The full-length of KRAS NRAS BRAF PIK3CA PTEN RRAS and MEK1 including UTR, exons, and introns were sequenced. These specific genomic regions were PCR-amplified many times to create enough material for sequencing. All reads derived from the same original fragment are **technical duplicates** (PCR duplicates). Therefore, a very high duplication rate is inherent to the technique. It does not reflect poor quality; it reflects the method.
      - Regarding the **high GC content**: Usually, a bimodal curve often suggests contamination (e.g., bacteria in a human sample) or a mixed sample. However, this is explained by the experimental design.
      You are not sequencing the whole human genome (which has a relatively uniform ~41% GC). You are sequencing a panel of specific amplicons. Different genes have different base compositions. The "camel" shape strongly suggests your targeted panel contains two distinct classes of amplicons:
        1. First Hump (Peak ~35% GC): Likely corresponds to a subset of your amplicons that are GC-poor.
        2. Second Hump (Peak ~62% GC): Likely corresponds to a subset of amplicons that are GC-rich (common in many coding regions).
      The red warning is because FastQC compares your distribution to a unimodal model based on a normal genome. Your amplicon-based distribution violates this model, so it gets flagged. This is not a problem for your data.





With these documentation and the information in **Table 1** 