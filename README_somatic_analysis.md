**WORK IN PROGRESS**

#  Part II – Somatic analysis (Bash pipelines)

## Introduction

Building on [Part I – Preparation & setup](README_setup.md), this tutorial demonstrates how to analyze DNA NGS data on a **local workstation or laptop**.
The workflow is designed to run under **limited computational resources**, and has been tested on macOS (Intel-based MacBook Pro). With minor adaptations, the same approach can be applied to Linux systems.

In **Part I**, steps I–IV were completed.

### Completed steps

I.	Create conda environment

II.	Create a folder structure

III.	Find & download small-sized FASTQ datasets

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

   - Understanding tumor biology: Subtyping the tumor more precisely, identifying driver vs passenger mutations and detecting resistance mutations
    
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
FastQC + MultiQC (QC)
 ↓
Trimming - Adapter/UMI trimming
 ↓
FastQC + MultiQC (QC)
 ↓
BWA-MEM (with Read Groups) → SAM		          ← 🔴 CRITICAL (GATK requires read groups.). Alignment step
 ↓
Convert SAM → BAM
 ↓
Sort BAM (samtools / Picard)			            ← 🔴 CRITICAL (Sorting before MarkDuplicates. Picard requires coordinate-sorted BAMs.)
 ↓
MarkDuplicates (Picard) (or UMI collapsing)
 ↓
Index BAM					                            ← 🔴 REQUIRED before GATK
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

## V. Bash scripting

