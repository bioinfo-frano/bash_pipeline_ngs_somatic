⚠️ **WORK IN PROGRESS**

#  Part II – Somatic analysis (Bash pipelines)

## Introduction

Building on [Part I – Preparation & setup](README_setup.md), this tutorial demonstrates how to analyze DNA NGS data on a **local workstation or laptop**.
The workflow is designed to run under **limited computational resources**, and has been tested on macOS (Intel-based MacBook Pro). With minor adaptations, the same approach can be applied to Linux systems.

In **Part I**, steps I–IV were completed. In this part, we focus on **bash scripting** to implement a complete somatic DNA-NGS analysis pipeline.

### Completed steps

I.	Create conda environment

II.	Create a folder structure

III.	Find and download small-sized FASTQ datasets

IV. Download a reference human genome and indexes

V.	Bash shell scripting (**covered in this part**)

---

## Types of DNA-NGS analysis

When investigating genomic alterations, two main categories of DNA-NGS analysis are commonly performed:

### 1. Somatic analysis 

Somatic analysis focuses on mutations acquired by tumor cells. These variants are **not inherited** and are typically absent from normal tissue. 

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


In many clinical settings, tumor (somatic) and matched normal (germline) samples are sequenced together. This matched analysis provides the most complete, accurate, and clinically actionable genomic portrait for a cancer patient and it is considered the best practice in modern genomic and precision oncology, **although** cost, sample availability, and access may limit its use.

---

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



In this **Part II**, bash scripting is used to perform a **somatic DNA-NGS analysis** on the small FASTQ dataset introduced in [Part I – Preparation & setup](README_setup.md).


## Bioinformatics overview: step-by-step somatic DNA-NGS pipeline

### Preprocessing

```bash
Raw FASTQ
 ↓
Quality control (FastQC + MultiQC)
 ↓
Trimming (adapters, low-quality bases, poly-A tails)      # FastQC + MultiQC
 ↓
Quality control (FastQC + MultiQC)                        # FastQC + MultiQC
 ↓
Alignment (BWA-MEM with Read Groups) → SAM                # 🔴 GATK requires read groups
 ↓
Conversion SAM → BAM
 ↓
Sort BAM (samtools)                                      # 🔴 REQUIRED before MarkDuplicates. Picard requires coordinate-sorted BAMs.)
 ↓
MarkDuplicates (Picard) or UMI collapsing
 ↓
MD/NM Tags (samtools)                                    # 🔴 REQUIRED before GATK
 ↓
Index BAM (samtools)
```

>**Note**: Base Quality Score Recalibration (BQSR) is often omitted for small targeted panels, UMI-based datasets or targeted amplicon sequencing and is therefore not included in this tutorial.
According to [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR) "The main case figure where you really might need to skip BQSR is when you have too little data (some small gene panels have that problem), or you're working with a really weird organism that displays insane amounts of variation."

### Somatic variant analysis pipeline 

```bash
Indexed BAM
 ↓
Mutect2 (somatic variant calling)
 ↓
FilterMutectCalls
 ↓
Variant Annotation (VEP, ClinVar, COSMIC, SnpEff)
```
---

## V. Bash scripting

According to [freeCodeCamp](https://www.freecodecamp.org/news/bash-scripting-tutorial-linux-shell-script-and-command-line-for-beginners/), a bash script "is a file containing a sequence of commands that are executed by the bash program line by line. It allows you to perform a series of actions, such as navigating to a specific directory, creating a folder, and launching a process using the command line. By saving these commands in a script, you can repeat the same sequence of steps multiple times and execute them by running the script". If you are new to bash scripting, the following [YouTube](https://www.youtube.com/watch?v=tK9Oc6AEnR4&t=614s) video tutorial, also from **freeCodeCamp**, provides a clear step-by-step introduction.

Each step in the somatic DNA-NGS analysis is implemented as a **separate bash script**. You do not need to write these scripts from scratch; they are provided and can be downloaded and executed directly on your computer.

### Quality control (QC) 👉 [01_qc.sh](bash_scripts/01_qc.sh)

#### FastQC & MultiQC

To learn how to interpret **FastQC** reports, check the official documentation from the [Babraham Bioinformatics](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- **Input FASTQ files**: `~/Genomics_cancer/data/SRR30536566/raw_fastq`
- **QC output path**: `~/Genomics_cancer/data/SRR30536566/qc`
- **QC reports**: 

  - `SRR30536566_1_fastqc.html` (R1)
      
  - `SRR30536566_2_fastqc.html` (R2)
      
  - `multiqc_report.html`
      
**Observations**: 

- A high content of duplicated reads (~60% and ~66% in R1 and R2, respectively) and high GC content.
      
- Overall high base quality across reads.
      
- Slightly reduced quality in the first and last ~5 bp, which will be trimmed (**Figure. 1**).
  
  
**Figure 1.** Per-base sequence quality plots for R1 and R2 show high-quality base calls across most read positions.
  
  
  ![Figure 1: Per-base sequence quality](images/FastQ_Per_base_sequence_quality.png)
  
  
      
- Regarding the **high content of duplicated reads** (see **Figure. 2**): The NCBI-SRA metadata for the dataset [RUN: SRR30536566](https://www.ncbi.nlm.nih.gov/sra/SRX25960056[accn]) reports: 
  
  - **Strategy**: `AMPLICON`
    
  - **Selection**: `PCR`

This means that the dataset was generated using **targeted amplicon sequencing** of full-length genes including KRAS NRAS BRAF PIK3CA PTEN RRAS and MEK1 including UTR, exons, and introns. These specific genomic regions were PCR-amplified many times to create enough material for sequencing. All reads derived from the same original fragment are **technical duplicates** (PCR duplicates), not biological duplicates. Therefore, a very high duplication rate is inherent to the technique (see **Figure.3**). It does not reflect poor quality; it reflects the method.

**Figure 2.** Sequence duplication levels of reads
  
  
  ![Figure 2: Sequence Duplication Levels](images/FastQ_Sequence_duplication_levels.png)
  

**Figure 3.** Proportion of Unique and duplicated reads 
  
  
  ![Figure 3: Sequence Duplication Levels](images/FastQ_Sequence_Counts_Duplicated_Reads.png)


- Regarding the **high GC content** (see **Figure 4**): Usually, a bimodal curve often suggests contamination (e.g., bacteria in a human sample) or a mixed sample that would raise concerns for whole-genome sequencing. However, this is expected for targeted amplicon panels. You are not sequencing the whole human genome (which has a relatively uniform ~41% GC). You are sequencing a panel of specific amplicons. Different genes have different base compositions, resulting in multiple GC peaks. The "camel" shape strongly suggests your targeted panel contains two distinct classes of amplicons:

  1. First Hump (Peak ~35% GC): Likely corresponds to a subset of your amplicons that are GC-poor.
  
  2. Second Hump (Peak ~62% GC): Likely corresponds to a subset of amplicons that are GC-rich (common in many coding regions).
  
The red warning is because FastQC compares your distribution to a unimodal model based on a normal genome. Your amplicon-based distribution violates this model, so it gets flagged. This is not a problem for your data.

**Figure 4.** Per sequence GC content showing a bimodal shape.
  
  
  ![Figure 4: Per sequence GC content](images/FastQ_Per_Sequence_GC_content.png)

---

### Trimming + QC 👉 [02_trim.sh](bash_scripts/02_trim.sh)

#### Cutadapt

Documentation for **Cutadapt** can be found in [Cutadapt](https://github.com/marcelm/cutadapt/blob/main/doc/guide.rst). Relevant sections include "***Trimming paired-end reads***" and "***Cutadapt's output***".

- **Trimmed FASTQ files output**: `~/Genomics_cancer/data/SRR30536566/trimmed`
- **Cutadapt log report**: `~/Genomics_cancer/logs/cutadapt_SRR30536566.log`
- **QC reports**: 
  -`SRR30536566_R1.trimmed_fastqc.html`
  -`SRR30536566_R2.trimmed_fastqc.html`
  -`multiqc_report_1.html`

**Trimming strategy**
  - Fixed trimming of the first and last 5 bp based on FastQC. 
  - Quality trimming of low-confidence bases (Phred < 20)
  - Removal of short reads (<30 bp).
  - Trimming of rare terminal poly-A stretches (≥10 A’s)
  
This trimming strategy is conservative and well-suited for targeted amplicon sequencing data, without affecting genuine sequence content.
  
**Cutadapt parameters**

```bash
-u 5 -u -5      → trim 5 bp from 5′ and 3′ of R1
-U 5 -U -5      → trim 5 bp from 5′ and 3′ of R2
-q 20,20        → trim low-quality bases (Phred < 20) from both ends
-m 30           → discard reads shorter than 30 bp
-a A{10}        → trim ≥10 A’s from 3′ end of R1
-A A{10}        → trim ≥10 A’s from 3′ end of R2
```

>**Note**: Fixed trimming (-u, -U) removes a constant number of bases and is therefore not reported separately in the Cutadapt log. Its effect is included in the difference between total processed and total written bases.

**Cutadapt summary interpretation** (`cutadapt_SRR30536566.log`):

After trimming and filtering:

- 97.2% of read pairs were retained

- 86.4% of total bases were preserved

- Quality trimming (`-q 20,20`) affected only 0.9% of bases

- Pairs that were too short (`-m 30`): 2.8%

- Poly-A contamination was minor (<5%)

Why is retained 86.4% of bases?

| Cause                         | Effect                         |
| ----------------------------- | ------------------------------ |
| Fixed trimming (`-u 5 -u -5`) | −10 bp per read                |
| Quality trimming (`-q 20,20`) | −0.9%                          |
| Poly-A trimming (`A{10}`)     | Variable tail removal          |
| Read removal (`-m 30`)        | 2.8% of reads dropped entirely |

✔ Excellent base quality

✔ Expected amplicon duplication

✔ Minor poly-A contamination

✔ Minimal read loss

This indicates **high-quality data** with minimal loss of informative reads.

---

### Alignment and BAM preprocessing workflow 👉 [03_align_&_bam_preprocessing.sh](bash_scripts/03_align_&_bam_preprocessing.sh)

**Documentation**: 

- **BWA**: <https://bio-bwa.sourceforge.net/bwa.shtml>
- **samtools** (v. 1.23): <https://www.htslib.org/doc/samtools.html>
- **Picard**: <https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard>

**Table 2**: Alignment and BAM preprocessing workflow (Script: 03_align_conv_sort_mdup_index.sh)

| Step | Step name                  | INPUT                                                                                                     | OUTPUT                                                                | Tool                    | Function / Role                                                                                                                                                                                                 |
| ---- | -------------------------- | --------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- | ----------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|  1 | Alignment                  | `SRR30536566_R1.trimmed.fastq.gz`<br>`SRR30536566_R2.trimmed.fastq.gz`<br>`Homo_sapiens_assembly38.fasta` | `SRR30536566.sam`                                                     | `bwa mem`               | Aligns paired-end reads to the reference genome. Adds **Read Group (RG)** information required by GATK and downstream tools. Output is an unsorted SAM alignment file.                                          |
|  2 | Conversion SAM → BAM       | `SRR30536566.sam`                                                                                         | `SRR30536566.bam`                                                     | `samtools view`         | Converts human-readable SAM into compressed binary BAM format for efficiency and downstream processing.                                                                                                         |
|  3 | Sort BAM (coordinate sort) | `SRR30536566.bam`                                                                                         | `SRR30536566.sorted.bam`                                              | `samtools sort`         | Sorts alignments by genomic coordinates (chromosome and position). **Required** for duplicate marking, indexing, and variant calling.                                                                           |
|  4 | MarkDuplicates    | `SRR30536566.sorted.bam`                                                                                  | `SRR30536566.sorted.markdup.bam`<br>`SRR30536566.markdup.metrics.txt` |`picard MarkDuplicates` | Identifies PCR/optical duplicates and **marks them in the BAM (FLAG + tags)** without removing reads. Duplicate sets are tagged (amplicon-aware), enabling variant callers to down-weight or ignore duplicates. |
|  5 | Add MD/NM tags             | `SRR30536566.sorted.markdup.bam`<br>`Homo_sapiens_assembly38.fasta`                                       | `SRR30536566.sorted.markdup.md.bam`                                   |`samtools calmd`        | Recalculates and adds **MD** (mismatch positions) and **NM** (Number of mismatches) tags. Improves robustness and compatibility with GATK and somatic variant callers.                                                 |
|  6 | Index final BAM            | `SRR30536566.sorted.markdup.md.bam`                                                                       | `SRR30536566.sorted.markdup.md.bam.bai`                               |`samtools index`        | Creates a BAM index enabling **random genomic access**. Required for variant calling (e.g. Mutect2), visualization (IGV), and QC tools.                                                                         |



**Table 3**: GATK-relevant read group fields

| Field | Meaning       | Recommendation      | Script                   | Interpretation
| ----- | ------------- | ------------------- |--------------------------|-----------------
| `ID`  | Read group ID | Unique per lane/run | `RG_ID="SRR30536566"`    | Patient ID
| `SM`  | Sample        | Biological sample   | `RG_SM="DMBEL-EIDR-071"` | Reflects the biological sample
| `LB`  | Library       | Library prep        | `RG_LB="AMPLICON"`       | Reflects library strategy
| `PL`  | Platform      | Sequencing platform | `RG_PL="ILLUMINA"`       | Mandatory
| `PU`  | Platform unit | Flowcell + lane     | `RG_PU="HiSeq4000"`      | Uniquely identifies the sequencing unit




### Samtools – Troubleshooting

When following this tutorial, it is very important to be aware of which version of samtools is installed in the conda environment.
In this course, the environment provides samtools v0.1.19 (released in 2013), which is much older than current samtools versions (≥1.x).

Because of this, some command-line options shown in modern tutorials and online documentation do NOT exist in this older version. Using newer options with an old version will lead to errors, even if the command looks perfectly correct.

You can always check the installed version with:

```bash
samtools
```

### 🔹 A common student concern: “Is the old version worse?”

During this tutorial, students often feel that using an old samtools version means lower-quality results.
**This is not true**.

✔️ The core algorithms used by samtools—such as:

- coordinate sorting of BAM files

- calculation of MD and NM tags

- interpretation of alignments

are conceptually identical across versions.

What has mainly improved in newer versions is:

- memory efficiency

- speed

- multi-threading support

- cleaner command-line syntax

The result of the analysis and biological meaning (sorted BAM order, MD/NM tag values, alignment interpretation) remain the same.


### 🔹 Key version-related differences in samtools versions

1. `samtools sort`

Old version (v0.1.19) does not support `-o`
```bash
THREADS=4
samtools sort -@ "$THREADS" \
  "$ALIGN_DIR/${SAMPLE}.bam" \
  "$ALIGN_DIR/${SAMPLE}.sorted"
```

Current version (≥1.x)
```bash
THREADS=4
samtools sort -@ "$THREADS" \
  -o "$ALIGN_DIR/${SAMPLE}.sorted.bam" \
  "$ALIGN_DIR/${SAMPLE}.bam"
```

Output:
```bash
SRR30536566.sorted.bam
```

2. `samtools calmd`

Old version (v0.1.19) does not support `-@`, `--threads`

```bash
FINAL_BAM="$ALIGN_DIR/${SAMPLE}.sorted.markdup.md.bam"
samtools calmd -b \
  "$ALIGN_DIR/${SAMPLE}.sorted.markdup.bam" \
  "$REF_FASTA" \
  > "$FINAL_BAM"
```

Current version (≥1.x)

```bash
FINAL_BAM="$ALIGN_DIR/${SAMPLE}.sorted.markdup.md.bam"
samtools calmd -@ "$THREADS" -b \
  "$ALIGN_DIR/${SAMPLE}.sorted.markdup.bam" \
  "$REF_FASTA" \
  > "$FINAL_BAM"
```

Output:
```bash
SRR30536566.sorted.markdup.md.bam.bai
```

### Folder structure: From QC → Trimming/Filtering → Alignment + BAM preprocessing.

```bash
Genomics_cancer/
├── reference/                 
│   └── GRCh38/
│       ├── fasta/
│       │   ├── Homo_sapiens_assembly38.fasta
│       │   ├── Homo_sapiens_assembly38.fasta.fai
│       │   └── Homo_sapiens_assembly38.dict
│       │   └── Homo_sapiens_assembly38.fasta.64.amb     
│       │   └── Homo_sapiens_assembly38.fasta.64.ann     
│       │   └── Homo_sapiens_assembly38.fasta.64.bwt     
│       │   └── Homo_sapiens_assembly38.fasta.64.pac    
│       │   └── Homo_sapiens_assembly38.fasta.64.sa     
│       │   └── Homo_sapiens_assembly38.fasta.64.alt                  
│       └── known_sites/       
│       └── bed/               
│       └── somatic_resources/ 
├── data/
│   └── SRR30536566/                
│       ├── raw_fastq/
│       │   ├── SRR30536566_1.fastq.gz
│       │   └── SRR30536566_2.fastq.gz
│       ├── qc/
│           └── multiqc_report.html                               # QC
│           └── multiqc_report_1.html                             # Trimm
│           └── multiqc_report_2.html                             # Align+BAM_prepross
│           └── SRR30536566_1_fastqc.html 
│           └── SRR30536566_R1.trimmed_fastqc.html 
│           └── SRR30536566_2_fastqc.html 
│           └── SRR30536566_R2.trimmed_fastqc.html                  
│       ├── trimmed/
│           └── SRR30536566_R1.trimmed.fastq.gz 
│           └── SRR30536566_R2.trimmed.fastq.gz          
│       ├── aligned/
│           └── SRR30536566.sam                                 # Removed
│           └── SRR30536566.bam                                 # Removed
│           └── SRR30536566.sorted.bam                          # Removed
│           └── SRR30536566.sorted.markdup.bam                  # Removed
│           └── SRR30536566.markdup.metrics.txt
│           └── SRR30536566.sorted.markdup.md.bam               # 
│           └── SRR30536566.sorted.markdup.md.bam.bai
│           
│       ├── variants/          
│       └── annotation/        
├── scripts/
│       └── 01_qc.sh
│       └── 02_trim.sh
│       └── 03_align_&bam_preprocessing.sh                  
└── logs/
        └── cutadapt_SRR30536566.log
        └── bwa_mem.log
        └── markduplicates.log
        └── SRR30536566.flagstat.txt                     
```
