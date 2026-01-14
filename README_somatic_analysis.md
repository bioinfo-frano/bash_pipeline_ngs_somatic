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
Learn read-orientation bias (LearnReadOrientationModel) 
 ↓
FilterMutectCalls (variant filtering)
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
---

### Alignment and BAM preprocessing workflow 👉 [03_align_&_bam_preprocess.sh](bash_scripts/03_align_&_bam_preprocess.sh)

**Documentation**: 

- **BWA**: <https://bio-bwa.sourceforge.net/bwa.shtml>
- **samtools** (v. 1.23): <https://www.htslib.org/doc/samtools.html>
- **Picard**: <https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard>


**Table 2**: Alignment and BAM preprocessing workflow 
**Script**: `03_align_&_bam_preprocess.sh`

| Step | Step name                  | INPUT                                                                                                     | OUTPUT                                                                | Tool                    | Function / Role                                                                                                                                                                                                 |
| ---- | -------------------------- | --------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- | ----------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|  1 | Alignment                  | `SRR30536566_R1.trimmed.fastq.gz`<br>`SRR30536566_R2.trimmed.fastq.gz`<br>`Homo_sapiens_assembly38.fasta` | `SRR30536566.sam`                                                     | `bwa mem`               | Aligns paired-end reads to the reference genome. Adds **Read Group (RG)** information required by GATK and downstream tools. Output is an **unsorted SAM** alignment file.                                          |
|  2 | Conversion SAM → BAM       | `SRR30536566.sam`                                                                                         | `SRR30536566.bam`                                                     | `samtools view`         | Converts human-readable SAM into compressed binary BAM format for efficiency and downstream processing.                                                                                                         |
|  3 | Sort BAM (coordinate sort) | `SRR30536566.bam`                                                                                         | `SRR30536566.sorted.bam`                                              | `samtools sort`         | Sorts alignments by genomic coordinates (chromosome and position). **Required** for duplicate marking, indexing, and variant calling.                                                                           |
|  4 | MarkDuplicates    | `SRR30536566.sorted.bam`                                                                                  | `SRR30536566.sorted.markdup.bam`<br>`SRR30536566.markdup.metrics.txt` |`picard MarkDuplicates` | Identifies PCR/optical duplicates and **marks them in the BAM (FLAG + tags)** without removing reads. Duplicate sets are tagged (amplicon-aware), enabling variant callers to down-weight or ignore duplicates. |
|  5 | Add MD/NM tags             | `SRR30536566.sorted.markdup.bam`<br>`Homo_sapiens_assembly38.fasta`                                       | `SRR30536566.sorted.markdup.md.bam`                                   |`samtools calmd`        | Recalculates and adds **MD** (mismatch positions) and **NM** (Number of mismatches) tags. Improves robustness and compatibility with GATK and somatic variant callers.                                                 |
|  6 | Index final BAM            | `SRR30536566.sorted.markdup.md.bam`                                                                       | `SRR30536566.sorted.markdup.md.bam.bai`                               |`samtools index`        | Creates a BAM index enabling **random genomic access**. Required for variant calling (e.g. Mutect2), visualization (IGV), and QC tools.                                                                         |


### Read Groups (RG): The "Birth Certificate" of Each Read

BWA-MEM assigns **Read Group (RG)** information to each read in the SAM/BAM file.
Read groups provide essential metadata that allows downstream tools—especially GATK—to correctly interpret the origin of each read. See **Table 3** for details.

**Table 3**: GATK-relevant read group fields

| Field | Meaning       | Recommendation      | Script                   | Interpretation
| ----- | ------------- | ------------------- |--------------------------|-----------------
| `ID`  | Read group ID | Unique per lane/run | `RG_ID="SRR30536566"`    | Technical identifier
| `SM`  | Sample        | Biological sample   | `RG_SM="DMBEL-EIDR-071"` | Reflects the biological sample
| `LB`  | Library       | Library preparation | `RG_LB="AMPLICON"`       | Reflects library strategy
| `PL`  | Platform      | Sequencing platform | `RG_PL="ILLUMINA"`       | **Mandatory for GATK**
| `PU`  | Platform unit | Flowcell + lane     | `RG_PU="HiSeq4000"`      | Uniquely identifies the sequencing unit


### Picard MarkDuplicates

After trimming and alignment of sample **SRR30536566**, a duplication level of **~54.6% of read pairs** was observed among ~7.7 million mapped reads (Figure 5).
This high duplication rate is **expected** because the data originate from **PCR-amplified amplicon regions** targeting specific genes in colon cancer cells.
Marked duplicates are **not removed**, but flagged so that downstream tools (e.g. Mutect2) can properly account for them.

**Figure 5.** Percentage of reads categorized by duplication states.

![Figure 5: Markduplicates_picard](images/Markduplicates_picard.png)

### Samtools – Troubleshooting

When following this tutorial, it is very important to be aware of which version of samtools is installed in the conda environment.
In this course, the environment provides **samtools v0.1.19 (released in 2013)**, which is much older than current samtools versions (≥1.x).

Because of this, some command-line options shown in modern tutorials and online documentation **do NOT exist** in this older version. Using newer options with an old version will lead to errors, even if the command looks perfectly correct.

You can always check the installed version with:

```bash
samtools
```

### 🔹 A common student concern: “Is the old version worse?”

During this tutorial, students often feel that using an old samtools version means lower-quality results.
**This is not true**.

✔️ The core algorithms used by samtools such as:

- coordinate sorting of BAM files

- calculation of MD and NM tags

- interpretation of alignments

are **conceptually identical across versions**.

What has mainly improved in newer versions is:

- memory efficiency

- speed

- multi-threading support

- cleaner command-line syntax

The **biological results and interpretation** (sorted BAM order, MD/NM values, variant evidence) remain the same.


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

>**Message about Samtools**:
This pipeline uses an older samtools version (0.1.19), which does not support multithreading (-@).
The MD and NM tags produced here are **biologically identical** to those generated by modern samtools versions; newer versions mainly improve speed and robustness, not the interpretation of variants.


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
│           └── multiqc_report_2.html                             # Align+BAM_preprocess
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
│           └── SRR30536566.sorted.markdup.md.bam               # Further analysis
│           └── SRR30536566.sorted.markdup.md.bam.bai           # Further analysis
│           
│       ├── variants/          
│       └── annotation/        
├── scripts/
│       └── 01_qc.sh
│       └── 02_trim.sh
│       └── 03_align_&_bam_preprocess.sh                  
└── logs/
        └── cutadapt_SRR30536566.log
        └── bwa_mem.log
        └── markduplicates.log
        └── SRR30536566.flagstat.txt                     
```

>**Key idea**:
>Each step of the [03_align_&_bam_preprocess.sh](bash_scripts/03_align_&_bam_preprocess.sh) script transforms the data into a format that is progressively **more structured, annotated, and analysis-ready**. 
By the end of this script, the BAM-generated file is:
>- sorted
>- duplicate-aware
>- MD/NM-tagged
>- indexed 
>This is the **expected starting point** for somatic variant calling with GATK **Mutect2**.
The pipeline is fully GATK-compatible and intentionally uses a legacy samtools version to keep the computational environment simple and reproducible.

---
---

### Variant calling with Mutect2 👉 [04_mutect2.sh](bash_scripts/04_mutect2.sh)

Variant calling is a bioinformatics process that identifies differences (variants) between a sample's DNA sequence and a reference genome. 


**Types of genetic variants associated with disease**:

| Variant Class      | Example          | Biological Meaning                    | Phenotypic Effect              | Example                         |
| ------------------ | ---------------- | ------------------------------------- | ------------------------------ | ------------------------------- |
| **Synonymous SNP** | GAA → GAG        | Base change without amino acid change | Usually none (neutral)         | Common population polymorphisms |
| **Missense SNP**   | A → T            | Amino acid substitution               | Variable (mild → severe)       | Sickle cell anemia (HBB)        |
| **Nonsense SNP**   | C → T            | Premature stop codon                  | Loss of protein function       | Duchenne muscular dystrophy     |
| **Indel**          | −3 bp            | Small insertion/deletion              | Frameshift or in-frame         | Cystic fibrosis (ΔF508)         |
| **CNV**            | Exon duplication | Copy number change                    | Gene dosage imbalance          | Charcot–Marie–Tooth disease     |
| **SV**             | t(9;22)          | Large chromosomal rearrangement       | Fusion protein / misregulation | Chronic myeloid leukemia        |



**Genetic Variant Types and Associated Diseases**

| Variant Type                      | Example          | Size Range | Molecular Consequence    | Disease Example                 | Gene Involved |
| --------------------------------- | ---------------- | ---------- | ------------------------ | ------------------------------- | ------------- |
| **Synonymous SNP (silent)**       | GAA → GAG        | 1 bp       | No amino acid change     | *None (neutral variant)*        | HBB           |
| **Missense SNP (non-synonymous)** | GAG → GTG        | 1 bp       | Amino acid substitution  | **Sickle cell anemia**          | HBB           |
| **Nonsense SNP (non-synonymous)** | CGA → TGA        | 1 bp       | Premature stop codon     | **Duchenne muscular dystrophy** | DMD           |
| **Indel (frameshift)**            | CTT deletion     | 1–50 bp    | Reading-frame disruption | **Cystic fibrosis**             | CFTR          |
| **Repeat expansion (Indel)**      | CAG expansion    | Variable   | Toxic protein elongation | **Huntington’s disease**        | HTT           |
| **CNV (duplication/deletion)**    | Exon duplication | 50 bp–5 Mb | Gene dosage alteration   | **Williams syndrome**           | ELN           |
| **SV (translocation)**            | t(9;22)          | >50 bp     | Fusion protein formation | **Chronic myeloid leukemia**    | BCR–ABL1      |
| **SV (inversion)**                | F8 inversion     | Variable   | Gene disruption          | **Hemophilia A**                | F8            |


**Cancer-Related Genetic Variants**

| Variant Type            | Cancer Example            | Oncogenic Mechanism   | Biological Consequence       | Targeted Therapy     |
| ----------------------- | ------------------------- | --------------------- | ---------------------------- | -------------------- |
| **Missense SNP**        | **BRAF V600E**            | Activating mutation   | Constitutive MAPK signaling  | Vemurafenib          |
| **Indel**               | **EGFR exon 19 deletion** | Gain-of-function      | Persistent EGFR activation   | Erlotinib, Gefitinib |
| **CNV (Amplification)** | **HER2 amplification**    | Increased gene dosage | Receptor overexpression      | Trastuzumab          |
| **SV (Fusion)**         | **BCR–ABL1** t(9;22)      | Gene fusion           | Constitutive tyrosine kinase | Imatinib             |
| **Missense SNP**        | **TP53 R175H**            | Loss-of-function      | Impaired DNA damage response | No direct therapy    |


In **Cancer Genomics**, variant calling analysis helps to:

1. Identify driver mutations that cause cancer

2. Find therapeutic targets (like EGFR mutations for lung cancer)

3. Track tumor evolution over time

4. Guide precision medicine (matching drugs to mutations)

Other applications:

5. Diagnose genetic diseases (finding mutations causing rare diseases)

6. Pharmacogenomics (predicting drug response)

7. Population genetics (studying human evolution)

8. Forensics (DNA fingerprinting)


### Variant callers

| Tool                | Best suited for                                             | Strengths                                                                                                                                     | Limitations                                                          |
| ------------------- | ----------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------- |
| **Mutect2 (GATK)**  | **Somatic variants in cancer** (tumor-only or tumor–normal) | ✔ Excellent sensitivity at low VAF<br>✔ Sophisticated error modeling<br>✔ Uses PoN + population AFs<br>✔ Industry standard in cancer genomics | Requires multiple input resources; more complex                      |
| **HaplotypeCaller** | Germline variant discovery                                  | ✔ Excellent for inherited variants<br>✔ Accurate diploid genotyping                                                                           | ❌ Not designed for somatic variants<br>❌ Poor sensitivity at low VAF |
| **VarScan2**        | Simple germline/somatic calling                             | ✔ Works at low coverage<br>✔ Easy to run                                                                                                      | ❌ High false-positive rate<br>❌ Limited error modeling               |
| **Strelka2**        | Somatic variants (especially indels)                        | ✔ Very accurate indel calling<br>✔ Fast in tumor–normal mode                                                                                  | ❌ Less robust in tumor-only mode<br>❌ Less transparent filtering     |
| **FreeBayes**       | Germline / pooled samples                                   | ✔ Flexible calling models                                                                                                                     | ❌ Not optimized for cancer somatic variants                          |

### Mutect2 is the correct choice for this tutorial

✔ Designed specifically for somatic cancer mutations

✔ Handles tumor-only data using PoN + gnomAD

✔ Detects low-frequency variants typical of cfDNA and amplicon panels

✔ Integrates seamlessly with GATK filtering and annotation steps

✔ Widely used in research and clinical pipelines


### Essential Files for Somatic Variant Calling with GATK Mutect2

Mutect2 compares tumor sequencing data against a reference genome and multiple external resources to distinguish **true somatic mutations** from **technical artifacts** and **germline variation**.

| #     | File category                          | Conceptual view (what this represents)                                                                                             | Purpose in somatic variant calling                                                                        | Example file(s)                                                                | Required                |
| ----- | -------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------ | ----------------------- |
| **1** | **Tumor BAM (aligned reads)**          | The *digital representation of the patient’s DNA*, where each read is a fragment sequenced from the tumor and mapped to the genome | Provides the raw evidence used to detect mismatches, insertions, and deletions relative to the reference  | `SRR30536566.sorted.markdup.md.bam`<br>`SRR30536566.sorted.markdup.md.bam.bai` | ✅ Yes                   |
| **2** | **Reference genome**                   | The *coordinate system and “normal” genome* against which tumor reads are compared                                                 | Defines genomic positions, enables alignment interpretation, and serves as baseline for variant detection | `Homo_sapiens_assembly38.fasta`                                                | ✅ Yes                   |
|       | Reference index files                  | The *navigation tools* that allow fast access to the reference genome                                                              | Required for efficient random access and full GATK compatibility                                          | `.fai` (samtools)<br>`.dict` (Picard)                                          | ✅ Yes                   |
| **3** | **Germline variant resource (gnomAD)** | A *population-level catalogue of normal human variation*                                                                           | Helps Mutect2 estimate whether a variant is likely germline vs somatic by using allele frequencies        | `af-only-gnomad.hg38.vcf.gz`<br>`af-only-gnomad.hg38.vcf.gz.tbi`               | ✅ Yes (tumor-only mode) |
| **4** | **Panel of Normals (PoN)**             | A *background noise model* built from many normal samples                                                                          | Removes recurrent technical artifacts that appear across samples but are not true mutations               | `1000g_pon.hg38.vcf.gz`<br>`1000g_pon.hg38.vcf.gz.tbi`                         | ⚠️ Strongly recommended |
| **5** | **BED / interval file**                | The *map of genomic regions that were actually sequenced*                                                                          | Restricts variant calling to targeted regions, reducing false positives and runtime                       | `crc_panel_7genes_sorted.hg38.bed`                                                       | ⚠️ Recommended          |
| **6** | **dbSNP (optional)**                   | A *catalogue of known common polymorphisms*                                                                                        | Used mainly for annotation and interpretation, not required by Mutect2 itself                             | `dbsnp_146.hg38.vcf.gz`                                                        | ❌ Optional              |
| **7** | **COSMIC (optional)**                  | A *knowledge base of known cancer mutations*                                                                                       | Enables biological interpretation and clinical relevance assessment after variant calling                 | `CosmicMutantExport.*`                                                         | ❌ Optional              |
| **8** | **F1R2 artifact data**                 | A *model of strand-orientation sequencing artifacts*                                                                               | Required to filter FFPE and orientation bias artifacts in downstream filtering                            | `SRR30536566.f1r2.tar.gz`                                                      | ✅ Yes (for filtering)   |

**Documentation**: 

- **Mutect2**: <https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2>


### Verify and download the essential files to run Mutect2

1. **Check presence of these files**:

✔ /Genomics_cancer/data/SRR30536566/aligned/SRR30536566.sorted.markdup.md.bam

✔ /Genomics_cancer/data/SRR30536566/aligned/SRR30536566.sorted.markdup.md.bam.bai

✔ /Genomics_cancer/reference/GRCh38/fasta/Homo_sapiens_assembly38.fasta

✔ /Genomics_cancer/reference/GRCh38/fasta/Homo_sapiens_assembly38.dic

✔ /Genomics_cancer/reference/GRCh38/fasta/Homo_sapiens_assembly38.fai

❌  `SRR30536566.f1r2.tar.gz` is an output from variant calling. The file will be used for **learn read-orientation bias**, a common artifact in Illumina sequencing (especially strong in amplicon and cfDNA data).

2. **Create folders**

Go to: /Genomics_cfDNA_SRR15506490/Somatic_SRR15506490/reference/GRCh38

Create these folders "somatic_resources" and "intervals":

```bash
mkdir somatic_resources intervals
```

3. **Download these files to ~/somatic_resources**

```bash
# gnomAD AF-only VCF
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz

# gnomAD index
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi

# PoN
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz

# PoN index
wget https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi
```

4. **Generation of BED file (if authors provided no BED)**

  - 4.1. Go to /Genomics_cfDNA_SRR15506490/Somatic_SRR15506490/reference/GRCh38/intervals

  - 4.2. Link to download file with `wget`

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
```
>**Note**: Do not download the 'GGF3' (from <https://www.gencodegenes.org/human/release_38.html>) nor the `Homo_sapiens_assembly38.contam.bed` (from <https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0?>) as BED files.

>**Note2**: GTF file structure <http://www.ensembl.org/info/website/upload/gff.html>:
```bash
field 1  chrom
field 2  source - name of the program that generated this feature, or the data source (database or project name)
field 3  feature (gene, transcript, exon, …)
field 4  start - Start position* of the feature.
field 5  end - End position* of the feature.
field 6  score - A floating point value.
field 7  strand - defined as + (forward) or - (reverse).
field 8  frame
field 9  attributes - a single long string semicolon-separated.
```

  - 4.3. Decompress

```bash
gunzip gencode.v38.annotation.gtf.gz
```
  - 4.4. Verify the presence of gene targets KRAS, NRAS, BRAF, PIK3CA, PTEN, RRAS, and MAP2K1 (MEK1) in **.gtf** file

>**Note**: The authors of dataset "SRR15506490" in <https://www.ncbi.nlm.nih.gov/sra/SRX25960056> point out in the following: 
"**Design**: Targeted sequencing of full-length, KRAS NRAS BRAF PIK3CA PTEN RRAS and MEK1 including UTR, exons, and introns."

  Use this code to verify genes individually:

```bash
less gencode.v38.annotation.gtf.gz | grep -w "gene" | head -n 3         # Copy/Paste the gene target in "gene" and hit Enter
```

  Use this code to verify all genes with a for-loop:

```bash

```

  - 4.5. Create the .bed file and sorted .bed file by running 👉 [make_crc_7genes_bed.sh](make_crc_7genes_bed.sh) 
  The `crc_panel_7genes.hg38.bed` and `crc_panel_7genes_sorted.hg38.bed` are located in: /Genomics_cancer/reference/GRCh38/intervals
  
  - 4.6. Lastly, check whether the `SRR30536566.sorted.markdup.md.bam` has RG (Read Groups)
  
  Run 
  ```bash
  samtools view -H SRR30536566.sorted.markdup.md.bam | grep "@RG"
  ```
  Expected output:
  ```bash
  @RG	ID:SRR30536566	SM:DMBEL-EIDR-071	LB:AMPLICON	PL:ILLUMINA	PU:HiSeq4000
@PG	ID:bwa	PN:bwa	VN:0.7.19-r1273	CL:bwa mem -t 4 -R @RG\tID:SRR30536566\tSM:DMBEL-EIDR-071\tLB:AMPLICON\tPL:ILLUMINA\tPU:HiSeq4000
  ```



















1️⃣ How to obtain / generate SRR30536566.f1r2.tar.gz
✅ Your intuition is 100% correct

You do NOT obtain this file beforehand

It is generated during the Mutect2 run

It is not an input, but an intermediate output used for downstream filtering

What this file is (conceptually)

*.f1r2.tar.gz contains data used to learn read-orientation bias, a common artifact in Illumina sequencing (especially strong in amplicon and cfDNA data).

Conceptually:

“Are false variants appearing preferentially on reads in one orientation (F1R2 vs F2R1)?”

This information is learned from your sample, not from a database.

📌 Key point:
Mutect2 writes this file while scanning reads for candidate variants.

