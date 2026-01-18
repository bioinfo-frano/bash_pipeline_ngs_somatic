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

## Bioinformatic workflows/pipelines

DNA-NGS analysis can be automated using several workflow approaches, including:
- Bash scripting
- Nextflow
- Snakemake
- Galaxy
- Others (Cromwell/WDL, Toil, Airflow, Argo Workflows, Make, Ruffus, Bpipe)

>**Note**: Nextflow, Snakemake, Cromwell+WDL, CWL, and Galaxy are all considered workflow managers or workflow management systems. Bash scripting is not considered a workflow manager. It is a shell scripting language. 


**Table 1**: Advantages and disadvantages of selected workflow approaches.

| Tool | Advantages | Disadvantages |
|------|------------|---------------|
| **Bash scripting** | • Maximum transparency and control<br>• No new syntax to learn<br>• Simple to start<br>• Universally available<br>• Ideal for learning and prototyping | • Limited scalability<br>• Manual error handling<br>• No native workflow management |
| **Snakemake** | • Python-based and readable<br>• Explicit file dependencies<br>• Excellent for complex workflows | • Scaling requires configuration<br>• Less cloud-native by default |
| **Nextflow** | • Highly portable (“write once, run anywhere”)<br>• Built-in parallelization and resume<br>• Strong community (nf-core) | • Learning curve (Groovy DSL)<br>• Debugging channels can be challenging |

---

## Bioiformatics tools

**GATK** (Genome Analysis Toolkit) is a comprehensive **software package** developed by the Broad Institute for analyzing high-throughput sequencing data, with a primary **focus on variant discovery** (**SNPs and indels**) in human genetics and medical genomics. 
It provides a step-by-step framework (pipeline) from raw reads → variant calls. To achieve that, GATK has a collection of tools specialized for:

- Data preprocessing (BQSR, MarkDuplicates)

- Variant calling (HaplotypeCaller, Mutect2)

- Variant filtering and refinement (VQSR, CNN)

- Variant manipulation (SelectVariants, CombineVariants)


### Why GATK has more advantages for SNP/Indel analysis

1. **Comprehensive, Validated Best Practices**

    - End-to-end pipeline: GATK provides a complete, battle-tested workflow from raw FASTQ to final VCF

    - Constantly updated: Regular releases incorporating latest research (e.g., CNN filtering)

    - Extensive documentation: The "GATK Best Practices" guides are the de facto standard in clinical genomics

2. **Sophisticated Error Modeling & Correction**

    - Base Quality Score Recalibration (BQSR): Corrects systematic errors in base quality scores

    - Variant Quality Score Recalibration (VQSR): Machine learning approach to filter variants using known resources (HapMap, 1000 Genomes, etc.)

    - Local realignment around indels: Critical for accurate indel calling in repetitive regions

3. **Specialized Calling Algorithms**

    - HaplotypeCaller: Uses local de novo assembly of haplotypes, which is superior for:

         - Calling in difficult regions (low complexity, repeats)

         - Resolving phase (which variants are on same chromosome)

         - Accurate indel calling (by realigning reads to candidate haplotypes)

    - Mutect2: State-of-the-art somatic caller with advanced tumor heterogeneity modeling

4. **Extensive Quality Control & Filtering**

    - Multiple filtering tiers: Hard filters, VQSR, and CNN-based filters

    - Rich annotation: Adds >100 annotations to each variant (QD, FS, MQ, etc.) for informed filtering

    - Joint genotyping: Enables population-level analysis by calling variants across multiple samples simultaneously

5. **Industry Adoption & Validation**

    - Clinical validation: Used in thousands of clinical labs and research studies

    - Benchmarking: Regularly tops precisionFDA challenges and other benchmarks

    - Integration: Works seamlessly with other Broad tools (Cromwell, WDL, Terra)

6. **Active Community & Support**

    - Forum support: Active community with Broad engineers answering questions

    - Training: Regular workshops and detailed online materials

    - Interoperability: Compatible with major public databases (gnomAD, dbSNP, ClinVar)


**GATK dominates clinical and population-scale human genomics because it offers a rigorously validated, end-to-end solution with sophisticated error correction that's essential for reliable variant calling**. Its "Best Practices" represent 10+ years of accumulated knowledge about NGS artifacts and how to mitigate them.

Using GATK (particularly Mutect2 for somatic, HaplotypeCaller for germline) would be considered the gold standard in clinical genomics, providing the accuracy and reliability needed for treatment decisions. However, for non-human projects or quick analyses, simpler tools might be more appropriate.

While GATK excels at SNP and small indel analysis, **structural variant (SV)** analysis requires specialized tools that use different algorithms and signal detection methods.

### Structural Variant (SV) Analysis Tools

| Tool | Description | Primary Use / Strength |
|------|-------------|-------------------------|
| **Manta** | Fast SV and indel caller for mapped sequencing reads | Excellent for germline and somatic SV detection in paired tumor-normal samples |
| **Delly** | Integrated SV caller for paired-end, split-read, and read-pair signals | Good for both germline and somatic SV across multiple SV types (DEL, DUP, INV, BND) |
| **Lumpy** | Probabilistic SV detection using multiple signals (split-read, read-pair) | High sensitivity for breakpoint detection, often used in research |
| **GRIDSS** | Genome Rearrangement Identification Software Suite | Somatic SV calling with high precision, excellent for complex rearrangements |
| **SvABA** | Genome-wide SV and indel caller using local assembly | Optimized for tumor samples, can detect ultra-rare somatic SVs |
| **cn.MOPS** | Copy Number variation detection via Mixture Of PoissonS | CNV detection from NGS data, good for population-level analysis |
| **ERDS** | Estimation by Read Depth with SNV | CNV detection specifically for whole-genome sequencing data |
| **TARDIS** | Targeted Assembly for Rearrangement Detection in Sequencing | Local assembly-based SV calling, good for complex regions |
| **WHAM** | Whole-genome Association Method | Integrates multiple SV signals, good for association studies |
| **Sniffles2** | Fast SV caller for long-read sequencing (PacBio, ONT) | **For long-read data**, the gold standard for SV detection |

### Integrated SV Analysis Platforms
| Tool | Description |
|------|-------------|
| **SURVIVOR** | Toolset for simulating, comparing, and merging SV calls | Post-processing of multiple SV callers to create consensus calls |
| **Truvari** | SV benchmarking and comparison toolkit | Quality assessment and comparison of SV callers |
| **VariantAnnotation** (R/Bioconductor) | R package for variant annotation and interpretation | Annotating and filtering SV calls in R environment |

### When to use these toos instead of GATK:

  - Large deletions/duplications (>50 bp)

  - Chromosomal rearrangements (translocations, inversions)

  - Complex events (tandem duplications, breakend events)

  - Copy number variations (CNVs) across the genome

   - Long-read sequencing data (PacBio, Oxford Nanopore)


### Why GATK offers different tools for somatic and germline variant analysis?

The reason GATK has separate, specialized tools for somatic vs. germline analysis comes down to fundamental biological and technical differences that require different statistical models and filters.

### Core Biological Differences between these tools:

| Aspect | Germline Variants (HaplotypeCaller) | Somatic Variants (Mutect2) |
|--------|-----------------------------------|---------------------------|
| **Origin** | Inherited from parents, present in **every cell** of the body | Acquired during life, present **only in tumor/cancer cells** |
| **Allele Frequency** | ~50% (heterozygous) or ~100% (homozygous) in normal tissue | Can be **very low** (5-20%) due to tumor heterogeneity & purity |
| **"Normal" Sample** | No paired normal needed (but can use) | **REQUIRES matched normal** from same patient to subtract germline background |
| **Goal** | Find **all** inherited variants at appreciable frequency | Find **new, cancer-driving** mutations against patient's germline background |

### Why Different Algorithms Are Needed

**HaplotypeCaller (Germline)**

   - Assumes variants are at high frequency (~50% or 100%). 

   - Optimized for sensitivity to catch all potential Mendelian variants

   - Uses local de novo assembly to build haplotypes, excellent for tricky regions

   - Joint genotyping allows calling across multiple samples to improve accuracy

   - Filters common sequencing artifacts that appear at lower frequencies

**Mutect2 (Somatic)**

   - Designed to find needles in a haystack - very low allele frequency variants

   - Extensive error modeling to distinguish real somatic mutations from:

        - Sequencing errors

        - PCR artifacts

        - Mapping errors

       - Normal cell contamination in tumor sample

   - Panel of Normals (PON) - Uses data from many normal samples to identify and remove systematic artifacts

   - Tumor-specific filters that consider:

       - Tumor purity (what % of sample is actually cancer cells)

       - Ploidy (is the tumor diploid or aneuploid?)

       - Clonality (is the mutation present in all or just some tumor cells?)
        
**Key Technical Divergences**

   - Statistical Models:

        - HaplotypeCaller: Bayesian genotype likelihood model assuming diploidy

        - Mutect2: Somatic likelihood model that doesn't assume diploidy in tumor

   - Error Handling:

       - HaplotypeCaller: Filters out low-frequency noise (assumes real variants are high frequency)

       - Mutect2: Aggressively models and removes noise at ALL frequencies

   - Output Differences:

        - HaplotypeCaller: Outputs genotypes (0/0, 0/1, 1/1)

        - Mutect2: Outputs somatic status (REF, SOMATIC, GERMLINE) and tumor LOD scores

>**Bottom line**: The biological and technical differences between inherited and acquired mutations are so profound that specialized tools with different statistical foundations are necessary. GATK provides both because in clinical cancer genomics, you typically need to perform **both analyses** on the same patient—somatic for treatment targets, germline for inherited risk assessment.


In this **Part II**, bash scripting is used to perform a **somatic DNA-NGS analysis** on the small FASTQ dataset introduced in [Part I – Preparation & setup](README_setup.md) and find SNP and Indels using GATK toolkit.

---

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
  
The red warning is because FastQC compares your distribution to a unimodal model based on a normal genome. Your amplicon-based distribution violates this model, so it gets flagged. This is **NOT** a problem for your data.

**Figure 4.** Per sequence GC content showing a bimodal shape.
  
  
  ![Figure 4: Per sequence GC content](images/FastQ_Per_Sequence_GC_content.png)

---

### Trimming + QC 👉 [02_trim.sh](bash_scripts/02_trim.sh)

#### Cutadapt

**Documentation for Cutadapt**:

<https://github.com/marcelm/cutadapt/blob/main/doc/guide.rst> Relevant sections include "***Trimming paired-end reads***" and "***Cutadapt's output***".

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

The alignment steps and processing of BAM file are explained in **Table 2A** and the inputs/outputs files in **Table 2B**.

**Table 2A**: Alignment and BAM Processing Workflow Steps
| Step | Tool | Function / Role |
|------|------|-----------------|
| **1**<br>**Alignment** | `bwa mem` | Aligns paired-end reads to the reference genome.<br>Adds **Read Group (RG)** information required by GATK and downstream tools.<br> Output is an **unsorted SAM** alignment file. |
| **2**<br>**SAM → BAM Conversion** | `samtools view` | Converts human-readable SAM into compressed binary BAM format for efficiency and downstream processing.  |
| **3**<br>**Sort BAM** | `samtools sort` | Sorts alignments by genomic coordinates (chromosome and position).<br>**Required** for duplicate marking, indexing, and variant calling.  |
| **4**<br>**Mark Duplicates** | `picard MarkDuplicates` | Identifies PCR/optical duplicates and **marks them in the BAM (FLAG + tags)** without removing reads.<br>Duplicate sets are tagged (amplicon-aware), enabling variant callers to down-weight or ignore duplicates. |
| **5**<br>**Add MD/NM Tags** | `samtools calmd` | Recalculates and adds **MD** (mismatch positions) and **NM** (Number of mismatches) tags.<br>Improves robustness and compatibility with GATK and somatic variant callers. |
| **6**<br>**Index BAM** | `samtools index` | Creates a BAM index enabling **random genomic access**.<br>Required for variant calling (e.g. Mutect2), visualization (IGV), and QC tools. |

**Table 2B**: Input/Output Files
| Step | Input Files | Output Files |
|------|-------------|--------------|
| **1** | `*_R1.trimmed.fastq.gz`, `*_R2.trimmed.fastq.gz`, `Homo_sapiens_assembly38.fasta` | `SRR30536566.sam` |
| **2** | `SRR30536566.sam` | `SRR30536566.bam` |
| **3** | `SRR30536566.bam` | `SRR30536566.sorted.bam` |
| **4** | `SRR30536566.sorted.bam` | `SRR30536566.sorted.markdup.bam`, `*.markdup.metrics.txt` |
| **5** | `SRR30536566.sorted.markdup.bam`, `Homo_sapiens_assembly38.fasta` | `SRR30536566.sorted.markdup.md.bam` |
| **6** | `SRR30536566.sorted.markdup.md.bam` | `SRR30536566.sorted.markdup.md.bam.bai` |


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
│       └── 0_wget_Hsapiens_assem38.sh
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
>Each step of the [alignment & BAM preprocessing](bash_scripts/03_align_&_bam_preprocess.sh) script transforms the data into a format that is progressively **more structured, annotated, and analysis-ready**. 
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

**Variant calling** is the bioinformatics process of identifying positions in a sequenced genome where the sample’s DNA differs from a reference genome. See examples of variant types related to various diseases and in cancer in **Table 4** and **Table 5**.

**Table 4**: Genetic Variant Types and Associated Diseases

| Variant Type                               | Example          | Size Range | Molecular Consequence    | Typical Biological Effect               | Disease Example                | Gene     |
| ------------------------------------------ | ---------------- | ---------- | ------------------------ | --------------------------------------- | ------------------------------ | -------- |
| **Synonymous SNP (silent)**                | GAA → GAG        | 1 bp       | No amino acid change     | Usually neutral                         | None (population polymorphism) | HBB      |
| **Missense SNP (non-synonymous)**          | GAG → GTG        | 1 bp       | Amino acid substitution  | Altered protein function (gain or loss) | Sickle cell anemia             | HBB      |
| **Nonsense SNP (non-synonymous)**          | CGA → TGA        | 1 bp       | Premature stop codon     | Truncated protein, loss-of-function     | Duchenne muscular dystrophy    | DMD      |
| **Indel (frameshift)**                     | CTT deletion     | 1–50 bp    | Reading-frame disruption | Severely altered protein                | Cystic fibrosis (ΔF508)        | CFTR     |
| **Indel (in-frame / repeat expansion)**    | CAG expansion    | Variable   | Protein elongation       | Toxic gain-of-function                  | Huntington’s disease           | HTT      |
| **Copy Number Variant (CNV)**<br>**(duplication/deletion)** | Exon duplication | 50 bp–5 Mb | Gene dosage alteration   | Overexpression or<br>haploinsufficiency    | Williams syndrome | ELN      |
| **Structural Variant**<br>**(SV: translocation)** | t(9;22)   | >50 bp     | Fusion gene formation    | Constitutive signaling                  | Chronic myeloid leukemia       | BCR–ABL1 |
| **Structural Variant**<br>**(SV: inversion)** | F8 inversion   | Variable   | Gene disruption          | Loss of gene function                   | Hemophilia A                   | F8       |


In **cancer genomics**, variant calling focuses on identifying somatic variants:

- Mutations present in tumor cells

- Absent (or rare) in the germline

- Often present at subclonal allele frequencies

Therefore, **cancer genomics** is a disciplie that helps to:

- Identify driver mutations that cause cancer

- Find therapeutic targets (like EGFR mutations for lung cancer)

- Reconstruction of tumor evolution and clonal architecture

- Support for precision oncology (matching drugs to mutations)

Additional applications include:

- Rare disease diagnosis (finding mutations causing rare diseases)

- Pharmacogenomics (predicting drug response)

- Population genetics (studying human evolution)

- Forensic genomics (DNA fingerprinting)

**Table 5**: Cancer-Related Genetic Variants

| Variant Type            | Cancer Example            | Oncogenic Mechanism   | Biological Consequence       | Targeted Therapy     |
| ----------------------- | ------------------------- | --------------------- | ---------------------------- | -------------------- |
| **Missense SNP**        | **BRAF V600E**            | Activating mutation   | Constitutive MAPK signaling  | Vemurafenib          |
| **Indel**               | **EGFR exon 19 deletion** | Gain-of-function      | Persistent EGFR activation   | Erlotinib, Gefitinib |
| **CNV (Amplification)** | **HER2 amplification**    | Increased gene dosage | Receptor overexpression      | Trastuzumab          |
| **SV (Fusion)**         | **BCR–ABL1** t(9;22)      | Gene fusion           | Constitutive tyrosine kinase | Imatinib             |
| **Missense SNP**        | **TP53 R175H**            | Loss-of-function      | Impaired DNA damage response | No direct therapy    |


### Variant callers

| Tool                | Best suited for                                             | Strengths                                                                                                                                     | Limitations                                                          |
| ------------------- | ----------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------- |
| **Mutect2 (GATK)**  | **Somatic variants in cancer** (tumor-only or tumor–normal) | ✔ Excellent sensitivity at low VAF<br>✔ Sophisticated error modeling<br>✔ Uses PoN + population AFs<br>✔ Clinical standard in cancer genomics | Requires multiple input resources; more complex                      |
| **HaplotypeCaller** | Germline variant discovery                                  | ✔ Excellent for inherited variants<br>✔ Accurate diploid genotyping                                                                           | ❌ Not designed for somatic variants<br>❌ Poor sensitivity at low VAF |
| **VarScan2**        | Simple germline/somatic calling                             | ✔ Works at low coverage<br>✔ Easy to run                                                                                                      | ❌ High false-positive rate<br>❌ Limited error modeling               |
| **Strelka2**        | Somatic variants (especially indels)<br>(tumor–normal)      | ✔ Very accurate indel calling<br>✔ Fast in tumor–normal mode                                                                                  | ❌ Less robust in tumor-only mode<br>❌ Less transparent filtering     |
| **FreeBayes**       | Germline / pooled samples                                   | ✔ Flexible calling models                                                                                                                     | ❌ Not optimized for cancer somatic variants                          |

### Mutect2 is the correct choice for this tutorial

✔ Designed specifically for somatic cancer mutations

✔ Supports tumor-only analysis using PoN + gnomAD

✔ Detects low-frequency variants (low-VAF) typical of cfDNA and amplicon panels

✔ Integrates downstream artifact and contamination filtering, seamlessly integrates with GATK filtering and annotation steps

✔ Widely adopted in research and clinical pipelines

> **Important point**
> Mutect2 is designed primarily to detect:
>  - SNPs
>  - small indels
> It is **not** intended for reliable detection of:
>  - large CNVs
>  - large structural variants
> Those require different tools (e.g., CNVkit, DELLY, Manta).


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

3. **Download Germline variant resource (gnomAD) and Panel of Normals (PoN) files to ~/somatic_resources**. 

    Run 👉 [0_wget_gnomad_PoN.sh](bash_scripts/0_wget_gnomad_PoN.sh) 

4. **Download GTF file to BED file generation (if authors provided no BED)**

  - 4.1. Go to /Genomics_cfDNA_SRR15506490/Somatic_SRR15506490/reference/GRCh38/intervals

  - 4.2. Download **GTF** file with `wget`

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

>"**Design**: Targeted sequencing of full-length, KRAS NRAS BRAF PIK3CA PTEN RRAS and MEK1 including UTR, exons, and introns."

  Use this code to verify genes individually:

```bash
less gencode.v38.annotation.gtf.gz | grep -w "gene" | head -n 3     # Copy/Paste the gene target in "gene" and hit Enter
```

  Use this code to verify all genes with a for-loop:

```bash

```

5. **Generation of BED file (if authors provided no BED)**

Create the **.bed** file and **sorted.bed** file by running 👉 [04_make_crc_7genes_bed.sh](bash_scripts/04_make_crc_7genes_bed.sh) 
  
The `crc_panel_7genes.hg38.bed` and `crc_panel_7genes_sorted.hg38.bed` are located in: /Genomics_cancer/reference/GRCh38/intervals
  

6. Lastly, check whether the `SRR30536566.sorted.markdup.md.bam` has RG (Read Groups)
  
  Run 
  ```bash
  samtools view -H SRR30536566.sorted.markdup.md.bam | grep "@RG"
  ```
  Expected output:
  ```bash
  @RG	ID:SRR30536566	SM:DMBEL-EIDR-071	LB:AMPLICON	PL:ILLUMINA	PU:HiSeq4000
@PG	ID:bwa	PN:bwa	VN:0.7.19-r1273	CL:bwa mem -t 4 -R @RG\tID:SRR30536566\tSM:DMBEL-EIDR-071\tLB:AMPLICON\tPL:ILLUMINA\tPU:HiSeq4000
  ```
7. Now, run 👉 [04_mutect2.sh](bash_scripts/04_mutect2.sh)


### Mutect2 output files (all expected, all correct)

| File                            | Status | Meaning                  |
| ------------------------------- | ------ | ------------------------ |
| `SRR30536566.unfiltered.vcf.gz` | ✅      | Raw somatic calls        |
| `.tbi`                          | ✅      | Indexed                  |
| `.stats`                        | ✅      | Internal GATK statistics |
| `SRR30536566.f1r2.tar.gz`       | ✅      | Strand artifact data     |
| Logs                            | ✅      | Clean                    |

### Mutect2 logs

a) mutect2.stderr.log

```bash
Processing 573443 bp from intervals
```

**Interpretation**: It matches the 7-gene CRC panel, calculated from `crc_panel_7genes_sorted.hg38.bed`
The "573443 bp" is the sum of all 7-gene lengths.

```bash
Final CRC 7-gene BED file contents:
Chr Gene_Start  Gene_End  Gene_Name
chr1	114704468	114716771	NRAS
chr10	87863624	87971930	PTEN
chr12	25205245	25250936	KRAS
chr15	66386836	66491544	MAP2K1
chr19	49635291	49640143	RRAS
chr3	179148113	179240093	PIK3CA
chr7	140719326	140924929	BRAF
```
```bash
Unable to find native library: native/libgkl_pairhmm_omp.dylib
Ignoring request for 4 threads
```
**Interpretation**: Expected on macOS; Does NOT affect results, only speed.

```bash
832807 read(s) filtered by: NotDuplicateReadFilter 
863479 total reads filtered out of 1479947 reads processed
```
**Interpretation**: Confirms MarkDuplicates worked; Mutect2 behaves as expected; 832807 duplicates removed.

b) mutect2.stdout.log

```bash
Tool returned:
SUCCESS
```
**Interpretation**: self explanatory.

**Verdict on Mutect2 variant calling step**: There are no errors, no conceptual problems, and no missing inputs.


>**Important:**
>
>Mutect2 does not “decide” whether a variant is cancer or not in isolation.
>
> It **models probabilities** using:
> - sequencing evidence
>
> - population allele frequencies
>
> - technical artifact profiles
>
> - contamination estimates
>
> Final biological interpretation happens after filtering and annotation.

### Folder structure: necessary and output files from Mutect2-Variant calling

```bash
Genomics_cancer/
├── reference/                 
│   └── GRCh38/
│       ├── fasta/
│       └── known_sites/       
│       └── intervals/                                           # BED files
│           └── gencode.v38.annotation.gtf
│           └── gencode.v38.annotation.gtf.gz
│           └── crc_panel_7genes.hg38.bed
│           └── crc_panel_7genes_sorted.hg38.bed                 # For Mutect2 analysis
│       └── somatic_resources/                                   # PoN & GNOMAD: For Mutect2 analysis
│           └── 1000g_pon.hg38.vcf.gz
│           └── 1000g_pon.hg38.vcf.gz.tbi
│           └── af-only-gnomad.hg38.vcf.gz
│           └── af-only-gnomad.hg38.vcf.gz.tbi
├── data/
│   └── SRR30536566/                
│       ├── raw_fastq/
│       ├── qc/
│       ├── trimmed/
│       ├── aligned/
│           └── SRR30536566.markdup.metrics.txt
│           └── SRR30536566.sorted.markdup.md.bam               # For Mutect2 analysis
│           └── SRR30536566.sorted.markdup.md.bam.bai           # For Mutect2 analysis
│           
│       ├── variants/                                           # Outputs
│           └── SRR30536566.f1r2.tar.gz          
│           └── SRR30536566.unfiltered.vcf.gz
│           └── SRR30536566.unfiltered.vcf.gz.stats 
│           └── SRR30536566.unfiltered.vcf.gz.tbi
│       └── annotation/        
├── scripts/
│       └── 0_wget_gnomad_PoN.sh
│       └── 0_wget_Hsapiens_assem38.sh
│       └── 01_qc.sh
│       └── 02_trim.sh
│       └── 03_align_&_bam_preprocess.sh
│       └── 04_make_crc_7genes_bed.sh
│       └── 04_mutect2.sh                                 
└── logs/
        └── cutadapt_SRR30536566.log
        └── bwa_mem.log
        └── markduplicates.log
        └── SRR30536566.flagstat.txt                     
        └── mutect2.stderr.log                                  # Output
        └── mutect2.stdout.log                                  # Output
```


### Orientation bias model 👉 [05_learn_read_orientation_model.sh](bash_scripts/05_learn_read_orientation_model.sh)

**Documentation**

- **LearnReadOrientationModel**: <https://gatk.broadinstitute.org/hc/en-us/articles/360051305331-LearnReadOrientationModel>

Learn the prior probability of read orientation artifact from the output of CollectF1R2Counts of Mutect2, which is `~/variants/SRR30536566.f1r2.tar.gz` file.
The `*.f1r2.tar.gz` contains data used to **learn read-orientation bias**, a common artifact in Illumina sequencing (especially strong in amplicon and cfDNA data).

**Conceptually**: “Are false variants appearing preferentially on reads in one orientation (F1R2 vs F2R1)?”

This information is learned from your sample, not from a database, enabling accurate read-orientation artifact filtering.

**Outputs**

a) From Terminal

```bash
Tool returned:
SUCCESS
```

b) **Output files**: `~/variants/SRR30536566.read-orientation-model.tar.gz`  (4.8 KB). Orientation models are always small.

c) Log file inspection (this is the most important part): `~/logs/learn_read_orientation_model.log`

- `IOUtils - Extracting data from archive: file:/variants/SRR30536566.f1r2.tar.gz`

- Sample name consistency

```bash
./DMBEL-EIDR-071.alt_histogram
./DMBEL-EIDR-071.ref_histogram
```
**Interpretation**
✔️ Sample name matches:

BAM @RG SM:DMBEL-EIDR-071

d) EM convergence (KEY QUALITY INDICATOR)

```bash
Context AAC: with 2316 ref and 459 alt examples, EM converged in 8 steps
Context ACG: with 301 ref and 142 alt examples, EM converged in 15 steps
Context AAA: with 7018 ref and 1195 alt examples, EM converged in 8 steps
...
```
✔️ Excellent signs:

- Hundreds to thousands of examples per context

- EM algorithm converged normally (5–15 iterations)

- No “failed to converge” warnings

**Interpretation**: Enough data for modeling; suitable for cfDNA + amplicon data

e) Orientation model contents

```bash
>SAMPLE=DMBEL-EIDR-071
context rev_comp f1r2_* f2r1_* ...
```

```bash
zless SRR30536566.read-orientation-model.tar.gz | head -n4
./DMBEL-EIDR-071.orientation_priors0100644 0000765 0000024 00000035607 15132156174 022020 0ustar00Franostaff0000000 0000000 15132156174 15132156174 #<METADATA>SAMPLE=DMBEL-EIDR-071
context	rev_comp	f1r2_a	f1r2_c	f1r2_g	f1r2_t	f2r1_a	f2r1_c	f2r1_g	f2r1_t	hom_ref	germline_het	somatic_het	hom_var	num_examples	num_alt_examples
ATT	AAT	2.761677961939111E-4	2.122446680952301E-5	2.5886223424394114E-5	0.0	2.8843202841039807E-5	2.650355750455237E-4	2.096714156143896E-5	0.0	0.999028937511455	2.205326254098435E-5	2.9108978973800586E-4	1.9795030390106504E-5	5051	761
CTT	AAG	4.629726714038336E-4	3.596216200562292E-5	3.370755788104045E-4	0.0	3.586133538480483E-4	5.915063999029915E-4	3.057171739011366E-5	0.0	0.9970301940849805	4.017032570311504E-5	5.413645457595304E-4	5.715691601958141E-4	3674	537
...
```
**Interpretation**:  👍
✔️ Correct metadata
✔️ Correct probability tables
The ^@ characters are normal binary padding from the tar format. They are not corruption.

>**IMPORTANT**: The model in file `SRR30536566.read-orientation-model.tar.gz` must later on be passed to `FilterMutectCalls`


### GetPileupSummaries 👉 [06a_get_pileup_summaries.sh](bash_scripts/06a_get_pileup_summaries.sh)

**Documentation**

- **GetPileupSummaries**: <https://gatk.broadinstitute.org/hc/en-us/articles/27007916224539-GetPileupSummaries>

Tabulates pileup metrics for inferring contamination.

`GetPileupSummaries` does **NOT** filter variants itself — it informs the next step (`CalculateContamination`)

**Outputs**

a) From Terminal

```bash
Tool returned:
SUCCESS
```

b) **Output files**: `~/variants/SRR30536566.pileups.table` (3.9 KB) and `~/logs/get_pileup_summaries.log`

c) Sample identity consistency

From `SRR30536566.pileups.table`:

```bash
#<METADATA>SAMPLE=DMBEL-EIDR-071
```
✔️ Confirms the **biological sample name (SM)**

✔️ Must match BAM `@RG SM` and Mutect2 `--tumor-sample`

This avoids one of the most common GATK contamination failures.


d) Read filtering statistics (NORMAL)

From `get_pileup_summaries.log`:

```bash
901738 total reads filtered out of 1464279 reads processed
```
This is expected and healthy for:

- Amplicon data

- MarkDuplicates enabled

- High mapping quality thresholds

```bash
797528 read(s) filtered by: NotDuplicateReadFilter 
```
✔️ Confirms duplicates were correctly marked earlier

✔️ Confirms Mutect2 and downstream tools are respecting duplicate flags

```bash
Processed 88955 total loci in 0.2 minutes.
```
✔️ Exactly what we expect for:

- gnomAD common sites

- restricted to ~573 kb panel space

This is a good number for robust contamination estimation.

e) `SRR30536566.pileups.table` content:

It is a **quality-control / modeling** table used **only** to estimate cross-sample contamination in tumor-only Mutect2 analyses.

GATK is asking:

>“At sites where the population allele frequency is known, does the observed data look like it came from one individual, or does it look like a mixture of multiple individuals?”

This table contains summary counts at common germline SNPs (from gnomAD), restricted to your panel intervals.

```bash
contig	position	ref_count	alt_count	other_alt_count	allele_frequency
chr1	114705427	809	0	0	0.034
```
As table:

|contig |	position  |	ref_count 	| alt_count	| other_alt_count 	| allele_frequency |
|-------|-----------|-------------|-----------|-------------------|------------------|
|chr1	  | 114705427 |	809         |	0         | 0             	  | 0.034            | 

```bash
contig  position  ref_count  alt_count  other_alt_count  allele_frequency
```
These columns are **per-SNP pileup summaries**.

**Interpretation**:

🔹 contig, position
```bash
chr1	114705427	809	0	0	0.034
```
This is a **known common SNP site** from gnomAD
It lies inside the 7-gene CRC panel

🔹 allele_frequency
```bash
0.034
```
This is the population allele frequency (from gnomAD):

> ~3.4% of people carry the ALT allele at this position

This is **not** your sample’s VAF (Variant Allele Frequency).

🔹 ref_count
```bash
809
```
Number of reads supporting the reference allele (after all filters: mapping quality, duplicates, etc.)

🔹 alt_count
```bash
0
```
Number of reads supporting the **ALT allele** (the gnomAD SNP allele)

🔹 other_alt_count
```bash
0
```
Reads supporting neither REF nor expected ALT (usually sequencing noise or errors)

At a SNP with population AF = 3.4%:

- ~96.6% of individuals are **homozygous REF**

- ~3.4% are **heterozygous**

Homozygous ALT is extremely rare

So in a **single uncontaminated individual**, most sites will show:

```bash
ref_count ≫ alt_count
```
That is **normal, healthy, and expected**. This is considered a **good signal**.

✅ The observed data matches the  **single-individual model**

Specifically:

- High coverage

- Very few unexpected ALT reads

- No systematic inflation of ALT counts

This tells GATK:

“This locus behaves like it comes from one genome, not a mixture.”

**What would bad signal look like?**. Examples of problematic pileups:

🚨 **Contamination**
```bash
ref_count = 500
alt_count = 40
allele_frequency = 0.03
```
ALT fraction (~7%) ≫ expected (~1–3%)

➡ suggests DNA from **another individual** mixed in

🚨 **Severe technical artifacts**
```bash
other_alt_count = 50
```

Why there are low ALT counts?

Even in tumor tissue, these effects still apply:

🔹 1. Population SNPs are mostly absent

Most gnomAD SNPs are not present in your patient.

🔹 2. Tumor purity < 100%

Normal stromal cells dilute tumor DNA.

🔹 3. Somatic events ≠ gnomAD SNPs

This table is not looking for tumor mutations.

🔹 4. Sequencing noise

Occasional ALT reads are normal.


**More interpretations**

Somatic variants are intentionally ignored in this step.

`GetPileupSummaries`:

- Only uses known germline SNPs

- Uses population AF (allele frequency)

- Does not attempt to detect cancer mutations

So:

- KRAS mutations?

- PIK3CA mutations?

- BRAF mutations?

➡ **Irrelevant for this table**


### Contamination estimation 👉 [06b_calculate_contamination.sh](bash_scripts/06b_calculate_contamination.sh)

**Documentation**

- **CalculateContamination**: <https://gatk.broadinstitute.org/hc/en-us/articles/360036888972-CalculateContamination>

Given pileup data from `GetPileupSummaries`, calculates the fraction of reads coming from cross-sample contamination. The resulting contamination table is used with `FilterMutectCalls`. 
This tool estimates contamination based on the signal from ref reads at hom alt sites.

**Outputs**

a) From Terminal:

```bash
Tool returned:
SUCCESS
```

b) Output files: 

`~variants/SRR30536566.contamination.table`

`~logs/calculate_contamination.log`

c) ⚠️ About the warnings in the log

```bash
WARN  KernelSegmenter - ... number of data points (1)
WARN  KernelSegmenter - No changepoint candidates were found
```
✔ These are expected and harmless
✔ They occur because you did not provide segmentation (by design)
✔ GATK attempted segmentation, realized there is only one data point, and correctly fell back to a single global estimate

👉 Nothing is wrong here

> **Note**: Tumor segmentation is useful when:
>
>Whole exome (WES)
>
>Whole genome (WGS)
>
>Many thousands of loci
>
>Copy-number variation affects allele fractions
>
>Your data:
>
>7-gene amplicon panel
>
>No meaningful CNV inference possible

d) **IMPORTNAT**: `SRR30536566.contamination.table`

```bash
sample              contamination          error
DMBEL-EIDR-071       0.016233858827731117   0.027470745560259448
```
Meaning of each column

| Column          | Meaning                                  |
| --------------- | ---------------------------------------- |
| `sample`        | Biological sample name (from `@RG:SM`)   |
| `contamination` | Estimated fraction of contaminating DNA  |
| `error`         | Statistical uncertainty of that estimate |

**Interpreting contamination** = 0.0162 ≈ 1.62%

This means:

**~1.6% of reads behave like they come from another genome**

This “**other genome**” is typically:

- Normal DNA (blood, stromal cells)

- Cross-sample contamination

- Index hopping (minor contributor)

**Is 1.6% contamination good or bad?**

Context-dependent answer:

| Context           | Interpretation                      |
| ----------------- | ----------------------------------- |
| WGS / WES         | Very good                           |
| Tumor-only panels | Acceptable                          |
| cfDNA             | Normal                              |
| Amplicon panels   | Slightly on the high side, but fine |

🔬 For **targeted cancer panels**, anything below:

- 2–3% → acceptable

- <1% → excellent

- >5% → concerning

👉 Your value is acceptable and usable

What does the `**error = 0.027**` mean?

It means:

- There is high uncertainty

- Due to limited loci

- Expected for:

  - Small panels

  - Few gnomAD SNPs

  - Targeted intervals

Important:

⚠️ Do **NOT** interpret the error literally as ±2.7%

Instead:

- It reflects low statistical power

GATK is warning you:

“I estimated contamination, but with limited confidence”

This is normal for 7 genes.


### Variant Filtering 👉 [07_filter_mutect_calls.sh](bash_scripts/07_filter_mutect_calls.sh)

**Documentation**

- **FilterMutectCalls**: <https://gatk.broadinstitute.org/hc/en-us/articles/360036888972-CalculateContamination>

Filter somatic SNVs and indels called by Mutect2.  To filter based on sequence context artifacts, specify the `--orientation-bias-artifact-priors`. This input is generated by `LearnReadOrientationModel`. 
If given a `--contamination-table` file, e.g. results from `CalculateContamination`, the tool will additionally filter variants due to contamination. This argument may be specified with a table for one or more tumor samples. 
`FilterMutectCalls` can also be given one or more `--tumor-segmentation` files, which are also output by `CalculateContamination`. 

**Outputs**

`~/variants/SRR30536566.filtered.vcf.gz`  47 KB

`~/variants/SRR30536566.filtered.vcf.gz.filteringStats.tsv`

`~/variants/SRR30536566.filtered.vcf.gz.tbi`

`~/logs/filter_mutect_calls.log`

1️⃣ What FilterMutectCalls actually did

From the log:
```bash
Processed 948 total variants
```
**Interpretation**

- Mutect2 found 948 candidate variants in your 7-gene CRC panel

- FilterMutectCalls evaluated each allele using:

  - strand bias

  - orientation bias

  - contamination

  - germline probability

  - weak evidence model

  - PON

  - base quality

  - haplotype consistency


2️⃣ Why almost everything is FILTERED (important!)

You will notice something striking: **Almost no variants are PASS**. This is expected and biologically correct for your dataset.

3️⃣ Understanding the filteringStats.tsv (very important)

This file answers:
“Which filters removed variants, and how strong are they?”

Key lines:
```bash
filter          FN     FNR
weak_evidence   1.32   0.24
contamination   1.32   0.24
orientation     1.39   0.26
germline        1.33   0.25
```
FN: false negatives

**Interpretation**

| Filter        | Meaning                                                  |
| ------------- | -------------------------------------------------------- |
| weak_evidence | Alt allele supported by too few reads given depth        |
| contamination | Alt allele frequency consistent with ~1.6% contamination |
| orientation   | F1R2 / F2R1 strand artifact                              |
| germline      | Looks like constitutional variant                        |

4️⃣ Let’s interpret real variants from your VCF (`SRR30536566.filtered.vcf.gz`)

Example 1 — High-AF “germline” variant (AF: allele fraction)
```bash
chr1 114705278 A>G FILTER=germline
AF=0.487
DP=1013
TLOD=1676.10
```
**Interpretation**

- ~50% allele fraction

- Very strong signal

- Extremely **high TLOD** (Tumor Log Odds score in somatic variant calling: It is a statistical score that represents the confidence that a detected variant signal in a tumor sample is a true somatic mutation rather than just background noise or a sequencing artifact.Variants with a TLOD score below a certain threshold are typically filtered out, as they are considered to have insufficient evidence of presence in the tumor.)

Example 2 — Low-AF somatic-looking but filtered
```bash
chr12 25227460 T>A
AF=0.020
DP=340
FILTER=base_qual;contamination;orientation;strand_bias;weak_evidence
```
**Interpretation**

- AF ~2%

- Only 7 alt reads

- Strong strand imbalance

- Low alt base quality

- Compatible with contamination + artifacts

➡️ Not reliable enough
Correctly filtered.

Example 3 — PON + germline + STR
```bash
chr15 66387311 T>TC
AF=0.549
FILTER=germline;panel_of_normals;STR
```
This has a high AF

➡️ Definitely not somatic

5️⃣ Why you see MANY germline calls

This is expected because:

- You used tumor-only mode

- There is no matched normal

- CRC tissue still contains germline DNA

Mutect2 therefore:

- Calls germline variants

- Labels them as germline

- Removes them in filtering

✔ This is exactly correct behavior.

6️⃣ Do you have any PASS variants?

Possible outcomes:

**Case A — No PASS variants**

This is plausible if:

- No hotspot somatic mutations in these 7 genes

- Or tumor purity is low

- Or variants are below ~2–3% VAF

**Case B — Few PASS variants**

These are your high-confidence somatic mutations

```bash
zless ..//data/SRR30536566/variants/SRR30536566.filtered.vcf.gz | grep "PASS"
```
Output: There are 4 PASS variants.

```bash
##FILTER=<ID=PASS,Description="Site contains at least one allele that passes filters">
##filtering_status=These calls have been filtered by FilterMutectCalls to label false positives with a list of failed filters and true positives with PASS.
zless ..//data/SRR30536566/variants/SRR30536566.filtered.vcf.gz | grep "chr*" | grep "PASS"
chr1	114713909	.	G	T	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=328,320|61,54;DP=817;ECNT=2;ECNTH=1;GERMQ=93;MBQ=41,41;MFRL=158,156;MMQ=60,60;MPOS=24;POPAF=5.60;ROQ=93;TLOD=323.24	GT:AD:AF:DP:F1R2:F2R1:FAD:SB	0/1:648,115:0.154:763:241,52:289,43:567,103:328,320,61,54
chr3	179210338	.	AGTAAGGTTTTTATTGTCATAAATTAGATATTTTTTATGGCAGTCAAACCTTCTCTCTTATGTATATATAATAGCTTTTCTTCCATCTCTTAG	A	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=1161,1314|78,88;DP=2641;ECNT=1;ECNTH=1;GERMQ=93;MBQ=41,41;MFRL=171,203;MMQ=60,60;MPOS=18;POPAF=5.60;ROQ=93;TLOD=106.19	GT:AD:AF:DP:F1R2:F2R1:FAD:SB	0/1:2475,166:0.019:2641:397,35:411,27:1915,146:1161,1314,78,88
chr3	179218294	.	G	A	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=454,461|162,185;DP=1324;ECNT=1;ECNTH=1;GERMQ=93;MBQ=41,41;MFRL=170,174;MMQ=60,60;MPOS=24;POPAF=5.60;ROQ=93;TLOD=1026.18	GT:AD:AF:DP:F1R2:F2R1:FAD:SB	0/1:915,347:0.277:1262:363,138:382,148:786,300:454,461,162,185
chr3	179226113	.	C	G	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=23,142|63,331;DP=589;ECNT=2;ECNTH=1;GERMQ=93;MBQ=41,41;MFRL=186,175;MMQ=60,60;MPOS=21;POPAF=1.34;ROQ=93;TLOD=1387.91	GT:AD:AF:DP:F1R2:F2R1:FAD:SB	0/1:165,394:0.698:559:68,145:75,182:151,350:23,142,63,331
```

7️⃣ Is this a failure biologically?

❌ No
✅ This is a clean, conservative result

In fact:

- Over-filtering is preferred in clinical pipelines

- **False positives are worse** than false negatives


### Post-filter (amplicon-specific) 👉 [08_postfilter.sh](bash_scripts/08_postfilter.sh)

**Documentation**

- **bcftools**: <https://gatk.broadinstitute.org/hc/en-us/articles/360036888972-CalculateContamination>

**Output files**

`~/variants/SRR30536566.postfilter_summary.txt`

`~/variants/SRR30536566.postfiltered.vcf.gz`

`~/variants/SRR30536566.postfiltered.vcf.gz.csi`

`~/logs/SRR30536566.postfilter.log`

**Output important information**

`SRR30536566.postfilter.log`
```bash
Variants retained after post-filtering: 3
```

`SRR30536566.postfiltered.vcf.gz`
```bash
zless ..//data/SRR30536566/variants/SRR30536566.postfiltered.vcf.gz | grep "chr*" | grep "PASS"
chr1	114713909	.	G	T	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=328,320|61,54;DP=817;ECNT=2;ECNTH=1;GERMQ=93;MBQ=41,41;MFRL=158,156;MMQ=60,60;MPOS=24;POPAF=5.6;ROQ=93;TLOD=323.24	GT:AD:AF:DP:F1R2:F2R1:FAD:SB	0/1:648,115:0.154:763:241,52:289,43:567,103:328,320,61,54
chr3	179218294	.	G	A	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=454,461|162,185;DP=1324;ECNT=1;ECNTH=1;GERMQ=93;MBQ=41,41;MFRL=170,174;MMQ=60,60;MPOS=24;POPAF=5.6;ROQ=93;TLOD=1026.18	GT:AD:AF:DP:F1R2:F2R1:FAD:SB	0/1:915,347:0.277:1262:363,138:382,148:786,300:454,461,162,185
chr3	179226113	.	C	G	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=23,142|63,331;DP=589;ECNT=2;ECNTH=1;GERMQ=93;MBQ=41,41;MFRL=186,175;MMQ=60,60;MPOS=21;POPAF=1.34;ROQ=93;TLOD=1387.91	GT:AD:AF:DP:F1R2:F2R1:FAD:SB	0/1:165,394:0.698:559:68,145:75,182:151,350:23,142,63,331
```

**Why 3 variants is the “right” answer**

You started with 4 PASS calls from `FilterMutectCalls`.

Your post-filter thresholds **removed** 1 borderline subclonal variant:

| Variant        | AF        | Result                |
| -------------- | --------- | --------------------- |
| chr3:179210338 | **0.019** | ❌ dropped (AF < 0.02) |


### Biological interpretation of the 3 retained variants

✅ chr1:114713909 (AF 15%)

```bash
DP=763
AD=648,115
AF=0.154
TLOD=323
```
**Interpretation**

- Strong allelic signal

- Balanced strand counts

- High TLOD

- Likely clonal or early driver

✅ chr3:179218294 (AF 27%)

```bash
DP=1262
AD=915,347
AF=0.277
TLOD=1026
```
**Interpretation**

- Very strong variant

- Likely founder clone

- Excellent depth and support

- Virtually impossible to be an artifact

✅ chr3:179226113 (AF 70%)
```bash
DP=559
AD=165,394
AF=0.698
TLOD=1387
```
**Interpretation**

- Near-homozygous or LOH-associated

- Possibly:

  - copy-number loss of REF allele

  - strong clonal expansion

- Extremely strong confidence


### Post-Filtering Results Summary

Your post-filtering successfully removed 1 variant (chr3:179210338) while keeping 3 variants:

    chr1:114713909 - G>T (VAF=15.4%, DP=763)

    chr3:179218294 - G>A (VAF=27.7%, DP=1262)

    chr3:179226113 - C>G (VAF=69.8%, DP=559)

The removed variant failed because its VAF (1.9%) was below your 2% threshold.

**VCF Format Field Explanations**
**GT:AD:AF:DP:F1R2:F2R1:FAD:SB**

This is the FORMAT field describing how the variant was called in your sample:

    GT (Genotype): 0/1 = heterozygous variant

    AD (Allele Depth): Reads supporting reference and alternate alleles

        Format: REF,ALT (e.g., 648,115 = 648 ref reads, 115 alt reads)

    AF (Allele Frequency): VAF = ALT/(REF+ALT)

    DP (Depth): Total reads at this position

    F1R2 & F2R1: Strand bias metrics

        F1R2: Forward reads supporting ref/alt (read1 in forward orientation)

        F2R1: Reverse reads supporting ref/alt (read2 in reverse orientation)

        Used to detect PCR artifacts/strand bias

    FAD (Filtered Allele Depth): AD after filtering low-quality reads

    SB (Strand Bias Table): ref_forward,ref_reverse,alt_forward,alt_reverse

**TLOD (Tumor LOD Score)**

    What it is: Log-odds score comparing variant vs. no variant hypotheses

    Interpretation:

        Higher TLOD = stronger evidence for variant

        TLOD > 6.3 is typical Mutect2 threshold for somatic variants

        Your variants have very high TLODs (323-1387), indicating excellent confidence

    Formula: TLOD = log₁₀[P(data|variant)/P(data|no variant)]

**Clinical Relevance for Colorectal Cancer**

Let me map these to your targeted genes:

    chr1:114713909 = NRAS codon 61 (likely Q61K/L/R)

        NRAS mutations in CRC (~3-5% of cases)

        Predictive: Anti-EGFR resistance (similar to KRAS)

    chr3:179218294 = PIK3CA exon 9 (likely E545K)

        Common in CRC (~15-20%)

        Associated with poor prognosis

        Emerging therapeutic target (PI3K inhibitors)

    chr3:179226113 = PIK3CA exon 20 (likely H1047R)

        Most common PIK3CA mutation in CRC

        Constitutively activates PI3K pathway

        High VAF (69.8%) suggests clonal/dominant mutation

**Quality Assessment**

All 3 variants look high-quality:

    ✅ High depth (>500x)

    ✅ Good VAF (15-70%)

    ✅ Strong strand balance (not strand-biased)

    ✅ Excellent TLOD scores

    ✅ Passed all Mutect2 filters

**Key Implications for Your Patient**

    NRAS mutation = Likely resistance to anti-EGFR therapies (cetuximab/panitumumab)

    Dual PIK3CA mutations = Strong PI3K pathway activation

        Consider PI3K/mTOR inhibitors in clinical trials

        Associated with poorer outcomes

    No KRAS/BRAF mutations detected = May still benefit from EGFR inhibitors if NRAS is wild-type (but you have NRAS mutation)


### Folder structure: after running post-filtering variant.

```bash
Genomics_cancer/
├── reference/                 
│   └── GRCh38/
│       ├── fasta/
│       └── known_sites/       
│       └── intervals/                                    
│       └── somatic_resources/                                   
├── data/
│   └── SRR30536566/                
│       ├── raw_fastq/
│       ├── qc/
│       ├── trimmed/
│       ├── aligned/
│       ├── variants/                                       
│           └── SRR30536566.f1r2.tar.gz          
│           └── SRR30536566.unfiltered.vcf.gz
│           └── SRR30536566.unfiltered.vcf.gz.stats 
│           └── SRR30536566.unfiltered.vcf.gz.tbi
│           └── SRR30536566.read-orientation-model.tar.gz
│           └── SRR30536566.pileups.table
│           └── SRR30536566.contamination.table
│           └── SRR30536566.filtered.vcf.gz.filteringStats.tsv
│           └── SRR30536566.filtered.vcf.gz
│           └── SRR30536566.filtered.vcf.gz.tbi
│           └── SRR30536566.postfiltered.vcf.gz
│           └── SRR30536566.postfiltered.vcf.gz.csi
│           └── SRR30536566.postfilter_summary.txt
│       └── annotation/        
├── scripts/
│       └── 0_wget_gnomad_PoN.sh
│       └── 0_wget_Hsapiens_assem38.sh
│       └── 01_qc.sh
│       └── 02_trim.sh
│       └── 03_align_&_bam_preprocess.sh
│       └── 04_make_crc_7genes_bed.sh
│       └── 04_mutect2.sh
│       └── 05_learn_read_orientation_model.sh          
│       └── 06a_get_pileup_summaries.sh            
│       └── 06b_calculate_contamination.sh             
│       └── 07_filter_mutect_calls.sh             
│       └── 08_postfilter.sh                                              
└── logs/
        └── cutadapt_SRR30536566.log
        └── bwa_mem.log
        └── markduplicates.log
        └── SRR30536566.flagstat.txt                     
        └── mutect2.stderr.log                                  
        └── mutect2.stdout.log                                  
        └── learn_read_orientation_model.log                                  
        └── get_pileup_summaries.log                                 
        └── calculate_contamination.log                                 
        └── filter_mutect_calls.log                             
        └── SRR30536566.postfilter.log                                 
