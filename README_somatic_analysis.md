**WORK IN PROGRESS**

#  Part II – Somatic analysis (Bash pipelines)

## Introduction

From [Part I – Preparation & setup](README_setup.md), the aim of this tutorial is to show how to analyse NGS data under restricted computational power and memory space.
Also, the tutorial is adapted to Macbook Pro laptops (Intel), therefore the tutorial should be adapted to other type operative systems.

### Types of DNA-NGS analysis

When finding out alterations of the genome, there are two types of DNA-NGS analysis:

1. **Somatic**: It focuses on mutations acquired by a tumor. These variants are not inherited and not present in normal tissue. Usually, somatic analyses are focused on

- Identifying actionable mutations for targeted therapies, e.g., PIK3CA, ERBB2, ESR1, AKT1 in breast cancer.

- Monitoring disease: Minimal residual disease (MRD) tracking by liquid biopsy follow-up, e.g., ctDNA.

- Clinical trials: Matching the patient to mutation-specific trials

- Understanding tumor biology: Subtyping the tumor more precisely, identifying driver vs passenger mutations and Detecting resistance mutations

2. **Germline**: It looks for inherited variants present in every cell of the body. These variants can predispose to cancer.

- Assessing hereditary cancer risk: Known familial mutation, strong family history, bilateral breast cancer, triple-negative breast cancer < 60 years

- Guiding treatment when germline mutations matter

- Informing family members

- Sample type: Blood or saliva (not tumor tissue)

**Somatic NGS Analysis**: For Understanding & Treating This Tumor
**Germline NGS Analysis**: For Assessing Inherited Risk

From one patient, both somatic (tumor) and germline (normal) DNA samples can be sequenced simultaneously and do somatic and germline analyses, respectively. This matched approach provides the most complete, accurate, and clinically actionable genomic portrait for a cancer patient and it is considered the best practice in modern genomic oncology, **though** cost and access can be factors.

### Bioinformatic DNA-NGS analysis

DNA-NGS analysis can be achieved by automation using these three main tools:
- Bash scripting
- Nextflow
- Snakemake
- Galaxy
- others (Cromwell/WDL, Toil, Airflow, Argo Workflows, Make, Ruffus, Bpipe)


**Table 1**: Advantages and disadvantages between 

| Tool | Advantages | Disadvantages |
|------|------------|---------------|
| **Bash Scripting** | • Maximum control & transparency<br>• Simple to start (no new syntax)<br>• Universally available | • Poor reproducibility & scalability<br>• Manual error handling & no built-in parallelization<br>• Difficult to maintain for complex pipelines |
| **Snakemake** | • Python-based, readable syntax<br>• Explicit file-based dependency graph<br>• Excellent for complex DAGs & Python integration | • Less built-in cloud portability<br>• "Pull" model (request final output) can be less intuitive<br>• Scaling requires profile configuration |
| **Nextflow** | • Extreme portability (write once, run anywhere)<br>• Implicit parallelism & powerful resume feature<br>• Strong community (nf-core) & production-ready | • Groovy DSL has a learning curve<br>• Debugging complex channels can be challenging<br>• Configuration can be multilayered |



- Simple, universal, no dependencies

- Perfect for small pipelines

- Transparent and easy to debug

- Works everywhere (HPC, cloud, local)


In this **PART II** tutorial, it will be used bash scripting to run somatic analysis of the small raw FASTQ dataset (read [Part I – Preparation & setup](README_setup.md)).

In Part I, the steps I to IV were achieved.

STEPS:

	I.	Create conda environment
	II.	Find & download small size FASTQ files
	III.	Find & download an indexed reference human genome
	IV.	Create folder structure (example)
	V.	Bash Scripting