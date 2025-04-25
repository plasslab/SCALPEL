<<<<<<< HEAD
# SCALPEL: A Nextflow-based Pipeline for Isoform Quantification at Single-Cell Resolution

<div align="center">
  <img src="https://data.cyverse.org/dav-anon/iplant/home/franzx5/SPERMATOGENESIS/SCALPEL_pipeline.png" alt="SCALPEL" width="800">
=======

# SCALPEL , a nextflow based tool for the quantification of isoforms at single-cell resolution

<!-- PROJECT LOGO -->
<!-- <br />
<div align="center">
  <a href="">
    <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/b/ba/Dessin_scalpel.svg/1200px-Dessin_scalpel.svg.png" alt="SCALPEL" width="300" height="300">
  </a>
</div>-->

<div align="right">
  <a href="">
    <img src="https://data.cyverse.org/dav-anon/iplant/home/franzx5/TUTO/PLOT0.png" alt="SCALPEL" >
  </a>
>>>>>>> main
</div>

## About the Project

**SCALPEL** is a robust pipeline designed for transcript isoform quantification and alternative polyadenylation (APA) characterization using 3'-tagged single-cell RNA-seq (scRNA-seq) data. Built with **Nextflow**, it integrates multiple processing steps including read quantification, APA annotation, and isoform usage analysis.

## Installation Options

SCALPEL can be installed and run using one of the following options:

### Option 1: Using Conda

<<<<<<< HEAD
1. Clone the repository
```bash
git clone https://github.com/p-CMRC-LAB/SCALPEL.git
cd SCALPEL
```

2. Create the Conda environment
```bash
conda env create -f requirements.yml
conda activate scalpel_conda
```

3. Run SCALPEL
```bash
nextflow run -resume main.nf \
  --sequencing chromium \
  --samplesheet path/to/samplesheet.csv \
  --transcriptome path/to/gencode.transcripts.fa \
  --gtf path/to/gencode.annotation.gtf \
  --ipdb path/to/mm10.polyA.track \
  --barcodes path/to/barcodes.csv \
  --clusters path/to/clusters.txt
```

---

### Option 2: Using Apptainer (Recommended for reproducibility)

You can download a prebuilt Apptainer container with all SCALPEL dependencies from the following link:  
[Download SCALPEL Container](https://data.cyverse.org/dav-anon/iplant/home/franzx5/SCALPEL.container.sif)

1. Download the container and clone the repository
```bash
wget https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/scalpel_container.sif](https://data.cyverse.org/dav-anon/iplant/home/franzx5/SCALPEL.container.sif
git clone https://github.com/p-CMRC-LAB/SCALPEL.git
cd SCALPEL
```

2. Run SCALPEL using the container
```bash
nextflow run -resume main.nf \
  -with-apptainer ./scalpel_container.sif \
  --sequencing chromium \
  --samplesheet path/to/samplesheet.csv \
  --transcriptome path/to/gencode.transcripts.fa \
  --gtf path/to/gencode.annotation.gtf \
  --ipdb path/to/mm10.polyA.track
```

---

## Required Input Files

| Parameter        | Description                                                               |
|------------------|---------------------------------------------------------------------------|
| `--samplesheet`  | CSV with sample name and paths to FASTQ/BAM/CellRanger outputs            |
| `--transcriptome`| FASTA of reference transcriptome                                          |
| `--gtf`          | GTF annotation file                                                       |
| `--ipdb`         | Internal priming annotation file                                          |
| `--barcodes`     | (Optional) Barcode whitelist per sample                                   |
| `--clusters`     | (Optional) Tab-delimited file with cell-to-cluster mappings               |
| `--sequencing`   | Must be `chromium` or `dropseq`                                           |

Reference files:

- [Mouse IP Annotation - mm10](https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/databases/mm10_polya.track.tar.gz)
- [Human IP Annotation - GRCh38](https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/databases/GRCh38_2020_A_polyA.track.tar.gz)
- [GENCODE Annotation (Mouse Example)](https://www.gencodegenes.org/mouse/)

---

## Output Files

After successful execution, SCALPEL generates a `results/` directory containing final outputs and intermediate files structured for downstream analysis.

### Overview of Key Output Files

| Filename pattern                       | Description |
|----------------------------------------|-------------|
| `*_filtered.bam`                       | Indexed BAM files with deduplicated reads that are not associated with internal priming positions. These are the reads used for final isoform quantification. |
| `*_filtered.bam.bai`                   | BAM index files corresponding to the filtered BAMs. |
| `*_APADGE.txt`                         | Isoform-level expression count matrix for each sample, reflecting APA-aware quantification. |
| `*_seurat.RDS`                         | Sample-specific Seurat object containing metadata, isoform expression data, and clustering results if applicable. |
| `iDGE_seurat.RDS`                      | Merged Seurat object across all samples for integrated downstream analysis of isoform-based expression. |
| `DIU_table.csv`                        | Table of Differential Isoform Usage (DIU) comparing conditions or clusters. Includes p-values, fold changes, and significance annotations. |
| `Runfiles/`                            | Execution metadata, logs, and intermediate outputs from Nextflow processes. May be useful for debugging or auditing. |

### Notes

- Output filenames are automatically prefixed with the corresponding sample name or labeled appropriately for merged/integrated results.
- Seurat objects (`.RDS`) can be loaded directly into R for visualization and further downstream analysis.
- Count matrices (`*_APADGE.txt`) can be used with external tools for statistical modeling, clustering, or plotting.
- Differential usage results (`DIU_table.csv`) highlight isoform switching events between groups or conditions of interest.

Make sure to clean up large temporary folders such as `work/` after completion, unless you plan to resume or reprocess.

---
=======
## Installation

1. Clone the repo and enter in the folder
```sh
git clone https://github.com/p-CMRC-LAB/SCALPEL.git
```
2. Install the required packages using the requirement.txt file in the SCALPEL folder
```sh
CONDA_CHANNEL_PRIORITY=flexible conda env create -f SCALPEL/requirements.yml
conda activate scalpelEnv
```
   
Another solution (if conda installation takes long) can be to create a Conda environment, install Mamba (faster implementation of Conda) and install the packages using mamba:
```sh
mamba env create --file SCALPEL/requirements.yml
mamba activate scalpelEnv
```

## SCALPEL usage

### Input files
For running, SCALPEL requires to provide specific input files path & parameters:
-   _SAMPLE path_  **[--samplesheet]**
-   _FASTA transcriptome reference path_  **[--transcriptome]**
-   _GTF annotation file path_  **[--gtf]**
-   _Internal priming annotation file for the organism_  **[--ipdb]**  (download link below)
-   _Sequencing type (chromium or dropseq)_  **[--sequencing]**

  1.  Required, **[\-\-samplesheet]**: Provide within a  **CSV**  (ex: samplesheet.csv) file the following paths : \
    -(In case of 10X based scRNA-seq sample [--sequencing  **chromium**]:
```sh
<SAMPLE_CELLRANGER_REPOSITORY_NAME>,<FASTQ1_FILE_PATH>,<FASTQ2_FILE_PATH>,<10X_CELLRANGER_REPOSITORY_PATH>
```
or DropSeq based scRNA-seq sample [--sequencing  **dropseq**]):
```sh
<SAMPLE_NAME>,<FASTQ1_PATH>,<FASTQ2_PATH>,<BAM_PATH>,<BAM_INDEX_PATH>,<DGE_PATH>
```
  2. Optional, **[\-\-barcodes]**: Provide within a **CSV** (ex: barcodes_whitelist.csv) for each input sample, a barcode whitelist file path:
```sh
<SAMPLE_NAME>,<BARCODE_WHITELIST_FILE_PATH>
```
  3. (Optional, **[\-\-clusters]**: To perform a differential analysis on defined cells cluster, provide within a tab-separated file (ex: clusters.txt) containing the cells annotation information. This file should contain 3 column associated to:
```sh
<SAMPLE_NAME>,<BARCODE_CELL_TAGS>,<CLUSTER_ANNOTATION>
```
>>>>>>> main

## Documentation and Examples

- [SCALPEL application on 10X scRNA-seq](https://github.com/p-CMRC-LAB/SCALPEL/wiki/SCALPEL-application-on-10X-scRNA%E2%80%90seq)
- [SCALPEL application on DropSeq scRNA-seq](https://github.com/p-CMRC-LAB/SCALPEL/wiki/SCALPEL-application-on-DropSeq-scRNA%E2%80%90seq)
- [Downstream analysis](https://github.com/p-CMRC-LAB/SCALPEL/wiki)

<<<<<<< HEAD
---
=======
- **Internal priming files which reference all the internal priming positions**
- [(Human) Internal priming annotation - GRCh38](https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/databases/GRCh38_2020_A_polyA.track.tar.gz) 
- [(Mouse) Internal priming annotation - mm10](https://data.cyverse.org/dav-anon/iplant/home/franzx5/Scalpel_docs/databases/mm10_polya.track.tar.gz) 


## SCALPEL execution
After activating the _scalpel_conda_ CONDA environment, SCALPEL can be executed using Nextflow: \
You can print the Scalpel help documentation by running the following command:
```
> nextflow run -resume SCALPEL/main.nf --help

	SCALPEL - NF  P I P E L I N E
	===============================

	Execution:
	Ex: nextflow run -resume scalpel.nf --sequencing <Sequencing type>
	 --samples <SAMPLE files folder path>
	 --reads <FASTQs files folder path> --transcriptome <FASTA transcriptome reference path>
	 --annot <GTF annotation file path> --ipdb <Internal priming annotation file> 
	
	Input files:
    - Annotation required files(required):
        - transcriptome reference [--transcriptome]
        - annotation GTF reference [--gtf]
        - internal priming annotation [--ipdb]
      
    - Reads processing files (required):
        - samples files [--samples]
        - fastqs files [--reads]
    
    - Params:
        Required:
        - sequencing type (required): ${params.sequencing}

        Optional:
        - transcriptomic distance threshold [--dt_threshold] (optional, default 600bp)
        - transcriptomic end distance threhsold [--dt_exon_end_threshold] (optional, default 30bp)
        - minimal distance of Ip from isoform 3'ends (optional, default 60bp)
        - params.threads [--threads] (default 30)
        - params.cpus [--cpus] (default 30)
```
- All the computational ressource required for the execution by Nextflow can be defined within the _**SCALPEL/nextflow.config**_ file:
```
/* Define Nextflow configuration settings for SCALPEL pipeline execution */
/* ===================================================================== */
>>>>>>> main

## Customizing Execution with `nextflow.config`

To modify resource usage and process settings, edit the `nextflow.config` file. For example:

```groovy
executor {
    name = 'slurm'               // Use 'local', 'slurm', etc.
    cpus = 64
}

process {
    withLabel: big_mem {
        cpus = 4
        memory = '8 GB'
    }
    withLabel: small_mem {
        cpus = 2
        memory = '2 GB'
    }
    // Additional process-specific settings...
}
```

<<<<<<< HEAD
### Apptainer Users: Configure `runOptions`

If using Apptainer, make sure to bind your local SCALPEL repository path inside the container by editing the following block in `nextflow.config`:

```groovy
apptainer {
    enabled = true
    autoMounts = true
    runOptions = "--bind /path/to/SCALPEL:/path/to/SCALPEL"
}
```

Adjust `/path/to/SCALPEL` to the full absolute path where the SCALPEL repository is located on your system.

---
=======
- **Execution**

1. Configurate a sample file **[--samplesheet]**:
```sh
> cat samplesheet.csv
SRR6129050,SRR6129050_S1_L001_R1_001.fastq.gz,SRR6129050_S1_L001_R2_001.fastq.gz,/data/fake_data/DATAS/GSE104556/SAMPLES/SRR6129050/
SRR6129051,SRR6129051_S1_L001_R1_001.fastq.gz,SRR6129051_S1_L001_R2_001.fastq.gz,/data/fake_data/DATAS/GSE104556/SAMPLES/SRR6129051/
```

2. (**Optional**), Configurate a barcode whitelist file **[--barcodes]**:
```sh
> cat barcodes_whitelist.csv
SRR6129050,/data/fake_data/NEW/SCALPEL/10X/SRR6129050_curatedBarcodes.txt
SRR6129051,/data/fake_data/NEW/SCALPEL/10X/SRR6129051_curatedBarcodes.txt

> head /data/fake_data/NEW/SCALPEL/10X/SRR6129050_curatedBarcodes.txt
AAACCTGAGCTTATCG-1
AAACCTGGTTGAGTTC-1
AAACCTGTCAACGAAA-1
AAACGGGCACAGGTTT-1
AAACGGGTCATTTGGG-1
...
```

3. (**Optional**), Configurate a cluster annotation file **[--barcodes]**:
```sh
> cat clusters.txt
SRR6129050  AAACCTGAGCTTATCG-1  RS1
SRR6129050  AAACCTGGTTGAGTTC-1  RS1
SRR6129050  AAACCTGTCAACGAAA-1  ES
SRR6129051  AAACGGGCACAGGTTT-1  RS1
SRR6129051  AAACGGGTCATTTGGG-1  ES
...
```

3. **Running**
```sh
nextflow run -resume SCALPEL/main.nf \
    --sequencing chromium \
    --samplesheet samplesheet.csv \
    --transcriptome gencode.vM10.transcripts.fa \
    --gtf gencode.vM21.annotation.gtf \
    --ipdb mm10.polyA.track \
    --barcodes barcodes_whitelist.csv \ (Optional)
    --clusters clusters.txt \ (Optional)
```


(**Optional**), In case of memory error issue in the execution of SCALPEL using SLURM executor on HPC, it can be run within an interactive session using srun with a 'local' executor:
```sh
srun --pty --mincpus XX --mem XX bash
```




## Results

During its execution, **SCALPEL** shows interactively in the Console all the information about the executed tasks:
Then, a _**./results**_ folder is generated encompassing all the final and intermediated files generated by **SCALPEL** that can be used for further downstream analysis:

_**./results/final_results**_ contains:
  - the final differential isoform usage table between the input samples
  - the Seurat object encompassing the iDGE of the input samples
  - the iDGE table of each input sample

```sh
> cat ./results/final_results
DIU_table.csv  iDGE_seurat.RDS.  SRR6129050/SRR6129050_APA_DGE.txt   SRR6129051/SRR6129051_APA_DGE.txt
```

Information about the collapsed isoforms by SCALPEL (see paper) can be found in _**./results/annotation_processing/isoform_processing**_

```sh
> cat ./results/annotation_processing/isoform_processing
chr10_collapsed_isoforms.txt  chr13_collapsed_isoforms.txt  chr16_collapsed_isoforms.txt ...
```

**NB:**
_**Important: The chromosome names must to be consistent between the BAM file and the annotation files. If the BAM file contains the '_chr_' character, the GTF and FASTA annotation files, and the internal priming reference annotation should contains the '_chr_' character, and inversely !**_

_**The internal priming reference files provided contains by default the '_chr_' character in the chromosome names !!**_

_**Be careful to delete the _work_ directory containing nextflow temporary files** when scalpel runs all its processs sucessfully and you don’t plan to relaunch scalpel with modified parameters. (This folder can fill an high memory physical space depending of the size of input files analyzed)_


## More about SCALPEL Usage

- [Example of SCALPEL application on 10X scRNA-seq](https://github.com/p-CMRC-LAB/SCALPEL/wiki/SCALPEL-application-on-10X-scRNA%E2%80%90seq)
- [Example of SCALPEL application on DropSeq scRNA-seq](https://github.com/p-CMRC-LAB/SCALPEL/wiki/SCALPEL-application-on-DropSeq-scRNA%E2%80%90seq)


## Downstream analysis

See [SCALPEL Wiki](https://github.com/p-CMRC-LAB/SCALPEL/wiki)


## Contact
>>>>>>> main

## Reference

Franz Ake, Sandra M. Fernández-Moya, Marcel Schilling, Akshay Jaya Ganesh, Ana Gutiérrez-Franco, Lei Li, Mireya Plass  
**Quantification of transcript isoforms at the single-cell level using SCALPEL**  
bioRxiv 2024.06.21.600022; doi: [10.1101/2024.06.21.600022](https://doi.org/10.1101/2024.06.21.600022)

---

## Contact

<<<<<<< HEAD
Franz AKE – [@aerodx5](https://twitter.com/aerodx5) – fake@idibell.cat  / aerod7710@gmail.com

GitHub: [https://github.com/p-CMRC-LAB/SCALPEL](https://github.com/p-CMRC-LAB/SCALPEL)

<p align="right"><a href="#top">Back to top</a></p>
=======
## Paper
[Access SCALPEL_PAPER](https://www.biorxiv.org/content/10.1101/2024.06.21.600022v1)

Quantification of transcript isoforms at the single-cell level using SCALPEL \
Franz Ake, Sandra M. Fernández-Moya, Marcel Schilling, Akshay Jaya Ganesh, Ana Gutiérrez-Franco, Lei Li, Mireya Plass \
bioRxiv 2024.06.21.600022; doi: https://doi.org/10.1101/2024.06.21.600022

>>>>>>> main
