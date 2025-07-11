# SCALPEL: A Nextflow-based Pipeline for Isoform Quantification at Single-Cell Resolution

<div align="center">
  <img src="docs/SCALPEL_diagram.png" alt="SCALPEL" width="800">
</div>

## About the Project

**SCALPEL** is a robust pipeline designed for transcript isoform quantification and alternative polyadenylation (APA) characterization using 3'-tagged single-cell RNA-seq (scRNA-seq) data. Built with **Nextflow**, it integrates multiple processing steps including read quantification, APA annotation, and isoform usage analysis.

## Installation Options

### Prerequisite:

**Nextflow>=v24.10.6** : [Official page](https://www.nextflow.io/docs/latest/install.html) / [CONDA](https://anaconda.org/bioconda/nextflow)


SCALPEL can be installed and run using one of the following options:

### Option 1: Using Conda environment

1. Clone the repository
```bash
> git clone https://github.com/p-CMRC-LAB/SCALPEL.git
```

2. Run SCALPEL using CONDA yml file
```bash
> nextflow run -resume SCALPEL/main.nf \
  --sequencing chromium \
  --samplesheet path/to/samplesheet.csv \
  --transcriptome path/to/gencode.transcripts.fa \
  --gtf path/to/gencode.annotation.gtf \
  --ipdb path/to/mm10.polyA.track \
  --cpus 40 \
  -with-conda SCALPEL/requirements.yml
```

or create CONDA environment and activate
```bash
> conda env create --file SCALPEL/requirements.yml
> conda activate scalpelEnv
> nextflow run -resume main.nf \
  --sequencing chromium \
  --samplesheet path/to/samplesheet.csv \
  --transcriptome path/to/gencode.transcripts.fa \
  --gtf path/to/gencode.annotation.gtf \
  --cpus 40 \
  --ipdb path/to/mm10.polyA.track
```

---

### Option 2: Using Apptainer container

You can download a prebuilt Apptainer container with all SCALPEL dependencies from the following link:  
[Download SCALPEL Container](https://zenodo.org/records/15717636/files/SCALPEL.container.sif?download=1)

1. Download the container and clone the repository
```bash
wget https://data.cyverse.org/dav-anon/iplant/home/franzx5/SCALPEL.container.sif
git clone https://github.com/p-CMRC-LAB/SCALPEL.git
```

2. Run SCALPEL using the container
```bash
nextflow run SCALPEL/main.nf \
  -with-apptainer SCALPEL.container.sif \
  --sequencing chromium \
  --samplesheet path/to/samplesheet.csv \
  --transcriptome path/to/gencode.transcripts.fa \
  --gtf path/to/gencode.annotation.gtf \
  --ipdb path/to/mm10.polyA.track
```

## Required Input Files

| Parameter         | Description                                                             |
|------------------|-------------------------------------------------------------------------|
| `--samplesheet`  | CSV with sample names and paths to FASTQ/BAM/CellRanger output          |
| `--transcriptome`| FASTA of reference transcriptome                                        |
| `--gtf`          | GTF annotation file                                                     |
| `--ipdb`         | Internal priming annotation file                                        |
| `--barcodes`     | (Optional) Barcode whitelist per sample                                 |
| `--clusters`     | (Optional) Tab-delimited file with cell-to-cluster mappings             |
| `--sequencing`   | Must be `chromium` or `dropseq`                                         |

Reference files:
- [Mouse IP Annotation - mm10](https://zenodo.org/records/15664563/files/mm10_polya.track.tar.gz?download=1)
- [Human IP Annotation - GRCh38](https://zenodo.org/records/15717592/files/hg38_ipriming_sites.bed.tar.gz?download=1)
- [GENCODE Annotations](https://www.gencodegenes.org/)

## Output Files and Execution Notes

After execution, SCALPEL generates a `results/` directory containing key outputs for downstream analysis.

### Key Output Files

| File / Pattern             | Description                                                                 |
|----------------------------|-----------------------------------------------------------------------------|
| `*_filtered.bam`           | BAM files with deduplicated reads excluding internal priming artifacts.     |
| `*_filtered.bam.bai`       | BAM index files.                                                           |
| `*_APADGE.txt`             | APA-aware isoform-level expression matrix per sample.                      |
| `*_seurat.RDS`             | Seurat object per sample.                                                  |
| `iDGE_seurat.RDS`          | Merged Seurat object across all samples.                                  |
| `DIU_table.csv`            | Differential isoform usage table.                                         |
| `Runfiles/`                | Execution logs and process metadata.                                       |

### Notes

- Output filenames are prefixed by the sample name.
- Seurat `.RDS` files are ready for downstream visualization and clustering in R.
- `*_APADGE.txt` matrices are compatible with other statistical environments.

For downstream analysis tutorials, visit:
- [Example of SCALPEL application on 10X scRNA-seq](https://github.com/p-CMRC-LAB/SCALPEL/wiki/SCALPEL-application-on-10X-scRNA%E2%80%90seq)
- [Example of SCALPEL application on DropSeq scRNA-seq](https://github.com/p-CMRC-LAB/SCALPEL/wiki/SCALPEL-application-on-DropSeq-scRNA%E2%80%90seq)
- [Downstream analysis Wiki](https://github.com/p-CMRC-LAB/SCALPEL/wiki)

---

## Customizing Execution with `nextflow.config`

To modify resource usage and process settings, optionally, you can edit the `nextflow.config` file. For example:

```groovy
/* Define Nextflow configuration settings for SCALPEL pipeline execution */
/* ===================================================================== */

/* -> Processes */

/* Enter here the desired Nextflow parameters for execution (see https://www.nextflow.io/docs/latest/index.html)*/

params {
    // General parameters
    workDir = './work'          // Directory for Nextflow work files
    process.executor = 'local'  // Default executor, can be overridden in the executor block
    process.queue = 'default'   // Default queue, can be overridden in the executor block

    //=====1. Specific parameters for SCALPEL pipeline====
    outputDir = "./results" // Directory for output files

    // 1a. Reads:
    samplesheet = null
    sequencing = null
    barcodes = null
    clusters = null

    // 1b. Annotation:
    transcriptome = null
    gtf = null
    ipdb = null

    // 1c. Thresholds:
    dt_threshold = 600
    de_threshold = 30
    ip_threshold = 60
    gene_fraction = "98%"
    binsize = 20
    subsample = 1
    help = null

    //=====2. Running parameters=====
    params.cpus = 16           /* Adjust as needed */
    params.memory = 40.GB     /* Adjust as needed */
}
```

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

## Reference

Franz Ake, Marcel Schilling, Sandra M. Fernández-Moya, Akshay Jaya Ganesh, Ana Gutiérrez-Franco, Lei Li, Mireya Plass  
**Quantification of transcript isoforms at the single-cell level using SCALPEL**  
*Nat Commun* 16, 6402 (2025). [doi://10.1038/s41467-025-61118-0](https://doi.org/10.1038/s41467-025-61118-0)

## Contact
Franz AKE – [@aerodx5](https://twitter.com/aerodx5) – franzarnold.ake@gmail.com  
GitHub: [https://github.com/p-CMRC-LAB/SCALPEL](https://github.com/p-CMRC-LAB/SCALPEL)

<p align="right"><a href="#top">Back to top</a></p>
