# SCALPEL: A Nextflow-based Pipeline for Isoform Quantification at Single-Cell Resolution

<div align="center">
  <img src="https://data.cyverse.org/dav-anon/iplant/home/franzx5/SPERMATOGENESIS/SCALPEL_pipeline.png" alt="SCALPEL" width="800">
</div>

## About the Project

**SCALPEL** is a robust pipeline designed for transcript isoform quantification and alternative polyadenylation (APA) characterization using 3'-tagged single-cell RNA-seq (scRNA-seq) data. Built with **Nextflow**, it integrates multiple processing steps including read quantification, APA annotation, and isoform usage analysis.

## Installation Options

SCALPEL can be installed and run using one of the following options:

### Option 1: Using Conda

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

## Documentation and Examples

- [SCALPEL application on 10X scRNA-seq](https://github.com/p-CMRC-LAB/SCALPEL/wiki/SCALPEL-application-on-10X-scRNA%E2%80%90seq)
- [SCALPEL application on DropSeq scRNA-seq](https://github.com/p-CMRC-LAB/SCALPEL/wiki/SCALPEL-application-on-DropSeq-scRNA%E2%80%90seq)
- [Downstream analysis](https://github.com/p-CMRC-LAB/SCALPEL/wiki)

---

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

Franz Ake, Sandra M. Fernández-Moya, Marcel Schilling, Akshay Jaya Ganesh, Ana Gutiérrez-Franco, Lei Li, Mireya Plass  
**Quantification of transcript isoforms at the single-cell level using SCALPEL**  
bioRxiv 2024.06.21.600022; doi: [10.1101/2024.06.21.600022](https://doi.org/10.1101/2024.06.21.600022)

---

## Contact

Franz AKE – [@aerodx5](https://twitter.com/aerodx5) – fake@idibell.cat  / aerod7710@gmail.com

GitHub: [https://github.com/p-CMRC-LAB/SCALPEL](https://github.com/p-CMRC-LAB/SCALPEL)

<p align="right"><a href="#top">Back to top</a></p>
