#!/usr/bin/env Rscript

# =============================================================================
# Script: quantification_processing.R
# Description: Merges Salmon quantifications and GTF for transcript abundance
# Inputs:
#   1. Genome GTF file path
#   2. Output file path for merged TPM table
# Outputs:
#   - Merged quantification table with TPM percentages (to args$output_path)
#   - Split GTFs per chromosome in the current directory
# =============================================================================


#libs
suppressWarnings({
  library(argparse, quietly=T, verbose=F)
  library(dplyr, quietly=T, verbose=F)
})

#- Argument parsing
parser <- ArgumentParser(description = 'Process Salmon quantifications and GTF')
parser$add_argument('genome', type = "character", help = 'GTF file path')
parser$add_argument('output_path', type = "character", help = 'Output path for merged TPM summary')
args <- parser$parse_args()

#- Opening
#gtf.file...
gtf = rtracklayer::import(args$genome)
gtf = data.frame(gtf[gtf$type=="exon",])
gtf.attributes = c("seqnames", "start", "end", "width", "strand", "gene_id", "gene_name", "transcript_id", "transcript_name")
gtf = gtf[, gtf.attributes]

#bulk quantification...
bulk.files = list.files(path = ".", pattern = "*.sf", full.names = T)

#- Processing bulk quantification sample files... and write merge_quant file
RECEPT = lapply(bulk.files, function(x){
    #reading bulk counts
    curr = data.table::fread(x, nThread=1)
    curr$sample = basename(x)
    return(curr)
}) %>% data.table::rbindlist() %>%
    #avoiding division by zero
    dplyr::mutate(TPM = TPM + 1) %>%
    #get MeanTPM & discard null transcrips in bulk
    dplyr::group_by(Name) %>%
    dplyr::reframe(meanTPM = mean(TPM)) %>%
    # dplyr::filter(meanTPM!=0) %>%
    dplyr::rename(transcript_id="Name") %>%
    #include gene metadata
    dplyr::left_join(distinct(gtf, gene_name, transcript_name, transcript_id), by="transcript_id") %>%
    #only considers all the transcript in the annotation GTF...
    na.omit() %>%
    #calculate Transcripts bulk % in gene
    dplyr::group_by(gene_name) %>%
    dplyr::reframe(transcript_name, bulk_TPMperc = meanTPM / sum(meanTPM)) %>%
    dplyr::arrange(gene_name,transcript_name) %>%
    data.table::data.table() %>%
    #write bulk quantification file
    data.table::fwrite(file=args$output_path, nThread=1, row.names=F, sep="\t")

#- Write GTF by chromosomes...
data.table::fwrite(gtf, file = "reduced_annotation.tsv", sep="\t", nThread=1, row.names=F, col.names=F)
