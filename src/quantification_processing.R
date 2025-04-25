#!/usr/bin/env Rscript

<<<<<<< HEAD
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
=======
#libs
suppressWarnings(library(argparse, quietly=T, verbose=F))
suppressWarnings(library(dplyr, quietly=T, verbose=F))

#Argument parser
parser = ArgumentParser(description='Processing of salmon quantification_file')
parser$add_argument('genome', type="character", help='path of gtf file')
parser$add_argument('output_path', type="character", help='output path')
args = parser$parse_args()

#Opening
#GTF
gtf = rtracklayer::import(args$genome) %>%
    data.frame() %>%
    dplyr::filter(type=="exon")

#check if exon entries into the GTF
if(nrow(gtf)==0){ stop("Error: No 'exon' feature into the GTF ! Please provide valid GTF") }

#-check presence of required attributes
if( FALSE %in% c(c("gene_id","transcript_id") %in% colnames(gtf)) ){ 
    stop("Error! The attributes 'gene_id', 'transcript_id' are expected in your GTF file") 
}

#-integrate gene_name & transcript_name features
if( FALSE %in% c(c("gene_name", "transcript_name") %in% colnames(gtf)) ){
    gtf = dplyr::mutate(gtf, gene_name=gene_id, transcript_name==transcript_id)
}

#-subset GTF file to specific columns
gtf = dplyr::distinct(gtf, seqnames, start, end, width, strand, gene_id, gene_name, transcript_id, transcript_name)

#-Processing of Salmon bulk quantification files
bulk_quants = lapply(list.files(path = ".", pattern = "*.sf", full.names = T), function(x){
    #--reading bulk counts
>>>>>>> main
    curr = data.table::fread(x, nThread=1)
    curr$sample = basename(x)
    return(curr)
}) %>% data.table::rbindlist() %>%
<<<<<<< HEAD
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
=======
    #get MeanTPM
    dplyr::group_by(Name) %>%
    #get the average TPM in the input samples
    dplyr::summarise(meanTPM = mean(TPM)) %>%
    dplyr::rename(transcript_id="Name") %>%
    ungroup() %>%
    distinct()

#-Integrate priors calculation to the gtf annotation
bulk_quants = dplyr::left_join(distinct(gtf,gene_name,transcript_name,transcript_id), bulk_quants, by=c("transcript_id")) %>%
    distinct() %>%
    dplyr::group_by(gene_name) %>%
    #calculate Transcripts bulk % in gene
    dplyr::mutate(transcript_name, bulk_TPMperc = meanTPM / sum(meanTPM)) %>%
    #change NA values resulting from null transcript mean calculation to 0
    dplyr::mutate(bulk_TPMperc=if_else(is.na(bulk_TPMperc),0,bulk_TPMperc)) %>%
    ungroup() %>%
    distinct(gene_name, transcript_name, meanTPM, bulk_TPMperc) %>% 
    arrange(gene_name, transcript_name)

#writing bulk quantification priors scores
bulk_quants %>% data.table::fwrite(file=args$output_path, nThread=1, row.names=F, sep="\t")

#split the gtf file by chromosome for processing
lapply(split(gtf, as.character(gtf$seqnames)), function(x){
    #writing
    x = distinct(x, seqnames,start,end,width,strand,gene_name,transcript_name)
    data.table::fwrite(x, file = paste0(as.character(x$seqnames)[1],".gtf"), sep="\t", nThread=1)
}) -> writeds
>>>>>>> main
