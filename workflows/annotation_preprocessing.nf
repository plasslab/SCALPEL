#!/usr/bin/env nextflow

/*
Loading of ANNOTATION data & preprocessing
==========================================
*/

process salmon_transcriptome_indexing {
    publishDir "${params.outputDir}/Runfiles/annotation_processing/salmon_bulk_quantification", overwrite: true
    cache true
    label 'big_rec'

    input:
        path(transcriptome_reference)
    output:
        path "transcriptome_index"
    script:
        """
        salmon index -t ${transcriptome_reference} -i transcriptome_index --gencode -p 4
        """
}


process salmon_bulk_quantification{
    tag "${pair_id}, ${fastq1}, ${fastq2}"
    publishDir "${params.outputDir}/Runfiles/annotation_processing/salmon_bulk_quantification", overwrite: true
    cache true
    label 'big_rec'
    
    input:
        path(transcriptome_index)
        tuple val(pair_id), path(fastq1), path(fastq2)
    output:
        path "${pair_id}.sf"
    script:
        """
        salmon quant -i ${transcriptome_index} -l ISR -1 ${fastq1} -2 ${fastq2} -o ${pair_id} --minScoreFraction 0.7 --vbPrior 5 --perNucleotidePrior -p 4
        mv ${pair_id}/quant.sf ${pair_id}.sf
        """
}


process tpm_counts_average{
    tag "${bulk_quants}"
    publishDir "${params.outputDir}/Runfiles/annotation_processing/salmon_bulk_quantification", overwrite: true
    cache true
    label "big_rec"

    input:
        path bulk_quants
        path gtf_annotation_reference
    output:
        tuple file("merge_quants.txt"), file("*.gtf")
    script:
        """
        Rscript ${baseDir}/src/quantification_processing.R ${gtf_annotation_reference} merge_quants.txt
        gawk '{print>\$1".gtf"}' reduced_annotation.tsv
        """
}


process isoform_selection_weighting{
    tag "${gtf.baseName}, ${merged_quants}"
    publishDir "${params.outputDir}/Runfiles/annotation_processing/isoform_processing", overwrite: true, mode: 'copy'
    cache true
    label "small_rec"

    input:
        tuple path(merged_quants), path(gtf)
        val(dt_thr)
        val(de_thr)
    output:
        tuple val("${gtf.baseName}"), path("${gtf.baseName}.exons"), path("${gtf.baseName}.exons_unique"), path("${gtf.baseName}_collapsed_isoforms.txt")
    script:
        """
        Rscript ${baseDir}/src/gtf_processing.R ${gtf} ${merged_quants} ${dt_thr} ${de_thr} ${gtf.baseName}.exons ${gtf.baseName}.exons_unique
        """
}
