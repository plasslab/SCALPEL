#!/usr/bin/env nextflow

/*
Loading of SAMPLES data & preprocessing according to the sequencing type
=====================================================================
*/

workflow samples_loading {
    /* Loading of samples */
    take:
        samples_paths
        selected_isoforms

    main:

        /* - In case of 10X file type, extract required files... */
        if ( "${params.sequencing}" == "chromium" ) {

            // extract sampleIDs and associated paths
            read_10x( samples_paths.map{ it=tuple(it[0], it[3]) }).set{ samples_selects }

        } else {

            samples_paths.map{ it = tuple( it[0], file(it[3]), file(it[4]), file(it[5]) ) }.set{ samples_selects }
        }

        if ( params.barcodes != null ) {

            ( Channel.fromPath(params.barcodes) | splitCsv(header:false) ).set{ barcodes_paths }
            (samples_selects.join(barcodes_paths, by:[0])).map{ it = tuple( it[0], it[1], it[2], it[5], it[4] )}.set{ samples_selects }

        }
        /* Integrate isoforms metadata */
        selected_isoforms.flatMap { it = it[0] }.combine(samples_selects).set{ samples }

        /* split bam files */
        bam_splitting( samples )

    emit:
        selected_bams = bam_splitting.out
        sample_files = samples
}



process read_10x {
    tag "${sample_id}, ${repo}"
    cache true
    label "small_mem"

    input:
        tuple val(sample_id), path(repo)
    output:
        tuple val(sample_id), path("${sample_id}.bam"), path("${sample_id}.bam.bai"), path("${sample_id}.barcodes"), path("${sample_id}.counts.txt")
    script:
    """
        Rscript ${baseDir}/src/read_10X.R ${repo}
        ln -s ${repo}/outs/possorted_genome_bam.bam ${sample_id}.bam
        ln -s ${repo}/outs/possorted_genome_bam.bam.bai ${sample_id}.bam.bai
        zcat ${repo}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > ${sample_id}.barcodes
    """
}


process bam_splitting {
    tag "${sample_id}, ${chr}"
    publishDir "${params.outputDir}/Runfiles/reads_processing/bam_splitting/${sample_id}"
    cache true
    label 'small_mem'

    input: 
        tuple val(chr), val(sample_id), path(bam), path(bai), val(bc_path), path(dge_matrix)
    output:
        tuple val(sample_id), val("${chr}"), path("${chr}.bam")
    script:
    if( params.sequencing == "dropseq" )
        if ( params.barcodes != null )
            """
            #Dropseq
            #Filter reads , Remove duuplicates and split by chromosome
            samtools view --subsample ${params.subsample} -b ${bam} ${chr} -D XC:${bc_path} --keep-tag "XC,XM" | samtools sort > tmp.bam
            #Remove all PCR duplicates ...
            samtools markdup tmp.bam ${chr}.bam -r --barcode-tag XC --barcode-tag XM
            rm tmp.bam
            """
        else
            """
            #Filter reads , Remove duuplicates and split by chromosome
            samtools view --subsample ${params.subsample} -b ${bam} ${chr} --keep-tag "XC,XM" | samtools sort > tmp.bam
            #Remove all PCR duplicates ...
            samtools markdup tmp.bam ${chr}.bam -r --barcode-tag XC --barcode-tag XM
            rm tmp.bam
            """
    else if( params.sequencing == "chromium" )
        """
        #Chromium_seq
        #Filter reads , Remove duuplicates and split by chromosome
        samtools view --subsample ${params.subsample} -b ${bam} ${chr} -D CB:${bc_path} --keep-tag "CB,UB" | samtools sort > tmp.bam
        #Remove all PCR duplicates ...
        samtools markdup tmp.bam ${chr}.bam -r --barcode-tag CB --barcode-tag UB
        rm tmp.bam
        """
}

    
process bedfile_conversion{
    tag "${sample_id}, ${chr}"
    publishDir "${params.outputDir}/Runfiles/reads_processing/bedfile_conversion/${sample_id}"
    cache true
    label "small_mem"

    input:
        tuple val(sample_id), val(chr), path(bam)
    output:
        tuple val(sample_id), val(chr), path("${chr}.bed")
    script:
    """
        #Convertion of bam files to bed files
        bam2bed --all-reads --split --do-not-sort < ${bam} | gawk -v OFS="\\t" '{print \$1,\$2,\$3,\$6,\$4,\$14"::"\$15}' > ${chr}.bed
    """
}


process reads_mapping_and_filtering {
    tag "${sample_id}, ${chr}, ${bed}, ${exons}"
    publishDir "${params.outputDir}/Runfiles/reads_processing/mapping_filtering/${sample_id}"
    cache true
    label "big_mem"

    input:
        tuple val(sample_id), val(chr), path(bed), path(exons), path(exons_unique), path(collapsed)
    output:
        tuple val(sample_id), val(chr), path("${chr}.reads"), path("${sample_id}_${chr}_read_distrib.txt")
    script:
    """
        Rscript ${baseDir}/src/mapping_filtering.R ${bed} ${exons} ${chr}.reads
        mv read_distance_distribution_on_transcriptomic_scope.txt ${sample_id}_${chr}_read_distrib.txt
    """
}


process ip_splitting {
    tag "${chr}"
    publishDir "${params.outputDir}/Runfiles/internalp_filtering/ipdb_splitted"
    cache true
    label "small_mem"

    input:
        tuple val(chr), path(ipref)
    output:
        tuple val(chr), path("${chr}_sp.ipdb")
    script:
    """
    #select chromosome
    gawk '{ if(\$1=="${chr}") print }' ${ipref} > ${chr}_sp.ipdb
    """
}


process ip_filtering {
    tag "${sample_id}, ${chr}, ${ipdb}"
    publishDir "${params.outputDir}/Runfiles/internalp_filtering/internalp_filtered/${sample_id}"
    cache true
    label "big_mem"

    input:
        tuple val(sample_id), val(chr), path(reads), path(ipdb)
        val(ip_thr)
    output:
        tuple val(sample_id), val(chr), path("${chr}_mapped.ipdb"), path("${chr}_unique.reads"), path("${chr}.readid"), path("${chr}_ipf.reads")
    script:
        """
        #filtering
        Rscript ${baseDir}/src/ip_filtering.R ${reads} ${ipdb} ${ip_thr} ${chr}.readid ${chr}_unique.reads ${chr}_mapped.ipdb ${chr}_ipf.reads
        """
}
