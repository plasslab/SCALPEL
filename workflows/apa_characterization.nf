

/* 
	Downstream analysis of samples (default params)
*/

process differential_isoform_usage{
	tag "${objs}"
	publishDir "${params.outputDir}/", overwrite: true, mode: 'copy'
    label 'big_rec'

	input:
		file(objs)

	output:
		tuple path("DIU_table.csv"), file("iDGE_seurat.RDS")
	
	script:
		if( params.clusters == null )
			"""
			cp ${baseDir}/src/scalpel_library.R .
			Rscript ${baseDir}/src/analysis_samples.R ./ NULL DIU_table.csv
			"""

		else if( params.clusters != null )
			"""
			cp ${baseDir}/src/scalpel_library.R .
			Rscript ${baseDir}/src/analysis_samples.R ./ ${params.clusters} DIU_table.csv
			"""
}


process generation_filtered_bams{
    tag "${sampleID}"
    publishDir "${params.outputDir}/", overwrite: true, mode: 'copy'
    label 'big_rec'

    input:
        tuple val(sampleID), file(bams)
        tuple val(sampleID), file(rids)
    output:
        file("${sampleID}_filtered.bam*")
    script:
        """
        cat *.readid > ALL_RIDS.txt
	    samtools merge -o tmp.bam ${bams} -@ ${task.cpus}
        samtools view -b -N ALL_RIDS.txt tmp.bam > ${sampleID}_filtered.bam 
        samtools index ${sampleID}_filtered.bam
        rm tmp.bam
        """
}
