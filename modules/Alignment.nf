process runAlignment {
    publishDir "${params.out_bam_dir}", mode: 'copy', pattern: "${bam}", enabled: ( params.save_all_BAM || params.last_BAM_output == "alignment" )
    publishDir "${params.out_bam_dir}/logs/", mode: 'move', pattern: "*.Log.final.out"
    publishDir "${params.out_bam_dir}/tabs/", mode: 'move', pattern: "*.tab"

    input:
    val ready  // for state dependency
    tuple path(fastq1), path(fastq2)

    output:
    path bam
    path '*.Log.final.out'
    path '*.tab'

    script:
    fastq_name = fastq1.toString().split("\\.")[0]
    bam = fastq_name + ".Aligned.bam"

    """
    STAR \
    --runMode alignReads \
    --readFilesCommand "gunzip -c" \
    --genomeDir $params.index_dir \
    --readFilesIn ${fastq1} ${fastq2} \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix "${fastq_name}." \
    --quantMode GeneCounts \
    --runThreadN $params.alignment_nt

    # just rename for nicer output
    mv "${core_name}.Aligned.out.bam" ${bam}
    """
}
