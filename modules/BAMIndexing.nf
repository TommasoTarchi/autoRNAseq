process runBAMIndexing {
    publishDir "${params.out_bam_dir}", mode: 'copy', pattern: "${bai}"

    input:
    path bam

    output:
    tuple path(bam), path(bai)

    script:
    bai = bam.toString() + ".bai"

    """
    samtools index -@ $params.BAM_indexing_nt "${bam}"
    """
}
