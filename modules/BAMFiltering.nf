process runBAMFiltering {
    publishDir "${params.out_bam_dir}", mode: 'copy', enabled: ( params.save_all_BAM || params.last_BAM_output == "filtering" )

    input:
    path bam

    output:
    path bam_filtered

    script:
    bam_filtered = ""
    if (params.first_BAM_output == "filtering") {
        bam_filtered = bam.toString().split("\\.")[0] + ".Aligned.filtered.bam"
    } else {
        bam_filtered = bam.toString()[0..-5] + ".filtered.bam"
    }

    """
    samtools view --threads $params.BAM_filtering_nt -b -q $params.BAM_quality_thres ${bam} > ${bam_filtered}
    """
}
