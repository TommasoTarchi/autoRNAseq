process runBAMStats {
    input:
    path bam

    output:
    val true  // for state depencency

    script:
    """
    bam_name=\$(basename "${bam}")

    samtools flagstat -@ $params.BAM_stats_nt ${bam} > "$params.out_bam_dir/stats/\${bam_name}.stats.txt"
    """
}
