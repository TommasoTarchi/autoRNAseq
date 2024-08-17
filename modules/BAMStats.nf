process runBAMStats {
    input:
    path bam

    output:
    val true  // for state depencency

    script:
    """
    # create needed subdirectory if not existing
    if [[ -d "$params.out_bam_dir/stats" ]]; then
        mkdir "$params.out_bam_dir/stats"
    fi

    # define BAM base name
    bam_name=\$(basename "${bam}")

    # run SAMtools
    samtools flagstat -@ $params.BAM_stats_nt ${bam} > "$params.out_bam_dir/stats/\${bam_name}.stats.txt"
    """
}
