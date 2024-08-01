process runSumResults {
    input:
    val ready1
    val ready2
    val ready3
    val ready4

    script:
    """
    # set available directories
    dirs=()
    if [[ -n "$params.index_dir" ]]; then
        dirs+=("$params.index_dir")
    fi
    if [[ -n "$params.trimmed_fastq_dir" ]]; then
        dirs+=("$params.trimmed_fastq_dir")
    fi
    if [[ -n "$params.out_bam_dir" ]]; then
        dirs+=("$params.out_bam_dir")
    fi
    if [[ -n "$params.gene_counts_dir" ]]; then
        dirs+=("$params.gene_counts_dir")
    fi
    if [[ -n "$params.splicing_dir" ]]; then
        dirs+=("$params.splicing_dir")
    fi

    # run multiQC if any directory available
    if ! [[ -z "\${dirs[@]}" ]]; then
        multiqc "\${dirs[@]}" \
        --force \
        --export \
        --outdir $params.report_dir
    fi
    """
}
