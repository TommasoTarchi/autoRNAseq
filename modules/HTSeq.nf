process runHTSeq {
    input:
    tuple path(bam), path(bai)

    output:
    val true  // for state depencency

    script:
    """
    # set strandedness parameter
    if [ "$params.gc_strandedness" -eq 0 ]; then
        strand="no"
    elif [ "$params.gc_strandedness" -eq 1 ]; then
        strand="yes"
    elif [ "$params.gc_strandedness" -eq 2 ]; then
        strand="reverse"
    fi

    bam_name=\$(basename "${bam}")
    core_name="\${bam_name%%.*}"

    htseq-count \
    -f bam \
    -r pos \
    -i gene_id \
    -t exon \
    -s "\${strand}" \
    ${bam} \
    $params.annotation_file \
    > "$params.gene_counts_dir/\${core_name}.counts.txt"
    """
}
