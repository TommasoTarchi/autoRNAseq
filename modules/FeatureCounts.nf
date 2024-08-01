process runFeatureCounts {
    input:
    tuple path(bam), path(bai)

    output:
    val true  // for state depencency

    script:
    """
    bam_name=\$(basename "${bam}")
    core_name="\${bam_name%%.*}"

    featureCounts \
    -a $params.annotation_file \
    -g gene_id \
    -t exon \
    -p \
    -s 2 \
    -o "$params.gene_counts_dir/\${core_name}.counts.txt" \
    -T $params.gene_counts_nt \
    ${bam}
    """
}
