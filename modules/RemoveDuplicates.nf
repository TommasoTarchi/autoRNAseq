process runRemoveDuplicates {
    publishDir "${params.out_bam_dir}", mode: 'copy', enabled: ( params.save_all_BAM || params.last_BAM_output == "duplicates" )

    input:
    path bam

    output:
    path bam_marked

    script:
    // define output BAM and metrics file name
    bam_marked = ""
    if (params.first_BAM_output == "duplicates") {
        bam_marked = bam.toString().split("\\.")[0] + ".Aligned.marked.bam"
    } else {
        bam_marked = bam.toString()[0..-5] + ".marked.bam"
    }
    metrics = bam.toString().split("\\.")[0] + ".dup_metrics.txt"

    """
    # create needed subdirectory if not existing
    if [[ ! -d "${params.out_bam_dir}/stats" ]]; then
        mkdir "${params.out_bam_dir}/stats"
    fi

    # run picard
    picard MarkDuplicates \
    --INPUT ${bam} \
    --OUTPUT ${bam_marked} \
    --REMOVE_SEQUENCING_DUPLICATES $params.remove_seq_duplicates \
    --METRICS_FILE "${params.out_bam_dir}/stats/${metrics}"
    """
}
