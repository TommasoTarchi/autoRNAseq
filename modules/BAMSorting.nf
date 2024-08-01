process runBAMSorting {
    publishDir "${params.out_bam_dir}", mode: 'copy', enabled: ( params.save_all_BAM || params.last_BAM_output == "sorting" )

    input:
    path bam

    output:
    path bam_sorted

    script:
    bam_sorted = ""
    if (params.first_BAM_output == "sorting") {
        bam_sorted = bam.toString().split("\\.")[0] + ".Aligned.sortedByCoord.bam"
    } else {
        bam_sorted = bam.toString()[0..-5] + ".sortedByCoord.bam"
    }

    """
    if samtools view -H ${bam} | grep -q '@HD.*SO:coordinate'; then
        mv ${bam} ${bam_sorted}
    else
        samtools sort -@ $params.BAM_sorting_nt -o ${bam_sorted} ${bam}
    fi
    """
}
