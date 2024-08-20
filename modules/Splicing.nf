import java.nio.file.Paths


process runSplicing {
    publishDir "${params.splicing_dir}", mode: 'move', pattern: "*.rmats"
    publishDir "${params.splicing_dir}", mode: 'move', pattern: "*_read_outcomes_by_bam.txt"

    input:
    path bam_list  // not single path but list
    path bai_list  // not single path but list
    val treatment  // first condition to compare
    val control  // second condition to compare
    val comparison  // comparison (i.e. name of the folder to store results)

    output:
    val true  // for state depencency
    path "*.rmats"  // stats files
    path "*_read_outcomes_by_bam.txt"  // report on used reads

    script:
    // redefine bamlist and path to avoid conflicts
    def bamlist = bam_list
    def paths = Paths

    // define strings listing BAMs with requested conditions
    def string_treatment = ""
    def string_control = ""

    // extract core names of original input files
    def core_names = []
    if (params.run_trimming || params.run_alignment) {
       for (file_pair in params.fastq_files) { 
            core_names << paths.get(file_pairs[0].toString()).getFileName().toString().split("\\.")[0]
       }
    } else {
        for (file_path in params.bam_files) {
            core_names << paths.get(file_path.toString()).getFileName().toString().split("\\.")[0]
        }
    }

    // write BAMs matching requested conditions to corresponding strings
    for (int i=0; i<bamlist.size(); i++) {  // iterate over BAM list
        def current_core_name = paths.get(bamlist[i].toString()).getFileName().toString().split("\\.")[0]

        for (int j=0; j<core_names.size(); j++) {  // iterate over original inputs

            if (current_core_name == core_names[j]) {  // look for matching pattern and condition
                if (params.conditions[j] == treatment.toString()) {
                    string_treatment = string_treatment + bamlist[i] + ","
                } else if (params.conditions[j] == control.toString()) {
                    string_control = string_control + bamlist[i] + ","
                }
            }
        }
    }

    // correct strings for commas
    string_treatment = string_treatment[0..-2]
    string_control = string_control[0..-2]

    // set options for rMATS-turbo (in case paired statistics and/or novel splice
    // sites detection are required)
    def rmats_options = ""
    if (params.use_paired_stats) {
        rmats_options = rmats_options + "--paired-stats"
    }
    if (params.detect_novel_splice) {
        rmats_options = rmats_options + " --novelSS --mil " + params.spl_min_intron_len.toString() + " --mel " + params.spl_max_exon_len.toString()
    }

    """
    # set strandedness parameter
    strand=""
    if [ "$params.spl_strandedness" -eq 0 ]; then
        strand="fr-unstranded"
    elif [ "$params.spl_strandedness" -eq 1 ]; then
        strand="fr-firststrand"
    elif [ "$params.spl_strandedness" -eq 2 ]; then
        strand="fr-secondstrand"
    fi

    # write paths to files matching conditions to files
    echo ${string_treatment} > list_treatment.txt
    echo ${string_control} > list_control.txt

    # create needed subdirectory if not existing
    if [[ ! -d "$params.splicing_dir/${comparison}" ]]; then
        mkdir "$params.splicing_dir/${comparison}"
    fi

    # run rMATS-turbo
    rmats.py \
    --task both \
    --b1 list_treatment.txt \
    --b2 list_control.txt \
    --gtf $params.annotation_file \
    -t paired \
    --libType "\${strand}" \
    --readLength $params.spl_read_len \
    --variable-read-length \
    --cstat $params.spl_cutoff_diff \
    --allow-clipping \
    --nthread $params.splicing_nt \
    --od "$params.splicing_dir/${comparison}" \
    --tmp . \
    $rmats_options \
    1> rmats.log

    # remove temporary files
    rm -r "$params.splicing_dir/${comparison}/tmp/"
    """
}
