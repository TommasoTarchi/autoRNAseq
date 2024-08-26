import java.nio.file.Paths


process runSplicing {
    input:
    path bam_list  // not single path but list
    path bai_list  // not single path but list
    tuple val(contrast), val(treatment), val(control)  // contrast to perform

    output:
    val true  // for state depencency

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
        strand="fr-secondstrand"
    elif [ "$params.spl_strandedness" -eq 2 ]; then
        strand="fr-firststrand"
    fi

    # write paths to files matching conditions to files
    echo ${string_treatment} > list_treatment.txt
    echo ${string_control} > list_control.txt

    # remove needed subdirectory if already existing (cleaning)
    if [[ -d "$params.splicing_dir/${contrast}" ]]; then
        rm -r "$params.splicing_dir/${contrast}"
    fi

    # create needed subdirectory
    mkdir "$params.splicing_dir/${contrast}"

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
    --od "$params.splicing_dir/${contrast}" \
    --tmp "$params.splicing_dir/${contrast}" \
    $rmats_options \
    1> rmats.log

    # move log file to log directory
    if [[ ! -d "$params.splicing_dir/${contrast}/logs" ]]; then
        mkdir "$params.splicing_dir/${contrast}/logs"
    fi
    mv rmats.log "$params.splicing_dir/${contrast}/logs/"
    """
}
