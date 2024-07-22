#! /usr/bin/env nextflow


// help message (can be displayed using --help option)
def helpMessage() {
    log.info """
    """
}

if (params.help) {
    helpMessage()
    exit 0
}


// check valid option for gene count algorithm and strandedness
def validCountAlgos = ['featureCounts', 'HTSeq']
def validStrandedness = [0, 1, 2]

if (!(params.count_algo in validCountAlgos)) {
    throw new IllegalArgumentException("Invalid value for 'count_algo'. Allowed values are: ${validCountAlgos.join(', ')}")
}
if (!(params.strandedness in validStrandedness)) {
    throw new IllegalArgumentException("Invalid value for 'strandedness'. Allowed values are: ${validStrandedness.join(', ')}")
}


process runGenomeIndexing {
    publishDir "${params.index_dir}", mode: 'move'
    
    output:
   val true  // for state dependency
    path "*"

    script:
    """
    STAR \
    --runMode genomeGenerate \
    --genomeDir . \
    --genomeFastaFiles $params.fasta_file \
    --sjdbGTFfile $params.annotation_file \
    --limitGenomeGenerateRAM $params.max_RAM_indexing \
    --runThreadN $params.genome_indexing_nt
    """
}

process runTrimming {
    publishDir "${params.trimmed_fastq_dir}", mode: 'copy', pattern: "*.fq.gz"
    publishDir "${params.trimmed_fastq_dir}/reports/", mode: 'move', pattern: "*.txt"
    publishDir "${params.trimmed_fastq_dir}/reports/", mode: 'move', pattern: "*.html"
    publishDir "${params.trimmed_fastq_dir}/reports/", mode: 'move', pattern: "*.zip"

    input:
    tuple path(input_fastq1), path(input_fastq2)

    output:
    tuple path(output_fastq1), path(output_fastq2)
    path '*.txt'
    path '*.html'
    path '*.zip'

    script:
    // extract core and output names
    core_fastq1 = ""
    core_fastq2 = ""
    if (input_fastq1.toString().endsWith('.fq.gz')) {
        core_fastq1 = input_fastq1.toString() - '.fq.gz'
    } else if (input_fastq1.toString().endsWith('.fastq.gz')) {
        core_fastq1 = input_fastq1.toString() - '.fastq.gz'
    }
    if (input_fastq2.toString().endsWith('.fq.gz')) {
        core_fastq2 = input_fastq2.toString() - '.fq.gz'
    } else if (input_fastq2.toString().endsWith('.fastq.gz')) {
        core_fastq2 = input_fastq2.toString() - '.fastq.gz'
    }
    output_fastq1 = core_fastq1 + '_val_1.fq.gz'
    output_fastq2 = core_fastq2 + '_val_2.fq.gz'
    
    // set number of threads
    n_threads = 1
    if (params.trimming_multithreaded) {
        n_threads = 4
    }

    """
    trim_galore \
    --quality $params.fq_quality_thres \
    --fastqc \
    --paired \
    --length $params.min_read_len \
    --cores ${n_threads} \
    ${input_fastq1} ${input_fastq2}
    """
}

process runAlignment {
    publishDir "${params.out_bam_dir}", mode: 'copy', pattern: "${bam}", enabled: ( params.save_all_BAM || params.last_BAM_output == "alignment" )
    publishDir "${params.out_bam_dir}/logs/", mode: 'move', pattern: "*.Log.final.out"
    publishDir "${params.out_bam_dir}/tabs/", mode: 'move', pattern: "*.tab"

    input:
    val ready  // for state dependency
    tuple path(fastq1), path(fastq2)

    output:
    path bam
    path '*.Log.final.out'
    path '*.tab'

    script:
    fastq_name = fastq1.toString().split("\\.")[0]
    bam = fastq_name + ".Aligned.bam"

    """
    STAR \
    --runMode alignReads \
    --readFilesCommand "gunzip -c" \
    --genomeDir $params.index_dir \
    --readFilesIn ${fastq1} ${fastq2} \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix "${fastq_name}." \
    --quantMode GeneCounts \
    --runThreadN $params.alignment_nt

    # just rename for nicer output
    mv "${core_name}.Aligned.out.bam" ${bam}
    """
}

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

process runRemoveDuplicates {
    publishDir "${params.out_bam_dir}", mode: 'copy', enabled: ( params.save_all_BAM || params.last_BAM_output == "duplicates" )

    input:
    path bam

    output:
    path bam_marked

    script:
    bam_marked = ""
    if (params.first_BAM_output == "duplicates") {
        bam_marked = bam.toString().split("\\.")[0] + ".Aligned.marked.bam"
    } else {
        bam_marked = bam.toString()[0..-5] + ".marked.bam"
    }
    metrics = bam.toString().split("\\.")[0] + ".dup_metrics.txt"

    """
    picard MarkDuplicates \
    --INPUT ${bam} \
    --OUTPUT ${bam_marked} \
    --REMOVE_SEQUENCING_DUPLICATES $params.remove_seq_duplicates \
    --METRICS_FILE "${params.out_bam_dir}/stats/${metrics}"
    """
}

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

process runBAMIndexing {
    publishDir "${params.out_bam_dir}", mode: 'move', pattern: "${bai}"

    input:
    path bam

    output:
    tuple path(bam), path(bai)

    script:
    bai = bam.toString() + ".bai"

    """
    samtools index -@ $params.BAM_indexing_nt "${bam}"
    """
}

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

process runHTSeq {
    input:
    tuple path(bam), path(bai)

    output:
    val true  // for state depencency

    script:
    """
    # set strandedness parameter
    if [ "$params.strandedness" -eq 0 ]; then
        strand="no"
    elif [ "$params.strandedness" -eq 1 ]; then
        strand="yes"
    elif [ "$params.strandedness" -eq 2 ]; then
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

process runSplicing {
    publishDir "${params.splicing_dir}", mode: 'move', pattern: '*.rmats'

    input:
    path bam_list
    path bai_list

    output:
    val true  // for state depencency

    script:
    // define files listing BAMs with requested conditions
    def file_condition1 = new File('file_condition1.txt')
    def file_condition2 = new File('file_condition2.txt')
    def string_condition1 = ""
    def string_condition2 = ""

    // write BAMs matching requested conditions to corresponding files
    for (int i=0; i<params.conditions.size(); i++) {
        if (params.conditions[i] == params.spl_condition1) {
            string_condition1 = string_condition1 + params.bam_files[i] + ","
        } else if (params.conditions[i] == params.spl_condition2) {
            string_condition2 = string_condition2 + params.bam_files[i] + ","
        }
    }
    file_condition1.write(string_condition1[0..-2] + "\n")
    file_condition2.write(string_condition2[0..-2] + "\n")

    """
    ./run_rmats \
    --b1 "file_condition1.txt" \
    --b2 "file_condition2.txt" \
    --gtf $params.annotation_file \
    -t paired \
    --readLength $params.read_length \
    --variable-read-length \
    --nthread $params.splicing_nt \
    --od $params.splicing_dir \
    --tmp .
    """
}

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


workflow {

    // run genome indexing
    def index_ready = true  // any value would be fine (just for state dependency)
    if (params.run_genome_indexing) {

        index_ready = runGenomeIndexing()[0]
    }


    // run trimming and fastQC reports
    def fastq_ch_trimmed = false
    def trimming_ready = false  // for multiQC
    if (params.run_trimming) {

        def fastq_ch = channel.from(params.fastq_files)

        fastq_ch_trimmed = runTrimming(fastq_ch)[0]

        // signal to multiQC that this process is done
        trimming_ready = fastq_ch_trimmed.collect()

    } else if (params.run_alignment) {

        fastq_ch_trimmed = channel.from(params.fastq_files)
    }


    // run reads alignment
    def bam_ch = false
    if (params.run_alignment) {

        bam_ch = runAlignment(index_ready, fastq_ch_trimmed)[0]

    } else if (params.run_BAM_sorting || params.run_remove_duplicates || params.run_BAM_filtering || params.run_BAM_indexing || params.run_BAM_stats || params.run_gene_counts){

        bam_ch = channel.fromPath(params.bam_files, checkIfExists: true)
    }


    // run BAM sorting
    def bam_ch_sorted = false
    if (params.run_BAM_sorting) {

        bam_ch_sorted = runBAMSorting(bam_ch)

    } else {

        bam_ch_sorted = bam_ch
    }


    // run remove duplicates of BAM
    def bam_ch_marked = false
    if (params.run_remove_duplicates) {

        bam_ch_marked = runRemoveDuplicates(bam_ch_sorted)

    } else {

        bam_ch_marked = bam_ch_sorted
    }


    // run alignment filtering
    def bam_ch_filtered = false
    if (params.run_BAM_filtering) {

        bam_ch_filtered = runBAMFiltering(bam_ch_marked)

    } else {

        bam_ch_filtered = bam_ch_marked
    }


    // run BAM files indexing
    def bam_ch_indexed = false
    if (params.run_BAM_indexing) {

        bam_ch_indexed = runBAMIndexing(bam_ch_filtered)
    }


    // run alignment stats summary
    def bam_stats_ch = false
    def bam_stats_ready = false
    if (params.run_BAM_stats) {

        bam_stats_ch = runBAMStats(bam_ch_filtered)

        bam_stats_ready = bam_stats_ch.collect()
    }


    // build BAM-BAI channels if BAM indexing was not run
    if (params.run_gene_counts || params.run_splicing) {

        // if indexing not run, extract BAM and corresponding index file and define channel
        if (!bam_ch_indexed) {
            def bam_bai_pairs = params.bam_files.collect{ path -> return (path.toString() + "{,.bai}") }
            bam_ch_indexed = channel.fromFilePairs(bam_bai_pairs, checkIfExists: true).map{baseName, fileList -> fileList}
        }
    }


    // run gene counts
    def counts_ready = false
    if (params.run_gene_counts) {

        // run featureCounts
        if (params.count_algo == 'featureCounts') {

            def counts_ready_par = runFeatureCounts(bam_ch_indexed)

            // make sure all counts have been performed
            counts_ready = counts_ready_par.collect()

        // run HTSeq
        } else if (params.count_algo == 'HTSeq') {

            def counts_ready_par = runHTSeq(bam_ch_indexed)

            // make sure all counts have been performed
            counts_ready = counts_ready_par.collect()
        }
    }


    // run splicing analysis
    def splicing_ready = false
    if (params.run_splicing) {

        // gather all BAM-BAI couples for running splicing
        def bam_bai_list = bam_ch_indexed.collect()
        def bam_list = []
        def bai_list = []
        bam_bai_list.each { pair ->
            bam_list << pair[0]
            bai_list << pair[1]
        }
        
        // sort paths to be able to retrieve conditions
        bam_list.sort()
        bai_list.sort()

        // run proper analysis
        splicing_ready = runSplicing(bam_list, bai_list)
    }


    // run results summary
    if (params.run_summarize_results) {

        // make sure all data have been processed before going on
        if (!trimming_ready) {
            trimming_ready = channel.of(1)  // "dummy" channel
        }

        if (!bam_stats_ready) {
            bam_stats_ready = channel.of(1)  // "dummy" channel
        }

        if (!counts_ready) {
            counts_ready = channel.of(1)  // "dummy" channel
        }

        if (!splicing_ready) {
            splicing_ready = channel.of(1)  // "dummy" channel
        }

        runSumResults(trimming_ready, bam_stats_ready, counts_ready, splicing_ready)
    }
}
