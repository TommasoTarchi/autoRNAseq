#! /usr/bin/env nextflow


// help message (can be displayed using --help option)
def helpMessage() {
    log.info """
        USAGE:


        This program can be used to run a RNAseq analysis pipeline, consisting of the following
        steps (also called "processes" in the following):

        1. Genome Indexing: preprocess the genome for alignment.
        2. Alignment: properly align reads to the reference genome.
        3. BAM Sorting: sort BAM files.
        4. Remove duplicates: remove duplicates in BAM files.
        5. BAM Filtering: quality filtering of aligned reads.
        6. BAM Indexing: index the alignment files.
        7. BAM Stats: generate a statistical summary of the alignment.
        8. Gene Counts: quantify gene expression.
        9. Results Summary: summarize the results.


        Any step of the pipeline can be run (multiple steps can be run together), see 'ARGUMENTS'
        and 'HOW TO RUN A PIPELINE' below for instructions.

        The program is designed to run on a cluster using SLURM as job scheduler: each step is run
        as an independent (and parallel, when no data dependency occurs) multithreaded job.

        Optionally, it can also be run locally. However, we strongly recommend to run it on a
        cluster, especially for more demanding pipelines (like those including genome indexing
        and/or alignment).

        ______________________________________________________________________________________


        ARGUMENTS:


        ALL variables must be set in the "config.json" file.

        "config.json" is organized in the following hierarchical manner:

        ---
        {
          "run_processes": {
            ... boolean variables indicating whether each process should be run or not
            ("all" is to run all pipeline from first to last step (BAM sorting excluded)) ...
          },
          "data_paths": {
            ... path variables to data ...
          },
          "processes": {
            ...
            process_name: {
              ... variables specific to process (both SLURM and function call variables) ...
            },
            ...
          },
          "run_locally": ...
        }
        ---


        The following describes path variables and lists all processes each variable
        is required for. If at least one step of your pipeline is included in one list,
        the corresponding path variable NEEDS to be specified.

        - "index_dir": path to directory for genome index files, required by:
        1. genome indexing
        2. alignment

        - "fasta_file": complete path to fasta file with reference genome, required by:
        1. genome indexing

        - "annotation_file": complete path to GTF/GFF files, required by:
        1. genome indexing
        8. gene counts

        - "fastq_files": list of complete paths to input read files, required by:
        2. alignment

        - "bam_dir": path to directory to store output alignment files, required by:
        2. alignment
        3. BAM sorting
        4. remove duplicates
        5. BAM filtering
        7. BAM stats
        8. gene counts

        - "bam_files": list of complete paths to input alignment files, required by:
        3. BAM sorting
        4. remove duplicates
        5. BAM filtering
        6. BAM indexing
        7. BAM stats
        8. gene counts
        (NOTICE: this variable is NEVER needed when your pipeline contains alignment step)

        - "gene_counts_dir": path to directory to store gene counts files, required by:
        8. gene counts

        - "report_dir": path to directory to store produced reports and plots, required by:
        9. summarize results

        _____________________________________________________________________________________


        HOW TO RUN A PIPELINE:


        1. Prepare a directory with Singularity images for required functions if not already available.

        2. Edit the "config.json" file as follows:

        2a. Set variables in "run_processes" section to true for the processes you wish to execute
            (if you set "all": true, then all steps will be run regardless of following variables'
            values).

        2b. Configure "data_paths" to specify paths to your data. Remember to use lists for "fastq_files"
            and "bam_files". Each path should be complete, in particular:
            - for "fastq_files", only the common prefix of reads pair should be passed, i.e. one
              full path without the "_R#_001.fastq.gz" suffix (glob patterns are allowed);
            - for "bam_files", include complete file paths including extensions (glob patterns are
              allowed);
            - if you don't need a path variable set it to an empty string/list.

        2c. Set "run_locally" variable to true if you want to run the pipeline on your local machine
            (not recommended for most applications).

        2d. Customize settings for each process under the "processes" section in "config.json". Refer to
            your cluster's specifications for SLURM settings, especially for the "queue" variable.

        3. Run the pipeline using the command:
           ---
           nextflow run main.nf
           ---

        Ensure all dependencies are properly configured and accessible before running the pipeline.
        See "NOTES" section below for details.

        Also, be sure to have "main.nf", "nextflow.config" and "config.json" in the same directory.

        _____________________________________________________________________________________


        NOTES:


        - The alignment step is designed for paired-end reads; however, subsequent steps can
          process both paired and single-end reads.

        - Ensure complete paths are provided for "data_paths" variables in "config.json".

        - Fastq files should be in ".gz" format, matching the pattern "*_R#_001.fastq.gz" or
          "*_R#_001.fq.gz".

        - Input alignment files must always be in BAM format.

        - "bam_files" is a list of input BAM files, "bam_dir" is the directory where all BAM
          file produced by the pipeline will be stored. Hence, files in "bam_files" do not need
          to be located in "bam_dir".

        - File names should contain relevant experimental information (e.g., cell line, sample
          number) only after dots. If not, change dots to dashes or other separators. All
          information after dots will be lost in subsequent files.
          Examples:
          - invalid: "COV362-TREATED-replica1.Tot_S11.Aligned.sortedByCoord.out.bam"
          - valid: "COV362-TREATED-replica1-Tot_S11.Aligned.sortedByCoord.out.bam"

        - "bam_dir" directory must contain the following subdirectories:
          - "logs/": For log files needed for quality control.
          - "stats/": For statistics summaries from SAMtools and metrics reports from Picard.
          - "tabs/": For tab files produced by STAR in alignment.
    """
}

if (params.help) {
    helpMessage()
    exit 0
}


// check valid option for gene count algorithm
def validCountAlgos = ['featureCounts', 'HTSeq']

if (!(params.count_algo in validCountAlgos)) {
    throw new IllegalArgumentException("Invalid value for 'count_algo'. Allowed values are: ${validCountAlgos.join(', ')}")
}


process runGenomeIndexing {
    output:
    val true  // for state dependency

    script:
    """
    STAR \
    --runMode genomeGenerate \
    --genomeDir $params.index_dir \
    --genomeFastaFiles $params.fasta_file \
    --sjdbGTFfile $params.annotation_file \
    --limitGenomeGenerateRAM $params.max_RAM_indexing \
    --runThreadN $params.genome_indexing_nt
    """
}

process runAlignment {
    if (params.save_all_BAM || params.last_BAM_output == "alignment") {
        publishDir "${params.bam_dir}", mode: 'copy', pattern: "${bam}"
    }
    publishDir "${params.bam_dir}/logs/", mode: 'copy', pattern: "*.Log.final.out"
    publishDir "${params.bam_dir}/tabs/", mode: 'copy', pattern: "*.tab"

    input:
    val ready  // for state dependency
    tuple path(fastq1), path(fastq2)

    output:
    path bam
    path '*.Log.final.out'
    path '*.tab'

    script:
    fastq_name = fastq1.toString().split("\\.")[0]
    core_name = fastq_name.substring(0, fastq_name.length() - 7)
    bam = core_name + ".Aligned.out.bam"

    """
    STAR \
    --runMode alignReads \
    --readFilesCommand "gunzip -c" \
    --genomeDir $params.index_dir \
    --readFilesIn ${fastq1} ${fastq2} \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix "${core_name}." \
    --quantMode GeneCounts \
    --runThreadN $params.alignment_nt
    """
}

process runBAMSorting {
    if (params.save_all_BAM || params.last_BAM_output == "sorting") {
        publishDir "${params.bam_dir}", mode: 'copy'
    }

    input:
    path bam

    output:
    path bam_sorted

    script:
    bam_sorted = bam.toString().split("\\.")[0] + ".Aligned.sortedByCoord.bam"

    """
    if samtools view -H ${bam} | grep -q '@HD.*SO:coordinate'; then
        cp ${bam} ${bam_sorted}
    else
        samtools sort -@ $params.BAM_sorting_nt -o ${bam_sorted} ${bam}
    fi
    """
}

process runRemoveDuplicates {
    if (params.save_all_BAM || params.last_BAM_output == "duplicates") {
        publishDir "${params.bam_dir}", mode: 'copy'
    }

    input:
    path bam

    output:
    path bam_marked

    script:
    bam_marked = bam.toString().split("\\.")[0] + ".Aligned.noDuplicates.bam"
    metrics = bam.toString().split("\\.")[0] + ".dup_metrics.txt"

    """
    picard MarkDuplicates \
    --INPUT ${bam} \
    --OUTPUT ${bam_marked} \
    --REMOVE_DUPLICATES true \
    --METRICS_FILE "${params.bam_dir}/stats/${metrics}"
    """
}

process runBAMFiltering {
    if (params.save_all_BAM || params.last_BAM_output == "filtering") {
        publishDir "${params.bam_dir}", mode: 'copy'
    }

    input:
    path bam

    output:
    path bam_filtered

    script:
    bam_filtered = bam.toString().split("\\.")[0] + ".Aligned.filtered.bam"

    """
    samtools view --threads $params.BAM_filtering_nt -b -q $params.quality_thres ${bam} > ${bam_filtered}
    """
}

process runBAMIndexing {
    publishDir "${params.bam_dir}", mode: 'copy'

    input:
    path bam

    output:
    val true  // for state depencency
    path bam_index

    script:
    bam_index = bam.toString() + ".bai"

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

    samtools flagstat -@ $params.BAM_stats_nt ${bam} > "$params.bam_dir/stats/\${bam_name}.stats.txt"
    """
}

process runFeatureCounts {
    input:
    path bam

    output:
    val true  // for state depencency

    script:
    """
    bam_name=\$(basename "${bam}")
    core_name="\${bam_name%%.*}"

    featureCounts \
    -a $params.annotation_file \
    -o "$params.gene_counts_dir/\${core_name}.counts.txt" \
    -T $params.gene_counts_nt \
    ${bam}
    """
}

process runHTSeq {
    input:
    path bam

    output:
    val true  // for state depencency

    script:
    """
    bam_name=\$(basename "${bam}")
    core_name="\${bam_name%%.*}"

    HTSeq-count \
    --format="bam" \
    ${bam} \
    $params.gene_counts_nt \
    > "$params.gene_counts_dir/\${core_name}.counts.txt"
    """
}

process runSumResults {
    input:
    val ready  // for state dependency

    script:
    """
    # set available directories
    dirs=()
    if [[ -n "$params.index_dir" ]]; then
        dirs+=("$params.index_dir")
    fi
    if [[ -n "$params.bam_dir" ]]; then
        dirs+=("$params.bam_dir")
    fi
    if [[ -n "$params.gene_counts_dir" ]]; then
        dirs+=("$params.gene_counts_dir")
    fi

    # run multiQC if any directory available
    if [[ -z "\${dirs[@]}" ]]; then
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
    if (params.run_genome_indexing || params.run_all) {

        index_ready = runGenomeIndexing()
    }


    // run reads alignment
    def bam_ch = false
    if (params.run_alignment || params.run_all) {

        // extract complete fastq files
        def fastq_files_complete = params.fastq_files.collect{ path -> return (path.toString() + "_R{1,2}_001.f*q.gz") }

        def fastq_ch = channel.fromFilePairs(fastq_files_complete, checkIfExists: true).map{baseName, fileList -> fileList}

        bam_ch = runAlignment(ready: index_ready, fastq_ch)[0]

    } else if (params.run_BAM_sorting || params.run_remove_duplicates || params.run_BAM_filtering || params.run_BAM_indexing || params.run_BAM_stats || params.run_gene_counts){

        bam_ch = channel.fromPath(params.bam_files, checkIfExists: true)
    }


    // run BAM sorting
    def bam_ch_sorted = false
    if (params.run_BAM_sorting || params.run_all) {

        bam_ch_sorted = runBAMSorting(bam_ch)

    } else {

        bam_ch_sorted = bam_ch
    }


    // run remove duplicates of BAM
    def bam_ch_marked = false
    if (params.run_remove_duplicates || params.run_all) {

        bam_ch_marked = runRemoveDuplicates(bam_ch)

    } else {

        bam_ch_marked = bam_ch
    }


    // run alignment filtering
    def bam_ch_filtered = false
    if (params.run_BAM_filtering || params.run_all) {

        bam_ch_filtered = runBAMFiltering(bam_ch_marked)

    } else {

        bam_ch_filtered = bam_ch_marked
    }


    // run BAM files indexing
    def BAMindex_ch = false
    def BAMindex_ready = false
    if (params.run_BAM_indexing || params.run_all) {

        BAMindex_ch = runBAMIndexing(bam_ch_filtered)[0]

        BAMindex_ready = BAMindex_ch.collect()
    }


    // run alignment stats summary
    def BAMstats_ch = false
    def BAMstats_ready = false
    if (params.run_BAM_stats || params.run_all) {

        BAMstats_ch = runBAMStats(bam_ch_filtered)

        BAMstats_ready = BAMstats_ch.collect()
    }


    // run gene counts
    def counts_ready = false
    if (params.run_gene_counts || params.run_all) {

        // run featureCounts
        if (params.count_algo == 'featureCounts') {

            counts_ready = runFeatureCounts(bam_ch_filtered)

        // run HTSeq
        } else if (params.count_algo == 'HTSeq') {

            counts_ready = runHTSeq(bam_ch_filtered)
        }
    }


    // run results summary
    if (params.run_summarize_results || params.run_all) {

        // make sure all data have been processed before going on
        if (!BAMindex_ready) {
            BAMindex_ready = channel.of(1, 2)  // "dummy" channel
        }

        if (!BAMstats_ready) {
            BAMstats_ready = channel.of(1, 2)  // "dummy" channel
        }

        if (!counts_ready) {
            counts_ready = channel.of(1, 2)  // "dummy" channel
        }

        def ready_for_sum = Channel.mix(BAMindex_ready, BAMstats_ready, counts_ready)

        runSumResults(ready: ready_for_sum)
    }
}
