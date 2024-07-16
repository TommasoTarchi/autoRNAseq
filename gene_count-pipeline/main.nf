#! /usr/bin/env nextflow


// help message (can be displayed using --help option)
def helpMessage() {
    log.info """
        USAGE:
           This pipeline processes paired-end FastQ files to produce aligned BAM files and gene expression counts. Below is a comprehensive guide to set up and run the pipeline.
           Any combination of its steps can be run.

        PIPELINE STEPS:
           This pipeline performs the following steps:
           - Genome Indexing: Preprocess the genome for alignment (using STAR).
           - FastQ Trimming: Trim reads based on quality scores and adapters (using Trim Galore!).
           - Alignment: Align reads to the reference genome (using STAR).
           - BAM Sorting: Sort alignment files by coordinates (using SAMtools).
           - Remove Duplicates: Remove or mark duplicate reads in BAM files (using Picard).
           - BAM Filtering: Filter aligned reads by MAPQ score (using SAMtools).
           - BAM Indexing: Index BAM files for fast retrieval (using SAMtools).
           - BAM Stats: Generate statistics from BAM files (using SAMtools).
           - Gene Counts: Quantify gene expression from BAM files (using either featureCounts or HTSeq).
           - Results Summary: Generate a summary report of pipeline results (using multiQC).

        REQUIREMENTS:
           - Nextflow and Singularity must be installed on your machine. Refer to their documentation for installation instructions.
           - Download or build the required container images for each pipeline step:
             * STAR v2.7.11b
             * Trim Galore! v0.6.7
             * SAMtools v1.3.1
             * Picard v3.1.1
             * featureCounts v2.0.6
             * HTSeq v2.0.2
             * multiQC v1.18

           Use wget or curl to download the container images if operating from the command line:
           ```
           wget <url_to_container_image> -O /path/to/your/container/image
           ```

        PARAMETERS:
           - All parameters are set in the `config.json` file. Avoid modifying `main.nf` or `nextflow.config`.
           - `config.json` is organized into:
             * run_processes: Boolean flags to enable or disable each process.
             * data_paths: Paths to input and output data.
             * processes: Specific parameters for each process, including SLURM and function call variables.
             * run_locally: Boolean to specify if the pipeline runs locally.
             * save_all_bams: Boolean to save intermediate BAM files.

            1. Data Paths:
               - Set the following paths in `data_paths` according to the steps you intend to run:
                 * index_dir: Directory for genome index files (required for genome indexing and alignment).
                 * fasta_file: Path to the reference genome fasta file (required for genome indexing).
                 * annotation_file: Path to the GTF/GFF file (required for genome indexing and gene counts).
                 * fastq_files: List of paths to input FastQ files (required for FastQ trimming and alignment).
                 * trimmed_fastq_dir: Directory for trimmed FastQ files (required for FastQ trimming).
                 * bam_dir: Directory for output BAM files (required for alignment, BAM sorting, removing duplicates, filtering, stats, and gene counts).
                 * bam_files: List of paths to input BAM files (required for BAM sorting, removing duplicates, filtering, indexing, stats, and gene counts).
                 * gene_counts_dir: Directory for gene count files (required for gene counts).
                 * report_dir: Directory for summary reports (required for results summary).

            2. Process Specific Parameters:
               - Each process has its own parameters:
                 * Common parameters: queue, time, memory, container_path, num_threads.
                 * Specific parameters:
                   - genome_indexing: max_RAM
                   - fastq_trimming: quality_thres, min_length, multithreaded
                   - remove_duplicates: remove_seq_duplicates
                   - BAM_filtering: quality_thres
                   - gene_counts: algo, strandedness

            3. Example of Input FastQ Files:
               - Ensure paired-end FastQ files are named appropriately (_R1_001.fastq.gz, _R2_001.fastq.gz).
               - Example:
                 ```
                 TREATED-replica1-S11_R1_001.fastq.gz
                 TREATED-replica1-S11_R2_001.fastq.gz
                 ```
               - Set fastq_files in config.json to include the common prefix of read pairs without the suffix.

        OUTPUT FILES:
           - Each step produces specific output files:
             * Genome Indexing: Indexed genome files in index_dir.
             * FastQ Trimming: Trimmed FastQ files in trimmed_fastq_dir.
             * Alignment: BAM files in bam_dir with suffix .Aligned.out.bam.
             * BAM Sorting: Sorted BAM files with suffix .Aligned.sortedByCoord.bam.
             * Remove Duplicates: BAM files with suffix .Aligned.marked.bam and duplicate metrics report.
             * BAM Filtering: Filtered BAM files with suffix .Aligned.filtered.bam.
             * BAM Indexing: Index files (.bai) in bam_dir.
             * BAM Stats: Statistics summary in bam_dir/stats/.
             * Gene Counts: Gene expression counts in gene_counts_dir.
             * Results Summary: HTML reports in report_dir.

        8. How to Run Your Pipeline:
           1. Clone the repository:
              ```
              git clone git@github.com:TommasoTarchi/autoRNAseq.git
              ```
           2. Navigate to the gene_count-pipeline directory.
           3. Edit `config.json`:
              - Set `run_processes` to true for desired steps.
              - Configure `data_paths` with correct paths.
              - Customize process parameters under `processes`.
              - Adjust `run_locally` and `save_all_bams`.
           4. Run the pipeline:
              ```
              nextflow run main.nf
              ```

        DOCS:
           For detailed instructions and troubleshooting, refer to the full README file provided with this pipeline.
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

process runTrimming {
    publishDir "${params.trimmed_fastq_dir}", mode: 'copy', pattern: "*.fq.gz"
    publishDir "${params.trimmed_fastq_dir}/reports/", mode: 'copy', pattern: "*.txt"
    publishDir "${params.trimmed_fastq_dir}/reports/", mode: 'copy', pattern: "*.html"
    publishDir "${params.trimmed_fastq_dir}/reports/", mode: 'copy', pattern: "*.zip"

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
    publishDir "${params.bam_dir}", mode: 'copy', pattern: "${bam}", enabled: {params.save_all_BAM || params.last_BAM_output == "alignment"}
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
    publishDir "${params.bam_dir}", mode: 'copy', enabled: {params.save_all_BAM || params.last_BAM_output == "sorting"}

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
    publishDir "${params.bam_dir}", mode: 'copy', enabled: {params.save_all_BAM || params.last_BAM_output == "duplicates"}

    input:
    path bam

    output:
    path bam_marked

    script:
    bam_marked = bam.toString().split("\\.")[0] + ".Aligned.marked.bam"
    metrics = bam.toString().split("\\.")[0] + ".dup_metrics.txt"

    """
    picard MarkDuplicates \
    --INPUT ${bam} \
    --OUTPUT ${bam_marked} \
    --REMOVE_SEQUENCING_DUPLICATES $params.remove_seq_duplicates \
    --METRICS_FILE "${params.bam_dir}/stats/${metrics}"
    """
}

process runBAMFiltering {
    publishDir "${params.bam_dir}", mode: 'copy', enabled: {params.save_all_BAM || params.last_BAM_output == "filtering"}

    input:
    path bam

    output:
    path bam_filtered

    script:
    bam_filtered = bam.toString().split("\\.")[0] + ".Aligned.filtered.bam"

    """
    samtools view --threads $params.BAM_filtering_nt -b -q $params.BAM_quality_thres ${bam} > ${bam_filtered}
    """
}

process runBAMIndexing {
    publishDir "${params.bam_dir}", mode: 'copy', pattern: "${bai}"

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

    samtools flagstat -@ $params.BAM_stats_nt ${bam} > "$params.bam_dir/stats/\${bam_name}.stats.txt"
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

process runSumResults {
    input:
    val ready1
    val ready2
    val ready3

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
    if [[ -n "$params.bam_dir" ]]; then
        dirs+=("$params.bam_dir")
    fi
    if [[ -n "$params.gene_counts_dir" ]]; then
        dirs+=("$params.gene_counts_dir")
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
    if (params.run_genome_indexing || params.run_all) {

        index_ready = runGenomeIndexing()
    }


    // run trimming and fastQC reports
    def fastq_ch_trimmed = false
    def trimming_ready = false  // for multiQC
    if (params.run_trimming || params.run_all) {

        // extract complete fastq files and define channel
        def fastq_files_complete = params.fastq_files.collect{ path -> return (path.toString() + "_R{1,2}_001.f*q.gz") }
        def fastq_ch = channel.fromFilePairs(fastq_files_complete, checkIfExists: true).map{baseName, fileList -> fileList}

        fastq_ch_trimmed = runTrimming(fastq_ch)[0]

	// signal to multiQC that this process is done
	trimming_ready = fastq_ch_trimmed.collect()

    } else if (params.run_alignment) {

        // extract complete fastq files and define channel
        def fastq_files_complete = params.fastq_files.collect{ path -> return (path.toString() + "_R{1,2}_001.f*q.gz") }
        fastq_ch_trimmed = channel.fromFilePairs(fastq_files_complete, checkIfExists: true).map{baseName, fileList -> fileList}
    }


    // run reads alignment
    def bam_ch = false
    if (params.run_alignment || params.run_all) {

        bam_ch = runAlignment(index_ready, fastq_ch_trimmed)[0]

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

        bam_ch_marked = runRemoveDuplicates(bam_ch_sorted)

    } else {

        bam_ch_marked = bam_ch_sorted
    }


    // run alignment filtering
    def bam_ch_filtered = false
    if (params.run_BAM_filtering || params.run_all) {

        bam_ch_filtered = runBAMFiltering(bam_ch_marked)

    } else {

        bam_ch_filtered = bam_ch_marked
    }


    // run BAM files indexing
    def bam_ch_indexed = false
    if (params.run_BAM_indexing || params.run_all) {

        bam_ch_indexed = runBAMIndexing(bam_ch_filtered)
    }


    // run alignment stats summary
    def bam_stats_ch = false
    def bam_stats_ready = false
    if (params.run_BAM_stats || params.run_all) {

        bam_stats_ch = runBAMStats(bam_ch_filtered)

        bam_stats_ready = bam_stats_ch.collect()
    }


    // run gene counts
    def counts_ready = false
    if (params.run_gene_counts || params.run_all) {

        // if indexing not run, extract BAM and corresponding index file and define channel
        if (!bam_ch_indexed) {
            def bam_bai_pairs = params.bam_files.collect{ path -> return (path.toString() + "{,.bai}") }
            bam_ch_indexed = channel.fromFilePairs(bam_bai_pairs, checkIfExists: true).map{baseName, fileList -> fileList}
        }

        // run featureCounts
        if (params.count_algo == 'featureCounts') {

            counts_ready = runFeatureCounts(bam_ch_indexed)

        // run HTSeq
        } else if (params.count_algo == 'HTSeq') {

            counts_ready = runHTSeq(bam_ch_indexed)
        }
    }


    // run results summary
    if (params.run_summarize_results || params.run_all) {

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

        runSumResults(trimming_ready, bam_stats_ready, counts_ready)
    }
}
