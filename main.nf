#! /usr/bin/env nextflow


// help message (can be displayed using --help option)
def helpMessage() {
    log.info """
        USAGE:


	This program can be used to run a RNAseq analysis pipeline, consisting of the following
	steps (also called "processes" in the following):

	1. Genome Indexing: Preprocess the genome for alignment.
	2. Alignment: Properly align reads to the reference genome.
	3. BAM Sorting: Sort BAM files (redundant if the previous step was run).
	4. Duplicates Marking: Mark, but DO NOT remove, duplicates in BAM files.
	5. BAM Filtering: Quality filtering of aligned reads.
	6. BAM Indexing: Index the alignment files.
	7. BAM Stats: Generate a statistical summary of the alignment.
	8. Gene Counts: Count the genes.
	9. Results Summary: Summarize the results.


	Any step of the pipeline can be run (multiple steps can be run together), see
	'ARGUMENTS' and 'HOW TO RUN A PIPELINE' below for instructions.


	The program is designed to bu run on a cluster using SLURM as job scheduler: each
	step is run as an independent (and parallel, when no data dependency occurs)
	multithreaded job.

	Optionally, it can also be run locally. However, for most demanding pipelines (in
	particular for pipelines including genome indexing and alignment steps) we suggest to
	use a cluster).

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
	  "additional_bindings": ...,
	  "run_locally": ...,
	  "processes": {
	    ...
	    process_name: {
	      ... variables specific to process (both SLURM and function call variables) ...
	    },
	    ...
	  }
	}
	---


	The following describes path variables and lists all processes each variable
	is required for. If at least one step of your pipeline is included in one list,
	the corresponding path variable NEED to be specified.

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
	  4. mark duplicates
	  5. BAM filtering
	  7. BAM stats
	  8. gene counts

	- "bam_files": list of complete paths to input alignment files, required by:
	  3. BAM sorting
	  4. mark duplicates
	  5. BAM filtering
	  6. BAM indexing
	  7. BAM stats
	  8. gene counts

	- "gene_counts_dir": path to directory to store gene counts files, required by:
	  8. gene counts

	- "data_dirs": list of paths to directories containg data for reports, required by:
	  9. summarize results

	- "report_dir": path to directory to store produced reports and plots, required by:
	  9. summarize results

	_____________________________________________________________________________________


	HOW TO RUN A PIPELINE:


	1. Prepare Singularity images for required functions if not already available.

	2. Edit the "config.json" file as follows:

	   - Set "run_processes" to true for the processes you wish to execute.

	   - Configure "data_paths" to specify paths to your data. Remember:
	     - Use lists for "fastq_files" and "bam_files". Each path should be complete and include
	       necessary glob patterns.
	     - For "fastq_files", provide one path per fastq pair without the "*_R#_001.fastq.gz" suffix.
	     - For "bam_files", include complete file paths including extensions.

	   - If files are located outside the current directory or its subdirectories, specify them in
	     the "additional_bindings" variable as a comma-separated list of full paths (for instance:
	     "additional_bindings": "/path/to/dir1,/path/to/dir2")

	   - Set "run_locally" variable to true if you want to run the pipeline on your local machine

	   - Customize settings for each process under the "processes" section in "config.json". Refer to
	     your cluster's specifications for SLURM settings, especially for the "queue" variable.

	4. Run the pipeline on your cluster using the command:
	   ---
	   nextflow run main.nf
	   ---

	Ensure all dependencies are properly configured and accessible before running the pipeline.

        _____________________________________________________________________________________


	NOTES:
	

	- Ensure complete paths are provided for "data_paths" variables in "config.json".

	- "bam_files" is a list of input BAM files, "bam_dir" is the directory where all BAM
	  file produced by the pipeline will be stored. Hence, files in "bam_files" do not need
	  to be located in "bam_dir".
	
	- File names should contain relevant experimental information (e.g., cell line, sample
	  number) only after dots. If not, change dots to dashes or other separators. All
	  information after dots will be lost in subsequent files.
	  Examples:
	  - invalid: "COV362-TREATED-replica1.Tot_S11.Aligned.sortedByCoord.out.bam"
	  - valid: "COV362-TREATED-replica1-Tot_S11.Aligned.sortedByCoord.out.bam"

	- The alignment step is designed for paired-end reads; however, subsequent steps can
	   process both paired and single-end reads.

	- Fastq files should be in ".gz" format, matching the pattern "*_R#_001.fastq.gz" or
	  "*_R#_001.fq.gz".

	- Alignment files must be in BAM format if the alignment step is skipped.

	- Mark duplicates step overwrites BAM files without data loss, if "bam_files" are contained
	  "bam_dir".

	- The BAM directory must contain the following subdirectories:
	  - "logs/": For log files needed for quality control.
	  - "stats/": For statistics summaries from SAMtools and metrics reports from Picard.
	  - "tabs/": For tab files produced by STAR in alignment.
	
	- If running the "gene counts" step alone, ensure only one BAM file per sample is present.
	  Remove duplicate BAM files (e.g., filtered and unfiltered reads from the same sample).
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
    publishDir "${params.bam_dir}", mode: 'copy', overwrite: false, pattern: "${bam}"
    publishDir "${params.bam_dir}/logs/", mode: 'copy', overwrite: false, pattern: "*.Log.final.out"
    publishDir "${params.bam_dir}/tabs/", mode: 'copy', overwrite: false, pattern: "*.tab"

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
    bam = core_name + ".Aligned.sortedByCoord.out.bam"

    """
    STAR \
    --runMode alignReads \
    --readFilesCommand "gunzip -c" \
    --genomeDir $params.index_dir \
    --readFilesIn ${fastq1} ${fastq2} \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix "${core_name}." \
    --quantMode GeneCounts \
    --runThreadN $params.alignment_nt
    """
}

process runBAMSorting {
    publishDir "${params.bam_dir}", mode: 'copy', overwrite: false

    input:
    path bam

    output:
    path bam_sorted

    script:
    bam_sorted = bam.toString().split("\\.")[0] + ".Aligned.sortedByCoord.out.bam"
 
    """
    samtools sort -@ $params.BAM_sorting_nt -o ${bam_sorted} ${bam}
    """
}

process runMarkDuplicates {
    publishDir "${params.bam_dir}", mode: 'copy', overwrite: false

    input:
    path bam

    output:
    path bam_marked

    script:
    bam_marked = bam.toString()  // just copy the name (we overwrite the BAM)
    metrics = bam.toString().split("\\.")[0] + ".dup_metrics.txt"

    """
    picard MarkDuplicates \
    I=${bam} \
    O=temp.bam \
    M="${params.bam_dir}/stats/${metrics}"

    mv temp.bam ${bam_marked}
    """
}

process runBAMFiltering {
    publishDir "${params.bam_dir}", mode: 'copy', overwrite: false

    input:
    path bam

    output:
    path bam_filtered

    script:
    bam_filtered = bam.toString().split("\\.")[0] + ".Aligned.sortedByCoord.filtered.bam"

    """
    samtools view --threads $params.BAM_filtering_nt -b -q $params.quality_thres ${bam} > ${bam_filtered}
    """
}

process runBAMIndexing {
    input:
    path bam

    output:
    val true  // for state depencency

    script:
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
    multiqc $params.data_dirs \
    --force \
    --export \
    --outdir $params.report_dir
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

        def fastq_ch = channel.fromFilePairs("$params.fastq_files_R{1,2}_001.f*q.gz", checkIfExists: true).map{baseName, fileList -> fileList}

	bam_ch = runAlignment(ready: index_ready, fastq_ch)[0]

    } else {

	bam_ch = channel.fromPath(params.bam_files, checkIfExists: true)
    }


    // run BAM sorting
    def bam_ch_sorted = false
    if (params.run_BAM_sorting) {

	bam_ch_sorted = runBAMSorting(bam_ch)
    
    } else {

	bam_ch_sorted = bam_ch
    }


    // run mark duplicates of BAM
    def bam_ch_marked = false
    if (params.run_mark_duplicates || params.run_all) {

	bam_ch_marked = runMarkDuplicates(bam_ch)

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

	BAMindex_ch = runBAMIndexing(bam_ch_filtered)

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
    if (params.run_sum_results || params.run_all) {

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
