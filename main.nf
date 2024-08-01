#! /usr/bin/env nextflow


import java.io.File


//
// import processes from definition files
//
include { runGenomeIndexing } from './modules/GenomeIndexing.nf'
include { runTrimming } from './modules/Trimming.nf'
include { runAlignment } from './modules/Alignment.nf'
include { runBAMSorting } from './modules/BAMSorting.nf'
include { runRemoveDuplicates } from './modules/RemoveDuplicates.nf'
include { runBAMFiltering } from './modules/BAMFiltering.nf'
include { runBAMIndexing } from './modules/BAMIndexing.nf'
include { runBAMStats } from './modules/BAMStats.nf'
include { runFeatureCounts } from './modules/FeatureCounts.nf'
include { runHTSeq } from './modules/HTSeq.nf'
include { runSplicing } from './modules/Splicing.nf'
include { runSumResults } from './modules/SumResults.nf'


//
// help message (can be displayed using --help option)
//
def helpMessage() {
    log.info """
    """
}

if (params.help) {
    helpMessage()
    exit 0
}


//
// check valid options for function call parameters
//
def validCountAlgos = ['featureCounts', 'HTSeq']
def validStrandedness = [0, 1, 2]

if (!(params.count_algo in validCountAlgos)) {
    throw new IllegalArgumentException("Invalid value for 'count_algo'. Allowed values are: ${validCountAlgos.join(', ')}")
}
if (!(params.gc_strandedness in validStrandedness)) {
    throw new IllegalArgumentException("Invalid value for 'strandedness' in gene counts. Allowed values are: ${validStrandedness.join(', ')}")
}
if (!(params.spl_strandedness in validStrandedness)) {
    throw new IllegalArgumentException("Invalid value for 'strandedness' in splicing analysis. Allowed values are: ${validStrandedness.join(', ')}")
}


//
// check that needed data path variables were given and are valid
//
def checkPath(path, file_name) {
    def input = new File(path)
    
    if (!input.exists()) {
        println "ERROR: data_path variable '${file_name}' does not exist or was not provided"
	    System.exit(1)
    }
}

if (params.run_genome_indexing) {
    checkPath(params.index_dir, "index_dir")
    checkPath(params.fasta_file, "fasta_file")
    checkPath(params.annotation_file, "annotation_file")
}
if (params.run_trimming){
    checkPath(params.trimmed_fastq_dir, "trimmed_fastq_dir")
}
if (params.run_alignment){
    checkPath(params.index_dir, "index_dir")
    checkPath(params.out_bam_dir, "out_bam_dir")
}
if (params.run_BAM_sorting){
    checkPath(params.out_bam_dir, "out_bam_dir")
}
if (params.run_remove_duplicates){
    checkPath(params.out_bam_dir, "out_bam_dir")
}
if (params.run_BAM_filtering){
    checkPath(params.out_bam_dir, "out_bam_dir")
}
if (params.run_BAM_indexing){
    checkPath(params.out_bam_dir, "out_bam_dir")
}
if (params.run_BAM_stats){
    checkPath(params.out_bam_dir, "out_bam_dir")
}
if (params.run_gene_counts){
    checkPath(params.annotation_file, "annotation_file")
    checkPath(params.gene_counts_dir, "gene_counts_dir")
}
if (params.run_splicing){
    checkPath(params.annotation_file, "annotation_file")
    checkPath(params.splicing_dir, "splicing_dir")
}
if (params.run_summarize_results){
    checkPath(params.report_dir, "report_dir")
}


//
// pipeline
//
workflow {

    //
    // run genome indexing
    //
    def index_ready = true  // any value would be fine (just for state dependency)
    if (params.run_genome_indexing) {

        index_ready = runGenomeIndexing()[0]
    }


    //
    // run trimming and fastQC reports
    //
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


    //
    // run reads alignment
    //
    def bam_ch = false
    if (params.run_alignment) {

        bam_ch = runAlignment(index_ready, fastq_ch_trimmed)[0]

    } else if (params.run_BAM_sorting || params.run_remove_duplicates || params.run_BAM_filtering || params.run_BAM_indexing || params.run_BAM_stats || params.run_gene_counts){

        bam_ch = channel.fromPath(params.bam_files, checkIfExists: true)
    }


    //
    // run BAM sorting
    //
    def bam_ch_sorted = false
    if (params.run_BAM_sorting) {

        bam_ch_sorted = runBAMSorting(bam_ch)

    } else {

        bam_ch_sorted = bam_ch
    }


    //
    // run remove duplicates of BAM
    //
    def bam_ch_marked = false
    if (params.run_remove_duplicates) {

        bam_ch_marked = runRemoveDuplicates(bam_ch_sorted)

    } else {

        bam_ch_marked = bam_ch_sorted
    }


    //
    // run alignment filtering
    //
    def bam_ch_filtered = false
    if (params.run_BAM_filtering) {

        bam_ch_filtered = runBAMFiltering(bam_ch_marked)

    } else {

        bam_ch_filtered = bam_ch_marked
    }


    //
    // run BAM files indexing
    //
    def bam_ch_indexed = false
    if (params.run_BAM_indexing) {

        bam_ch_indexed = runBAMIndexing(bam_ch_filtered)
    }


    //
    // run alignment stats summary
    //
    def bam_stats_ch = false
    def bam_stats_ready = false
    if (params.run_BAM_stats) {

        bam_stats_ch = runBAMStats(bam_ch_filtered)

        bam_stats_ready = bam_stats_ch.collect()
    }


    //
    // build BAM-BAI channels if BAM indexing was not run
    //
    if (params.run_gene_counts || params.run_splicing) {

        // if indexing not run, extract BAM and corresponding index file and define channel
        if (!bam_ch_indexed) {
            def bam_bai_pairs = params.bam_files.collect{ path -> return (path.toString() + "{,.bai}") }
            bam_ch_indexed = channel.fromFilePairs(bam_bai_pairs, checkIfExists: true).map{baseName, fileList -> fileList}
        }
    }


    //
    // run gene counts
    //
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


    //
    // run splicing analysis
    //
    def splicing_ready = false
    if (params.run_splicing) {

        // gather all BAM-BAI couples for running splicing
        bam_ch_indexed.multiMap { pair ->
            bam: pair[0]
            bai: pair[1]
        }.set { bam_bai_list }

        def bam_list = bam_bai_list.bam.collect()
        def bai_list = bam_bai_list.bai.collect()

        // run proper analysis
        splicing_ready = runSplicing(bam_list, bai_list)[0]
    }


    //
    // run results summary
    //
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
