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
