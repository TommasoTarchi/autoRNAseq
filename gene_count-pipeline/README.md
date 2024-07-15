# Gene expression count pipeline

This pipeline can be used to produce compressed alignment files (BAM) and gene expression count files from "raw" paired-end read files (FastQ).


## Outline

- [Pipeline steps](#pipeline-steps)
- [Requirements](#requirements)
- [Parameters description](#parameters-description)
  - [Data paths](#data-paths)
  - [Process specific parameters](#process-specific-parameters)
  - [Output files](#output-files)
  - [Example of input FastQ files](#example-of-input-fastq-files)
- [How to run your pipeline](#how-to-run-your-pipeline)


## Pipeline steps

This pipeline implements the following steps (between parentheses you have the functions used
for implementation):

1. **Genome Indexing**: preprocess the genome for alignment ([*STAR*][STAR]).
2. **FastQ trimming**: trim reads based on length, Phred score and adapters, and quality control ([*Trim Galore!*][trim_galore]).
3. **Alignment**: properly align reads to the reference genome ([*STAR*][STAR]).
4. **BAM Sorting**: sort alignment files ([*SAMtools*][SAMtools]).
5. **Remove duplicates**: remove (or mark only) duplicates in alignment files ([*picard*][picard]).
6. **BAM Filtering**: filter aligned reads by MAPQ quality score ([*SAMtools*][SAMtools]).
7. **BAM Indexing**: index the alignment files ([*SAMtools*][SAMtools]).
8. **BAM Stats**: generate a statistical summary of the alignment ([*SAMtools*][SAMtools]).
9. **Gene Counts**: quantify gene expression ([*featureCounts*][featureCounts] or [*HTSeq*][HTSeq]).
10. **Results Summary**: summarize the results ([*multiQC*][multiQC]).


## Requirements

- You need to have Nextflow and Singularity installed on your machine. You can look at the
related documentation [here][nextflow] and [here][singularity] for instructions.

- You need to have the containers on which the steps of the pipeline will run (**remember**
  that you only need the container images related to the steps you want to run). You can either
  build the containers by yourself or download them from the assets of this program's release.
  The following is a list of the available containers in assets:
  - STAR v2.7.11b ([https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/STAR-v2.7.11b.sif](https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/STAR-v2.7.11b.sif))
  - Trim Galore! v0.6.7 ([https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/trim_galore-v0.6.7.sif](https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/trim_galore-v0.6.7.sif))
  - SAMtools v1.3.1 ([https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/SAMtools-v1.3.1.sif](https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/SAMtools-v1.3.1.sif))
  - picard v3.1.1 ([https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/picard-v3.1.1.sif](https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/picard-v3.1.1.sif))
  - featureCounts v2.0.6 ([https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/featureCounts-v2.0.6.sif](https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/featureCounts-v2.0.6.sif))
  - HTSeq v2.0.2 ([https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/HTSeq-v2.0.2.sif](https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/HTSeq-v2.0.2.sif))
  - multiQC v1.18 ([https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/multiQC-v1.18.sif](https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/multiQC-v1.18.sif))
  
  If you are operating from command line, you can use [*wget*][wget] (or [*curl*][curl]) to download the images:
  ````
  $ wget <url_to_container_image> -O /path/to/your/container/image
  ````

- If your pipeline uses FastQ files, please make sure they are **paired-end**, **zipped**, and
ending in `_R#_001.fastq.gz` or `_R#_001.fq.gz`, with `#` equal to 1 for first read and 2 for second
read.

- Make sure that in all input files **all relevant information is placed after dots**. If this is
not the case, you can replace these dots with other seprators.
Example: `<info1>.<info2>.Aligned.out.bam` should be renamed something like
`<info1>_<info2>.Aligned.out.bam` if you don't want to lose `<info2>` from the final
output file name.

- If you want to use the pipeline only to remove or mark duplicates (without running previous steps), please
make sure that input BAM files are sorted (you can include the BAM sorting step if you are not sure).

- If your pipeline contains the gene count step, please include the BAM indexing step as well. If you
only want to run the gene count step, make sure your input BAM files are indexed.

- If you want to run the gene count step, you need to know the *strandedness* of your data. If you don't
know it, you can infer it from BAM files using [RSeQC infer_experiment.py](http://rseqc.sourceforge.net/#infer-experiment-py)
or from FastQ files using [how_are_we_stranded_here](https://github.com/signalbash/how_are_we_stranded_here).
**Notice** that the second option may not work because of some bug in output file processing. In this case
you can just manually check the output file `strandedness_check.txt` and interpret its content, remembering:
  - `1++,1--,2+-,2-+` refers to forward strandedness;
  - `1+-,1-+,2++,2--` refers to reverse strandedness.


## Parameters description

All parameters can be set from the `config.json` file. Please, do not modify neither
`main.nf` nor `nextflow.config`, unless you are familiar with nextflow.

`config.json` is organized according to the following hierarchical structure:

````
{
  "run_processes": {
    ... boolean variables indicating whether each process should be run or not
        ("all" is to run all pipeline from first to last step) ...
  },
  
  "data_paths": {
    ... path variables to input and output data ...
  },
  
  "processes": {
    
    "process_1": {
      ... variables specific to process 1 (both SLURM and function call variables) ...
    },
    
    "process_2": {
      ... variables specific to process 2 (both SLURM and function call variables) ...
    },
    
    ...
  
  },
  
  "run_locally" -> boolean: whether the pipeline should be run locally

  "save_all_bams" -> boolean: whether output BAM files should be saved at each step
}
````


### Data paths

Particular care must be taken in setting the `data_path` variables.

If you run only some of the pipeline steps (as it is usually the case), you will only need some of
these variables. The following list shows for each data path variable which steps of the pipeline need
it to be set. If **at least one** of the steps you intend to run is listed for a variable, then you
need to set that variable.

- `index_dir`: path to directory for genome index files. Required by: 1. genome indexing, 3. alignment; **optional**: 10. summarize results.

- `fasta_file`: complete path to fasta file with reference genome. Required by: 1. genome indexing.

- `annotation_file`: complete path to GTF/GFF file. Required by: 1. genome indexing, 9. gene counts.

- `fastq_files`: list of complete paths to input (**zipped**) read files. Required by: 2. FastQ trimmming, 3. alignment.

- `trimmed_fastq_dir`: path to directory to store trimmed read files. Required by: 2. FastQ trimmming; **optional**: 10. summarize results.

- `bam_dir`: path to directory to store output alignment files. Required by: 3. alignment, 4. BAM sorting,
  5. remove duplicates, 6. BAM filtering, 8. BAM stats, 9. gene counts; **optional**: 10. summarize results.

- `bam_files`: list of complete paths to input alignment files. Required by: 4. BAM sorting, 5. remove duplicates,
  6. BAM filtering, 7. BAM indexing, 8. BAM stats, 9. gene counts (**Notice**: this variable is **never** needed when
  your pipeline contains the alignment step, i.e. number 3.).

- `gene_counts_dir`: path to directory to store gene counts files. Required by: 9. gene counts; **optional**: 10. summarize results.

- `report_dir`: path to directory to store produced reports and plots. Required by: 10. summarize results.


### Process specific parameters

As mentioned previously, inside the `processes` scope, each process has its own scope for parameter setting.

Common to all processes are the following variables:

````
"queue" -> string: name of cluster partition to run job on
"time" -> string: maximum time for job (example: "2h")
"memory" -> string: RAM used for the job (example: "2GB")
"container_path" -> string: full path to singularity image for the process
"num_threads" -> integer: number of threads (excluded "FastQ trimming", "alignment" when "algo": "HTSeq", "remove_duplicates" and "summarize_results")
````

For `container_path` setting refer to the [requirements](#requirements).

Other process-specific parameteres are:

````
"genome_indexing": {
  "max_RAM" -> string: maximum RAM for STAR indexing **in bytes** (should be same amount as in "memory")
}

"fastq_trimming": {
  "quality_thres" -> integer: Phred threshold for quality filtering
  "min_length" -> integer: length threshold (in bp) for reads length filtering (a pair of reads is dropped if
                           at least one of them is below the threshold)
  "multithreaded" -> boolean: whether Trim Galore! should be multithreaded (the number of threads is fixed to
                              4 to make sure the function works correctly)
}

"remove_duplicates": {
  "remove_seq_duplicates" -> boolean: whether duplicates likely caused by sequencing process should be removed
                                      (if false, duplicates are only marked and not removed)
}

"BAM_filtering": {
  "quality_thres" -> integer: MAPQ threshold for quality filtering
}

"gene_counts": {
  "algo" -> string: algorithm for gene expression quantification (allowed options: "featureCounts","HTSeq")
  "strandedness" -> integer: 0 for non-stranded, 1 for forward-stranded, 2 for reverse-stranded
}
````

**All** process-specific parameters of the processes you intend to run must be set to some value.

Notice that in `config.json` some of the parameters are set to a default value. However, the value set is **not**
guaranteed to work.


### Output files

The following is a list of output files of each step of the pipeline.

If not specified, the output file is always saved, independently of which steps are run after. Exeption are
BAM files: by default only the last version is saved, but they can optionally be saved at each step by setting
the related variable in `config.json`.

1. Genome Indexing:
    - Various files with indexed and preprocessed genome, saved into `inded_dir`.

2. FastQ trimming:
    - Trimmed FastQ files, called as the original FastQ but with suffix `_val_#`, with `#` equal to
      1 and 2.
    - Trimming reports, saved as text files into "`trimmed_fastq_dir`/reports/".
    - *FastQC* reports, saved as html files into "`trimmed_fastq_dir`/reports/".

3. Alignment:
    - BAM files with alignment (one per paired-end pair of fastq file), **unsorted** and saved
    into `bam_dir`. The name of the file will have all the relevant information contained in the fastq
    names followed by the suffix: `.Aligned.out.bam`.
    - Alignment log files, saved into "`bam_dir`/logs/".
    - Alignment tab files, saved into "`bam_dir`/tabs/".

4. BAM Sorting:
    - BAM files sorted by coordinates, saved into `bam_dir` with extention `.Aligned.sortedByCoord.bam`.

5. Remove duplicates:
    - BAM files with duplicates removed, saved into `bam_dir` with extention `.Aligned.marked.bam`.
    - Duplicate metrics report, saved in "`bam_dir`/stats/" with extention `.dup_metrics.txt`.

6. BAM Filtering:
    - BAM files filtered according to some threshold on quality score, saved into `bam_dir` with extention
    `.Aligned.filtered.bam`.

7. BAM Indexing:
    - Index files of input BAM files (one for each BAM), saved into `bam_dir` with same name as input BAM plus
    `.bai`.

8. BAM Stats:
    - Statistics summary of input BAM files (one for each BAM), saved into "`bam_dir`/stats/" with extention
    `.stats.txt`.

9. Gene Counts:
    - Gene expression count files (one for each input BAM), saved into `gene_counts_dir` with extention
    `.counts.txt`.

10. Results Summary:
    - html reports of all steps run, saved into `report_dir`.


### Example of input FastQ files

Suppose you want to run reads alignment and suppose you have a directory containing the following files:

- `TREATED-replica1-S11_R1_001.fastq.gz`
- `TREATED-replica1-S11_R2_001.fastq.gz`
- `TREATED-replica2-S22_R1_001.fastq.gz`
- `TREATED-replica3-S34_R1_001.fq.gz`
- `TREATED-replica3-S34_R2_001.fq.gz`
- `TREATED-replica4-S25_R1_001.fastq`
- `TREATED-replica4-S25_R2_001.fastq`
- `RES_PT-replica1-S12_R1_001.fq.gz`
- `RES_PT-replica1-S12_R2_001.fq.gz`
- `RES_PT-replica2-S24_R1_001.fq.gz`
- `RES_PT-replica2-S24_R2_001.fq.gz`
- `RES_PT-replica3-S24_R1_001.fq.gz`
- `RES_PT-replica3-S24_R2_001.fq.gz`

Now, for instance if you set `fastq_files` to the list: [`TREATED-replica*`, `RES_PT-replica1-S12`,
`RES_PT-replica2-S24_R?_001.fq.gz`], the files in the directory will be treated in the following way:

- `TREATED-replica1-S11_R1_001.fastq.gz` and `TREATED-replica1-S11_R2_001.fastq.gz`: processed
- `TREATED-replica2-S22_R1_001.fastq.gz`: not processed, since it does not have a corresponding paired read file
(notice that the file matches one of the glob pattern passed).
- `TREATED-replica3-S34_R1_001.fq.gz` and `TREATED-replica3-S34_R2_001.fq.gz`: processed
- `TREATED-replica4-S25_R1_001.fastq` and `TREATED-replica4-S25_R2_001.fastq`: not processed, since
they do not match the expected format
- `RES_PT-replica1-S12_R1_001.fq.gz` and `RES_PT-replica1-S12_R2_001.fq.gz`: processed
- `RES_PT-replica2-S24_R1_001.fq.gz` and `RES_PT-replica2-S24_R2_001.fq.gz`: not processed, since the
provided pattern wrongly includes the suffix `_R#_001.fq.gz`
- `RES_PT-replica3-S24_R1_001.fq.gz` and `RES_PT-replica3-S24_R2_001.fq.gz`: not processed, since it does
not match any glob pattern passed.

The output of alignment will therefore be:

- `TREATED-replica1-S11.Aligned.out.bam`
- `TREATED-replica3-S34.Aligned.out.bam`
- `RES_PT-replica1-S12.Aligned.out.bam`


## How to run your pipeline

1. Make sure you satisfy all requirements listed in [this section](#requirements).

2. Clone this repository, using:
   ````
   $ git clone git@github.com:TommasoTarchi/autoRNAseq.git
   ````

3. Navigate to the `gene_count-pipeline` directory and edit the `config.json` file as follows:

    - Set variables in `run_processes` section to true for the processes you wish to execute
      (if you set `all` to true, then all steps will be run regardless of the values of the following
      variables).

    - Configure `data_paths` to specify paths to your data following the descriptions in [this section](#data-paths).
      Remember to use lists for `fastq_files` and `bam_files`. Each path should be complete, in particular:
      - for `fastq_files`, only the common prefix of reads pair should be passed, i.e. one full path per read pair
        **without** the `_R#_001.fastq.gz` suffix (**glob patterns are allowed** - see [here](#example-of-input-fastq-files)
        for more details);
      - `trimmed_fastq_dir` should contain a subdirectory called "reports/";
      - for `bam_files`, include full paths (**glob patterns are allowed**);
      Also:
      - `bam_dir` should contain three subdirectories called "logs/", "stats/" and "tabs/";
      - if you don't need a path variable set it to an empty string/list.
      
    - Customize settings for each process you are running under the corresponding `processes` section in `config.json`. Refer
      to your cluster's specifications for SLURM settings. See [here](#process-specific-parameters) for details. Parameters of
      processes you are not running will be ignored (you can leave them as they are).
    
    - Change `run_locally` variable to true if you want to run the pipeline on your local machine (**not recommended** for
      most applications). If you want to run it on a cluster, leave it to false.

    - Change `save_all_bams` variable to true if you want to keep all BAM files produced along the pipeline (usually **not
      recommended**, especially when working with many files). If you only want the final result, leave it to false.

5. Run the pipeline using:
   ````
   $ nextflow run main.nf
   ````







[nextflow]: https://www.nextflow.io/docs/latest/index.html
[singularity]: https://apptainer.org/user-docs/master/
[nextflow_and_singularity]: https://nextflow.io/docs/edge/container.html#singularity
[build_containers]: https://github.com/fburic/notes/blob/master/singularity_conda.md
[STAR]: https://docs.csc.fi/apps/star/
[trim_galore]: https://github.com/FelixKrueger/TrimGalore/tree/master/Docs
[SAMtools]: https://www.htslib.org/doc/samtools.html
[picard]: https://broadinstitute.github.io/picard/
[featureCounts]: https://subread.sourceforge.net/featureCounts.html
[HTSeq]: https://htseq.readthedocs.io/en/master/overview.html
[multiQC]: https://multiqc.info/docs/ 
[wget]: https://www.gnu.org/software/wget/manual/wget.html
[curl]: https://curl.se/docs/
