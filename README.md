# AutoRNAseq: an automated pipeline for paired-end RNAseq data analysis

This pipeline can be used to perform alignment, FastQ and BAM processing, gene expression count
and splicing analysis, for **paired-end reads**.

It is composed of several steps. The user can choose to run any combination of these steps with
restrictions defined in [Requirements](#requirements), and given that the correct arguments and containers
are provided.

The implementation uses *Nextflow* for automation and *Singularity* for containerization. If you
are new to these tools, you can find the documentation, respectively, [here][nextflow] and
[here][singularity]. In particular, it can be useful to consult [this section][nextflow_and_singularity]
for using Nextflow and Singularity together.

Each step of each pipeline runs on a different container, so that you can just download the
containers needed by the pipeline steps you intend to run. Container images ready to use can be found among
the assets of the [program release][assets].

Pipelines support both *SLURM* (default) and local run. However, the second one is not recommended for
most demanding pipelines. When run with SLURM, each process is launched as a different job on the cluster.
Resources parameters can be adjusted differently for each step of the pipeline.


## Table of contents

- [Pipeline steps](#pipeline-steps)
- [What you will find in this repository](#what-you-will-find-in-this-repository)
- [Requirements](#requirements)
- [Parameters description](#parameters-description)
  - [Input files](#input-files)
  - [Data paths](#data-paths)
  - [Process specific parameters](#process-specific-parameters)
  - [Output files](#output-files)
- [How to run your pipeline](#how-to-run-your-pipeline)


## Pipeline steps

This pipeline implements the following steps (between parentheses you have the functions used
for implementation):

1. **Genome Indexing**: preprocess the genome for alignment ([*STAR*][STAR]).
2. **FastQ trimming**: trim reads based on length, Phred score and adapters, and run quality control ([*Trim Galore!*][trim_galore]).
3. **Alignment**: properly align reads to the reference genome ([*STAR*][STAR]).
4. **BAM Sorting**: sort alignment files ([*SAMtools*][SAMtools]).
5. **Remove duplicates**: remove (or mark only) duplicates in alignment files ([*picard*][picard]).
6. **BAM Filtering**: filter aligned reads by MAPQ quality score ([*SAMtools*][SAMtools]).
7. **BAM Indexing**: index the alignment files ([*SAMtools*][SAMtools]).
8. **BAM Stats**: generate a statistical summary of the alignment ([*SAMtools*][SAMtools]).
9. **Gene Counts**: quantify gene expression ([*featureCounts*][featureCounts] or [*HTSeq*][HTSeq]).
10. **Splicing Analysis**: comparative splicing analysis between given conditions ([*rMATS-turbo*][rMATS-turbo]).
11. **Results Summary**: summarize the results ([*multiQC*][multiQC]).


## What you will find in this repository

- This README file: description of the pipeline and instructions to run it
- [`main.nf`](./main.nf): nextflow main
- [`nextflow.config`](./nextflow.config): nextflow configuration file for parameters
- [`modules/`](./modules/): directory containing definition files of all processes
- [`config.json`](./config.json): configuration file for user


## Requirements

- You need to have Nextflow and Singularity installed on your machine. You can look at the
related documentation [here][nextflow] and [here][singularity] for instructions.

- You need to have the container images on which the steps of the pipeline will run (**remember** that you only need
  the containers for the steps you want to run). The images can be downloaded from the assets of this program's release
  (make sure you place **all container images in the same directory**):
  - STAR v2.7.11b ([https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/STAR-v2.7.11b.sif](https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/STAR-v2.7.11b.sif))
  - Trim Galore! v0.6.7 ([https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/trim_galore-v0.6.7.sif](https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/trim_galore-v0.6.7.sif))
  - SAMtools v1.3.1 ([https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/SAMtools-v1.3.1.sif](https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/SAMtools-v1.3.1.sif))
  - picard v3.1.1 ([https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/picard-v3.1.1.sif](https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/picard-v3.1.1.sif))
  - featureCounts v2.0.6 ([https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/featureCounts-v2.0.6.sif](https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/featureCounts-v2.0.6.sif))
  - HTSeq v2.0.2 ([https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/HTSeq-v2.0.2.sif](https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/HTSeq-v2.0.2.sif))
  - rMATS v4.3.0 ([https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/rMATS-v4.3.0.sif](https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/rMATS-v4.3.0.sif))
  - multiQC v1.18 ([https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/multiQC-v1.18.sif](https://github.com/TommasoTarchi/autoRNAseq/releases/download/v0.1.0-alpha/multiQC-v1.18.sif))
  
  If you are operating from command line, you can use [*wget*][wget] (or [*curl*][curl]) to download the images:
  ````
  $ wget <url_to_container_image> -O /path/to/your/container/image
  ````

- If your pipeline uses FastQ files, please make sure they are **paired-end** and **zipped**.

- Make sure that in all input files **all relevant information is placed before dots**. If this is
not the case, you can replace these dots with other seprators.
Example: `<info1>.<info2>.Aligned.out.bam` should be renamed something like
`<info1>_<info2>.Aligned.out.bam` if you don't want to lose `<info2>` from the final
output file name. Also make sure that all input files have different names.

- If you want to use the pipeline only to remove or mark duplicates (without running previous steps), please
make sure that input BAM files are sorted (you can include the BAM sorting step if you are not sure).

- If your pipeline contains the gene count and/or the splicing analysis steps, please include the BAM indexing step
as well. If you only want to run the gene count step, make sure your input BAM files are indexed and that each index
file is placed in the same directory as the corresponding BAM file.

- At the moment,the pipeline only supports splicing analysis from BAM files. If you only have FastQ files and you want to
perform splicing analysis, please include (at least) steps 3,4 and 7 in the pipeline.

- If you want to run the gene count and/or splicing analysis steps, you need to know the *strandedness* of your data. If you don't
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
    ... boolean variables indicating whether each process should be run or not ...
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

  "container_dir" -> string: path to directory with container images

  "nf_work_dir" -> string: path to work directory for pipeline (default: "./work/")
  
  "run_locally" -> boolean: whether the pipeline should be run locally

  "save_all_bams" -> boolean: whether output BAM files should be saved at each step
}
````


### Input files

Input files must be passed through a txt file.

There are three possible scenarios:
1. If your pipeline contains only steps 1 and/or 11, then you don't need any input files.
2. If your pipeline starts at step 2 or 3, then your txt file should look like:
   ````
   /path/to/read1_of_sample1,/path/to/read2_of_sample1,condition_of_sample1
   /path/to/read1_of_sample2,/path/to/read2_of_sample2,condition_of_sample2
   ...
   ````
   i.e. each line should contain complete path to the FastQ with first read, complete path to the FastQ with
   second read and condition of the sample, in this order and comma-separated.
3. If your pipeline starts at step 4, 5, 6, 7, 8, 9, or 10 then your txt file should look like:
   ````
   /path/to/bam1,condition1
   /path/to/bam2,condition2
   ...
   ````
   i.e. each line should contain complete path to the BAM file and condition of the sample, in this order and comma-separated.

**Please check** that your txt file does not contain any empty lines, as they would most likely produce an error.


### Data paths

Particular care must be taken in setting the `data_path` variables.

If you run only some of the pipeline steps (as it is usually the case), you will only need some of
these variables. The following list shows for each data path variable which steps of the pipeline need
it to be set (for step reference numbers see [this section](#pipeline-steps)). If **at least one** of the steps you intend to run is listed for a variable, then you
need to set that variable.

- `input_list`: complete path to txt file containing input files described in the previous section.
  Required by steps: 2, 3, 4, 5, 6, 7, 8, 9, 10.

- `index_dir`: path to directory for genome index files. Required by steps: 1, 3.

- `fasta_file`: complete path to fasta file with reference genome. Required by steps: 1.

- `annotation_file`: complete path to GTF/GFF file. Required by steps: 1, 9, 10.

- `trimmed_fastq_dir`: path to directory to store trimmed read files. Required by steps: 2.

- `out_bam_dir`: path to directory to store output alignment files. Required by steps: 3, 4, 5, 6, 7, 8, 9.

- `gene_counts_dir`: path to directory to store gene counts files. Required by steps: 9.

- `splicing_dir`: path to directory to store results from splicing analysis. Required by steps: 10.

- `report_dir`: path to directory to store produced reports and plots. Required by steps: 11.

Please, make sure to use **absolute paths**.


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

"splicing_analysis": {
  "condition1" -> string: first condition to be compared (MUST correspond to at least one input BAM file)
  "condition2" -> string: second condition to be compared (MUST correspond to at least one input BAM file)
  "strandedness" -> integer: 0 for non-stranded, 1 for forward-stranded, 2 for reverse-stranded
  "read_length" -> integer: reads length (not all reads have to be of this length - rMATS-turbo is set to
                            to handle varying length reads; in this case a reasonable approach is to use
                            the length of reads before trimming)
  "use_paired_stats" -> boolean: whether to use paired stats model
  "cutoff_diff" -> float: cutoff difference used in hypothesis test for differential alternative splicing
                          (ignored if "use_paired_stats": false); example: 0.0001 for 0.01% difference
  "detect_novel_splice" -> boolean: whether to detect unannotated splice sites
  "min_intron_len" -> integer: minimum intron length (ignored if "detect_novel_splice": false)
  "max_exon_len" -> integer: maximum exon length (ignored if "detect_novel_splice": false)
}
````

**All** process-specific parameters of the processes you intend to run must be set to some value.

Notice that in `config.json` some of the parameters are set to a default value. However, the value set is **not**
guaranteed to work.


### Output files

The following is a list of the output files produced at each step of the pipeline.

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
    - BAM files with alignment (one per paired-end pair of fastq file), **unsorted** and saved into `bam_dir`.
      The name of the file will have all the information contained before dots in the name of the corresponding
      FastQ file with first reads followed by the suffix: `.Aligned.bam`.
    - Alignment log files, saved into "`bam_dir`/logs/".
    - Alignment tab files, saved into "`bam_dir`/tabs/".

4. BAM Sorting:
    - BAM files sorted by coordinates, saved into `bam_dir`. All files that have undergone this step will
      contain `sortedByCoord` in their name's suffix.

5. Remove duplicates:
    - BAM files with duplicates removed (or just marked), saved into `bam_dir`. All files that have undergone
      this step will contain `marked` in their name's suffix.
    - Duplicate metrics report, saved in "`bam_dir`/stats/" with extention `.dup_metrics.txt`.

6. BAM Filtering:
    - BAM files filtered according to some threshold on quality score, saved into `bam_dir`. All files that have
      undergone this step will contain `filtered` in their name's suffix.

7. BAM Indexing:
    - Index files of input BAM files (one for each BAM), saved into `bam_dir` with same name as input BAM plus
    `.bai`.

8. BAM Stats:
    - Statistics summary of input BAM files (one for each BAM), saved into "`bam_dir`/stats/" with extention
    `.stats.txt`.

9. Gene Counts:
    - Gene expression count files (one for each input BAM), saved into `gene_counts_dir` with extention
    `.counts.txt`.

10. Splicing analysis:
    - Files with differential splicing data, saved into `splicing_dir` with extention `.txt`.
    - `summary.txt` containing summary of all differential splicing events detected, saved into `splicing_dir`.
    - `.rmats` files with summary of BAM processing, saved into `splicing_dir`.

11. Results Summary:
    - html reports of all steps run, saved into `report_dir`.


## How to run your pipeline

1. Make sure you satisfy all requirements listed in [this section](#requirements).

2. Clone this repository, using:
   ````
   $ git clone git@github.com:TommasoTarchi/autoRNAseq.git
   ````

3. If your pipeline does not contain exclusively steps 1 and/or 11, produce a txt file listing input files, as described in
   [this section](#input-files).

4. Edit the `config.json` file as follows:

    - Set variables in `run_processes` section to true for the processes you wish to execute.

    - Configure `data_paths` to specify paths to your data following the descriptions in [this section](#data-paths).
      Make sure that:
      - `trimmed_fastq_dir` contains a subdirectory called "reports/";
      - `out_bam_dir` contains three subdirectories called "logs/", "stats/" and "tabs/";
      - if you don't need a path variable set it to an empty string or leave it as it is.
      
    - Customize settings for each process you are running under the corresponding `processes` section in `config.json`, see
      [here](#process-specific-parameters) for details. Refer to your cluster's specifications for SLURM settings. Parameters of
      processes you are not running will be ignored (you can leave them as they are).

    - Set `container_dir` to the path to directory with container images.

    - Set `nf_work_dir` to the path to work directory of your choice (choose a directory with sufficient disk available, if running
      on a cluster, we suggest to place this directory in the `scratch/`, if available). If you leave it empty, `./work/` (nextflow
      default) will be passed.
    
    - Change `run_locally` variable to true if you want to run the pipeline on your local machine (**not recommended** for
      most applications). If you want to run it on a cluster, leave it to false.

    - Change `save_all_bams` variable to true if you want to keep all BAM files produced along the pipeline (usually **not
      recommended**, especially when working with many files). If you only want the final result, leave it to false.

6. Run the pipeline using:
   ````
   $ nextflow run main.nf
   ````

7. (**optional**) If your pipeline was run successfully and you think you will not need any of the temporary files (i.e
   those not included among the outputs), we strongly suggest to clean `nf_work_dir`. The program is optimized to use
   less disk possible, however temporary files could still occupy a lot of disk space.








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
[rMATS-turbo]: https://github.com/Xinglab/rmats-turbo
[multiQC]: https://multiqc.info/docs/ 
[wget]: https://www.gnu.org/software/wget/manual/wget.html
[curl]: https://curl.se/docs/
