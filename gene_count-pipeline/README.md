# Gene expression count pipeline

This pipeline can be used to produce compressed alignment files (*BAM*) from "raw" paired-end
read files (*FastQ*).


## Outline

In this file you will find:

- [Pipeline steps](#pipeline-steps)
- [Requirements](#requirements)
- [Parameters description](#parameters-description)
  - [Data paths](#data-paths)
  - [Process specific parameters](#process-specific-parameters)
  - [Output files](#output-files)
- [How to run your pipeline](#how-to-run-your-pipeline)
  - [Example of input FastQ files](#example-of-input-fastq-files)


## Pipeline steps

1. Genome Indexing: preprocess the genome for alignment.
2. Alignment: properly align reads to the reference genome.
3. BAM Sorting: sort alignment files (redundant if the previous step was run).
4. Remove duplicates: remove duplicates in alignment files.
5. BAM Filtering: quality filtering of aligned reads.
6. BAM Indexing: index the alignment files.
7. BAM Stats: generate a statistical summary of the alignment.
8. Gene Counts: quantify gene expression.
9. Results Summary: summarize the results.


## Requirements

- You need to have Nextflow and Singularity installed on your machine. You can look at the
related documentation [here][nextflow] and [here][singularity] for instructions.

- You need to build, or download, containers with the following functions installed (each point
of the list describes the function needed by the corresponding step in the pipeline):

1. *STAR* ([docs][STAR])
2. *STAR* ([docs][STAR])
3. *SAMtools* ([docs][SAMtools])
4. *picard* ([docs][picard])
5. *SAMtools* ([docs][SAMtools])
6. *SAMtools* ([docs][SAMtools])
7. *SAMtools* ([docs][SAMtools])
8. *featureCounts* ([docs][featureCounts]) **or** *HTSeq* ([docs][HTSeq])
9. *multiQC* ([docs][multiQC])

(Remember that you need **only** the containers corresponding to the steps you want to run).


## Parameters description

All parameters can be set from the `config.json` file. We strongly suggest to **not modify** either
`main.nf` nor `nextflow.config`.

`config.json` is organized according to the following hierarchical structure:

````
{
  "run_processes": {
    ... boolean variables indicating whether each process should be run or not
    ("all" is to run all pipeline from first to last step (BAM sorting excluded)) ...
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
  
  "run_locally": boolean (default "false")
}
````


### Data paths

Particular care must be taken in setting the `data_path` variables.

If you run only some of the pipeline steps (as it is usually the case), you will only need some of
these variables. The following list shows for each data path variable which steps of the pipeline need
it to be set. If **at least one** of the steps you intend to run is listed under a variable, then you
need to set that variable.

- `index_dir`: path to directory for genome index files, required by: 1. genome indexing, 2. alignment.

- `fasta_file`: complete path to fasta file with reference genome, required by: 1. genome indexing.

- `annotation_file`: complete path to GTF/GFF files, required by: 1. genome indexing, 8. gene counts.

- `fastq_files`: list of complete paths to input (**zipped**) read files, required by: 2. alignment.

- `bam_dir`: path to directory to store output alignment files, required by: 2. alignment, 3. BAM sorting,
  4. remove duplicates, 5. BAM filtering, 7. BAM stats, 8. gene counts.

- `bam_files`: list of complete paths to input alignment files, required by: 3. BAM sorting, 4. remove duplicates,
  5. BAM filtering, 6. BAM indexing, 7. BAM stats, 8. gene counts (**Notice**: this variable is **never** needed when
  your pipeline contains the alignment step, i.e. number 2.).

- `gene_counts_dir`: path to directory to store gene counts files, required by: 8. gene counts.

- `report_dir`: path to directory to store produced reports and plots, required by: 9. summarize results.

`fastq_files` should be passed in couples of paired-end zipped read files, and should end in `_R#_001.fastq.gz` or
`_R#_001.fq.gz`, with `#` equal to 1 and 2 (de facto standard).

`bam_dir` should contain three subdirectories, called "logs/", "stats/" and "tabs/".


### Process specific parameters

As mentioned previously, inside the `processes` scope, each process has its own scope for parameter setting.

Common to all processes are the following variables:

````
"queue" -> string: name of cluster partition to run job on
"time" -> string: maximum time for job (example: "2h")
"memory" -> string: RAM used for the job (example: "2GB")
"container_path" -> string: full path to singularity image for the process
"num_threads" -> integer: number of threads (not supported for "remove_duplicates" and "summarize_results")
````

Other process-specific parameteres are:

````
"genome_indexing": {
  "max_RAM" -> string: maximum RAM for STAR indexing (should be equal to "memory")
}

"BAM_filtering": {
  "quality_thres" -> integer: threshold for quality filtering
}

"gene_counts": {
  "algo" -> string: algorithm for gene expression quantification (allowed options: "featureCounts","HTSeq")
}
````

Notice that in `config.json` some of the parameters are set to a default value. However, the value set is **not**
guaranteed to work.

**All** process-specific parameters of te processes you intend to run must be set to some value.


### Output files

The following is a list of output files of each step of the pipeline. If not specified, the output
file is always saved, independently of which steps are run after.

1. Genome Indexing:
    - Various files with indexed and preprocessed genome, saved into `inded_dir`.

2. Alignment:
    - BAM files with alignment (one per paired-end pair of fastq file), sorted by coordinates and saved
    into `bam_dir`. The name of the file will have all the relevant information contained in the fastq
    names followed by the suffix: `.Aligned.sortedByCoord.out.bam`.
    - Alignment log files, saved into "`bam_dir`/logs/".
    - Alignment tab files, saved into "`bam_dir`/tabs/".

3. BAM Sorting:
    - BAM files sorted by coordinates, saved into `bam_dir` with extention `.Aligned.sortedByCoord.out.bam`.

4. Remove duplicates:
    - BAM files with duplicates removed, saved into `bam_dir` with same name as the input BAM. Notice that
    the input BAM will be overwritten by this step if contained in `bam_dir`.
    - Duplicate metrics report, saved in "`bam_dir`/stats/" with extention `.dup_metrics.txt`.

5. BAM Filtering:
    - BAM files filtered according to some threshold on quality score, saved into `bam_dir` with extention
    `.Aligned.sortedByCoord.filtered.bam`.

6. BAM Indexing:
    - Index files of input BAM files (one for each BAM), saved into `bam_dir` with same name as input BAM plus
    `.bai`.

7. BAM Stats:
    - Statistics summary of input BAM files (one for each BAM), saved into "`bam_dir`/stats/" with extention
    `.stats.txt`.

8. Gene Counts:
    - Gene expression count files (one for each input BAM), saved into `gene_counts_dir` with extention
    `.counts.txt`.

9. Results Summary:
    - html reports of all steps run, saved into `report_dir`.


### Example of input FastQ files

Suppose you want to run reads alignment and suppose you have a directory containing the following files:

- `TREATED-replica1-Tot_S11_R1_001.fastq.gz`
- `TREATED-replica1-Tot_S11_R2_001.fastq.gz`
- `TREATED-replica2-Tot_S22_R1_001.fastq.gz`
- `TREATED-replica3-Tot_S34_R1_001.fq.gz`
- `TREATED-replica3-Tot_S34_R2_001.fq.gz`
- `TREATED-replica4-Tot_S25_R1_001.fastq`
- `TREATED-replica4-Tot_S25_R2_001.fastq`
- `RES_PT-replica1-Tot_S12_R1_001.fq.gz`
- `RES_PT-replica1-Tot_S12_R2_001.fq.gz`
- `RES_PT-replica2-Tot_S24_R1_001.fq.gz`
- `RES_PT-replica2-Tot_S24_R2_001.fq.gz`

Now, for instance if you set `fastq_files` to the list: [`TREATED-replica*`, `RES_PT-replica1-Tot_S12`,
`RES_PT-replica2-Tot_S24_R?_001.fq.gz`], the files in the directory will be treated in the following way:

- `TREATED-replica1-Tot_S11_R1_001.fastq.gz` and `TREATED-replica1-Tot_S11_R2_001.fastq.gz` (processed)
- `TREATED-replica2-Tot_S22_R1_001.fastq.gz` (not processed, since it does not have a corresponding paired read file)
- `TREATED-replica3-Tot_S34_R1_001.fq.gz` and `TREATED-replica3-Tot_S34_R2_001.fq.gz` (processed)
- `TREATED-replica4-Tot_S25_R1_001.fastq` and `TREATED-replica4-Tot_S25_R2_001.fastq` (not processed, since
they do not match the expected format)
- `RES_PT-replica1-Tot_S12_R1_001.fq.gz` and `RES_PT-replica1-Tot_S12_R2_001.fq.gz` (processed)
- `RES_PT-replica2-Tot_S24_R1_001.fq.gz` and `RES_PT-replica2-Tot_S24_R2_001.fq.gz` (not processed, since the
provided pattern wrongly includes the suffix `_R#_001.fq.gz`)

The output of alignment will therefore be:

- `TREATED-replica1-Tot_S11.Aligned.sortedByCoord.out.bam`
- `TREATED-replica3-Tot_S34.Aligned.sortedByCoord.out.bam`
- `RES_PT-replica1-Tot_S12.Aligned.sortedByCoord.out.bam`


## How to run your pipeline

1. Make sure you satisfy all requirements listed in [this section](#requirements).

2. Clone this repository, using:
   ````
   $ git clone git@github.com:TommasoTarchi/autoRNAseq.git
   ````

3. Navigate to the `gene_count-pipeline` directory and edit the `config.json` file as follows:

    3a. Set variables in "run_processes" section to true for the processes you wish to execute
    (if you set `all` to true, then all steps will be run regardless of the values of the following
    variables).

    3b. Configure `data_paths` to specify paths to your data. Remember to use lists for `fastq_files`
    and `bam_files`. Each path should be complete, in particular:

   - for `fastq_files`, only the common prefix of reads pair should be passed, i.e. one
        full path **without** the `_R#_001.fastq.gz` suffix (**glob patterns are allowed** - see
        [here](#example-of-input-fastq-files) for more details);
    
   - for `bam_files`, include full paths (**glob patterns are allowed**);
   
   - if you don't need a path variable set it to an empty string/list.

    3c. Customize settings for each process under the `processes` section in `config.json`. Refer to
    your cluster's specifications for SLURM settings. See [here](#process-specific-parameters) for details.

    3d. Set `run_locally` variable to true if you want to run the pipeline on your local machine
    (not recommended for most applications).

4. Run the pipeline using:
   ````
   $ nextflow run main.nf
   ````

**Notice** that all the information contained in the name of input FastQ and BAM files after the first dot
will be lost. If some relevant information is placed after dots, please change these dots with other separators.
Examples:
- invalid file name: `COV362-TREATED-replica1.Tot_S11.Aligned.sortedByCoord.out.bam`;
- valid file name: `COV362-TREATED-replica1-Tot_S11.Aligned.sortedByCoord.out.bam`.








[nextflow]: https://www.nextflow.io/docs/latest/index.html
[singularity]: https://apptainer.org/user-docs/master/
[nextflow_and_singularity]: https://nextflow.io/docs/edge/container.html#singularity
[build_containers]: https://github.com/fburic/notes/blob/master/singularity_conda.md
[STAR]: https://docs.csc.fi/apps/star/
[SAMtools]:https://www.htslib.org/doc/samtools.html
[picard]: https://broadinstitute.github.io/picard/
[featureCounts]: https://subread.sourceforge.net/featureCounts.html
[HTSeq]: https://htseq.readthedocs.io/en/master/overview.html
[multiQC]: https://multiqc.info/docs/
