# AutoRNAseq: an automated pipeline for paired-end RNAseq data analysis

In this repository you can find the following automated pipelines:

1. Gene expression count: from paired-end reads to gene expression.

All pipelines use *Nextflow* as automation tool and *Singularity* for containerization.
If you are new to these tools, you can find the documentation, respectively, [here][nextflow]
and [here][singularity]. In particular, it can be useful to consult [this section][nextflow_and_singularity]
for using Nextflow and Singularity together.

If you want to build the containers by yourself and you are familiar with conda, you can follow
the simple steps described in [this tutorial][build_containers].

All pipelines are composed of a certain number of steps. The user can choose to run **any
combination** of these steps, given that the correct arguments and containers are provided.

Each step of each pipeline runs on a different container, so that you can just install the
functions needed by the pipeline steps you intend to run.

Pipelines support both *SLURM* (default) and local run. However, the second one is not recommended
for most demanding pipelines.

When run with SLURM, each process is launched as a different job on the cluster. Resources
parameters can be adjusted differently for each step of the pipeline.


## Outline

bla bla bla


## What you will find in this repository

bla bla bla


## Gene expression count

### Pipeline steps

1. Genome Indexing: preprocess the genome for alignment.
2. Alignment: properly align reads to the reference genome.
3. BAM Sorting: sort alignment files (redundant if the previous step was run).
4. Remove duplicates: remove duplicates in alignment files.
5. BAM Filtering: quality filtering of aligned reads.
6. BAM Indexing: index the alignment files.
7. BAM Stats: generate a statistical summary of the alignment.
8. Gene Counts: quantify gene expression.
9. Results Summary: summarize the results.

### Requirements

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

Remember that you need **only** the containers corresponding to the steps you want to run.

- You need to download the three files in the [pipeline directory](./gene_expression_count/).
Make sure to have **all three files in the same directory**. However, notice that the location of
the directory containing these files is not relevant, as the pipeline can work on data located in
any directory of the filesystem, assumed needed permissions for the user.

### Parameters

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

Variables you don't need should be ignored by the program. However, to avoid any possible bug we
suggest to set them to an empty string (or an empty list if the variable is of type list).

`bam_dir` should contain three subdirectories, called "logs/", "stats/" and "tabs/".

For more details on supported file names see [this section](#how-to-run-your-pipeline).

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

### How to run your pipeline

1. Make sure you satisfy all requirements listed in [this section](#requirements).

2. Edit the "config.json" file as follows:

    2a. Set variables in "run_processes" section to true for the processes you wish to execute
    (if you set `all` to true, then all steps will be run regardless of the values of the following
    variables).

    2b. Configure `data_paths` to specify paths to your data. Remember to use lists for `fastq_files`
    and `bam_files`. Each path should be complete, in particular:

   - for `fastq_files`, only the common prefix of reads pair should be passed, i.e. one
        full path **without** the `_R#_001.fastq.gz` suffix (**glob patterns are allowed** - see
        [here](#example-of-input-fastq-files) for more details);
    
   - for `bam_files`, include full paths (**glob patterns are allowed**);
   
   - if you don't need a path variable set it to an empty string/list.

    2c. Customize settings for each process under the `processes` section in `config.json`. Refer to
    your cluster's specifications for SLURM settings, especially for the `queue` variable.

    2d. Set `run_locally` variable to true if you want to run the pipeline on your local machine
    (not recommended for most applications).

4. Run the pipeline using (make sure you have `main.nf`, `nextflow.config` and `config.json` in the same
directory):

````
$ nextflow run main.nf
````

**Notice** that all the information contained in the name of input FastQ and BAM files after the first dot
will be lost. If some relevant information is placed after dots, please change these dots with other separators.
Examples:
- invalid file name: `COV362-TREATED-replica1.Tot_S11.Aligned.sortedByCoord.out.bam`;
- valid file name: `COV362-TREATED-replica1-Tot_S11.Aligned.sortedByCoord.out.bam`.

#### Example of input FastQ files

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
- `TREATED-replica2-Tot_S22_R1_001.fastq.gz` (not processed, since it does not have a pair read file)
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
