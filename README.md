# AutoRNAseq: automated pipelines for paired-end RNAseq data analysis

In this repository you can find the following automated pipelines:

1. [Gene expression count](./gene_count-pipeline): from paired-end reads to gene expression counts.
2. Differential expression analysis: *soon available...*

All pipelines are composed of a certain number of steps. The user can choose to run **any
combination** of these steps, given that the correct arguments and containers are provided.

All pipelines use *Nextflow* as automation tool and *Singularity* for containerization.
If you are new to these tools, you can find the documentation, respectively, [here][nextflow]
and [here][singularity]. In particular, it can be useful to consult [this section][nextflow_and_singularity]
for using Nextflow and Singularity together.

Each step of each pipeline runs on a different container, so that you can just download the
containers needed by the pipeline steps you intend to run. Container images ready to use can be found among
the assets of the [program release][assets].

If you want to build the containers by yourself and you are familiar with conda, you can follow
the simple steps described in [this tutorial][build_containers].

Pipelines support both *SLURM* (default) and local run. However, the second one is not recommended
for most demanding pipelines.

When run with SLURM, each process is launched as a different job on the cluster. Resources
parameters can be adjusted differently for each step of the pipeline.


## What you will find in this repository

Content of the repository (each pipeline has a dedicated directory):

- This README file: general description of the repo.
- [`gene_count-pipeline/`](./gene_count-pipeline/): gene expression count pipeline, containing:
  - source code for pipeline (`main.nf`, `nextflow.config`, `config.json`);
  - README.md, with instruction for running it.


## Gene expression count pipeline

This pipeline can be used to produce compressed alignment files (*BAM*) and gene expression count files
from **paired-end** read files (*FastQ*).

It is composed of the following steps:

1. **Genome Indexing**: preprocess the genome for alignment.
2. **FastQ trimming**: trim reads based on length, Phred score and adapters, and quality control.
3. **Alignment**: properly align reads to the reference genome.
4. **BAM Sorting**: sort alignment files.
5. **Remove duplicates**: remove (or mark only) duplicates in alignment files.
6. **BAM Filtering**: filter aligned reads by MAPQ quality score.
7. **BAM Indexing**: index the alignment files.
8. **BAM Stats**: generate a statistical summary of the alignment.
9. **Gene Counts**: quantify gene expression.
10. **Results Summary**: summarize the results.

A more detailed description of the pipeline and instruction on how to run it can be found in the
[related README.md](./gene_count-pipeline/README.md).







[nextflow]: https://www.nextflow.io/docs/latest/index.html
[singularity]: https://apptainer.org/user-docs/master/
[nextflow_and_singularity]: https://nextflow.io/docs/edge/container.html#singularity
[assets]: https://github.com/TommasoTarchi/autoRNAseq/releases
[build_containers]: https://github.com/fburic/notes/blob/master/singularity_conda.md
