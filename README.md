# rnaseq_star

A simple Nextflow (DSL2) pipeline to quantify paired-end RNA-seq data against an existing index with `STAR`,
followed by creating a matrix of counts with `featureCounts`.

![CI](https://github.com/ATpoint/rnaseq_star/actions/workflows/ci.yml/badge.svg)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

**USAGE**

On the HPC this will run STAR with 16 cores and 30GB of RAM:

```bash

nextflow run main.nf \
    -profile singularity,slurm \
    --fastq 'path/to/fastq/*_{1,2}.fq.gz' \
    --index 'path/to/index/folder/' \
    --gtf 'path/to/gtf'

```