# Psolaren cross-linking to study Zika-host interactions

## Overview

This repository provides code for the paper 
**COMRADES determines in vivo RNA structures and interactions**
by [Ziv _et al._ (2018)](https://doi.org/10.1038/s41592-018-0121-0).

## Repository structure

The `scripts` subdirectory contains a number of common scripts to process the sequencing data.
These will be automatically called from elsewhere and do not need to be modified except for debugging purposes.

Each other subdirectory refers to a batch of data collected throughout the course of this study.
The data presented in the paper comes from `real_20170807`.
Other batches involve a variety of different protocol conditions and biological samples, and can be largely ignored.

## Repeating the read alignment 

Follow the instructions in `genomes/README.md` to create the genome annotation.
This assumes that:

- You are on a LSF system.
- You have [_STAR_](https://github.com/alexdobin/STAR) installed.

Download the FASTQ files from ArrayExpress (E-MTAB-6427) into `real_20170807/fastq`.
Then, enter `real_20170807` and execute the `run_me.sh` script.
This further assumes:

- You have Python installed with an up-to-date version of _pysam_.
- You have [_umi-tools_](https://github.com/CGATOxford/UMI-tools) installed.
- You have R installed with an up-to-date version of [_diffHic_](https://bioconductor.org/packages/diffHic).

This script will perform mapping, duplicate removal and construction of HDF5 files for differential analyses in _diffHic_.

## Repeating the differential analyses

The `real_20170807/analysis` subdirectory contains:

- `qc.Rmd` and `supp_qc.Rmd`, for quality control.
The `supp_qc.Rmd` should be run first to be included in the report upon compilation of `qc.Rmd`.
- `differential.Rmd`, to test for differential intensity between COMRADES (nicknamed "livefire") and the reverse control.
This is done for interactions between human genes and the ZIKA genome.
- `self.Rmd`, to test for differential intensity between COMRADES and the reverse control for ZIKA-ZIKA interactions.

