#!/bin/bash
#
#path to binaries
CONDA_BINS=/nfs/users/nfs_a/ac40/software/miniconda3/bin
#
#Variables:
adapter=TGGAATTCTCGGGTGCCAAGG
sample_no=$LSB_JOBINDEX
#
#make directory to store output
if [ ! -d output_cutadapt ]; then
		mkdir output_cutadapt
	fi
#
#Run cutadapt
$CONDA_BINS/cutadapt \
-a $adapter \
--minimum-length 15 \
--discard-untrimmed \
-o output_cutadapt/F${sample_no}.fastq.gz \
20200123.B-F${sample_no}_R1.fastq.gz
#
exit $?
