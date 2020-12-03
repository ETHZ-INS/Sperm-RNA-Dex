#!/bin/bash
# 
#path to binaries
CONDA_BINS=/nfs/users/nfs_a/ac40/software/miniconda3/bin
#
#cd into the directory containing the adapter-trimmed files
#
cd output_cutadapt/output_trimming
#
#Files to process
file_list=`ls *.fastq.gz`
#
#make directory to store fastqc output
	if [ ! -d output_fastqc ]; then
			mkdir output_fastqc
	fi
#
#Run fastqc
$CONDA_BINS/fastqc -o output_fastqc $file_list
#
exit $?
