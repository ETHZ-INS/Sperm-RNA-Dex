#!/bin/bash
#
#path to binaries
CONDA_BINS=/nfs/users/nfs_a/ac40/software/miniconda3/bin
#
#Variables:
	base_name=KG-7-
adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTC
sample_no=$LSB_JOBINDEX
#Add a 0 before the sample name if it is <10
if [ $sample_no -lt "10" ]; then
			sample_no=`echo 0$sample_no`
	fi
#
#
#cd into directory that contains the files merged in previous step
cd output_merge_2lanes
#make directory to store output
if [ ! -d output_cutadapt ]; then
		mkdir output_cutadapt
	fi
#
#Run cutadapt
$CONDA_BINS/cutadapt \
-a $adapter \
--minimum-length 15 \
-o output_cutadapt/cutadapt_${base_name}${sample_no}.fastq \
${base_name}${sample_no}.fastq.gz
#
exit $?
