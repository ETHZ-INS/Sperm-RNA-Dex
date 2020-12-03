#!/bin/bash
#
#path to binaries
CONDA_BINS=/nfs/users/nfs_a/ac40/software/miniconda3/bin
#
#Variables:
adapter=CCA$
sample_no=$LSB_JOBINDEX
#
#cd into directory that contains the files from previous cutadapt step
cd output_cutadapt/output_deduplicate
#make directory to store output
if [ ! -d output_cutadaptCCA ]; then
		mkdir output_cutadaptCCA
	fi
#
#Run cutadapt
$CONDA_BINS/cutadapt \
-a $adapter \
--minimum-length 16 \
-o output_cutadaptCCA/F${sample_no}_dedup_noCCA.fastq.gz \
F${sample_no}_dedup.fastq.gz
#
exit $?
