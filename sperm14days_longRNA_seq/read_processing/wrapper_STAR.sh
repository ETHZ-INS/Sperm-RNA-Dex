#!/bin/bash
#
#Variables:
base_name=cutadapt_KG-7-
sample_no=$LSB_JOBINDEX
#Add a 0 before the sample name if it is <10
if [ $sample_no -lt "10" ]; then
			sample_no=`echo 0$sample_no`
		fi
#
# 
#path to binaries
CONDA_BINS=/nfs/users/nfs_a/ac40/software/miniconda3/bin
#
#path to genome database
genome=/lustre/scratch115/teams/miska/users/ac40/databases/STARgenomes/ENSEMBL.mus_musculus.release-75
#
#cd into the directory containing the adapter-trimmed files
cd output_merge_2lanes/output_cutadapt
#
#make directory to store star output
	if [ ! -d output_STAR ]; then
			mkdir output_STAR 
	fi
#
#Run STAR
$CONDA_BINS/STAR \
--outFileNamePrefix output_STAR/${base_name}${sample_no}.STAR. \
--alignEndsType EndToEnd \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 12 \
--genomeDir $genome \
--readFilesIn ${base_name}${sample_no}.fastq
#
exit $?
