#!/bin/bash
#
#Variables:
	base_name=KG-7-
sample_no=$LSB_JOBINDEX
#
#make directory to store output
	if [ ! -d output_merge_2lanes ]; then
			mkdir output_merge_2lanes
	fi
#
#Add a 0 before the sample name if it is <10
if [ $sample_no -lt "10" ]; then
			sample_no=`echo 0$sample_no`
	fi
#
#merge 2 files
	gunzip -c \
${base_name}${sample_no}_S*_L001_R1_001.fastq.gz \
${base_name}${sample_no}_S*_L002_R1_001.fastq.gz \
> output_merge_2lanes/${base_name}${sample_no}.fastq 
#Compress output file
	 gzip output_merge_2lanes/${base_name}${sample_no}.fastq 
#
exit $?
