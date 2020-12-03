#!/bin/bash
# 
#path to binaries
CONDA_BINS=/nfs/users/nfs_a/ac40/software/miniconda3/bin
#
#path to gtrna database
genome=/lustre/scratch115/teams/miska/users/ac40/databases/gtrnadb/mm10-tRNAs_from_Wayo.gff
#
#cd into the directory containing the aligned .bam files
cd output_cutadapt/output_deduplicate/output_cutadaptCCA/output_tRNA_STAR
#
#Files to process
file_list=`ls *.out.bam`
#
#make directory to store output
	if [ ! -d output_tRNA_featureCounts ]; then
			mkdir output_tRNA_featureCounts 
	fi
#
#Run featureCounts
$CONDA_BINS/featureCounts \
-a $genome \
-t exon \
-M \
-o output_tRNA_featureCounts/tRNA_featurecounts.txt \
--fraction \
--minOverlap 16 \
-T 6 \
$file_list
#
# $? is a shell variable that is set to the exit code of the last command
exit $?
