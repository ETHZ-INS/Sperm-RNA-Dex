#!/bin/bash
# 
#path to binaries
CONDA_BINS=/nfs/users/nfs_a/ac40/software/miniconda3/bin
#
#path to gtrna database
genome=/lustre/scratch115/teams/miska/users/ac40/databases/STARgenomes/ENSEMBL.mus_musculus.release-75/Mus_musculus.GRCm38.75.gtf
#
#cd into the directory containing the aligned .bam files
cd output_cutadapt/output_deduplicate/output_STAR
#
#Files to process
file_list=`ls *.out.bam`
#
#make directory to store output
	if [ ! -d output_biotypes_plot_featureCounts ]; then
			mkdir output_biotypes_plot_featureCounts 
	fi
#
#Run featureCounts
$CONDA_BINS/featureCounts \
-a $genome \
-g gene_id \
-t exon \
-M \
-o output_biotypes_plot_featureCounts/exon_featurecounts.txt \
--fraction \
--minOverlap 15 \
-T 6 \
$file_list
#
# $? is a shell variable that is set to the exit code of the last command
exit $?
