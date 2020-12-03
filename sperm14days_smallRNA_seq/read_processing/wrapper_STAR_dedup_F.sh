#!/bin/bash
#
#Variables:
sample_no=$LSB_JOBINDEX
#
# 
#path to binaries
CONDA_BINS=/nfs/users/nfs_a/ac40/software/miniconda3/bin
#
#path to genome database
genome=/lustre/scratch115/teams/miska/users/ac40/databases/STARgenomes/ENSEMBL.mus_musculus.release-75
#
#cd into the directory containing the adapter-trimmed files
cd output_cutadapt/output_deduplicate
#
#make directory to store star output
	if [ ! -d output_STAR ]; then
			mkdir output_STAR 
	fi
#
#Run STAR
$CONDA_BINS/STAR \
--outFileNamePrefix output_STAR/F${sample_no}_dedup.STAR. \
--outFilterMultimapNmax 5000 \
--winAnchorMultimapNmax 5000 \
--outFilterMismatchNmax 0 \
--alignIntronMax 1 \
--alignEndsType EndToEnd \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 12 \
--genomeDir $genome \
--readFilesCommand gunzip -c \
--readFilesIn F${sample_no}_dedup.fastq.gz
#
exit $?
