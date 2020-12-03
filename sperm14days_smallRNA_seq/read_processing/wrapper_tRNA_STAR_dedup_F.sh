#!/bin/bash
#
#Variables:
sample_no=$LSB_JOBINDEX
# 
#path to binaries
CONDA_BINS=/nfs/users/nfs_a/ac40/software/miniconda3/bin
#
#path to genome database
genome=/lustre/scratch115/teams/miska/users/ac40/databases/STARgenomes/ENSEMBL.mus_musculus.release-75
#
#cd into the directory containing the adapter-trimmed files
cd output_cutadapt/output_deduplicate/output_cutadaptCCA
#
#make directory to store star output
	if [ ! -d output_tRNA_STAR ]; then
			mkdir output_tRNA_STAR 
	fi
#
#Run STAR
$CONDA_BINS/STAR \
--outFileNamePrefix output_tRNA_STAR/F${sample_no}_dedup_noCCA.STAR. \
--outFilterMultimapNmax 5000 \
--winAnchorMultimapNmax 5000 \
--outFilterMismatchNmax 0 \
--alignIntronMax 2 \
--alignEndsType EndToEnd \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 12 \
--genomeDir $genome \
--readFilesCommand gunzip -c \
--readFilesIn F${sample_no}_dedup_noCCA.fastq.gz
#
exit $?
