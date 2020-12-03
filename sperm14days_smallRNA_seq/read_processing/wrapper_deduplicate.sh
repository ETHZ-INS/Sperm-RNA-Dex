#!/bin/bash
#
#Deduplication and trimming of 4 random nucleotides at both ends of the reads
#
#cd into the directory containing the adapter-trimmed files
#
cd output_cutadapt
#
#make directory to store output
if [ ! -d output_deduplicate ]; then
		mkdir output_deduplicate
fi
#
for f in *fastq.gz; do
	echo "Deduplicating $f..."
	zcat $f | paste - - - - | \
	sort -u --parallel 4 --buffer-size 2G --stable -t $'\t' -k2,2 | \
	awk -F"\t" '{
		sl=length($2) 
		print $1"\n"substr($2,5,sl-8)"\n+\n"substr($4,5,sl-8) 
		}' | gzip > output_deduplicate/`basename $f fastq.gz`dedup.fastq.gz
done
#
exit $?
