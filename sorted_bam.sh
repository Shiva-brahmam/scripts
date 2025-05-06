#!/bin/bash

#give the input and output diretories
in_dir="/home/rnaseq_analysis/SAM_files"
out_dir="/home/rnaseq_analysis/SAM_files/sort_bam"

#no.of threads 
threads=32

#read all the files in the drectory
for f_bam in "$in_dir"/*_R.fastp_Q.bam; do
	filename=$(basename "$f_bam" _R.fastp_Q.bam)

	s_bam="$out_dir/$filename.sorted.bam"

	samtools sort -@ $threads "$f_bam" -o "$s_bam"
done
