#!/bin/bash

SECONDS=0
#give the input and output directories
in_dir="/home/angad/Vikash/rnaseq_analysis/fastp_QC"
out_dir="/home/angad/Vikash/rnaseq_analysis/SAM_files"
sum_dir="/home/angad/Vikash/rnaseq_analysis/SAM_files/fastp_sum" #summary of alignment

#index of reference genome
index="/home/angad/Vikash/rnaseq_analysis/genomes_oryza/oryza_indica_index/oryza_indica.toplevel.index"

#specify the number of threads
threads=32 # it can be changed based on the capability of system
p=20
#loop for all forward files in input directory
for f_file in "$in_dir"/*1.fastq.gz; do
        echo "input read1:"$f_file
        #extracting the base of file 
        filename=$(basename "$f_file" 1.fastq.gz)

        #reverse file name
        r_file="$in_dir/$filename"2.fastq.gz
        echo "input read2:"$r_file
        #check reverse file exist or not
        if [ ! -f "$r_file" ]; then
                echo "reverse reads file not found"
                continue
        fi
        #summary file name
        sum_name="$sum_dir/$filename.txt"
        echo "out summary:"$sum_name
        #output file name
        out_sam="$out_dir/$filename.sam"
        echo "out SAM file:"$out_sam 
	#convrsion of sam to bam
	o_bam="$out_dir/$filename.bam"
	echo "output bam:"$o_bam
	#conversion of bam to sorted bam
	s_bam="$out_dir/$filename.sorted.bam"
	echo "output sorted bam:"$s_bam

        #run hisa2
        hisat2 -p $threads --phred33 --summary-file "$sum_name" -x "$index" -1 "$f_file" -2 "$r_file" -S "$out_sam"
	#convert sam to bam 
	samtools view -S -b "$out_sam" > "$o_bam"
	#remove unwanted sam file
	rm "$out_sam"

	#convert bam to sorted bam
	#samtools sort -@ $P "$O_bam" -o "$s_bam" 

done

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds "

# f_file -> forward read file
# r_file -> reverrse read file 
# after execution of this file first four lines will input file names and out file names

