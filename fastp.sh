
#bin/bash
#input and output directories 
input_dir="/home/RNA-seq"
output_dir="/home/rnaseq_analysis/fastp_QC"
sum_out="/home/rnaseq_analysis/fastp_QC"
#threads=16
#loop to read the forward files
for f_file in "$input_dir"/*1_001.fastq.gz; do
	filename=$(basename "$f_file" 1_001.fastq.gz)
	echo $f_file
	r_file="$input_dir/$filename"2_001.fastq.gz
	echo $r_file
	#output file names
	out_f="$output_dir/$filename.fastp_Q1.fastq.gz"
	out_r="$output_dir/$filename.fastp_Q2.fastq.gz"
	out_sum="$sum_out/$filename.txt"
	out_html="$filename".html
	
	#execution of the code
	fastp -i "$f_file" -I "$r_file" -o "$out_f" -O "$out_r" -h "$out_html" -R "\"$out_sum\""
done
