#!/bin/bash
trimmomatic PE \
	-phred33 \
	-summary /home/angad/Vikash/rnaseq_analysis/trimmomatic/summary_files/SRR1_S206_L003_001_summary.txt \
	-threds 32 \ 
	./home/angad/Vikash/RNA-seq/GAR1_merged_R1_001.fastq.gz /home/angad/Vikash/RNA-seq/GAR1_merged_R2_001.fastq.gz \
	./home/angad/Vikash/rnaseq_analysis/trimmomatic/GAR1_merged_R1_001_Q.fastq.gz /home/angad/Vikash/rnaseq_analysis/trimmomatic/UQ/GAR1_merged_R1_001_uq.fastq.gz \
	./home/angad/Vikash/rnaseq_analysis/trimmomatic/GAR1_merged_R2_001_Q.fastq.gz /home/angad/Vikash/rnaseq_analysis/trimmomatic/UQ/GAR1_merged_R2_001_uq.fastq.gz \
	ILLUMINACLIP:/home/angad/miniconda3/share/trimmomatic-0.39-2/adapters/commanAdap.fa:2:30:10 \
	HEADCROP:15 \
	LEADING:30 \
	TRAILING:30 \
	MINLEN:100
