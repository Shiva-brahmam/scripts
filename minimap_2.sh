#!/bin/bash
# Define paths to folders
ASSEMBLY_DIR="/mnt/d/Documents/demo/depth/assembly"
READS_DIR="/mnt/d/Documents/demo/depth/reads"
OUTPUT_DIR="/mnt/d/Documents/demo/depth/minimap"

# Create output directories
mkdir -p $OUTPUT_DIR/index
mkdir -p $OUTPUT_DIR/sam
mkdir -p $OUTPUT_DIR/mapped
mkdir -p $OUTPUT_DIR/sorted_bam
mkdir -p $OUTPUT_DIR/indexed_bam
mkdir -p $OUTPUT_DIR/depth
mkdir -p $OUTPUT_DIR/bins

# Loop through each assembly file
for ASSEMBLY in $ASSEMBLY_DIR/*.fasta; do
    SAMPLE=$(basename "$ASSEMBLY" .fasta)  # Extract sample name (e.g., barcode55)

    echo "Processing $SAMPLE..."

    # 1. Index the assembly and save in separate folder
    echo "Indexing assembly for $SAMPLE..."
    minimap2 -d $OUTPUT_DIR/index/${SAMPLE}.mmi $ASSEMBLY

    # 2. Map reads and save SAM file
    echo "Mapping reads and saving SAM file for $SAMPLE..."
    minimap2 -ax map-ont $OUTPUT_DIR/index/${SAMPLE}.mmi $READS_DIR/${SAMPLE}.fastq.gz > $OUTPUT_DIR/sam/${SAMPLE}.sam

    # 3. Convert SAM to BAM
    echo "Converting SAM to BAM for $SAMPLE..."
    samtools view -bS $OUTPUT_DIR/sam/${SAMPLE}.sam > $OUTPUT_DIR/mapped/${SAMPLE}.bam

    # 4. Sort BAM file
    echo "Sorting BAM file for $SAMPLE..."
    samtools sort -o $OUTPUT_DIR/sorted_bam/${SAMPLE}_sorted.bam $OUTPUT_DIR/mapped/${SAMPLE}.bam

    # 5. Index BAM file
    echo "Indexing BAM file for $SAMPLE..."
    samtools index $OUTPUT_DIR/sorted_bam/${SAMPLE}_sorted.bam
    mv $OUTPUT_DIR/sorted_bam/${SAMPLE}_sorted.bam.bai $OUTPUT_DIR/indexed_bam/${SAMPLE}_sorted.bai

    # 6. Compute depth for MetaBAT2
    echo "Computing depth for $SAMPLE..."
    jgi_summarize_bam_contig_depths --outputDepth $OUTPUT_DIR/depth/${SAMPLE}_depth.txt $OUTPUT_DIR/sorted_bam/${SAMPLE}_sorted.bam

    # 7. Run MetaBAT2
    echo "Running MetaBAT2 for $SAMPLE..."
    mkdir -p $OUTPUT_DIR/bins/$SAMPLE
    metabat2 -i $ASSEMBLY -a $OUTPUT_DIR/depth/${SAMPLE}_depth.txt -o $OUTPUT_DIR/bins/$SAMPLE/bin --minContig 2500

    echo "Completed processing for $SAMPLE."
    echo "------------------------------"
done

echo "All samples processed successfully!"
