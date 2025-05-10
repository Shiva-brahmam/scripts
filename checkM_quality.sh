#!/bin/bash

# Define the parent directory containing barcode folders
BASE_DIR="/path/to/parent_directory"  # Change this to the actual path

# Define a single output directory for CheckM results
OUTPUT_DIR="/path/to/checkm_results"  # Change this to output path
mkdir -p "$OUTPUT_DIR"

# Loop through each barcode directory
for barcode in "$BASE_DIR"/barcode*; do
    if [ -d "$barcode" ]; then
        # Extract barcode name (e.g., barcode55)
        BARCODE_NAME=$(basename "$barcode")
        echo "Processing bins in $BARCODE_NAME"

        # Define output filenames
        LINEAGE_FILE="$OUTPUT_DIR/${BARCODE_NAME}.ms"
        SUMMARY_FILE="$OUTPUT_DIR/${BARCODE_NAME}_summary.txt"

        # Run CheckM lineage workflow and save lineage.ms as barcodeXX.ms
        checkm lineage_wf -x fa -t 36 "$barcode" "$OUTPUT_DIR" --file "$LINEAGE_FILE"

        # Generate summary as barcodeXX_summary.txt
        checkm qa "$LINEAGE_FILE" "$OUTPUT_DIR" -o 2 > "$SUMMARY_FILE"

        echo "CheckM analysis completed for $BARCODE_NAME."
        echo "Results saved: $LINEAGE_FILE and $SUMMARY_FILE"
    fi
done

echo "CheckM quality check completed for all bins."

