#!/bin/bash

# Define the parent directory containing sample folders
BASE_DIR="/path/to/parent_directory"  # Change this to the actual path

# Define a single output directory for GTDB-Tk results
OUTPUT_DIR="/path/to/gtdbtk_results"  # Change this to your desired output path
mkdir -p "$OUTPUT_DIR"

# Path to GTDB-Tk database (Change this to your database location)
GTDBTK_DB="/path/to/gtdbtk_db"

# Number of CPUs to use (adjust based on your system)
CPUS=36

# File extension for bins
BIN_EXTENSION="fa"

# Log file to record progress and errors
LOG_FILE="$OUTPUT_DIR/gtdbtk_analysis.log"

# Function to log messages
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOG_FILE"
}

# Loop through each barcode directory
for barcode in "$BASE_DIR"/barcode*; do
    if [ -d "$barcode" ]; then
        # Extract barcode name (e.g., barcode55)
        BARCODE_NAME=$(basename "$barcode")
        log_message "Running GTDB-Tk on bins in $BARCODE_NAME"

        # Define output directory for this barcode
        GTDBTK_OUTPUT="$OUTPUT_DIR/$BARCODE_NAME"
        mkdir -p "$GTDBTK_OUTPUT"

        # Run GTDB-Tk classification
        if gtdbtk classify_wf --genome_dir "$barcode" \
                              --out_dir "$GTDBTK_OUTPUT" \
                              --extension "$BIN_EXTENSION" \
                              --cpus "$CPUS" \
                              --prefix "$BARCODE_NAME" \
                              --scratch_dir "$GTDBTK_OUTPUT" \
                              --db_dir "$GTDBTK_DB" >> "$LOG_FILE" 2>&1; then
            log_message "GTDB-Tk classification completed for $BARCODE_NAME. Results saved in $GTDBTK_OUTPUT"
        else
            log_message "ERROR: GTDB-Tk classification failed for $BARCODE_NAME. Check $LOG_FILE for details."
        fi
    else
        log_message "WARNING: $barcode is not a directory. Skipping."
    fi
done

log_message "GTDB-Tk classification completed for all bins."

echo "Completed anlaysis for all the files"
