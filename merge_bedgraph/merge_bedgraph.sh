#!/bin/bash

# This script merges multiple Repli-seq-derived bedGraph files
# into a unified matrix by first creating standardized genomic intervals
# and then mapping values to those intervals by chromosome.

# Step 0: Check for the input file list
file_list="list.txt"

if [ ! -f "$file_list" ]; then
    echo "Error: Input file list '$file_list' not found." >&2
    exit 1
fi

# Read file list into array
files=()
while IFS= read -r line; do
    files+=("$line")
done < "$file_list"

# Step 1: Generate standardized merged genomic ranges
echo "[Step 1] Generating temp_merged_ranges.tsv..."
python scripts/generate_merged_ranges.py

if [ $? -ne 0 ]; then
    echo "Error: Failed to generate merged ranges." >&2
    exit 1
fi
echo "Merged ranges created successfully."

# Step 2: Process each chromosome in parallel
chromosomes=(chr{1..22} chrX chrY)

echo "[Step 2] Processing each chromosome in parallel..."
for chr in "${chromosomes[@]}"; do
    python scripts/map_values_to_ranges.py "$chr" "${files[@]}" &
done

wait
echo "Chromosome-wise processing complete."

# Step 3: Merge per-chromosome results into one file
echo "[Step 3] Merging final result..."
merged_file="final_merged_result.tsv"
head -n 1 "chr1_merged_result.tsv" > "$merged_file"

for chr in "${chromosomes[@]}"; do
    tail -n +2 "${chr}_merged_result.tsv" >> "$merged_file"
done

echo "Final result saved to '$merged_file'."
