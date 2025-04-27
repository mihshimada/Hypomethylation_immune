import pandas as pd
import numpy as np
import sys
from tqdm import tqdm

tqdm.pandas()

# Command-line arguments
if len(sys.argv) < 3:
    print("Error: Please specify the chromosome name and at least one input file.")
    sys.exit(1)

chromosome = sys.argv[1]
files = sys.argv[2:]

# Load input files and unified ranges
columns = ["Chr", "Start", "End"]
try:
    merged_ranges_all = pd.read_csv("temp_merged_ranges.tsv", sep="\t")
except FileNotFoundError:
    print("Error: 'temp_merged_ranges.tsv' not found. Please generate the unified ranges first.")
    sys.exit(1)

merged_ranges = merged_ranges_all[merged_ranges_all["Chr"] == chromosome]

# Process data for each chromosome
print(f"Processing chromosome {chromosome}...")
intermediate_results = []

for file in tqdm(files, desc="Processing files"):
    # Read the input file
    df = pd.read_csv(file, sep="\t", header=None, names=columns + ["Value"])
    df = df[df["Chr"] == chromosome]

    # Map values based on the unified ranges
    result = []
    df_idx = 0
    num_rows = len(df)

    for _, row in merged_ranges.iterrows():
        chr_, start, end = row["Chr"], row["Start"], row["End"]
        value = 0

        # Check input data sequentially
        while df_idx < num_rows:
            df_row = df.iloc[df_idx]

            if df_row["Chr"] != chr_ or df_row["End"] <= start:
                # If not matching, move to the next row
                df_idx += 1
                continue

            if df_row["Start"] >= end:
                # If the range exceeds, stop
                break

            # If the range matches, retain the value
            if df_row["Start"] <= start and df_row["End"] >= end:
                value = df_row["Value"]
                break

        result.append([chr_, start, end, value])

    # Convert the result into a DataFrame
    result_df = pd.DataFrame(result, columns=["Chr", "Start", "End", file.split(".")[0]])
    intermediate_results.append(result_df)

# Merge intermediate results
final_result = merged_ranges.copy()

for result_df in intermediate_results:
    final_result = pd.merge(final_result, result_df, on=["Chr", "Start", "End"], how="left")

# Fill NaN values with 0
final_result = final_result.fillna(0)

# Save the result
output_file = f"{chromosome}_merged_result.tsv"
final_result.to_csv(output_file, sep="\t", index=False)
print(f"Results saved to '{output_file}'.")
