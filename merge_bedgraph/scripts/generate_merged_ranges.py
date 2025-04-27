import pandas as pd
import numpy as np
from tqdm import tqdm
import re

# Enable progress bars for pandas operations
tqdm.pandas()

# Read the list of input files
print("Reading input file list...")
with open("list.txt", "r", encoding="utf-8") as f:
    files = [line.strip().split('\t')[0] for line in f]

# Define the columns for the bedGraph files
columns = ["Chr", "Start", "End"]

# Initialize a list to collect all ranges from the bedGraph files
all_ranges = []

# Read and process each input bedGraph file
dataframes = {}
print("Loading files...")
for file in tqdm(files, desc="Processing files"):
    name = file.split(".")[0]
    df = pd.read_csv(file, sep="\t", header=None, names=columns + [name])
    df = df.sort_values(by=["Chr", "Start", "End"])
    dataframes[name] = df
    all_ranges.append(df[["Chr", "Start", "End"]])

# Function to convert chromosome names to numeric values for sorting
def extract_chr_number(chr_value):
    chr_value = chr_value.replace("chr", "").replace("Chr", "")
    if chr_value == 'X':
        return 23
    elif chr_value == 'Y':
        return 24
    elif chr_value == 'M':
        return 25
    else:
        match = re.search(r'\d+', chr_value)
        return int(match.group()) if match else float('inf')

# Combine all ranges and remove duplicates
print("Consolidating unique ranges...")
all_ranges = pd.concat(all_ranges).drop_duplicates().reset_index(drop=True)
all_ranges["Chr_num"] = all_ranges["Chr"].apply(extract_chr_number)

# Generate unique keys for each range based on chromosome number and start position
SCALE_FACTOR = 10**12
all_ranges["UniqueKeyStart"] = all_ranges["Chr_num"] * SCALE_FACTOR + all_ranges["Start"]

# Function to merge overlapping ranges for each chromosome
def create_merged_ranges(chr_group):
    breakpoints = np.unique(chr_group[['Start', 'End']].values.flatten())
    merged_ranges = pd.DataFrame({
        "Chr": chr_group.iloc[0]["Chr"],
        "Start": breakpoints[:-1],
        "End": breakpoints[1:]
    })
    merged_ranges["Chr_num"] = chr_group.iloc[0]["Chr_num"]
    merged_ranges["UniqueKeyStart"] = chr_group.iloc[0]["Chr_num"] * SCALE_FACTOR + merged_ranges["Start"]
    return merged_ranges

# Create merged ranges for each chromosome
print("Creating merged ranges for each chromosome...")
merged_ranges_all = (
    all_ranges.groupby("Chr", group_keys=False)
    .progress_apply(create_merged_ranges)
    .reset_index(drop=True)
)

# Save the temporary merged ranges to a file
merged_ranges_all.to_csv("temp_merged_ranges.tsv", sep="\t", index=False)
print("Temporary merged ranges saved to 'temp_merged_ranges.tsv'.")
