import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor

# File paths
listA_file = "EPIC_context.txt"  # anonymized for public release
annotation_file = "EPIC_annotation.txt"
list_file = "list.txt"
output_dir = "intermediate_files"
os.makedirs(output_dir, exist_ok=True)

# Create dict_pos mapping ProbeID to chromosomal positions
dict_pos = {}
with open(annotation_file, "r", encoding="utf-8") as f_1:
    for line in f_1:
        cols = line.strip().split('\t')
        if len(cols) >= 4:
            dict_pos[cols[0]] = f"chr{cols[2]}:{cols[3]}"

# Split listA_file by chromosome (now only reading ProbeID)
listA = pd.read_csv(
    listA_file,
    sep="\t",
    header=None,
    usecols=[0],  # Only read the ProbeID column
    dtype=str
)
listA.columns = ["ProbeID"]  # Rename the single column to "ProbeID"

# Add chromosome information based on ProbeID
listA["Chr"] = listA["ProbeID"].map(lambda x: dict_pos[x].split(":")[0] if x in dict_pos else None)
for chr_name, group in listA.groupby("Chr"):
    group.to_csv(f"{output_dir}/listA_{chr_name}.txt", sep="\t", index=False)

# Load and split listB_file by chromosome
import sys
listB_file = sys.argv[1]
listB = pd.read_csv(
    listB_file,
    sep=",",
    header=0,
    names=["Chr", "Start", "End", "Early", "Late", "Log2_Ratio", "Smoothed_Log2_Ratio"],
    dtype=str
)
for chr_name, group in listB.groupby("Chr"):
    group.to_csv(f"{output_dir}/listB_{chr_name}.txt", sep="\t", index=False)

# Process listA and listB for each chromosome
def process_chromosome(chr_name):
    listA_path = f"{output_dir}/listA_{chr_name}.txt"
    listB_path = f"{output_dir}/listB_{chr_name}.txt"
    if not os.path.exists(listB_path):
        return  # skip if no data for this chromosome

    listA = pd.read_csv(listA_path, sep="\t", dtype=str)
    listB = pd.read_csv(listB_path, sep="\t", dtype=str)
    listB["Start"] = pd.to_numeric(listB["Start"], errors="coerce")
    listB["End"] = pd.to_numeric(listB["End"], errors="coerce")
    listB = listB.dropna(subset=["Start", "End"])
    listB["Start"] = listB["Start"].astype(int)
    listB["End"] = listB["End"].astype(int)

    def is_in_region(probe_id):
        if probe_id in dict_pos:
            chr_info, pos = dict_pos[probe_id].split(":")
            pos = int(pos)
            matches = listB[
                (listB["Chr"] == chr_info) &
                (listB["Start"] <= pos) &
                (listB["End"] >= pos)
            ]
            return matches["Smoothed_Log2_Ratio"].iloc[0] if not matches.empty else "No"
        else:
            return "No_Info"

    listA["Result"] = listA["ProbeID"].apply(is_in_region)
    listA.to_csv(f"{output_dir}/result_{chr_name}.txt", sep="\t", index=False)

if __name__ == "__main__":
    # Parallel processing per chromosome
    chromosomes = [f"chr{i}" for i in range(1, 25)]
    with ProcessPoolExecutor() as executor:
        executor.map(process_chromosome, chromosomes)

    # Concatenate results
    result_files = [f"{output_dir}/result_{chr_name}.txt" for chr_name in chromosomes if os.path.exists(f"{output_dir}/result_{chr_name}.txt")]
    if result_files:
        final_output_file = "final_output.txt"
        pd.concat([pd.read_csv(f, sep="\t") for f in result_files]).to_csv(final_output_file, sep="\t", index=False)
        print(f"Final results saved to {final_output_file}")
    else:
        print("No results to concatenate.")
