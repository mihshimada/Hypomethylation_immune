import pyranges as pr
import itertools
import glob
import os

def load_regions(file):
    """
    Load genomic regions from a file into a PyRanges object.
    Assumes a tab-delimited format with columns: chromosome, start, end, and other metadata.
    """
    regions = []
    with open(file, 'r') as f:
        for line in f:
            if line.strip():  # Skip empty lines
                parts = line.split("\t")
                chrom, start, end = parts[0], int(parts[1]), int(parts[2])
                regions.append([chrom, start, end])
    return pr.from_dict({"Chromosome": [r[0] for r in regions],
                         "Start": [r[1] for r in regions],
                         "End": [r[2] for r in regions]})

def calculate_shared_regions(file1, file2):
    """
    Calculate the shared overlap (in Mb) between two genomic region files.
    """
    # Load regions from files
    ranges1 = load_regions(file1)
    ranges2 = load_regions(file2)

    # Find shared regions
    shared = ranges1.intersect(ranges2)

    # If there are no shared regions, the resulting PyRanges object is empty
    if shared.df.empty:  # Check if DataFrame inside PyRanges is empty
        total_shared_length = 0
    else:
        # Calculate total shared length manually
        shared_df = shared.df
        total_shared_length = (shared_df["End"] - shared_df["Start"]).sum()

    # Convert to Mb and return
    return total_shared_length / 1e6

def compare_all_files(files):
    """
    Perform pairwise comparison of all files and calculate shared regions.
    """
    results = []
    for file1, file2 in itertools.combinations(files, 2):
        shared_length_mb = calculate_shared_regions(file1, file2)
        results.append((os.path.basename(file1), os.path.basename(file2), shared_length_mb))
    return results

# Main execution
if __name__ == "__main__":
    # Replace with the path to your files
    file_pattern = "/Users/mihoko/Desktop/Methylation_otherDisease/NT1_2410/PMD_all/Overlapped/*.txt"  # Adjust the path to your input files
    files = glob.glob(file_pattern)

    if len(files) < 2:
        print("At least two files are required for comparison.")
    else:
        # Perform pairwise comparison
        results = compare_all_files(files)

        # Write results to a file
        output_file = "result.txt"
        with open(output_file, "w") as f:
            f.write("Results of pairwise comparisons:\n")
            for file1, file2, shared_length in results:
                line = f"{file1} and {file2}: {shared_length:.2f} Mb\n"
                f.write(line)
                print(line)  # Also print to console
        print(f"\nResults saved to {output_file}")
