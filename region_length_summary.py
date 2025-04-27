import os
import glob

def calculate_region_length(file_path):
    """
    #Reads a file and calculates the total length of regions in megabases (Mb).
    """
    regions = []

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                start = int(parts[1])
                end = int(parts[2])
                region_length = (end - start) / 1e6  # Convert to megabases
                regions.append(region_length)

    return sum(regions)


def process_files(directory_path):
    """
    Processes all .txt files in the specified directory and calculates the total region length (Mb) per sample.
    """
    file_paths = glob.glob(os.path.join(directory_path, "*.txt"))
    results = {}

    for file_path in file_paths:
        sample_name = os.path.basename(file_path)  # Use file name as sample ID
        total_length_mb = calculate_region_length(file_path)
        results[sample_name] = total_length_mb

    return results


def print_results(results):
    """
    Prints the results in a formatted table.
    """
    print(f"{'Sample ID':<50} {'Total Length (Mb)'}")
    print("-" * 70)

    for sample, total_length in results.items():
        print(f"{sample:<50} {total_length:.2f} Mb")


# === User input ===
# Specify the path to the directory containing the input .txt files
directory_path = '/path/to/your/directory'  # <- Replace this with your actual path

# Process files and display results
results = process_files(directory_path)
print_results(results)
