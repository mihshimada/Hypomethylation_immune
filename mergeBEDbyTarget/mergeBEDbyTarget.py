import os

# Load mapping between file ID and target
file_to_target = {}
with open("files_assay_target.txt", "r", encoding="utf-8") as f:
    for line in f:
        parts = line.strip().split('\t')
        file_to_target[parts[0]] = parts[1]

# Get the list of unique targets
targets = set(file_to_target.values())
print(targets)

# Process each target
for target_item in targets:
    file_ids = [fid for fid, t in file_to_target.items() if t == target_item]
    print(target_item)
    print(len(file_ids))

    output_filename = f"{target_item}_merged.txt"
    with open(output_filename, "w", encoding="utf-8") as out_file:

        # Initialize peak lists per chromosome (chr1 to chr22)
        peak_list_all = [[] for _ in range(22)]

        # Read each file and collect peak regions
        for file_id in file_ids:
            filename = file_id + ".bed"
            with open(filename, "r", encoding="utf-8") as bed_file:
                for line in bed_file:
                    parts = line.strip().split()
                    chr_name = parts[0].replace('chr', '')
                    if chr_name.isdigit():
                        chr_idx = int(chr_name) - 1
                        peak_list_all[chr_idx].append([int(parts[1]), int(parts[2])])

        # Merge overlapping peaks (without counting overlaps)
        merged_peak_list_all = []
        for chrom_idx in range(22):
            chr_peaks = peak_list_all[chrom_idx]
            chr_peaks.sort()

            merged_chr_peaks = []
            if not chr_peaks:
                merged_peak_list_all.append([])
                continue

            start, end = chr_peaks[0]
            for i in chr_peaks[1:]:
                if i[0] > end:
                    merged_chr_peaks.append([start, end])
                    start, end = i
                else:
                    end = max(end, i[1])
            merged_chr_peaks.append([start, end])
            merged_peak_list_all.append(merged_chr_peaks)

        # Write merged peaks to output file
        for chrom_idx in range(22):
            chr_number = chrom_idx + 1
            chr_merged = merged_peak_list_all[chrom_idx]
            for peak in chr_merged:
                out_file.write(f"{chr_number}\t{peak[0]}\t{peak[1]}\n")
