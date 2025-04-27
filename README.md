# Hypomethylation_immune
This repository contains the Python scripts used to perform the analyses conducted in our study (doi: ).

The scripts compiled here were utilized during the research process and are intended for formatting and analyzing data specific to this project.
Although the programs are not designed as a general-purpose analysis pipeline, we are happy to respond to any inquiries and are open to modifications to support your usage.

## Python Requirements
Python >= 3.7

## How to use this repository
Please choose to either clone the entire directory or, if your focus is on a specific script, to download that particular script onto your local drive.

## Documentation
Please see [wiki](https://github.com/mihshimada/Hypomethylation_immune/wiki) for detailed explanation of scripts.

The list of programs is as follows:

**Count Flanking CpG Sites and Sequence Context**
* CpGContextCounter.py

**Calculate the total length of detected DMPs (in Mb)**
* region_length_summary.py
  
**Calculate the overlapping regions between PMDs**
* pairwise_overlap_length.py

**Align and merge multiple bedGraph files into unified genomic intervals**
* merge_bedgraph.sh
* generate_merged_ranges.py
* map_values_to_ranges.py

**Probe-wise Region Annotation Based on Chromosomal Coordinates**
* annotate_probes.sh
* CpGRegionAnnotator.py


**Merge datasets downloaded from the ENCODE project based on specified criteria**
* mergeBEDbyTarget.py

**Feature Region Overlap Checker**
* GenomicFeatureLocator.py


