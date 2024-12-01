#!/bin/bash
# census run batch script

# NCBI retrieval run iteratively using entrez
find ../non-redundant_mapping_tables/census_subtype_after_filtering -name "*unmatched.csv" -type "f" | xargs -n 1 -I {} -P 1 python3.7 ncbi_annotation_retrieval.py {} ../contig_environment_metadata/all_ncbi_annotations.csv {}_ncbi_retrieved_af.csv >> {}_ncbi_retrieval_readout.csv
