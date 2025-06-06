#!/bin/bash
# census run batch script

# NCBI retrieval run below


find ../non-redundant_mapping_tables/census_subtype_after_filtering -name "*_ncbi_retrieved_af.csv" -type "f" | xargs -n 1 -I {} -P 1 python3.7 krona_table_generation_ncbi_whole.py {} {}_krona_categories.csv
