#!/bin/bash
# census run batch script

# INPUT: NCBI tables from each subtype containing specific species annotations in csv format. change path to folder for different tables (in find)
# OUTPUT: tables of contig annotations usuable as input to Krona.html

find ../non-redundant_mapping_tables/census_subtype_after_filtering -name "*_ncbi_retrieved_af.csv" -type "f" | xargs -n 1 -I {} -P 1 python3.7 krona_table_generation_ncbi_whole.py {} {}_krona_categories.csv
