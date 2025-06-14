#!/bin/bash
# census run batch script

# JGI retrieval run below
# INPUT:vJGI tables from each subtype containing specific species annotations in csv format.
# OUTPUT: tables of contig annotations usuable as input to Krona.html


find ../non-redundant_mapping_tables/census_subtype_after_filtering -name "*_jgi_annotations_after_filtering.csv" -type "f" | xargs -n 1 -I {} -P 1 python3.7 krona_table_generation_jgi.py {} {}_krona_categories.csv
