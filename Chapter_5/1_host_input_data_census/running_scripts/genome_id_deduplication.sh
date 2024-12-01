module load python3/3.10.4

# deduplicate genome_ids for census among all subtypes

find input_files/ -name "*_filtered.csv" -type "f" -exec basename {} \; | xargs -n 1 -I {} -P 1 python3 subtype_id_names_as_set.py input_files/{} mapped_spacer_genome_ids/{}"_host_ids.csv"