#!/bin/bash

module load python3/3.10.4

# script to generate plasmid + genome type cluster conservation matrix for host-encoded contigs.

python3 /g/data/va71/crispr_pipeline_annotation/host_phage_interaction_final_7-5-2024/source_code/vCONTACT2_clade_component_annotation_virsort_plasmid_only.py host_only_graphs/casIIIB/host/c1.clusters host_only_graphs/casIIIB/host/final_table.csv host_only_graphs/casIIIB/host/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_p_ovrlps.fa_p.csv host_only_graphs/casIIIB/host/casIIIB_host_gold_annotations.csv host_only_graphs/casIIIB/host/casIIIB_host_ncbi_annotations_corrected.csv host_only_components_vCONTACT/casIIIB_host_only_components.csv host_only_matrix_generation/casIIIB_host_only_virsort_plasmid_matrix.csv
