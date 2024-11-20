#!/bin/bash


block_dr_name=$1
working_dir=$2
offset=$3
spacer_only_prefix=${block_dr_name%*_drs_*_reconciled_spacers*}
spacer_only_suffix=${block_dr_name##*_drs_*_reconciled_spacers}
spacer_only_name=${spacer_only_prefix}"_reconciled_spacers"${spacer_only_suffix}

python3 /g/data/va71/crispr_pipeline_annotation/pipe_line_source_files/spacer_mapping_hit_filtration_tony_by_block.py ${working_dir}"co_occurrance_bug_diagnostics_results/"${spacer_only_name} ${working_dir}"co_occurrance_bug_diagnostics_results2/"$block_dr_name ${working_dir}"co_occurrance_bug_diagnostics_filtered_results/"${spacer_only_name} ${offset}

