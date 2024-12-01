#!/bin/bash
# plasme_running_helper_script


export PATH=$PATH:/g/data/va71/crispr_pipeline_annotation/prodigal/prodigal
export PATH=$PATH:/g/data/va71/crispr_pipeline_annotation/plasme/PLASMe/diamond

python3 /g/data/va71/crispr_pipeline_annotation/plasme/PLASMe/PLASMe.py $1 $2 -c 0.6 -i 0.6 -p 0.5 -t $3 --temp $4
python3 /g/data/va71/crispr_pipeline_annotation/plasme_table_parser.py $1 $2 $2"_p.csv" 