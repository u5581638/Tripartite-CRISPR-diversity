#!/bin/bash

#parallel bash script runner
sequence=($( basename $1 ))
database=($( basename $2 ))
perc_id=$3
db_directory=$4
out_dir=$5
query_cover=$6
formatting=$7
# need to change the query
blastn -query $db_directory"spacer_distribution_analysis/"$sequence -db "/g/data/va71/labelled_genomes/"$database  -outfmt $7 -perc_identity $3 -max_target_seqs 100000000 -max_hsps 1 -qcov_hsp_perc $6  >> $db_directory$out_dir$sequence$database"_hits.csv" # need to optimise for spacers!!
