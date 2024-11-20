#!/bin/bash

#parallel bash script runner
sequence=($( basename $1 ))
database=($( basename $2 ))
# need to change the query
tblastn -query "../"$sequence -db "/g/data/va71/labelled_genomes/"$database  -outfmt 10 -evalue 0.0000001 -max_target_seqs 100000000 -max_hsps 10 >> ../co_occurrance_bug_diagnostics_results/$sequence$database"_hits.csv" # need to check redirect to stdout as csv works!
