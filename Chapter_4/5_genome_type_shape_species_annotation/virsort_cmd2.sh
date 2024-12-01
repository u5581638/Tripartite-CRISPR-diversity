#!/bin/bash

# helper script to run virsorter2.
virsorter run -i $1 -w $2 --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae -j $3 --tmpdir $2/tmp --rm-tmpdir # shouldn't need to add a db. Does virsorter make it's own folder?
# /g/data/va71/my_conda/new_conda/bin/conda deactivate