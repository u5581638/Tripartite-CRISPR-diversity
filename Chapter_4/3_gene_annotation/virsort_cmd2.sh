#!/bin/bash

# my_dir=($( pwd ))
# cd /g/data/va71/my_conda/new_conda/etc/profile.d
# source conda.sh
# /g/data/va71/my_conda/new_conda/bin/conda activate vs2
# cd $my_dir
# echo "Yes!!"
virsorter run -i $1 -w $2 --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae # shouldn't need to add a db. Does virsorter make it's own folder?
# /g/data/va71/my_conda/new_conda/bin/conda deactivate