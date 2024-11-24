#!/bin/bash

#PBS -P va71
#PBS -l walltime=48:00:00
#PBS -l mem=250GB
#PBS -l ncpus=52
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71
#PBS -q normalsr

module load blast
conda init
source activate vContact2
conda activate vContact2
export PATH=$PATH:/g/data/va71/Prodigal-GoogleImport
export PYTHONPATH=$PYTHONPATH:/g/data/va71/vConTACT2/MAVERICLab-vcontact2-c0413a6c92e8/build/lib

input_dir=/g/data/va71/vConTACT2/MAVERICLab-vcontact2-c0413a6c92e8/bin/cas13a/phage/
input_protein=cas13a.fasta_all_hits.csv_genomes.fasta_crisprs.lst_full_real_arr_positions.csv_all_hits_blast_filtered_hitmap.csv_standardised.csv_faidx_bp_window.fasta_nr.fa_aa_raw.fasta
output_protein=cas13a_phage_g2g.csv
output_dir=cas13a_phage_network
python3 vcontact2_gene2genome -p $input_dir$input_protein -o $input_dir$output_protein -s Prodigal-FAA
python3 vcontact2 --raw-proteins $input_dir$input_protein --db None --rel-mode Diamond --pcs-mode MCL --pc-inflation 1.2 --pc-evalue 0.001 --proteins-fp $input_dir$output_protein --output-dir $input_dir$output_dir --c1-bin /g/data/va71/vConTACT2/MAVERICLab-vcontact2-c0413a6c92e8/bin/cluster_one-1.0.jar -t 52 -e "cytoscape"