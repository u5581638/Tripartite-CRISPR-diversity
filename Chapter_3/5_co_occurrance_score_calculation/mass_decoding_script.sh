#!/bin/bash

#PBS -P va71
#PBS -q normalsr
#PBS -l walltime=48:00:00
#PBS -l mem=250GB
#PBS -l ncpus=8
#PBS -l wd
#PBS -l storage=scratch/xi88+gdata/xi88+gdata/va71


module load python3/3.10.4
cd run_0/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_2/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_3/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_4/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_5/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_6/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_7/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_8/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_9/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_10/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_11/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_12/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_13/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_14/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_15/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_16/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_17/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_18/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_19/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_20/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_21/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_22/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_23/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_24/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_25/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_26/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_27/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_28/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_29/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_30/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_31/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_32/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_33/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_34/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_35/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_36/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_37/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_38/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_39/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_40/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_41/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_42/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_43/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_44/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_45/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_46/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_47/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_48/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_49/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_50/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_51/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_52/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_53/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_54/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_55/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_57/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_58/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_59/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_60/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_61/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_62/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_63/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_64/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_65/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_66/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_67/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_68/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_69/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_70/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_71/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_72/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_73/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_74/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_75/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_76/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_77/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_78/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_80/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_81/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_82/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_83/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_84/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_85/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_86/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_87/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_88/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_89/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_90/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_91/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_92/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_93/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_94/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_95/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_96/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_97/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_98/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_99/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_100/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_101/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_102/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_103/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_104/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_105/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_106/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_107/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_108/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_109/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_110/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_111/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_112/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_113/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_114/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_115/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_116/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_117/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_118/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_119/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_120/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_121/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_122/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_123/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_124/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_125/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_126/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_127/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_128/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_129/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_130/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_131/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_132/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_133/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_134/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_135/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_136/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_137/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_138/
python3 ../sequence_decoder.py all_global_hits.csv
cd ../run_139/
python3 ../sequence_decoder.py all_global_hits.csv


