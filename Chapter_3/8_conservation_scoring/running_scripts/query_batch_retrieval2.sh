#!/bin/bash

python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12b.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12h.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12h.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12i.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12i.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12j.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12j.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12k.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12k.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13d.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas13d.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13a.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/spcas9.fasta__annotations_all.csv_redund.csv aa_raw_proteins/spcas9.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIF_csy1.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIF_csy1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIA_cas10.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIIIA_cas10.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeVU1.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeVU1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta CARF3..CARF
mkdir conserved_hit_phylo/CARF3..CARF 
mv aa_non_redundant/*_phylo_seqs.fa conserved_hit_phylo/CARF3..CARF/ 

#!/bin/bash

python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12b.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12h.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12h.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12i.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12i.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12j.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12j.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12k.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12k.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13d.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas13d.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13a.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/spcas9.fasta__annotations_all.csv_redund.csv aa_raw_proteins/spcas9.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIF_csy1.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIF_csy1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIA_cas10.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIIIA_cas10.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeVU1.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeVU1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta DrmB
mkdir conserved_hit_phylo/DrmB 
mv aa_non_redundant/*_phylo_seqs.fa conserved_hit_phylo/DrmB/ 

#!/bin/bash

python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12b.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12h.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12h.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12i.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12i.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12j.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12j.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12k.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12k.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13d.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas13d.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13a.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/spcas9.fasta__annotations_all.csv_redund.csv aa_raw_proteins/spcas9.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIF_csy1.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIF_csy1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIA_cas10.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIIIA_cas10.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeVU1.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeVU1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta retron_type_VI_hth_cluster2_1
mkdir conserved_hit_phylo/retron_type_VI_hth_cluster2_1 
mv aa_non_redundant/*_phylo_seqs.fa conserved_hit_phylo/retron_type_VI_hth_cluster2_1/ 

#!/bin/bash

python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12b.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12h.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12h.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12i.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12i.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12j.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12j.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12k.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12k.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13d.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas13d.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13a.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/spcas9.fasta__annotations_all.csv_redund.csv aa_raw_proteins/spcas9.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIF_csy1.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIF_csy1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIA_cas10.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIIIA_cas10.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeVU1.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeVU1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta IetA
mkdir conserved_hit_phylo/IetA 
mv aa_non_redundant/*_phylo_seqs.fa conserved_hit_phylo/IetA/ 

#!/bin/bash

python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12b.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12h.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12h.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12i.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12i.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12j.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12j.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12k.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12k.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13d.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas13d.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13a.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/spcas9.fasta__annotations_all.csv_redund.csv aa_raw_proteins/spcas9.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIF_csy1.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIF_csy1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIA_cas10.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIIIA_cas10.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeVU1.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeVU1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta SspH_00003
mkdir conserved_hit_phylo/SspH_00003 
mv aa_non_redundant/*_phylo_seqs.fa conserved_hit_phylo/SspH_00003/ 

#!/bin/bash

python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12b.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12h.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12h.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12i.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12i.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12j.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12j.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12k.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12k.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13d.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas13d.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13a.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/spcas9.fasta__annotations_all.csv_redund.csv aa_raw_proteins/spcas9.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIF_csy1.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIF_csy1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIA_cas10.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIIIA_cas10.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeVU1.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeVU1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Mtase_I
mkdir conserved_hit_phylo/Mtase_I 
mv aa_non_redundant/*_phylo_seqs.fa conserved_hit_phylo/Mtase_I/ 

#!/bin/bash

python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12b.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12h.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12h.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12i.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12i.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12j.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12j.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12k.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12k.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13d.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas13d.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13a.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/spcas9.fasta__annotations_all.csv_redund.csv aa_raw_proteins/spcas9.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIF_csy1.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIF_csy1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIA_cas10.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIIIA_cas10.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeVU1.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeVU1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta Rease_I
mkdir conserved_hit_phylo/Rease_I 
mv aa_non_redundant/*_phylo_seqs.fa conserved_hit_phylo/Rease_I/ 

#!/bin/bash

python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12b.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12h.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12h.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12i.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12i.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12j.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12j.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12k.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12k.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13d.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas13d.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13a.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/spcas9.fasta__annotations_all.csv_redund.csv aa_raw_proteins/spcas9.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIF_csy1.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIF_csy1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIA_cas10.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIIIA_cas10.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeVU1.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeVU1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEi
mkdir conserved_hit_phylo/AbiEi 
mv aa_non_redundant/*_phylo_seqs.fa conserved_hit_phylo/AbiEi/ 

#!/bin/bash

python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12b.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas12b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12h.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12h.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12i.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12i.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12j.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12j.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12k.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12k.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas12g.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas12g.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13d.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/cas13d.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13a.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/cas13b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/cas13b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/spcas9.fasta__annotations_all.csv_redund.csv aa_raw_proteins/spcas9.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIA_cas7a.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeIB_cas8b.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIF_csy1.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIF_csy1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIA_cas10.fasta__annotations.csv_redund.csv aa_raw_proteins/typeIIIA_cas10.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_annotations.csv_redund.csv aa_raw_proteins/typeIIIB_cmr2.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
python3.7 ../../subtype_tree_proteins_sum.py aa_non_redundant/typeVU1.fasta_all_hits.csv_genomes.fasta_annotations_all.csv_redund.csv aa_raw_proteins/typeVU1.fasta_all_hits.csv_genomes.fasta_aa_raw.fasta AbiEii
mkdir conserved_hit_phylo/AbiEii 
mv aa_non_redundant/*_phylo_seqs.fa conserved_hit_phylo/AbiEii/ 
