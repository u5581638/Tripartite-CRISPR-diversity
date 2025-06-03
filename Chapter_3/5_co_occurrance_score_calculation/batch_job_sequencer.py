# proccess sequences in batches of 30-50 for pipeline running
# according to previous KSU usage, 40 sequences should be optimal!!

from Bio import SeqIO
import sys

# INPUT: file containing input representative protein sequences to co-occurrance calculation (in FASTA format)
# i.e. representative_distance_5000_10000avg_dist_position_n_distance_relabelled_grt_3_grt_3_cluster_declustered_one_member_gt_300_renumbered_sorted_descending_deduplicated_flattened_all_sequences.fasta_p_formatted.fastanot_in_crisprs.fasta_labelled.fasta
# OUTPUT: FASTA file comprising 80 sequences labelled by run number (for running tBLASTN to compute co-occurrance/CRISPRicity)
# i.e. 80_sequence_runs/run_1/run_1_run_80_sequence_5_10kb.fasta 
# SHELL: python3 representative_distance_5000_10000avg_dist_position_n_distance_relabelled_grt_3_grt_3_cluster_declustered_one_member_gt_300_renumbered_sorted_descending_deduplicated_flattened_all_sequences.fasta_p_formatted.fastanot_in_crisprs.fasta_labelled.fasta
# Note: the directories containing the runs must be created prior to runnning this program.
sequences = list(SeqIO.parse(sys.argv[1], "fasta"))
i=0
run_number = 0
while (i < len(sequences)):
	k = 0 
	ret_list = []
	while (k < 80):
#		print(k)
		if (i + k == len(sequences)):
			break
		ret_list.append(sequences[k + i])
		k += 1
	SeqIO.write(ret_list, "80_sequence_runs/" + "run_" + str(run_number) + "/" + "run_" + str(run_number) + "_run_" +"80_sequence_5_10kb.fasta", "fasta")
	run_number += 1	
	i += 80
