# proccess sequences in batches of 30-50 for pipeline running
# according to previous KSU usage, 40 sequences should be optimal!!

from Bio import SeqIO
import sys
 
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
