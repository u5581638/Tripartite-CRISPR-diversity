from Bio import SeqIO
import sys
import csv


def rejoiner (split_ids):
	i = 0 
	ret_splits = ""
	ret_splits = "" + str(split_ids[1])
	while (i < len(split_ids)):
		ret_splits = ret_splits + " " + split_ids[i]
		i += 1
#	print("split_id_is:")
#	print(split_ids)
		
	return ret_splits	

# function to parse output from hhblits to table (first used on pdb70).
def pdb70_parse(sequence_url):

	sequences = SeqIO.parse(sequence_url, "fasta")

	outfile = open(sequence_url + "hits.csv", "a")
	spam_writer = csv.writer(outfile)

	ret_list = [] 
	for sequ in sequences:
		my_sequ_id = sequ.id
		my_sequ_description = sequ.description
		my_sequ_description = my_sequ_description.split(";")
		my_sequ_description = my_sequ_description[0]
		my_sequ_description = rejoiner(my_sequ_description.split(" "))
	#	print("The sequence is:")
		my_sequ_sequence = str(sequ.seq)
	#	print(my_sequ_sequence)
		probability = my_sequ_sequence.split("E-value=")
		probability = probability[0]
		probability = probability.split("Probab=")
		probability = probability[1]
		e_value = my_sequ_sequence.split("Score=")
		e_value = e_value[0]
		e_value = e_value.split("E-value=")
		e_value = e_value[1]
		score = my_sequ_sequence.split("Aligned_cols=")
		score = score[0]
		score = score.split("Score=")
		score = score[1]
		similarity = my_sequ_sequence.split("Sum_probs=")
		similarity = similarity[0]
		similarity = similarity.split("Similarity=")
		similarity = similarity[1]
	#	print(probability)
	#	print(score)
	#	print(similarity)
		ret_list.append([my_sequ_id, my_sequ_description, probability, e_value, score, similarity])
		spam_writer.writerow([my_sequ_id, my_sequ_description, probability, e_value, score, similarity])
	outfile.close()	
	return	ret_list
def pdb70_parse_pfam(sequence_url):

	sequences = SeqIO.parse(sequence_url, "fasta")

	outfile = open(sequence_url + "hits.csv", "a")
	spam_writer = csv.writer(outfile)

	ret_list = [] 
	for sequ in sequences:
		my_sequ_id = sequ.id
		my_sequ_description = sequ.description
		my_sequ_description = my_sequ_description.split(" ; ")
		my_sequ_description = rejoiner(my_sequ_description)
	#	print("The sequence is:")
	#	print(my_sequ_description)
		my_sequ_sequence = str(sequ.seq)
	#	print(my_sequ_sequence)
		probability = my_sequ_sequence.split("E-value=")
		probability = probability[0]
		probability = probability.split("Probab=")
		probability = probability[1]
		e_value = my_sequ_sequence.split("Score=")
		e_value = e_value[0]
		e_value = e_value.split("E-value=")
		e_value = e_value[1]
		score = my_sequ_sequence.split("Aligned_cols=")
		score = score[0]
		score = score.split("Score=")
		score = score[1]
		similarity = my_sequ_sequence.split("Sum_probs=")
		similarity = similarity[0]
		similarity = similarity.split("Similarity=")
		similarity = similarity[1]
	#	print(probability)
	#	print(score)
	#	print(similarity)
		ret_list.append([my_sequ_id, my_sequ_description, probability, e_value, score, similarity])
		spam_writer.writerow([my_sequ_id, my_sequ_description, probability, e_value, score, similarity])
	outfile.close()	
	return	ret_list	