# protein_reconciliation

# this file contains two functions to change the formatting of MetaGeneMark and GeneMarkS2 to be the same as prodigal.
# Note, these functions did not end up being used as all protein prediction was performed using prodigal.
import sys
from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord

def genemarkS2_reformatting (genemark2_orf_url):
	
	predicted_sequences = SeqIO.parse(genemark2_orf_url, "fasta")
	previous_id = ""
	protein_number = 1
	gene_number = 1
	ret_proteins = []
	for protein in predicted_sequences:
		protein_identifier = protein.description.split(" ")
		index_number = protein_identifier[0]
		protein_id = protein_identifier[1]
		protein_start = protein_identifier[2]
		protein_end = protein_identifier[3]
		protein_sense = protein_identifier[4]
		new_sequence = Seq(str(protein.seq))

		prodigal_seqRecord = SeqRecord(new_sequence, id=protein_id + "_" + str(protein_number),description="# " + protein_start + " # " + protein_end + " # " + protein_sense + " # " + "ID=" + str(gene_number) + "_" + str(protein_number) + ";partial=NA;start_type=NA;rbs_motif=AGG;rbs_spacer=NA;gc_cont=NA",name="")   
		ret_proteins.append(prodigal_seqRecord)

		if (protein_id == previous_id):
			protein_number += 1
		else:
			protein_number = 0	
			gene_number += 1
		previous_id = protein_id
	
	SeqIO.write(ret_proteins, genemark2_orf_url + "_geneS2prod_reformatted.fasta", "fasta")
	return ret_proteins 

def metagenemark_reformatting (metagenemark_orf_url):

	predicted_sequences = SeqIO.parse(metagenemark_orf_url)
	previous_id = ""
	protein_number = 1
	gene_number = 1
	for protein in predicted_sequences:
		fore_details = protein.id.split("|")  
		protein_id = protein.description.split(" \t>") [1].split(" ") [0] # Might be something incompatible at the end of the file
		protein_start = fore_details[4]
		protein_end = fore_details[5]
		protein_sense = fore_details[3]
		new_sequence = Seq(str(protein.seq))
		prodigal_seqRecord = SeqRecord(new_sequence, id=protein_id + "_" + str(protein_number),description=" # " + protein_start + " # " + protein_end + " # " + protein_sense + " # " + "ID=" + str(gene_number) + "_" + str(protein_number) + ";partial=NA;start_type=NA;rbs_motif=AGG;rbs_spacer=NA;gc_cont=NA",name="")   
		ret_proteins.append(prodigal_seqRecord)

		if (protein_id == previous_id):
			protein_number += 1
		else:
			protein_number = 0	
			gene_number += 1
		previous_id = protein_id
	ret_proteins = []
	SeqIO.write(ret_proteins, genemark2_orf_url + "_metagen2prod_reformatted.fasta", "fasta")

	return ret_proteins

# genemarkS2_reformatting(sys.argv[1])	