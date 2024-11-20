# genome identifier extractor

from Bio import SeqIO
import sys


def genomes_from_id (protein_ids, genome_url):
	genome_ids = set()
	for protein_id in protein_ids:
		print(protein_id)
	#	prot_id_only = protein_id.split(" # ")
	#	prot_id_only = prot_id_only[0]
	#	prot_id_only = protein_id.split("|")
	#	prot_id_only = prot_id_only[1]
		prot_id_only = protein_id
		prot_id_only = prot_id_only.split("::")
		prot_id_only_num = prot_id_only[1]
		prot_id_only_num = prot_id_only_num.split("_")
		prot_id_only_num = prot_id_only_num[0]
		prot_id_only = prot_id_only[0] + "::" + prot_id_only_num 
		genome_ids.add(prot_id_only)
	
	genomes = SeqIO.parse(genome_url, "fasta")
	ret_genomes = []
	for genome in genomes:
		if (genome.id in genome_ids):
			ret_genomes.append(genome)
	return ret_genomes

def protein_ids_from_id (protein_ids, protein_url):
	protein_ids = {}
	ret_list = []
	proteins = SeqIO.parse(protein_url, "fasta")
	for prot in proteins:
		if (prot.id in protein_ids):
			ret_list.append(prot.description)
	return ret_list		

def proteins_from_id (protein_ids, protein_url):
	protein_ids = {}
	ret_list = []
	proteins = SeqIO.parse(protein_url, "fasta")
	for prot in proteins:
		if (prot.id in protein_ids):
			print("great!!")
			ret_list.append(prot)
	return ret_list