# subprogram to cluster sequences in a representative number (assuming partial homology)

from Bio import SeqIO
import sys
import csv
import genemark_sequence_reformatting
import prodigal_genemark_reconciliation
import protein_hit_of_origin 
import genome_retriever_anno_prot 
import rep_cluster_extractor
import genome_identifier_extractor
import genome_retriever_target_protein_adder_rep
import table_identifier_extractor
import subprocess

# This may be mutually exclusive with phage clustering due to the BLAST query requirements!!
def rep_cluster (b, min_seq_id, db_directory_path, a, merge_protein=0,large_dataset=0,phylogeny_switch=0, is_phage=0):
	subprocess.run(["/g/data/va71/crispr_pipeline_annotation/prodigal/prodigal", "-i", b,"-o", b + "_prot.txt","-a" , b + "_trans_aa.fa"  ])
	if (merge_protein == 1):
		subprocess.run(["/g/data/va71/crispr_pipeline_annotation/annotation_database/annotation_upgrades/protein_reconciliation/genemark_reinstall/gms2_linux_64/gms2.pl", "--seq", b ,"--genome-type", "auto", "--out", b + "_trans_aa_gms2.out", "--faa", b + "_trans_aa_gms2.faa"])
		genemark_sequence_reformatting.genemarkS2_reformatting(b + "_trans_aa_gms2.faa")
		prodigal_genemark_reconciliation.union(b + "_trans_aa_gms2.faa" + "_geneS2prod_reformatted.fasta", b + "_trans_aa.fa", b + "_trans_aa.fa") # This will hopefully override the original protein file
	subprocess.run(["makeblastdb", "-in",b + "_trans_aa.fa", "-dbtype", "prot"] )
	subprocess.run(["blastp", "-query", file_path , "-db", b + "_trans_aa.fa", "-evalue", "0.00001", "-max_target_seqs", "1000000", "-outfmt", "10", "-out", b + "_trans_aa.fa" + "_prot.hits.csv", "-max_hsps", "1"])
	anno_prot_list_url = protein_hit_of_origin.anno_list_from_hits(b + "_trans_aa.fa", b + "_trans_aa.fa" + "_prot.hits.csv")
		
	homologs = genome_retriever_anno_prot.genome_retriever(b + "_trans_aa.fa", b + "_trans_aa.fa" + "_prot.hits.csv")
	SeqIO.write(homologs, b + "_trans_aa.fa" + "_cluster_specific.fasta", "fasta")
	if (large_dataset == 1):
		subprocess.run(["/g/data/va71/mmseqs/bin/mmseqs", "createdb", b + "_trans_aa.fa" + "_cluster_specific.fasta", b + "clusterdb" ])
			# identity cut-off must be low enough to minimise the number of sequences used for downstream MDS
		subprocess.run(["/g/data/va71/mmseqs/bin/mmseqs", "cluster", b + "clusterdb", b + "clusterdb_cluster", db_directory_path + "tmp", "--cov-mode", "0","-c" ,"0.6", "--min-seq-id", min_seq_id, "--cluster-mode", "1"])
		subprocess.run(["/g/data/va71/mmseqs/bin/mmseqs", "createsubdb", b + "clusterdb_cluster", b + "clusterdb", b + "clusterdb_cluster_rep"])
		subprocess.run(["/g/data/va71/mmseqs/bin/mmseqs", "convert2fasta", b + "clusterdb_cluster_rep", b + "_trans_aa.fa" + "_cluster_specific.fasta"]) # need to make sure this line works. Otherwise the testing input will be overwritten

	subprocess.run(["/g/data/va71/crispr_pipeline_annotation/clustalo", "-i", b + "_trans_aa.fa" + "_cluster_specific.fasta", "-o", b + "_trans_aa.fa" + "_cluster_specific.fasta" + "_aligned.fa"])
	subprocess.run(["Rscript", "mds_calculator.r", b + "_trans_aa.fa" + "_cluster_specific.fasta" + "_aligned.fa"]) # need to have bios2mds (or equivalent) installed. May be better/more compatible to use an alternative package?
	#	may as well plot the mds coordinates in a post script file. Add a cmd to do this here!!	
	if (phylogeny_switch == 1):
		phylo_renamer.rename(b + "_trans_aa.fa" + "_cluster_specific.fasta" + "_aligned.fa")
		subprocess.run(["/g/data/va71/FastTree " + b + "_trans_aa.fa" + "_cluster_specific.fasta" + "_aligned.fa" + "_phylo_relabelled.fasta " + " > " + b + "_trans_aa.fa" + "_cluster_specific.fasta_tree.newick"  ], shell=True)
		# need to export table containinf MDS values
	representative_cluster_ids = rep_cluster_extractor.rep_sequence_determinant(b + "_trans_aa.fa" + "_cluster_specific.fasta" + "_aligned.fa" + "_x_y_coords.csv")
		# need to use the protein ids to lookup the proteins then find just the subset of genomes corresponding to these proteins!
	representative_genomes = genome_identifier_extractor.genomes_from_id(representative_cluster_ids, b)
	anno_genomes = genome_retriever_target_protein_adder_rep.protein_adder(anno_prot_list_url, representative_genomes)		
	SeqIO.write(anno_genomes, b + "_rep_genomes.fasta", "fasta")

	# want to append to a specialised rep-genomes_cluster. This is used for annotation. However, this will cause problems!! Really want to only take the cases where spacer mapping occurs

	#	with open(file_path + "_all_hits.csv", "r") as csvfile:
				# believe block_dict is actually a duplication of these steps
				# tblastn_hits = list(csv.reader(csvfile))
	representative_table = table_identifier_extractor.hits_from_genomes(representative_genomes, a, False)
	
	return representative_table # ? not sure, may be best to omitt open statement
