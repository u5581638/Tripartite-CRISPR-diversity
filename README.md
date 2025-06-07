# Tripartite-CRISPR-diversity
This repository contains a copy of the workflow used for surveying the intra-subtype and spacer acquisition diversity across CRISPR-Cas subtypes. The workflow functions as 3 seperate yet inter-dependent layers. These layers comprise the following overall steps:

**Layer 1:**
1. Identify CRISPR-arrays/ or other distinguishing motifs, in a 10TB data block comprised of assembled metagenomic and NGS sequencing data.
2. Extract DNA upstream and downstream of the arrays. Predict Open Reading Frames (ORFs). Cluster into putative families
3. Verify that co-encoded genes are associated with CRISPR-arrays, by computing CRISPRicity, abundance and average distance from the CRISPR-array for a representative from each cluster.
4. Survey the space of co-associated genes using HMM based homology to determined the fraction of uncategorised genes.

**Layer 2:**

5. Retrieve CRISPR-Cas operon encoding contigs for each CRISPR-Cas subtype (base on the presence or absence of a signature identifying gene).
6. Expand and verify the contig encoding CRISPR-arrays, and map the spacers form these contigs onto the 10TB data block.
7. Use mapper-spacing data to analyse biases in the distribution of spacers during spaer acquisition.

**Layer 3:**

8. Retrieve the spacer-mapped target sequences. These are usually Phage or plasmid Mobile Genetic Elements (MGEs). These may be prophage encoded.
9. Retrieve and prediction contig DNA annotations describing whether contigs are linear/circular and host/phage encoded. Also performs ORF annotation using PFAM and a customised anti-phage defence database, referred to as "DEFLOC".
10. Use host and mapped MGE contigs to generate host and MGE gene cluster networks. These networks are then combined to generate host-phage interaction networks. Concurrently, heatmap describing the composition of each community (local cluster) are also generated.

![Overall workflow no numbering](https://github.com/user-attachments/assets/299c98a4-9d6c-4bb2-8fbd-5bb6d43f4239)


# Installation

The code requires the following dependencies to be installed. :

**Standalone Programs:**

python3		3.10.4

prodigal	 2.6.3

GeneMarkS	 1.14	(optional-unused)

PILER-CR	 1.06

CRISPRdetect	 2,2

CRISPRleader	 1.0.3

mmseqs2

BLAST  2.6.0

hhblits	 3.3.0

HMMER	 3.3.2

virsorter	 2.2.4

PlasME	 1.0

samtools	 1.10

seqkit	 2.4.0

eggnog	 2.1.10	(optional-unused)

viennafold	 2.5.1	(optional-unused)

vmatch	 2.3.0 (optional-unused)

IQtree2	 2.2.2.6

DIAMOND  2.0.11

vConTACT2	 0.11.3

clustal omega	 1.2.4

**Python modules:**

Biopython	 1.78	

matplotlib  3.5.3

plotly	 5.18.0

pandas	 2.1.4

**R packages:**

ggplot2	  3.5.1

pheatmap	 1.0.12

ggtree	 3.13.1

igraph	 2.0.3

r2r	 0.1.1

stringr	 1.5.1

randomcoloR	 1.1.0.1

dplyr	 1.1.4

hrbrthemes	 0.8.7

EMT	 1.3.1

ggpubr  0.6.0

ggprism  1.0.5




# Running Workflows

Within the repository, the files which were used to run the workflow are organised into 3 chapters. These have a relation to, but not a one-to-one correspondance to the layers described above. While the files stored in the folder named "Chapter 3" correspond to the 1st layer of the proccess, the files stored in "Chapter 5" refer exclusively the step (7.), which performs an analysis of the aggregate distribution of mapped spacers in cases where two or more spacers are mapped  to the same target sequence. The files stored in "Chapter 4" are used for running all of the remaining steps of the second and third layers described above. There is a significant degree of interdependency terms of input and output immediate files required, between the general steps described both within, and between each layer.

To run the workflow, each of the dependencies described above must be pre-installed. For each Chapter, a specification file detailing the exact steps run is given. This file includes a description of each step, as well as the script in the repository that was used to run. Each python, BASH and R file within the repository mentioned within the specification also includes a description of the input and output files required. For most files an example input and output file, which was used to run the workflow, is also given. Many of the exact filenames used will change as specified by the user, but the addition of these examples is useful, as the file extensions can be used to trace many of the intermediate files generated while running the workflow.

Because many of the scripts and code were originally run from a single fixed directory on a computing cluster, some of the directory paths specified within each script may need to be changed to run correctly for an external user. 

This code was originally designed to be run as a single integrated workflow, however this approach has not yet been implemented. A executable version of this workflow, which allows all the required file paths to be specified at input, is in development.
