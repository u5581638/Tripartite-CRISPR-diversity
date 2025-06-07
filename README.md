# Tripartite-CRISPR-diversity
This repository contains a copy of the workflow used for surveying the intra-subtype and spacer acquisition diversity across CRISPR-Cas subtypes. The workflow functions as 3 seperate yet inter-dependent layers. In the top layer, a computational pipeline mines and extracts DNA upstream and downstream of CRISPR-arrays, using assembled prokaryotic metagenomic and NGS sequencing data as the main feedstock. Genes were then predicted from this data, clustered into putative families, and assessed for their degree of co-association with the arrays. The categorical conservation of different CRISPR-associated genes, based on annotation by domain homology was surveyed to identify the fraction of associated genes whose function remains unknown. This was also performed at the CRISPR-Cas subtype level, where the identities of known and unknown co-associated genes, were examined. 

Note: Many of the scripts and code were originally run from a single fixed directory (/g/data/va71/crispr_pipeline_annotation) on a computing cluster. These have not yet been adapted for standalone use. As such running from each chapter specific folder may produce dependency and pathing errors. I will work to fix these issues progressively in a separate branch.
![image](https://github.com/user-attachments/assets/fb1c04dc-6cf1-4c80-9c6c-347d942e59e4)

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

This code was originally designed to be run as a single integrated workflow, however this approach has not yet been implemented. An executable version of this workflow, which allows all the required file paths to be specified at input, is in development.
