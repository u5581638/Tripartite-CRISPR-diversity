# Tripartite-CRISPR-diversity
Code for surveying the intra-subtype and spacer acquisition diversity across CRISPR-Cas subtypes. This code was originally designed to be run as a single integrated workflow, however this approach was not implemented and most code was instead run seperately. As a result, some scripts contain functions and unused sections of code which were not utilised. These will be progressively deleted.

Note: Many of the scripts and code were originally run from a single fixed directory (/g/data/va71/crispr_pipeline_annotation) on a computing cluster. These have not yet been adapted for standalone use. As such running from each chapter specific folder may produce dependency and pathing errors. I will work to fix these issues progressively in a seperate branch.

The code utilised the following dependencies:

Standalone Programs:

python3		3.10.4

prodigal	 2.6.3

GeneMarkS	 1.14	(optional-unused)

PILER-CR	 1.06

CRISPRdetect	 2,2

CRISPRleader	 1.0.3

mmseqs2	 b06bee91ea69c3e41f8479f69651e5c2068f3fe7©

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

Python modules:

Biopython	 1.78	

matplotlib  3.5.3

plotly	 5.18.0

pandas	 2.1.4

R packages:

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

