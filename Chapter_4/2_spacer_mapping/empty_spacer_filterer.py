# empty_spacer_filterer.py

import csv
import sys
from Bio import SeqIO

# remove any entries without spacers
def remove_blanks(input_url):
	with open(input_url, "r") as csvfile:
		hit_table = list(csv.reader(csvfile))

	csvfile2 = open(input_url + "_no_blanks.csv","w")
	spam_writer = csv.writer(csvfile2)

	for hit in hit_table:
		if (hit[14] != ''):
			spam_writer.writerow(hit)
	csvfile2.close()
	return 0