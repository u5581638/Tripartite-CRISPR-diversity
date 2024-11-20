# progrma to deduplicate identical hitmap hit.
import csv
import sys

def deduplicate(spacer_filtered_url, output_url):
	hit_table_url = open(spacer_filtered_url, "r")
	hit_table = csv.reader(hit_table_url)
	hit_dict = {}
	# skip the header
	next(hit_table)
	for row in hit_table:
		if (((row[1] + row[8] + row[9]) not in hit_dict) or ((row[1] + row[9] + row[8]) not in hit_dict)):
			hit_dict[(row[1] + row[8] + row[9])] = row		
	hit_table_url.close()
	out_table = open(output_url, "w")
	spam_writer = csv.writer(out_table)
	for row in hit_dict.values():
		spam_writer.writerow(row)
	out_table.close()	
	return 0
