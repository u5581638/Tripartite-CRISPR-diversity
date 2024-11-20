 # virsorter_parser.py

import sys
import csv
import pandas

def parse(virsorter_table1_url, virsorter_table2_url, virsorter_table3_url):
	boundary_table = pandas.read_csv(virsorter_table1_url, delimiter='\t')
	boundary_table = boundary_table.sort_values(by=['seqname'])
	extra_info = pandas.read_csv(virsorter_table2_url, delimiter='\t')
	extra_info = extra_info.sort_values(by=['seqname'])
	extra_columns = extra_info[['dsDNAphage','ssDNAphage','NCLDV','RNA','lavidaviridae']]
	boundary_table = boundary_table.join(extra_columns)
	boundary_table.to_csv(virsorter_table3_url)

