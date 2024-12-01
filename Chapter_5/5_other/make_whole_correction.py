# readjust the rows of the table

import sys
import csv

ret_in = open(sys.argv[1],"r")
csvfile = csv.reader(ret_in)
ret_out = open(sys.argv[2],"w")
spamwriter = csv.writer(ret_out)

for row in csvfile:
	my_row = row 
	if (len(my_row) > 30):
		my_row.pop(29)
	spamwriter.writerow(my_row)
ret_in.close()
ret_out.close()	
