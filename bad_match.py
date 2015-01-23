##########################################################################
#	Scans the alignment, pulls out sequences with more than 1/3 clipped. 	#
##########################################################################

import csv
import sys
import re

SAM_in = sys.argv[1]
SAM_out = sys.argv[2]
thresh = 1.0/3


phial = open(SAM_out, 'w')

with open(SAM_in, 'rb') as csvfile:
		spamreader = csv.reader(csvfile, delimiter='\t')
		for row in spamreader:
			if row[0][0] == '@':	#	rewrite header lines
				phial.write('%s\n'% tuple(['%s\t'*len(row)%tuple(row)]))	
			else:
				seq_len = float(len(row[9]))
				cigar = row[5]
				clips = sum(map(int, re.findall(r"([0-9]+)[SH]", cigar)))
				if clips/seq_len > thresh:
					phial.write('%s\n'% tuple(['%s\t'*len(row)%tuple(row)]))

phial.close()

