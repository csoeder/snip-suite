##########################################################################
#	Scans the alignment, pulls out sequences with more than 1/3 clipped. 	#
##########################################################################

import csv
import sys
import re

BED_in = sys.argv[1]
BED_out = sys.argv[2]
thresh = 0.25


phial = open(BED_out, 'w')

with open(BED_in, 'rb') as csvfile:
		spamreader = csv.reader(csvfile, delimiter='\t')
		for row in spamreader:
			cigar = row[6]
			seq_len = float(sum(map(int, re.findall(r"([0-9]+)[MISX=]", cigar))))
			clips = sum(map(int, re.findall(r"([0-9]+)[SH]", cigar)))
			if clips/seq_len > thresh:
				phial.write('%s\n'% tuple(['%s\t'*len(row)%tuple(row)]))

phial.close()

