

#	$1	Parent .SNPs 
#	$2
#	$3	output

import sys
import csv
import os
from subprocess import call, check_output
from random import randint as rn




parent1 = sys.argv[1]
parent2 = sys.argv[2]
bedmask = sys.argv[3]
fileout = sys.argv[4]

windowsize=10**4	# shuffle a 10kb window
L2=23011544	#Size of the 2L chromosome
numwins = 5

phial = open(bedmask, 'w')

for i in range(0, numwins):
	start = rn(1, L2-windowsize-1)
	phial.write('chr2L\t%s\t%s\n'%tuple([start, start+windowsize]))
phial.close()

call('bedtools intersect -wa %s -b %s > %s.unsort'%tuple([parent1, bedmask, fileout]), shell=True)
call('bedtools intersect -wa %s -b %s >> %s.unsort'%tuple([parent2, bedmask, fileout]), shell=True)
call('bedtools sort -sizeA -i %s.unsort > %s'%tuple([fileout, fileout]),shell=True)
call('rm %s.unsort'%tuple([fileout]))






