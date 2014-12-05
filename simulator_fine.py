

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
fileout = sys.argv[3]

windowsize=10**4	# shuffle a 10kb window
L2=23011544	#Size of the 2L chromosome
numwins = 5

phial = open('bedmask.bed', 'w')

for i in range(0, numwins):
	start = rn(1, L2-windowsize-1)
	phial.write('2L\t%s\t%s\n'%tuple([start, start+windowsize]))
phial.close()

call('bedtools intersect -abam %s -b bedmask.bed | samtools view -h - > %s.sam '%tuple([parent1, fileout]), shell=True)
call('bedtools intersect -v -abam %s -b bedmask.bed | samtools view - >> %s.sam '%tuple([parent2, fileout]), shell=True)
print "assembled"
call('samtools view -Shb %s.sam > %s.bam'%tuple([fileout, fileout]), shell=True)
print "converted"
call('samtools sort %s.bam %s.sort'%tuple([fileout, fileout]), shell=True)
print "sorted"
call('rm %s.sam'%fileout)
call('rm %s.bam'%fileout)






