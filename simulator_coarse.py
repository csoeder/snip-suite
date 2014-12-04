

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

breakpoint = rn(1, L2)

call('samtools view -h %s 2L:1-%s > %s.sam'%tuple([parent1, breakpoint, fileout]), shell=True)
call('samtools view %s 2L:%s-%s >> %s.sam'%tuple([parent2, L2, breakpoint, fileout]), shell=True)
print "print assembled"
call('samtools view -Shb %s.sam > %s.bam'%tuple([fileout, fileout]), shell=True)
print "converted"
call('samtools sort %s.bam %s.sort'%tuple([fileout, fileout]), shell=True)
print "sorted"
call('rm %s.sam'%fileout)
call('rm %s.bam'%fileout)




