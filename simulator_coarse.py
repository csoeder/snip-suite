

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

breakpoint = rn(1, L2)

call('samtools view -h %s chr2L:1-%s > %s.sam'%tuple([parent1, breakpoint, fileout]), shell=True)
call('samtools view %s chr2L:%s-%s >> %s.sam'%tuple([parent2, L2, breakpoint, fileout]), shell=True)
call('samtools view -Shb %s.sam > %s.bam'%tuple([fileout, fileout]), shell=True)
call('samtools sort -f %s.bam %s.sort.bam'%tuple([fileout, fileout]), shell=True)
call('rm %s.sam'%fileout)
call('rm %s.bam'%fileout)




