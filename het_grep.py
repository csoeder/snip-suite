'''

A script to look for heterozygoti (ie, incomplete paternal degradation) in crown-sekelsky flies
'''
import matplotlib 
matplotlib.use('agg')
import vcf
import sys
from matplotlib import pyplot as plt
import numpy as np



try:
	cov_thresh = float(sys.argv[1])	#set the threshold for coverage
except IndexError:
	cov_thresh = 10

het_dict = {}
agree_dict = {}

for i in range(1,10):
	het_dict[i] = []
	agree_dict[i] = []
	parser = vcf.Reader(open('REC%s/REC%s.vcf'%tuple([i,i]), 'r'))
	for rec in parser:
		if rec.is_snp and float(rec.INFO['DP']) > cov_thresh:
			agree_dict[i].append(float(sum(rec.INFO['DP4'][-2:]))/sum(rec.INFO['DP4']))
			het_dict[i].append(rec.heterozygosity)


plt.hist(agree_dict.values(), histtype='barstacked', color=['r', 'g', 'b', 'm', 'c', 'y', 'k', 'r', 'g'])
plt.show()




