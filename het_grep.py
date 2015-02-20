'''

A script to look for heterozygoti (ie, incomplete paternal degradation) in crown-sekelsky flies
'''
import matplotlib 
matplotlib.use('agg')
import vcf
import sys
from matplotlib import pyplot as plt
import numpy as np


file_in = sys.argv[1]
file_out = sys.argv[2]
lower_bound = float(sys.argv[3])
upper_bound = float(sys.argv[4])

try:
	cov_thresh = float(sys.argv[5])	#set the threshold for coverage
except IndexError:
	cov_thresh = 10

print ""

phial = open(file_out, 'w')
agreement = []
locations = []
total_num=0
parser = vcf.Reader(open(file_in, 'r'))
for rec in parser:	#	For each record in the VCF...
	if rec.is_snp and float(sum(rec.INFO['DP4'])) > cov_thresh:	#	if it records a SNP and is of sufficient depth...
		total_num += 1
		fraction_nonref = float(sum(rec.INFO['DP4'][-2:]))/sum(rec.INFO['DP4'])	#	calculate the fraction of reads supporting the non-ref alleles
		if fraction_nonref > lower_bound and fraction_nonref < upper_bound :
			phial.write('%s\t%s\t%s\n'%tuple([rec.CHROM, rec.POS, rec.POS+1]))
			agreement.append(fraction_nonref)
			locations.append(rec.CHROM)
phial.close()
		
locations = np.array(locations)
loc_dict={}
for loc in np.unique(locations):
	loc_dict[loc] = len(locations[locations == loc])


plt.hist(agreement, bins=18)
plt.subplot(211)
plt.title('Non-Homozygous Non-Reference Allele Frequencies\n(%s of %s Total Called SNPs)'%tuple([len(agreement), total_num]))
plt.xlabel('Fraction of Quality Reads Supporting a Non-Ref')
plt.ylabel('Number of SNPs with Measured Frequency')
plt.subplot(212)
plt.pie(loc_dict.values(), label=loc_dict.keys())
plt.legend()
plt.title('Non-Homozygous Non-Reference Alleles by Chromosome')
plt.savefig('heterozygous_stats.png')



# het_dict = {}
# agree_dict = {}

# for i in range(1,10):
# 	het_dict[i] = []
# 	agree_dict[i] = []
# 	parser = vcf.Reader(open('REC%s/REC%s.vcf'%tuple([i,i]), 'r'))
# 	for rec in parser:
# 		if rec.is_snp and float(rec.INFO['DP']) > cov_thresh:
# 			agree_dict[i].append(float(sum(rec.INFO['DP4'][-2:]))/sum(rec.INFO['DP4']))
# 			het_dict[i].append(rec.heterozygosity)


# plt.hist(agree_dict.values(), histtype='barstacked', color=['r', 'g', 'b', 'm', 'c', 'y', 'k', 'r', 'g'])
# plt.show()




