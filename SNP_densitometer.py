import csv
import sys
import numpy as np
import matplotlib 
matplotlib.use('agg')
from matplotlib import pyplot as plt
from Bio import SeqIO



Parent1 = sys.argv[1]
Parent2 = sys.argv[2]
ref_dict = sys.argv[3]

window_size = 10**3	#looking at a kbp window at a time
critical_thresh = 4	#we want 


def histogrammatical(input_SNPs, name, color ):

	snp_dict = {}
	with open(input_SNPs, 'rb') as csvfile:
		spamreader = csv.reader(csvfile, delimiter='\t')
		for row in spamreader:
			if row[0] in snp_dict.keys():
				snp_dict[row[0]].append(int(row[1]))
			else:
				snp_dict[row[0]] = [int(row[1])]

	chrom_lens = {}
	for record in SeqIO.parse(ref_dict, "fasta"):
		if record.id in snp_dict.keys():
			chrom_lens[record.id] = len(record)


	snp_rho = {}
	mux = 0
	for chro in snp_dict.keys():
		snp_rho[chro] = []
		start = 0
		snps = np.array(snp_dict[chro])
		while start < chrom_lens[chro]:
			snp_rho[chro].append(len(snps[snps>start][snps[snps>start]<start+window_size]))
			start += 100

		num, bins, patches = plt.hist(snp_rho[chro], bins=25, histtype='step', color=color)
		mux = max(mux, max(bins))

	plt.plot([],[], color, label=name)
	return mux


m1 = histogrammatical(Parent1, Parent1.split('/')[-1].split('.SNPS')[0], 'b')
m2 = histogrammatical(Parent2, Parent2.split('/')[-1].split('.SNPS')[0], 'r')


plt.vlines(critical_thresh, 0, max(m1,m2), 'k', label='threshold')
plt.legend()
plt.show()
plt.savefig('Parental_SNP_Densities.png')








