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
	print name
	for chro in snp_dict.keys():
		snp_rho[chro] = []
		start = 0
		snps = np.array(snp_dict[chro])
		while start < chrom_lens[chro]:
			snp_rho[chro].append(len(snps[snps>start][snps[snps>start]<start+window_size]))
			start += 100

		num, bins, patches = plt.hist(snp_rho[chro], bins=49, histtype='step', color=color)
		mux = max(mux, max(num))
		tot = float(sum(num))
		adequate = float(sum(num[bins[:-1]+0.5*np.diff(bins)>critical_thresh]))
#		print adequate
#		print len(num), num
#		print len(np.diff(bins)), np.diff(bins)
#		print chro, adequate
		print "%s:	%.2f%%	adequate;	%.2f%% sparse\n"%tuple([chro, 100*adequate/tot, 100-100*adequate/tot])


	plt.plot([],[], color, label=name)
	print
	return mux


m1 = histogrammatical(Parent1, Parent1.split('/')[-1].split('.SNPS')[0], 'b')
m2 = histogrammatical(Parent2, Parent2.split('/')[-1].split('.SNPS')[0], 'r')


plt.vlines(critical_thresh, 0, max(m1,m2)/2, 'k', label='threshold')
plt.legend()
plt.title('Diagnostic SNPs per %0.1e Base Pair Window'%tuple([window_size]))
plt.xlabel('# SNPs per Window')
plt.ylabel('# Windows')
plt.show()
plt.savefig('Parental_SNP_Densities.png')








