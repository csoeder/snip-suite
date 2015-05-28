import csv
import sys
import numpy as np
#import matplotlib 
#matplotlib.use('agg')
#from matplotlib import pyplot as plt
#from Bio import SeqIO


SNP_file = sys.argv[1]
color = sys.argv[2]


box_size=10**3
chrom="2L"
#chroms=["2L","2R","3L","3R","4","X", "Y"]


def archimedes(points):

	points = np.array(points)
	points.sort()

	print len(points)
	start, last = 0, points[-1]

	coord, density = [], []

	while start < last:
		coord.append(int(start+0.5*box_size))
		clip1 = points[points<start+box_size]
		clip2 = clip1[clip1>start]
		dens = float(len(clip2))/box_size
		density.append(dens)
		#print dens
		start += box_size	
	return coord, density


def write_to_varwig(coords, density, phial, colour, name):
        vial = open(phial, 'w')
        vial.write('browser position chr%s\n'%chrom)
        vial.write('browser hide all\n')
        vial.write('track type=wiggle_0 name="%s" description="varStep format" visibility=full autoScale=off viewLimits=0:%s color=%s graphType=points priority=20\n'%tuple([name, max(density), colour]))
        vial.write('variableStep chrom=chr%s\n'%tuple([chrom]))
        for pair in zip(coords, density):
                vial.write('%s\t%s\n'%tuple(pair))
        vial.close()



SNP_list=[]

with open("%s.SNPS"%SNP_file, 'rb') as csvfile:

	spamreader = csv.reader(csvfile, delimiter='\t')
	for row in spamreader:
		SNP_list.append(int(row[1]))


coord, rho  = archimedes(SNP_list)
write_to_varwig(coord, rho, "%s.density.wig"%SNP_file, color, "%s density"%SNP_file)




# Parent1 = sys.argv[1]
# Parent2 = sys.argv[2]
# ref_dict = sys.argv[3]

# window_size = 10**3	#looking at a kbp window at a time
# critical_thresh = 4	#we want 


# def histogrammatical(input_SNPs, name, color ):

# 	snp_dict = {}
# 	with open(input_SNPs, 'rb') as csvfile:
# 		spamreader = csv.reader(csvfile, delimiter='\t')
# 		for row in spamreader:
# 			if row[0] in snp_dict.keys():
# 				snp_dict[row[0]].append(int(row[1]))
# 			else:
# 				snp_dict[row[0]] = [int(row[1])]

# 	chrom_lens = {}
# 	for record in SeqIO.parse(ref_dict, "fasta"):
# 		if record.id in snp_dict.keys():
# 			chrom_lens[record.id] = len(record)


# 	snp_rho = {}
# 	mux = 0
# 	print name
# 	for chro in snp_dict.keys():
# 		snp_rho[chro] = []
# 		start = 0
# 		snps = np.array(snp_dict[chro])
# 		while start < chrom_lens[chro]:
# 			snp_rho[chro].append(len(snps[snps>start][snps[snps>start]<start+window_size]))
# 			start += 100

# 		num, bins, patches = plt.hist(snp_rho[chro], bins=49, histtype='step', color=color)
# 		mux = max(mux, max(num))
# 		tot = float(sum(num))
# 		adequate = float(sum(num[bins[:-1]+0.5*np.diff(bins)>critical_thresh]))
# #		print adequate
# #		print len(num), num
# #		print len(np.diff(bins)), np.diff(bins)
# #		print chro, adequate
# 		print "%s:	%.2f%%	adequate;	%.2f%% sparse\n"%tuple([chro, 100*adequate/tot, 100-100*adequate/tot])


# 	plt.plot([],[], color, label=name)
# 	print
# 	return mux


# m1 = histogrammatical(Parent1, Parent1.split('/')[-1].split('.SNPS')[0], 'b')
# m2 = histogrammatical(Parent2, Parent2.split('/')[-1].split('.SNPS')[0], 'r')


# plt.vlines(critical_thresh, 0, max(m1,m2)/2, 'k', label='threshold')
# plt.legend()
# plt.title('Diagnostic SNPs per %0.1e Base Pair Window'%tuple([window_size]))
# plt.xlabel('# SNPs per Window')
# plt.ylabel('# Windows')
# plt.show()
# plt.savefig('Parental_SNP_Densities.png')








