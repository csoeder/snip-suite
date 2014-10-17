"""
Given two lists of SNPs from SNP_call:
		1 	partition them into shared SNPs, SNPs from list 1 only, and SNPs from list 2 only
		2 	

"""


import matplotlib 
matplotlib.use('agg') 
import csv
import sys
import numpy as np
import matplotlib.pyplot as plt


parent1_SNPS_file = sys.argv[1]
title1 = sys.argv[2]
parent2_SNPS_file = sys.argv[3]
title2 = sys.argv[4]
hybrid_SNPS_file = sys.argv[5]
titleHyb = sys.argv[6]



chrom='2L'
box_size=10**2


def pool_snps(parent1, parent2):
	shared = {}
	disjoint1 = {}
	disjoint2 = {}

	with open(parent1, 'rb') as csvfile:
		spamreader = csv.reader(csvfile, delimiter='\t')
		for row in spamreader:
			if row[0] == chrom:
				disjoint1[int(row[1])] = row[3]###		Collect all the SNPs from parent 1



	with open(parent2, 'rb') as csvfile:
		spamreader = csv.reader(csvfile, delimiter='\t')
		for row in spamreader:
			if row[0] == chrom:
				if int(row[1]) in disjoint1.keys() and disjoint1[int(row[1])] == row[3]:#	If this site has a SNP in parent1
#																					#	and the polymorphism is the same as parent1
					shared[int(row[1])] = disjoint1.pop(int(row[1]))#				Then it's a shared SNP and toss it in the shared dict
				else:
					disjoint2[int(row[1])] = row[3]		#Else, toss it in with parent2

	return shared, disjoint1, disjoint2


def snp_grep(parent1, parent2, hybrid):
	parent1_snps = []
	parent2_snps = []


	with open(hybrid, 'rb') as csvfile:
		spamreader = csv.reader(csvfile, delimiter='\t')
		for row in spamreader:
			if row[0] == chrom:
				pos = int(row[1])
				var = row[3]

				if pos in parent1.keys() and parent1[pos] == var:
					parent1_snps.append(pos)
				elif pos in parent2.keys() and parent2[pos] == var:
					parent2_snps.append(pos)
	return parent1_snps, parent2_snps



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





shared_SNPs, disjoint1_SNPs, disjoint2_SNPs = pool_snps(parent1_SNPS_file, parent2_SNPS_file)
introgress1, introgress2 = snp_grep(disjoint1_SNPs, disjoint2_SNPs, hybrid_SNPS_file)

coord, dens = archimedes(disjoint1_SNPs.keys())
write_to_varwig(coord, dens, '%s_disjoint.wig'%title1, '255,0,0', '%s disjoint SNP density'%title1)
coord, dens = archimedes(disjoint2_SNPs.keys())
write_to_varwig(coord, dens, '%s_disjoint.wig'%title2, '0,0,255', '%s disjoint SNP density'%title1)


coord, dens = archimedes(introgress1)
write_to_varwig(coord, dens, '%s_SNPs_in_%s.wig'%tuple([title1, titleHyb]), '255,0,0', '%s-specific SNP density in %s'%tuple([title1, titleHyb]))
coord, dens = archimedes(introgress2)
write_to_varwig(coord, dens, '%s_SNPs_in_%s.wig'%tuple([title2, titleHyb]), '0,0,255', '%s-specific SNP density in %s'%tuple([title2, titleHyb]))



