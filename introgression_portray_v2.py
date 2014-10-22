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
from subprocess import call, check_output


parent1_SNPS_file = sys.argv[1]
title1 = sys.argv[2]
parent2_SNPS_file = sys.argv[3]
title2 = sys.argv[4]
hybrid_SNPS_file = sys.argv[5]
titleHyb = sys.argv[6]



chrom='2L'
box_size=10**3
missing_SNP_threshold = 10	#hybrid must have at least this coverage to declare that it is missing a parental SNP
#Super awesome

def pool_snps(parent1, parent2):
	""" given two .SNPS files, this loads them into python dict objects, then outputs
	which SNPs are disjoint and which are shared """
	shared = {}		#shared SNPs go here
	disjoint1 = {}	#those only in parent1
	disjoint2 = {}	#those only in parent2
	with open(parent1, 'rb') as csvfile:
		spamreader = csv.reader(csvfile, delimiter='\t')
		for row in spamreader:
			if row[0] == chrom:
				disjoint1[int(row[1])] = row[3]###		Collect all the SNPs from parent 1;
				###										load them into the parent1 dict in the form POSITION:ALT_ALLELE
	with open(parent2, 'rb') as csvfile:
		spamreader = csv.reader(csvfile, delimiter='\t')
		for row in spamreader:
			if row[0] == chrom:
				if int(row[1]) in disjoint1.keys() and disjoint1[int(row[1])] == row[3]:#	If this site has a SNP in parent1
#																					#	and the polymorphism is the same as parent1
					shared[int(row[1])] = disjoint1.pop(int(row[1]))#				Then it's a shared SNP. Remove and toss it in the shared dict
				else:
					disjoint2[int(row[1])] = row[3]		#Else, toss it in with parent2
	return shared, disjoint1, disjoint2


def snp_grep(parent1, parent2, hybrid):
	"""
	Given the two dicts of disjoint parental SNPs, load the hybrid .SNPS file, and 
	look at each SNP site in the two dicts. Classify each as present or absent 
	in the hybrid genome
	"""
	parent1_present = []#	SNPs from p1 which are present
	parent2_present = []#	SNPs from p2 which are present
	parent1_absent = []#	SNPs from p1 which are absent
	parent2_absent = []#	SNPs from p2 which are absent
	new_snps = {}#			SNPs which are unique to hybric
	hypervars = {}#			SNPs which are a third polymorph of a site variable between p1 and p2
	hybrid_snps={}#			SNPs which are present in the hybrid

	def grep_coverage(chro, point):
		f = open('test.bed', 'w')
		f.write('%s\t%s\t%s\n'%tuple([chro, point, point+1]))
		f.close()
		site_cov = int(check_output(['bedtools', 'intersect', '-a', '%s/%scoverage.bedgraph'%tuple([titleHyb, titleHyb]), '-b', 'test.bed' ]).split('\t')[-1])
		call(['rm', 'test.bed'])		
		return site_cov

	with open(hybrid, 'rb') as csvfile:
		spamreader = csv.reader(csvfile, delimiter='\t')
		for row in spamreader:
			if row[0] == chrom:
				hybrid_snps[int(row[1])] = row[3]###		Collect all the SNPs from hybrid
				###										load them into the parent1 dict in the form POSITION:ALT_ALLELE	
	for key in parent1.keys():	#		For each SNP site in P1...
		if key in hybrid_snps.keys():		#	if the site is variable in Hybrid...
			if parent1[key] == 	hybrid_snps[key]:#		and if the variation is the same...
				parent1_present.append(key)			#		then this snp is present int he hybrid
			elif (key in parent2.keys() and parent2[key] == hybrid_snps[key]):
				#	If site is also variable in P2, and the variation matches, then the P2 SNP is here, and this site is not 'absent' 
				parent2_present.append(key)
			else:	#	If the site doesn't have a recognized SNP, that SNP has gone missing. Where is it???
				if grep_coverage(chrom, key, key+1) > missing_SNP_threshold:#	If we have sufficient coverage...
					parent1_absent.append(key)#		Declare this SNP absent
		else:	#	if the key is not listed as a hybrid SNP, the site is missing! Where is it?
			if grep_coverage(chrom, key, key+1) > missing_SNP_threshold:#	If we have sufficient coverage...
				parent1_absent.append(key)#		Declare this SNP absent
	for key in parent2.keys():	#		For each SNP site in P1...
		if key in hybrid_snps.keys():		#	if the site is variable in Hybrid...
			if parent2[key] == 	hybrid_snps[key]:#		and if the variation is the same...
				parent2_present.append(key)			#		then this snp is present int he hybrid
			elif (key in parent1.keys() and parent1[key] == hybrid_snps[key]):
				#	If site is also variable in P2, and the variation matches, then the P2 SNP is here, and this site is not 'absent' 
				parent1_present.append(key)
			else:	#	If the site doesn't have a recognized SNP, that SNP has gone missing. Where is it???
				if grep_coverage(chrom, key, key+1) > missing_SNP_threshold:#	If we have sufficient coverage...
					parent2_absent.append(key)#		Declare it absent
		else:	#	if the key is not listed as a hybrid SNP, the site is definitely missing!
			if grep_coverage(chrom, key, key+1) > missing_SNP_threshold:#	If we have sufficient coverage...
				parent2_absent.append(key)#	Declare it absent

	for key in hybrid_snps.keys():	#	collect SNPs unique it the hybrid
		if not (key in parent1.keys() or key in parent2.keys()):	#If the SNP is at a unique site...
			new_snps[key] = hybrid_snps[key]
		if (key in parent1.keys() and hybrid_snps[key] not in parent1[key]):
			new_snps[key] = hybrid_snps[key] 
			hypervars[key] = hybrid_snps[key] 
		if (key in parent2.keys() and hybrid_snps[key] not in parent2[key]):
			new_snps[key] = hybrid_snps[key] 
			hypervars[key] = hybrid_snps[key] 

	parent1_present = list(set(parent1_present))
	parent2_present = list(set(parent2_present))
	parent1_absent = list(set(parent1_absent))
	parent2_absent = list(set(parent2_absent))

	return parent1_present, parent2_present, parent1_absent, parent2_absent, new_snps, hypervars


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

def write_to_bed(points, phial, colour, name):
	vial = open(phial, 'w')
	vial.write('browser hide all')
	vial.write('track name="%s" description="%s" visibility=1 itemRgb="On"\n'%tuple([name,name]))
	for point in points:
		vial.write('%s\t%s\t%s\therpderp\t0\t+\t%s\t%s\t%s\n'%tuple([chrom, point, point+1, point, point+1, colour]))
	vial.close()



shared_SNPs, disjoint1_SNPs, disjoint2_SNPs = pool_snps(parent1_SNPS_file, parent2_SNPS_file)
present1, present2, absent1, absent2, n00bs, hypervars = snp_grep(disjoint1_SNPs, disjoint2_SNPs, hybrid_SNPS_file)

print "\t\tREPORT:\t%s\t\t"%chrom
print "Between %s and %s, %s SNP variants were logged."%tuple([title1, title2, len(shared_SNPs)+len(disjoint1_SNPs)+len(disjoint2_SNPs)])
print "%s SNPs were identified unique to %s"%tuple([len(disjoint1_SNPs), title1])
print "%s SNPs were identified unique to %s"%tuple([len(disjoint2_SNPs), title2])
print "%s SNPs were found shared between %s and %s"%tuple([len(shared_SNPs), title1, title2])
print "%s contained %s SNPs unique to %s"%tuple([titleHyb, len(present1), title1])
print "%s was missing %s SNPs unique to %s"%tuple([titleHyb, len(absent1), title1])
print "%s contained %s SNP unique to %s"%tuple([titleHyb, len(present2), title2])
print "%s was missing %s SNPs unique to %s"%tuple([titleHyb, len(absent2), title2])
print "%s contained %s SNPs unseen in either %s or %s"%tuple([titleHyb, len(n00bs), title1, title2])
print "\t\t including %s sites representing a third polymorphism"%tuple([len(hypervars)])
print "\t\t~~~END REPORT~~~\t\t"



write_to_bed(present1, '%s_SNPs_present_in_%s.bed'%tuple([title1, titleHyb]), '255,0,0', '%s SNPs in %s'%tuple([title1, titleHyb]))
write_to_bed(present2, '%s_SNPs_present_in_%s.bed'%tuple([title2, titleHyb]), '0,0,255', '%s SNPs in %s'%tuple([title2, titleHyb]))
write_to_bed(absent1, '%s_SNPs_absent_in_%s.bed'%tuple([title1, titleHyb]), '255,0,255', '%s SNPs missing from %s'%tuple([title1, titleHyb]))
write_to_bed(absent2, '%s_SNPs_absent_in_%s.bed'%tuple([title2, titleHyb]), '0,255,0', '%s SNPs missing from %s'%tuple([title2, titleHyb]))

coord, dens = archimedes(disjoint1_SNPs.keys())
write_to_varwig(coord, dens, '%s_disjoint.wig'%title1, '255,0,0', '%s disjoint SNP density'%title1)
coord, dens = archimedes(disjoint2_SNPs.keys())
write_to_varwig(coord, dens, '%s_disjoint.wig'%title2, '0,0,255', '%s disjoint SNP density'%title1)

coord, dens = archimedes(present1)
write_to_varwig(coord, dens, '%s_SNPs_in_%s.wig'%tuple([title1, titleHyb]), '255,0,0', '%s-specific SNP density in %s'%tuple([title1, titleHyb]))
coord, dens = archimedes(present2)
write_to_varwig(coord, dens, '%s_SNPs_in_%s.wig'%tuple([title2, titleHyb]), '0,0,255', '%s-specific SNP density in %s'%tuple([title2, titleHyb]))





