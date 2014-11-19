"""
Given two lists of SNPs from SNP_call:
		1 	partition them into shared SNPs, SNPs from list 1 only, and SNPs from list 2 only
		2 	

"""

#		THIS VERSION OPTIMIZED FOR THE CROWN-SEKELSKY FLIES


import matplotlib 
matplotlib.use('agg') 
import csv
import sys
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call, check_output
import pickle


parent1_SNPS_file = sys.argv[1]
title1 = sys.argv[2]
parent2_SNPS_file = sys.argv[3]
title2 = sys.argv[4]
hybrid_SNPS_file = sys.argv[5]
titleHyb = sys.argv[6]
hybrid_Align=sys.argv[7]	#sorted BAM file of hybrid alignment

###yesyes


#chroms=["2L","2R","3L","3R","4","X", "YHet", "2RHet"]
chroms=["2L"]
box_size=10**3
missing_SNP_threshold = 10	#hybrid must have at least this coverage to declare that it is missing a parental SNP
#Super awesome

def pool_snps(parent1, parent2):
	""" given two .SNPS files, this loads them into python dict objects, then outputs
	which SNPs are disjoint and which are shared """
	p1=dict.fromkeys(chroms,{})
	p2=dict.fromkeys(chroms,{})

	with open(parent1, 'rb') as csvfile:
		spamreader = csv.reader(csvfile, delimiter='\t')
		for row in spamreader:
			try:
				p1[row[0]][int(row[1])] = row[3]###		Collect all the SNPs from parent 1;
					###										load them into the parent1 dict in the form {CHR:{POSITION:ALT_ALLELE}}
			except KeyError:	#	if the chromosome isn't in chroms
				pass

	with open(parent2, 'rb') as csvfile:
		spamreader = csv.reader(csvfile, delimiter='\t')
		for row in spamreader:
			try:
				p2[row[0]][int(row[1])] = row[3]###		Collect all the SNPs from parent 1;
					###										load them into the parent1 dict in the form {CHR:{POSITION:ALT_ALLELE}}
			except KeyError:	#	if the chromosome isn't in chroms
				pass

	shared = dict.fromkeys(chroms,{})		#shared SNPs go here
	disjoint1 = dict.fromkeys(chroms,{})	#those only in parent1
	disjoint2 = dict.fromkeys(chroms,{})	#those only in parent2

	for chro in chroms:						#For each chromosome...
		disj1 = list(set(p1[chro].keys()).difference(set(p2[chro].keys())))	#sites unique to p1
		disj2 = list(set(p2[chro].keys()).difference(set(p1[chro].keys())))	#sites unique to p2
		inter = list(set(p1[chro].keys()).intersection(set(p2[chro].keys())))#shared sites

		for i in inter:						#for each shared site...
			if p1[chro][i] != p2[chro][i]:	#if the phenotype at the site is actually different...
				inter.pop(inter.index(i))			#remove it from the shared SNPS
				disj1.append(i)			#insert each site 
				disj2.append(i)			#				in the appropriate list
		for j in disj1:						#
			disjoint1[chro][j]=p1[chro][j]	#Now take each site and load it into the master dicts for output
		for j in disj2:
			disjoint2[chro][j]=p2[chro][j]
		for j in inter:
			shared[chro][j]=p1[chro][j]

	return shared, disjoint1, disjoint2


def cov_grep(snp_list, bam_file):
	phial=open('%s_sites.bed'%titleHyb, 'w')
	for chro in chroms:
		for site in snp_list[chro]:
			phial.write('%s\t%s\t%s\n'%tuple([chro, site, site+1]))
	phial.close()
	call('bedtools coverage -abam %s -b %s_sites.bed > %s.cov'%tuple([hybrid_Align, titleHyb, titleHyb]), shell=True)
	coverage=dict.fromkeys(chroms,{})
	with open('%s.cov'%titleHyb, 'rb') as csvfile:
		spamreader=csv.reader(csvfile, delimiter='\t')
		for row in spamreader:
			coverage[row[0]][int(row[1])]=int(row[3])
	return coverage




def snp_grep(parent1, parent2, hybrid, hyb_cov):
	"""
	Given the two dicts of disjoint parental SNPs, load the hybrid .SNPS file, and 
	look at each SNP site in the two dicts. Classify each as present or absent 
	in the hybrid genome
	"""

	hybrid_snps=dict.fromkeys(chroms,{})
	parent_snps=dict.fromkeys(chroms, [])
	for chro in chroms:
		parent_snps[chro].extend(parent1[chro])
		parent_snps[chro].extend(parent2[chro])

	with open(hybrid, 'rb') as csvfile:
		spamreader = csv.reader(csvfile, delimiter='\t')
		for row in spamreader:
			try:
				hybrid_snps[row[0]][int(row[1])] = row[3]###		Collect all the SNPs from hybrid
				###										load them into the parent1 dict in the form POSITION:ALT_ALLELE	
			except KeyError:
				pass
	print "begin coverage grep"
	coverage = cov_grep(parent_snps, hyb_cov)
	print "done grepping "



	parent1_present = dict.fromkeys(chroms,[])	#	SNPs from p1 which are present
	parent2_present = dict.fromkeys(chroms,[])	#	SNPs from p2 which are present
	#parents 1 and 2 absent pending coverage grepping
	parent1_absent = dict.fromkeys(chroms,[])	#	SNPs from p1 which are present
	parent2_absent = dict.fromkeys(chroms,[])	#	SNPs from p2 which are present

	gnu_vars = dict.fromkeys(chroms,[])

	for chro in chroms:
		for site in hybrid_snps[chro].keys():
			if site in set(parent1[chro].keys()).intersection(set(parent2[chro].keys())):
				if hybrid_snps[chro][site] in parent1[chro][site]:
					parent1_present[chro].append(int(site))
				elif hybrid_snps[chro][site] in parent2[chro][site]:
					parent2_present[chro].append(int(site))
				else:#Record this as an additional SNP
					gnu_vars[chro].append(int(site))
			elif site in parent1[chro].keys():
				if hybrid_snps[chro][site] in parent1[chro][site]:
					parent1_present[chro].append(int(site))
			elif site in parent2[chro].keys():
				if hybrid_snps[chro][site] in parent2[chro][site]:
					parent2_present[chro].append(int(site))
			else:
				gnu_vars[chro].append(int(site))


		for site in parent1[chro].keys():
			if site not in parent1_present[chro] and site not in parent2_present[chro] and site not in gnu_vars[chro]:#make sure the site isn't recorded ANYwhere....
				if coverage[chro][site] >= missing_SNP_threshold:
					parent1_absent[chro].append(int(site))
		for site in parent2[chro].keys():
			if site not in parent1_present[chro] and site not in parent2_present[chro] and site not in gnu_vars[chro]:
				if coverage[chro][site] >= missing_SNP_threshold:
					parent2_absent[chro].append(int(site))

	return parent1_present, parent2_present, gnu_vars, parent1_absent, parent2_absent#, new_snps, hypervars




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
	for chro in points.keys():
		for site in points[chro]:
			vial.write('%s\t%s\t%s\therpderp\t0\t+\t%s\t%s\t%s\n'%tuple([chro, site, site+1, site, site+1, colour]))
	vial.close()



shared_SNPs, disjoint1_SNPs, disjoint2_SNPs = pool_snps(parent1_SNPS_file, parent2_SNPS_file)
present1, present2, n00bs, absent1, absent2= snp_grep(disjoint1_SNPs, disjoint2_SNPs, hybrid_SNPS_file, hybrid_Align)

pickle.dump( [shared_SNPs, disjoint1_SNPs, disjoint2_SNPs, present1, present2, n00bs, absent1, absent2], open('%s.pickle'%titleHyb, 'wb') )

###http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
print "\t\tREPORT:\t\t\t"
print "Between %s and %s, %s SNP variants were logged."%tuple([title1, title2, len( [item for sublist in shared_SNPs.values() for item in sublist] )+len( [item for sublist in disjoint1_SNPs.values() for item in sublist] )+len([item for sublist in disjoint2_SNPs.values() for item in sublist])])
print "%s SNPs were identified unique to %s"%tuple([len([item for sublist in disjoint1_SNPs.values() for item in sublist]), title1])
for chro in chroms:
	print "\t%s on %s\n"%tuple([len(disjoint1_SNPs[chro]), chro])
print "%s SNPs were identified unique to %s"%tuple([len([item for sublist in disjoint2_SNPs.values() for item in sublist]), title2])
for chro in chroms:
	print "\t%s on %s\n"%tuple([len(disjoint2_SNPs[chro]), chro])
print "%s SNPs were found shared between %s and %s"%tuple([len([item for sublist in shared_SNPs.values() for item in sublist]), title1, title2])
for chro in chroms:
	print "\t%s on %s\n"%tuple([len(shared_SNPs[chro]), chro])
print "%s contained %s SNPs unique to %s"%tuple([titleHyb, len([item for sublist in present1.values() for item in sublist]), title1])
for chro in chroms:
	print "\t%s on %s\n"%tuple([len(present1[chro]), chro])
print "%s was missing %s SNPs unique to %s"%tuple([titleHyb, len([item for sublist in absent1.values() for item in sublist]), title1])
for chro in chroms:
	print "\t%s on %s\n"%tuple([len(absent1[chro]), chro])
print "%s contained %s SNP unique to %s"%tuple([titleHyb, len([item for sublist in present2.values() for item in sublist]), title2])
for chro in chroms:
	print "\t%s on %s\n"%tuple([len(present2[chro]), chro])
print "%s was missing %s SNPs unique to %s"%tuple([titleHyb, len([item for sublist in absent2.values() for item in sublist]), title2])
for chro in chroms:
	print "\t%s on %s\n"%tuple([len(absent2[chro]), chro])
print "%s contained %s SNPs unseen in either %s or %s\n"%tuple([titleHyb, len([item for sublist in n00bs.values() for item in sublist]), title1, title2])
#print "\t\t including %s sites representing a third polymorphism"%tuple([len(hypervars)])
print "\t\t~~~END REPORT~~~\t\t"



write_to_bed(present1, '%s_SNPs_present_in_%s.bed'%tuple([title1, titleHyb]), '255,0,0', '%s SNPs in %s'%tuple([title1, titleHyb]))
write_to_bed(present2, '%s_SNPs_present_in_%s.bed'%tuple([title2, titleHyb]), '0,0,255', '%s SNPs in %s'%tuple([title2, titleHyb]))
write_to_bed(absent1, '%s_SNPs_absent_in_%s.bed'%tuple([title1, titleHyb]), '0,255,0', '%s SNPs missing from %s'%tuple([title1, titleHyb]))
write_to_bed(absent2, '%s_SNPs_absent_in_%s.bed'%tuple([title2, titleHyb]), '255,153,51', '%s SNPs missing from %s'%tuple([title2, titleHyb]))

#coord, dens = archimedes(disjoint1_SNPs.keys())
#write_to_varwig(coord, dens, '%s_disjoint.wig'%title1, '255,0,0', '%s disjoint SNP density'%title1)
#coord, dens = archimedes(disjoint2_SNPs.keys())
#write_to_varwig(coord, dens, '%s_disjoint.wig'%title2, '0,0,255', '%s disjoint SNP density'%title1)

#coord, dens = archimedes(present1)
#write_to_varwig(coord, dens, '%s_SNPs_in_%s.wig'%tuple([title1, titleHyb]), '255,0,0', '%s-specific SNP density in %s'%tuple([title1, titleHyb]))
#coord, dens = archimedes(present2)
#write_to_varwig(coord, dens, '%s_SNPs_in_%s.wig'%tuple([title2, titleHyb]), '0,0,255', '%s-specific SNP density in %s'%tuple([title2, titleHyb]))

#



