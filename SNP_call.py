"""
Given a VCF file and the BAM file it was based upon, a flat file of SNPS and coverages is produced. 
"""

import vcf
import sys
import os


parental_vcf = sys.argv[1]		#input variability data
parental_bam = sys.argv[2]		#for coverage data
SNP_log = sys.argv[3]			#output file for Good SNPs


parser = vcf.Reader(open(parental_vcf))
phial = open(SNP_log, 'w')

for rec in parser:
	if rec.is_snp and float(sum(rec.INFO['DP4'][-2:]))/sum(rec.INFO['DP4']) > 0 :#if it's a SNP and it's called by more than 75% of the reads.
		#																		Sekelsky data are monoploid, so no heterozygosity worries.


		eff = open('temp.bed', 'w')
		eff.write('%s\t%s\t%s\n'%tuple([rec.CHROM, int(rec.POS), int(rec.POS)+1]))
		eff.close()
		coverage = float(os.system("bedtools coverage -abam %s -b temp.bed | cut -f 5 "%tuple([parental_bam])))	#	Site coverage
		os.remove('temp.bed')	#	cleanup
		#					^^^ This pulls the coverage at the SNP site
		eff = open('temp.bed', 'w')
		eff.write('%s\t%s\t%s\n'%tuple([rec.CHROM, max(0, int(rec.POS)-1000), int(rec.POS)+1000]))
		eff.close()
		reg_cov = float(os.system("bedtools coverage -abam %s -b temp.bed | cut -f 5 | awk '{ total += $1; count++ } END { print total/count }' "%tuple([parental_bam])))	#	Regional coverage
		#					^^^	This pulls out average regional coverage to compare
		os.remove('temp.bed')	#	cleanup
		phial.write('%s\t%s\t%s\t%s\t%s\t%s\n'%tuple([rec.CHROM, rec.POS, rec.REF, '%s'*len(rec.ALT)%tuple(rec.ALT), coverage, reg_cov]))
		#			Writes a flatfile of approved SNPs:	Chromome	Position	Reference	Alternate alleles	Coverage at-site 	Covereage +/- 1kb



