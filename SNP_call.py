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
phial = open('%s.SNPS'%SNP_log, 'w')
print os.system('pwd')
print os.system('ls')

for rec in parser:
	if rec.is_snp and float(sum(rec.INFO['DP4'][-2:]))/sum(rec.INFO['DP4']) > 0 and float(rec.INFO['DP']) > 12 :
#																		#if it's a SNP and it's called by more than 75% of the reads, and there are 12+ reads covering.
		#																		Sekelsky data are monoploid, so no heterozygosity worries.

		eff = open('%s.bed'%SNP_log, 'w')
		eff.write('%s\t%s\t%s\n'%tuple([rec.CHROM, int(rec.POS), int(rec.POS)+1]))
		eff.close()
		coverage = float(os.system("bedtools coverage -abam %s -b %s.bed | cut -f 5 "%tuple([parental_bam, SNP_log])))	#	Site coverage
		os.remove('%s.bed'%SNP_log)	#	cleanup
		#					^^^ This pulls the coverage at the SNP site
		eff = open('%s.bed'%SNP_log, 'w')
		eff.write('%s\t%s\t%s\n'%tuple([rec.CHROM, max(0, int(rec.POS)-1000), int(rec.POS)+1000]))
		eff.close()
		reg_cov = float(os.system("bedtools coverage -abam %s -b %s.bed | cut -f 5 | awk '{ total += $1; count++ } END { print total/count }' "%tuple([parental_bam, SNP_log])))	#	Regional coverage
		#					^^^	This pulls out average regional coverage to compare
		os.remove('%s.bed'%SNP_log)	#	cleanup
		#coverage = 0
		#reg_cov = 0
		phial.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%tuple([rec.CHROM, rec.POS, rec.REF, '%s'*len(rec.ALT)%tuple(rec.ALT), coverage, reg_cov, rec.INFO['DP'], rec.QUAL]))
		#			Writes a flatfile of approved SNPs:	Chromome	Position	Reference	Alternate alleles	Coverage at-site 	Covereage +/- 1kb



