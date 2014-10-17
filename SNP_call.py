import vcf
import sys


parental_vcf = sys.argv[1]
SNP_log = sys.argv[2]


parser = vcf.Reader(open(parental_vcf))
phial = open(SNP_log, 'w')

for rec in parser:
	if rec.is_snp and float(sum(rec.INFO['DP4'][-2:]))/sum(rec.INFO['DP4']) > 0 :#if it's a SNP and it's called by more than 75% of the reads
		phial.write('%s\t%s\t%s\t%s\n'%tuple([rec.CHROM, rec.POS, rec.REF, '%s'*len(rec.ALT)%tuple(rec.ALT)]))

