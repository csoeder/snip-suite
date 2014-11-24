#!/bin/sh
#############################################
#	$1	the reference genome				#
#	$2	the FASTQ files of sequences to align
#	$3	job title							#
#		THIS VERSION OPTIMIZED FOR RICH M.'S SUGARFLIES
#############################################
#Load relevant modules 						#
. /nas02/apps/Modules/default/init/bash		#
module load bwa								#
module load samtools						#
#	Is the reference genome ready?			#
###	Align the reads! ########################################
bwa aln -O 10 -E 3 -M 2 -o 2 -e 10 -n 5 $1 $2 > $3.sai	#	Parameters from PSISeq procedure
bwa samse $1 $3.sai $2> $3.sam 			#
###	Cleanup the alignment	#################################
samtools view -bS $3.sam > $3.bam							#
samtools sort $3.bam $3.sort								#
samtools index $3.sort.bam									#
###	Generate the VCF ########################################
samtools mpileup -uf $1 $3.sort.bam | bcftools view -bvcg - > $3.bcf
bcftools view $3.bcf | vcfutils.pl varFilter -D100 > $3.vcf #	http://samtools.sourceforge.net/mpileup.shtml
###	Generate the .SNPS file #################################
#python SNP_call.py $4.vcf $4.sort.bam $4.SNPS				#
#############################################################