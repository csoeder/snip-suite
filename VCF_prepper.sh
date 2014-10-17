#!/bin/sh

#	$1	the reference genome
#	$2
#	$3	the FASTQ files of sequences to align
#	$4	the VCF output
#	$5	job title

#############################################
Load relevant modules 						#
. /nas02/apps/Modules/default/init/bash		#
module load bwa
module load samtools
#	Is the reference genome ready?

###	Align the reads!
bwa aln -O 10 -E 3 -M 2 -o 2 -e 10 -n 5 $1 $2 > $5_1.sai
bwa aln -O 10 -E 3 -M 2 -o 2 -e 10 -n 5 $1 $3 > $5_2.sai
bwa sampe $1 read1.sai read2.sai $2 $3 > $5.sam
###	Cleanup the alignment
samtools view -bS $5.sam > $5.bam
samtools sort $5.bam $5.sort
samtools index $5.sort.bam
###	Generate the VCF
samtools mpileup -uf $1 $5.sort.bam | bcftools view -bvcg - > $5.bcf
bcftools view $5.bcf | vcfutils.pl varFilter -D100 > $5.vcf
###	Generate the .SNPS file
python SNP_call.py $5.vcf $5.sort.bam $5.SNPS

