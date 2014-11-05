
#!/bin/sh
#############################################
#	$1	the reference genome				#
#	$2 										#
#	$3	the FASTQ files of sequences to align
#	$4	job title							#
#############################################
Load relevant modules 						#
. /nas02/apps/Modules/default/init/bash		#
module load bwa								#
module load samtools						#
#	Is the reference genome ready?			#
###	Align the reads! ########################################
bwa aln -O 10 -E 3 -M 2 -o 2 -e 10 -n 5 $1 $2 > $4_1.sai	#	Parameters from PSISeq procedure
bwa aln -O 10 -E 3 -M 2 -o 2 -e 10 -n 5 $1 $3 > $4_2.sai	#
bwa sampe $1 $4_1.sai $4_2.sai $2 $3 > $4.sam 			#
###	Cleanup the alignment	#################################
samtools view -bS $4.sam > $4.bam							#
samtools sort $4.bam $4.sort								#
samtools index $4.sort.bam									#
###	Generate the VCF ########################################
samtools mpileup -uf $1 $4.sort.bam | bcftools view -bvcg - > $4.bcf
bcftools view $4.bcf | vcfutils.pl varFilter -D100 > $4.vcf #	http://samtools.sourceforge.net/mpileup.shtml
###	Generate the .SNPS file #################################
#python SNP_call.py $4.vcf $4.sort.bam $4.SNPS				#
#############################################################