#!/bin/sh

ref_genome=$1
genome_index=$2 	#	cdjones_lab/csoeder/rich_introgression/ref_genomes/GMAP_ref/
single=$3
file_out=$4

#	re: building GMAP index:	gmap_build -d dmel-all-chromosome-r5.57 -D GMAP_ref dmel-all-chromosome-r5.57.fasta	
#	re: -d option ...	GSNAP seems to insist the the GMAP index be referenced relative to /proj/seq/data/gsnap/reference_genomes/ 

gsnap -A sam -d ../../../$genome_index $single > "$file_out".sam

samtools view -bS "$file_out".sam | samtools sort - "$file_out".sort								#
samtools index "$file_out".sort.bam	

samtools mpileup -uf $ref_genome "$file_out".sort.bam | bcftools view -bvcg - | bcftools view - | vcfutils.pl varFilter -D100 > "$file_out".vcf





