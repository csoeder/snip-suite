#!/bin/sh
#############################################
# $1	:	Parent1.SNPS
# $2	:	Parent2.SNPS
# $3	:	Hybrid.SNPS
# $4	:	Alignment.sort.bam
#############################################
#Load relevant modules 						#
. /nas02/apps/Modules/default/init/bash		#
module load bedtools
#############################################


title_1=$(echo $1 | cut -f 1 -d . | rev | cut -f 1 -d / | rev)
title_2=$(echo $2 | cut -f 1 -d . | rev | cut -f 1 -d / | rev)
title_hyb=$(echo $3 | cut -f 1 -d . | rev | cut -f 1 -d / | rev)

declare -a somes=("2L","2R","3L","3R","4","X","YHet","2RHet")


for chrom in "${somes[@]}"; 

	do echo $chrom;
	python snip-suite/Portrayal_allSNPs.py $1 $title_1 $2 $title_2 $3 $title_hyb $4 $chrom; 
	
done

echo track name=$title_hyb description=$title_hyb itemRgb="On" > $title_hyb.bed

sed -e 's/^/chr/' *in_$title_hyb* >> $title_hyb.pre
sort -k 1,1 -k2,2n $title_hyb.pre >> $title_hyb.bed
rm $title_hyb.pre
rm *in_$title_hyb*