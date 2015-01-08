#!/bin/sh

VCF_IN=$1
CLEAN_BED=$2



vcftools --vcf $VCF_IN --remove-indels --non-ref-af 0.25 --max-non-ref-af 0.75 --recode --stdout | vcf2bed > ambig.bed

BED_LEN=$(wc -l $CLEAN_BED | cut -f 1 -d ' ');

BOOKMARK=1;

echo $BED_LEN;
echo $BOOKMARK

while [[ $BOOKMARK -lt $BED_LEN ]]
do
        line1=$(head -n $BOOKMARK ambig.bed | tail -1);
        chr1=$(echo $line1 | tr ' ' '\t'| cut -f 1 );
        line2=$(head -n $((BOOKMARK+1)) ambig.bed | tail -1 | tr ' ' '\t');
        chr2=$(echo $line2 | tr ' ' '\t' | cut -f 1 );

        if [[ $chr1 == $chr2 ]];
                then
                        start=$(echo $line1 | tr ' ' '\t'| cut -f 2 );
                        stop=$(echo $line2 | tr ' ' '\t'| cut -f 2);
                        echo -e "$chr1\t$start\t$stop" >> bridged.bed;
                        ((BOOKMARK = BOOKMARK + 1));
                else
                        ((BOOKMARK = BOOKMARK + 1));
        fi
done

bedtools intersect -v -a bridged.bed -b $CLEAN_BED | bedtools intersect -a ambig.bed -b - > contiguous_middlers.bed
#rm bridged.bed
#rm ambig.bed































