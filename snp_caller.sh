#!/bin/sh
echo "Start the comparison!"
parent1=PARC1_GSNAP
parent2=PARG_GSNAP
hybrid=REC3_GSNAP

echo "	Filtering the variant sites...	"
vcftools --vcf "$parent1".vcf --stdout --remove-indels --non-ref-af 0.95 --minDP 10 --recode | vcf2bed > "$parent1".filtered.bed  
vcftools --vcf "$parent2".vcf --stdout --remove-indels --non-ref-af 0.95 --minDP 10 --recode | vcf2bed > "$parent2".filtered.bed  
vcftools --vcf "$hybrid".vcf --stdout --remove-indels --non-ref-af 0.95 --minDP 10 --recode | vcf2bed > "$hybrid".filtered.bed 
echo "	...done!	"

#	make 3-column .clipped.bed files to disregard all but position 
cut -f 1,2,3 "$parent1".filtered.bed > "$parent1".filtered.clipped.bed
cut -f 1,2,3 "$parent2".filtered.bed > "$parent2".filtered.clipped.bed

echo "	Collect all the locations which are variant in one but not the other parent...	"
bedtools intersect -v -wa -a "$parent1".filtered.clipped.bed -b "$parent2".filtered.clipped.bed | cut -f 1,2,3 > unique_in_"$parent1".bed
bedtools intersect -v -wa -a "$parent2".filtered.clipped.bed -b "$parent1".filtered.clipped.bed | cut -f 1,2,3 > unique_in_"$parent2".bed
echo "	...done!"
echo "	Collect all the locations which are variable in both parents..."
bedtools intersect -a "$parent1".filtered.clipped.bed -b "$parent2".filtered.clipped.bed > sites_in_common.bed
echo "	...	take note of the alleles at each of these sites...	"
bedtools intersect -wa -a "$parent1".filtered.bed -b sites_in_common.bed | cut -f 1,2,3,7 > "$parent1"_alleles_at_common_sites.bed
bedtools intersect -wa -a "$parent2".filtered.bed -b sites_in_common.bed | cut -f 1,2,3,7 > "$parent2"_alleles_at_common_sites.bed

#	Here's a recursive function which scans the *_alleles_at_common_sites.bed files, and collects any lines between the files which differ,
#		i.e., any lines which have differing alleles
check_the_common_sites() {
	report=$(cmp "$parent1"_alleles_at_common_sites.bed "$parent2"_alleles_at_common_sites.bed) 	#	Where do the files first differ?
	line_number=$(echo $report | cut -f 2 -d : | cut -f 3 -d e) 								#	Find the line number
	sed -n "$line_number"p "$parent1"_alleles_at_common_sites.bed >> diagnostic_sites_"$parent1".bed #	Push that line into a file just for parent1
	sed -i "$line_number"d "$parent1"_alleles_at_common_sites.bed
	sed -n "$line_number"p "$parent2"_alleles_at_common_sites.bed >> diagnostic_sites_"$parent2".bed #	And from the other file into just parent2
	sed -i "$line_number"d "$parent2"_alleles_at_common_sites.bed
	cmp --silent "$parent1"_alleles_at_common_sites.bed "$parent2"_alleles_at_common_sites.bed || check_the_common_sites #	Any more lines to go? 
}

echo "	... and collect all the shared sites with different alleles! (This may take a moment)."
cmp --silent "$parent1"_alleles_at_common_sites.bed "$parent2"_alleles_at_common_sites.bed || check_the_common_sites
echo "	Done!"

echo "	Collect the locations from $hybrid which are unique to the parents..."
bedtools intersect -wa -a "$hybrid".filtered.bed -b unique_in_"$parent1".bed | cut -f 1,2,3 > "$parent1"_unique_in_"$hybrid".bed
len=$(wc -l "$parent1"_unique_in_"$hybrid".bed | cut -f 1 -d ' ')
echo -e "$parent1\t666\t+" | perl -ne 'print "$_" x'$len'' > filler.bed.temp
echo "255,0,0" | perl -ne 'print "$_" x'$len'' > colour.bed.temp
cut -f 1,2,3 "$parent1"_unique_in_"$hybrid".bed| paste - filler.bed.temp > partial.bed.temp
cut -f 2,3 "$parent1"_unique_in_"$hybrid".bed | paste partial.bed.temp - colour.bed.temp > "$parent1"_in_$hybrid.bed
rm *.temp

bedtools intersect -wa -a "$hybrid".filtered.bed -b unique_in_"$parent2".bed | cut -f 1,2,3 > "$parent2"_unique_in_"$hybrid".bed
len=$(wc -l "$parent2"_unique_in_"$hybrid".bed | cut -f 1 -d ' ')
echo -e "$parent2\t666\t+" | perl -ne 'print "$_" x'$len'' > filler.bed.temp
echo "0,0,255" | perl -ne 'print "$_" x'$len'' > colour.bed.temp
cut -f 1,2,3 "$parent2"_unique_in_"$hybrid".bed | paste - filler.bed.temp > partial.bed.temp
cut -f 2,3 "$parent2"_unique_in_"$hybrid".bed | paste partial.bed.temp - colour.bed.temp > "$parent1"_in_$hybrid.bed
rm *.temp
echo "	...done!"

echo "	Mask the diagnostic sites, keeping only sites which are variants in the hybrid..."
bedtools intersect -wa -a "$hybrid".filtered.bed -b diagnostic_sites_"$parent1".bed| cut -f 1,2,3,7 > "$hybrid"_alleles_at_common_sites.bed
bedtools intersect -wa -a diagnostic_sites_"$parent1".bed -b "$hybrid"_alleles_at_common_sites.bed > "$parent1"_diagnostics_masked.bed
bedtools intersect -wa -a diagnostic_sites_"$parent2".bed -b "$hybrid"_alleles_at_common_sites.bed > "$parent2"_diagnostics_masked.bed
echo "	... done!"


check_the_diagnostic_sites() {
	report=$(cmp "$hybrid"_alleles_at_common_sites.bed "$parental"_diagnostics_masked.bed)	#	Any variant sites whose alleles disagree with the parent?
	line_number=$(echo $report | cut -f 2 -d : | cut -f 3 -d e) 	#	Find the disagreements
	sed -i "$line_number"d "$hybrid"_alleles_at_common_sites.bed 	#	and remove them
	sed -i "$line_number"d "$parental"_diagnostics_masked.bed
	cmp --silent "$hybrid"_alleles_at_common_sites.bed "$parental"_diagnostics_masked.bed || check_the_diagnostic_sites 	#	If any more disagreements, repeat!
}

parental=$parent1
echo "	Finding common sites in hybrid belonging to $parental ...	"
cmp --silent "$hybrid"_alleles_at_common_sites.bed "$parental"_diagnostics_masked.bed || check_the_diagnostic_sites
mv "$hybrid"_alleles_at_common_sites.bed "$parental"_diagnostic_sites_in_"hybrid".bed
len=$(wc -l "$parental"_diagnostic_sites_in_"hybrid".bed | cut -f 1 -d ' ')
#	Package everything, color-coded
echo -e "$parental\t666\t+" | perl -ne 'print "$_" x'$len'' > filler.bed.temp
echo "255,0,0" | perl -ne 'print "$_" x'$len'' > colour.bed.temp
cut -f 1,2,3 "$parental"_diagnostic_sites_in_"hybrid".bed | paste - filler.bed.temp > partial.bed.temp
cut -f 2,3 "$parental"_diagnostic_sites_in_"hybrid".bed | paste partial.bed.temp - colour.bed.temp >> "$parental"_in_"$hybrid".bed
rm *.temp
echo "	...done!	"

#	Rebuild the base file...	
bedtools intersect -wa -a "$hybrid".filtered.bed -b diagnostic_sites_"$parent1".bed| cut -f 1,2,3,7 > "$hybrid"_alleles_at_common_sites.bed

parental=$parent2
echo "	Finding common sites in hybrid belonging to $parental ...	"
cmp --silent "$hybrid"_alleles_at_common_sites.bed "$parental"_diagnostics_masked.bed || check_the_diagnostic_sites
mv "$hybrid"_alleles_at_common_sites.bed "$parental"_diagnostic_sites_in_"hybrid".bed
len=$(wc -l "$parental"_diagnostic_sites_in_"hybrid".bed | cut -f 1 -d ' ')
#	Package everything, color-coded
echo -e "$parental\t666\t+" | perl -ne 'print "$_" x'$len'' > filler.bed.temp
echo "255,0,0" | perl -ne 'print "$_" x'$len'' > colour.bed.temp
cut -f 1,2,3 "$parental"_diagnostic_sites_in_"hybrid".bed | paste - filler.bed.temp > partial.bed.temp
cut -f 2,3 "$parental"_diagnostic_sites_in_"hybrid".bed | paste partial.bed.temp - colour.bed.temp >> "$parental"_in_"$hybrid".bed
rm *.temp
echo "	...done!	"

echo 'track name=""$hybrid" ancestry" description=""$hybrid" ancestry" itemRgb="On"' > "$hybrid"_ancestry.bed
cat "$parent1"_in_"$hybrid".bed "$parent2"_in_"$hybrid".bed >> "$hybrid"_ancestry.bed


