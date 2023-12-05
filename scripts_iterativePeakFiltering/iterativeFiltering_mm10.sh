#!/bin/bash

set -e

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

module use.own /project/wherrylab/wlabtools/privatemodules/
module load bedtoolsv2.29.2
export PATH=$PATH:/project/wherrylab/wlabtools/scripts/iterative_peak_filtering/

peak_half_width=250
chrom_sizes_file='/project/wherrylab/wlabtools/Ref_files/mm10/mm10.chrom.sizes'
blacklist_file='/project/wherrylab/wlabtools/Ref_files/mm10/mm10.blacklist.bed'
genome_file='/project/wherrylab/wlabtools/genomes/mouse.mm10.genome'

##Step2 after macs2 call,use summit bed files
#Extend peak summmits to X bp

for i in `ls *_summits.bed`;
do
	sampleID="${i%.*}"
	echo "${sampleID}"
	calculate_center_peaks_and_extend.sh $i ${sampleID}_extended.bed $chrom_sizes_file $peak_half_width
	
done

echo $?
echo 'Summit bed files are extended'


##Step3 Remove peaks overlapping blacklist regions

for i in `ls *_extended.bed`;
do
	sampleID="${i%.*}"
	echo "${sampleID}"
	remove_peaks_overlapping_blacklisted_regions.sh $i ${sampleID}_debl.bed $chrom_sizes_file $blacklist_file
	
done

echo $?
echo 'Blacklisted regions are removed'

##Step4 Iteratively filter out less sig peaks that overlap with a more significant one

for i in `ls *_extended_debl.bed`;
do
	sampleID="${i%.*}"
	echo "${sampleID}"
	iterative_peak_filtering.sh $i ${sampleID}_filtered.bed $chrom_sizes_file
done

echo $?
echo 'Iterative filtering is done!'

##Step5 Normalize peaks by creating a score per million

for i in `ls *_extended_debl_filtered.bed`;
do
	sampleID="${i%.*}"
	echo "${sampleID}"
	normalize_macs2_peak_scores.sh $i ${sampleID}_norm.bed
	
done

echo $?
echo 'Peak Scores normalization is done!'

# Create 3 column bed file
for i in `ls *_extended_debl_filtered_norm.bed`;
do
	sampleID="${i%.*}"
	echo "${sampleID}"
	awk 'BEGIN {OFS="\t"} $1 ~ /chr./ {print $1, $2, $3, $5}' $i > ${sampleID}_norm2.bed	
done


## Combine peak files and rerun iterative filtering

cat *_extended_debl_filtered_norm.bed > combined_bed_cat.bed

iterative_peak_filtering.sh combined_bed_cat.bed combined_cat_filtered.bed $chrom_sizes_file






