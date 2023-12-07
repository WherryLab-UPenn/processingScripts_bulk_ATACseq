#!/bin/bash

if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi

module load python/3.9.1

source $HOME/my_python-3.9.1/bin/activate
export PYTHONPATH=$HOME/my_python-3.9.1/lib/python3.9/site-packages


######################################################################

### Make INDIVIDUAL coverage files


# Generate 3 column bed files from macs3 output called peaks file
for i in `ls *_peaks.narrowPeak`;
do
 sampleID="${i%.*}"
 awk '{print $1"\t"$2"\t"$3}' ${sampleID}.narrowPeak > ${sampleID}.bed
done
echo $?


# Calculate coverage of each bam compared to the sample's called peaks
for i in `ls *_shift.bam`;
do
 sampleID="${i%.*}"
 bedtools coverage -counts -b $i  -a ${sampleID}_peaks.bed > ${sampleID}_indPeaks.bedgraph 
done
echo $?


# Make file with the total number of reads in the individually called peaks --> use to calculate FRiPs

for i in `ls *_indPeaks.bedgraph`;
do
 sampleID="${i%_*_*_*_*_*_*.*}"
 awk -v sampleID="$sampleID" 'BEGIN {OFS="\t"} {s+=$4}END{print sampleID,s}' $i >> combTotalIndCounts.txt
done


# Get number of peaks called for each sample
for i in *_peaks.bed;
do
	sampleID="${i%_*_*_*_*_*_*.*}"
	echo "${sampleID}" >> ${sampleID}_indPeakNumber.txt
	wc -l $i | awk '{print $1}' >> ${sampleID}_indPeakNumber.txt
done

find . -name '*_indPeakNumber.txt' | sort | xargs paste | column -s $'\t' -t > combIndPeakNumber.txt

######################################################################

# This section will need to be repeated for every group comparison -- it is basically creating a consensus list of peak regions (like genes or transcripts for RNA-seq).

# Make 3 column bed file for coverage analysis
awk 'BEGIN {OFS="\t"} $1 ~ /chr./ {print $1, $2, $3}' combined_cat_filtered.bed > combined_cat_filtered_v2.bed


## Make counts bedgraph files with the 4th column being # of tags in each union peak:

for i in `ls *_shift.bam`;
do
 sampleID="${i%.*}"
 bedtools coverage -counts -b $i  -a combined_cat_filtered_v2.bed > ${sampleID}_counts.bedgraph 
done
echo $?
echo 'The union coverage files are ready!'


## Make file with the total number of reads in the union peak list
for i in `ls *_counts.bedgraph`;
do
 echo $i
 sampleID="${i%.*}"
 awk '{s+=$4}END{print FILENAME,s}' ${sampleID}.bedgraph >> combTotalCountUnionPeaks.txt
done

echo $?
echo 'File with total tags in union_merged list is ready!'

# Combine all count bedgraph files into one file with chr, start, end AND headers
for i in *_counts.bedgraph;
do
	sampleID="${i%_*_*_*_*_*_*.*}"
 	echo "${sampleID}"
	echo "${sampleID}" > ${i}_unionTags.txt
	awk '{print $4}' $i >> ${i}_unionTags.txt
done

find . -name '*_unionTags.txt' | sort | xargs paste |  paste -d ' ' > combUnionReads.txt


echo $'chr\tstart\tend' | cat - combined_cat_filtered_v2.bed > combined_cat_filtered_v3.bed

paste -d ' ' combined_cat_filtered_v3.bed combUnionReads.txt > combUnionReadsWithLabels.txt

# Get number of peaks in union&merged list
wc -l combined_cat_filtered_v2.bed | cut -d' ' -f1 >> combUnionPeakNumber.txt

rm *_unionTags.txt

