#!/bin/bash



######################################################################

### Organize output
# Make folders and sort files

mkdir output

mkdir int_files_toRm

mkdir int_files_toRm/sams
mkdir int_files_toRm/mapBams
mkdir int_files_toRm/dedupBams
mkdir int_files_toRm/deblBams
mkdir int_files_toRm/deMitoBams


mkdir output/fastqs
mkdir output/sortBams
mkdir output/shiftedBams

mkdir output/bws
mkdir output/histos
mkdir output/combFiles
mkdir output/indStats
mkdir output/insertMetrics
mkdir output/dedupMetrics

mkdir output/combFiles

mkdir output/peakFiles
mkdir output/peakFiles/iterativeFiltering_Ints
mkdir output/peakFiles/macs3_output

mkdir output/bedgraphs


# Move intermediate files that will be deleted eventually
mv *sam int_files_toRm/sams
mv *out int_files_toRm/sams
mv *map.bam *map.bam.bai int_files_toRm/mapBams
mv *dedup.bam *dedup.bam.bai int_files_toRm/dedupBams
mv *map_deMito.bam *map_deMito.bam.bai int_files_toRm/deMitoBams
mv *debl.bam *debl.bam.bai int_files_toRm/deblBams
rm *_tmp.bam


# Move files to keep
mv *fastq.gz output/fastqs
mv *sort.bam *sort.bam.bai output/sortBams
mv *shift.bam *shift.bam.bai output/shiftedBams

mv *.bw output/bws
mv *_log output/bws

mv *peaks.xls *lambda.bdg *peaks.narrowPeak *summits.bed *control_lambda.bdg *treat_pileup.bdg output/peakFiles/macs3_output
mv *extended.bed *extended_debl.bed *filtered.bed *_norm.bed *_norm2.bed  output/peakFiles/iterativeFiltering_Ints
mv *_peaks.bed output/peakFiles

mv *bedgraph output/bedgraphs

mv *_mapStats.txt output/indStats
mv *_chrStats.txt output/indStats
mv *_finalMapStats.txt output/indStats

mv *map_deMito_dedupMetrics.txt output/dedupMetrics
mv *map_deMito_dedup_insertMetrics output/insertMetrics

mv *indPeakNumber.txt output/indStats

mv *.pdf output/histos
mv output/histos output/combFiles

mv comb* output/combFiles
########################################################################

### Make folders and sort files

mkdir output

mkdir output/fastqs
mkdir output/sams
mkdir output/sortBams
mkdir output/mapBams
mkdir output/dedupBams
mkdir output/deblBams
mkdir output/deMitoBams
mkdir output/dedupBeds
mkdir output/deblBeds
mkdir output/peaks
mkdir output/unionBedgraphs
mkdir output/indBedgraphs
mkdir output/bws
mkdir output/histos
mkdir output/combFiles
mkdir output/indStats
mkdir output/insertMetrics
mkdir output/dedupMetrics

mv *fastq.gz output/fastqs
mv *sam output/sams
mv *out output/sams
mv *sort.bam *sort.bam.bai output/sortBams
mv *map.bam *map.bam.bai output/mapBams
mv *dedup.bam *dedup.bam.bai output/dedupBams
mv *map_deMito.bam *map_deMito.bam.bai output/deMitoBams
mv *map_deMito_dedup.bed output/dedupBeds
mv *debl.bed output/deblBeds
mv *debl.bam *debl.bam.bai output/deblBams
mv *_peaks.xls *_peaks.narrowPeak *_summits.bed *_peaks.bed output/peaks
mv *_counts.bedgraph output/unionBedgraphs
mv *indCounts.bedgraph output/indBedgraphs
mv *.bw output/bws
mv *_log output/bws

mv *_mapStats.txt output/indStats
mv *_chrStats.txt output/indStats
mv *_indPeakNumber.txt output/indStats

mv *map_deMito_dedupMetrics.txt output/dedupMetrics
mv *map_deMito_dedup_insertMetrics output/insertMetrics

mv *.pdf output/histos
mv output/histos output/combFiles
mv comb* *union* output/combFiles
