#!/bin/bash

if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi

# Load in required modules
module load python/3.9.1
module load bowtie2/2.1.0 
module load samtools-1.1
module add java/openjdk-1.7.0
module load picard-tools-1.141
module load FastQC-0.11.2


# Pyton environment
source $HOME/my_python-3.9.1/bin/activate
export PYTHONPATH=$HOME/my_python-3.9.1/lib/python3.9/site-packages

# Copy reference files into working directory
cp $HOME/ref_files/hg38_bt2index/hg38* .
cp $HOME/ref_files/hg38.blacklist.bed .
cp $HOME/ref_files/hg38.chrom.sizes .


######################################################################

### QC Fastq files

mkdir ./fastqc
fastqc -t 8 *fastq.gz -o ./fastqc

# Combine FastQC files into one html document
cd ./fastqc/
multiqc .

######################################################################

cd ../..

### Alignment
# Mapping atac paired end files to hg38 using bowtie2

for i in *_R1.fastq.gz;
do
 sampleID="${i%_*_*.*}"
 echo "${sampleID}"
 bowtie2 --very-sensitive -X 2000 -p16 -x hg38 -1 ${sampleID}_R1.fastq.gz -2 ${sampleID}_R2.fastq.gz -S ${sampleID}.sam 2> ${sampleID}.out
done

# Combine all of the .out files produced from the sam alignment

for i in *.out;
do
	sampleID="${i%_*_*.*}"
	echo "${sampleID}" > ${sampleID}_out.txt
	cat $i >> ${sampleID}_out.txt
done

find . -name '*_out.txt' | xargs cat >> combOutSams.txt
rm *_out.txt

# STOPPING POINT -- Alignment reports

######################################################################


### Make sorted bam
for i in `ls *.sam`;
do
	sampleID="${i%.*}"
 	echo "${sampleID}"
	samtools view -bS -@ 16 $i | samtools sort - -m 2G -@ 16 -o ${sampleID}_sort.bam
done
echo $?
echo 'You have sorted bam files :)!'




# Make flagstat file to check conversion to bam
for i in `ls *_sort.bam`;
do
	echo $i >> combFlagstat_bam.txt
	samtools flagstat $i >> combFlagstat_bam.txt
done

echo $?




# Make index for sorted bam files
for i in `ls *_sort.bam`;
do
 samtools index $i
done
echo $?
echo 'You have indexes for sorted bam files!'



# Make bam files with only mapped reads
for i in `ls *_sort.bam`;
do
	sampleID="${i%_*.*}"
 	echo "${sampleID}"
	samtools view -bS -f 2 -@ 16 $i > ${sampleID}_map.bam
	samtools index ${sampleID}_map.bam
done
echo $?
echo 'You have mapped-only bam files!'




# Make bam files without mitochondrial reads

for i in `ls *_map.bam`;
do
    sampleID="${i%*.*}"
    echo "${sampleID}"
	samtools view -@ 16 -h $i | grep -v chrM | samtools sort -@ 16 -O bam -o ${sampleID}_deMito.bam
	samtools index ${sampleID}_deMito.bam
done
echo $?
echo 'You have bams without mitocondrial reads!'





### Remove duplicates

# Change 'Xmx_g for amount of available ram

for i in `ls *_deMito.bam`;
do
	sampleID="${i%.*}"
 	echo "${sampleID}"
	java -Xmx48g -jar $PICARD_JAR MarkDuplicates INPUT=$i OUTPUT=${sampleID}_dedup.bam METRICS_FILE=${sampleID}_dedupMetrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true
done

echo $?



# Combine picard duplicate metrics
for i in `ls *_dedupMetrics.txt`;
do
	sampleID="${i%_*_*_*.*}"
	awk 'FNR==8 {print FILENAME, $0}' ${sampleID}_map_deMito_dedupMetrics.txt >> comb_picardDedupStats1.txt
done

awk '{print $1"\t"$9"\t"$10"\t"$11}' comb_picardDedupStats1.txt > comb_picardDedupStats2.txt


echo $'sampleID\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE' | cat - comb_picardDedupStats2.txt > comb_picardDedupStats.txt
rm comb_picardDedupStats1.txt
rm comb_picardDedupStats2.txt




# Make index for deduped bam files
for i in `ls *_dedup.bam`;
do
 samtools index $i
done
echo $?
echo 'You have indexes for sorted bam files!'




### Get insert metrics -- fragment size information

# Get insert metrics -- file and histogram
for i in `ls *_dedup.bam`;
do
	sampleID="${i%.*}"
 	echo "${sampleID}"
	java -Xmx48g -jar $PICARD_JAR CollectInsertSizeMetrics I=$i OUTPUT=${sampleID}_insertMetrics H=${sampleID}_histo.pdf
done

echo $?





### Remove blacklisted regions

for i in `ls *_dedup.bam`;
do
	sampleID="${i%.*}"
	echo "${sampleID}"
	
	bedtools intersect -abam ${sampleID}.bam -b hg38.blacklist.bed -v | samtools sort -o ${sampleID}_debl.bam -	
	samtools index ${sampleID}_debl.bam
done

echo $?
echo 'You have bed files without blacklist regions!'






# Make flagstat file to check blacklist removal
for i in `ls *_dedup.bam`;
do
	echo $i >> combFlagstat_dedup_bam.txt
	samtools flagstat $i >> combFlagstat_dedup_bam.txt
done
echo $?

for i in `ls *_debl.bam`;
do
	echo $i >> combFlagstat_debl_bam.txt
	samtools flagstat $i >> combFlagstat_debl_bam.txt
done
echo $?



### Shift reads for Tn5 insertion

export PYTHONPATH=/gpfs/fs02/apps/python/3.9.1/lib/python3.9/site-packages


for i in `ls *_debl.bam`;
do
	sampleID="${i%.*}"
	echo "${sampleID}"
	alignmentSieve --numberOfProcessors 16 --ATACshift --bam ${sampleID}.bam -o ${sampleID}_shift_tmp.bam
	samtools sort -@ 16 -O bam -o ${sampleID}_shift.bam ${sampleID}_shift_tmp.bam
	samtools index ${sampleID}_shift.bam
done
echo $?





### Make normalized bigwig files

for i in `ls *_shift.bam`;
do
 sampleID="${i%.*}"
 echo "${sampleID}"
	
 lines=$(samtools view -c ${sampleID}.bam);\
 bedtools genomecov -ibam ${sampleID}.bam -bg -scale $(echo "1000000 / ${lines} " | bc -l) -g hg38.chrom.sizes | \
 wigToBigWig -clip stdin hg38.chrom.sizes ${sampleID}.bw 2> ${sampleID}_log
done

echo $?




### Get map and chromosome numbers with sorted bam and deduped_debl_bams


### With sort_bam

## Create file with mapping numbers

# Create file with row names for mapStats
echo sampleID >> rowNames_mapStats.txt
echo totalReads >> rowNames_mapStats.txt
echo totalMappedReads >> rowNames_mapStats.txt
echo unmappedReads >> rowNames_mapStats.txt
echo pairedMappedReads >> rowNames_mapStats.txt

# Get total reads, total mapped reads, total unmapped reads, properly paired reads -- for each sample
for i in *_sort.bam;
do
	sampleID="${i%_*.*}"
 	echo "${sampleID}"
	echo "${sampleID}" >> ${sampleID}_mapStats.txt
	samtools view $i -@ 16 | wc -l >> ${sampleID}_mapStats.txt
	samtools view -F4 $i -@ 16 | wc -l >> ${sampleID}_mapStats.txt
	samtools view -f4 $i -@ 16 | wc -l >> ${sampleID}_mapStats.txt
	samtools view -f2 $i -@ 16 | wc -l >> ${sampleID}_mapStats.txt

done

# Combine individual mapStat files into one file
find . -name '*_mapStats.txt' | sort | xargs paste | column -s $'\t' -t > combMapStats.txt




### Create file with chromosome stats

# Create file with row names for chrStats# Get number of reads per chromosome

echo sampleID >> rowNames_chrStats.txt
echo chr1 >> rowNames_chrStats.txt
echo chr2 >> rowNames_chrStats.txt
echo chr3 >> rowNames_chrStats.txt
echo chr4 >> rowNames_chrStats.txt
echo chr5 >> rowNames_chrStats.txt
echo chr6 >> rowNames_chrStats.txt
echo chr7 >> rowNames_chrStats.txt
echo chr8 >> rowNames_chrStats.txt
echo chr9 >> rowNames_chrStats.txt
echo chr10 >> rowNames_chrStats.txt
echo chr11 >> rowNames_chrStats.txt
echo chr12 >> rowNames_chrStats.txt
echo chr13 >> rowNames_chrStats.txt
echo chr14 >> rowNames_chrStats.txt
echo chr15 >> rowNames_chrStats.txt
echo chr16 >> rowNames_chrStats.txt
echo chr17 >> rowNames_chrStats.txt
echo chr18 >> rowNames_chrStats.txt
echo chr19 >> rowNames_chrStats.txt
echo chr20 >> rowNames_chrStats.txt
echo chr21 >> rowNames_chrStats.txt
echo chr22 >> rowNames_chrStats.txt
echo chr23 >> rowNames_chrStats.txt
echo chrX >> rowNames_chrStats.txt
echo chrY >> rowNames_chrStats.txt
echo chrM >> rowNames_chrStats.txt

# Get number of reads per chromosome -- for each samples
for i in *_sort.bam;
do
	sampleID="${i%_*.*}"
 	echo "${sampleID}" >> ${sampleID}_chrStats.txt
	samtools view $i chr1 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr2 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr3 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr4 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr5 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr6 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr7 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr8 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr9 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr10 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr11 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr12 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr13 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr14 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr15 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr16 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr17 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr18 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr19 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr20 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr21 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr22 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chr23 -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chrX -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chrY -@ 16 | wc -l >> ${sampleID}_chrStats.txt
	samtools view $i chrM -@ 16 | wc -l >> ${sampleID}_chrStats.txt
done

# Combine individual chrStat files into one file
find . -name '*_chrStats.txt' | sort | xargs paste > combChrStats.txt



### With fully processed bam file

echo sampleID >> rowNames_finalMapStats.txt
echo totalReads >> rowNames_finalMapStats.txt


for i in *_map_deMito_dedup_debl.bam;
do
	sampleID="${i%_*_*_*_*.*}"
 	echo "${sampleID}"
	echo "${sampleID}" >> ${sampleID}_finalMapStats.txt
	samtools view $i | wc -l >> ${sampleID}_finalMapStats.txt
done

# Combine individual mapStat files into one file
find . -name '*_finalMapStats.txt' | sort | xargs paste | column -s $'\t' -t > combFinalMapStats.txt




# STOPPING POINT -- Check all raw metrics and bigwigs generated above -- samples look ok? May want to remove bad samples before proceeding.

######################################################################



