#!/bin/bash

if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi


module load python/3.9.1

source $HOME/my_python-3.9.1/bin/activate
export PYTHONPATH=$HOME/my_python-3.9.1/lib/python3.9/site-packages


###############################################################

# Call Peaks


for i in `ls *_shift.bam`;
do
	sampleID="${i%.*}"
	echo "${sampleID}"
	
	macs3 callpeak -f BAMPE -t  ${sampleID}.bam -g hs -n  ${sampleID} -B -q 0.01 --nomodel --shift -75 --extsize 150 --keep-dup all
	
done
echo $?

