#!/bin/bash

#PBS -N bedmap_bam    
#PBS -l nodes=1:ppn=6   
#PBS -q long            
#PBS -M abirami.soundiramourtty@univ-perp.fr  
#PBS -m e          

    
#software
source ~/.bashrc
module load mosdepth/0.3.3

####You need to change the variables
###################################

FQ1=$file
FQ2=$(echo $file | sed -e "s/_1/_2/")
out=$(basename $file | awk -F "_" '{print $2}')

DIR=/home/abirami/MASSANE
bed=/home/abirami/MASSANE/all_500ltr_activeTE_ins.bed
bam=/home/database/Beech_genomes/Massane/MAPPING_v3/BWA/$out-vs-Fagus_sylvatica_v3.sort.bam

########################
#######################
#work in local
mkdir /scratch/TRACKPOSON-$out
cd  /scratch/TRACKPOSON-$out

#copy data in local directecoty
cp $bed .
ref=$(basename $bed)

######coveragebed by 10kb windows
#bam2bed < $bam | sort-bed - > $out-bam.bed
#bedmap --echo --fraction-map 0.8 --skip-unmapped --mean $bed $out-bam.bed > $out-all_ins_TE.bedmap_meanreads.bed

mosdepth --no-per-base --use-median --by $bed $out-all_500ltr_activeTE_ins.cov $bam
zcat $out-all_500ltr_activeTE_ins.cov.regions.bed.gz > $out-all_500ltr_activeTE_ins.cov.regions.bed

######cleaning temporary files
cp $out-all_500ltr_activeTE_ins.cov.regions.bed $DIR

#rm *vs-MITE_ABI_per10kb.regions.bed.gz
rm *mosdepth.*.txt
rm *.gz.csi
rm *.bed


#cleaning
rm -r /scratch/TRACKPOSON-$out
