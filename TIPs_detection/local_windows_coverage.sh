# The aim of this script is to calculate coverage of all insertions containing windows to calculate heterozygoty score reflecting allelic frequency.
#!/bin/bash

#PBS -N coverage    
#PBS -l nodes=1:ppn=6   
#PBS -q long            
#PBS -M abirami.soundiramourtty@univ-perp.fr  
#PBS -m e          

    
# load tools
source ~/.bashrc
module load mosdepth/0.3.3

# Variables
FQ1=$file
FQ2=$(echo $file | sed -e "s/_1/_2/")
out=$(basename $file | awk -F "_" '{print $2}')
DIR=/home/abirami/MASSANE
bed=/home/abirami/MASSANE/all_LTR_ins.bed
bam=/home/database/Beech_genomes/Massane/MAPPING_v3/BWA/$out-vs-Fagus_sylvatica_v3.sort.bam

# Work in local
mkdir /scratch/TRACKPOSON-$out
cd  /scratch/TRACKPOSON-$out

# Copy data in local directecoty
cp $bed .
ref=$(basename $bed)

# Total read coverage by all windows containing TIPs
mosdepth --no-per-base --use-median --by $bed $out-all_10kb.cov $bam
zcat $out-all_10kb.cov.regions.bed.gz > $out-all_10kb.cov.regions.bed

# Cleaning temporary files
cp $out-all_10kb.cov.regions.bed $DIR
rm *mosdepth.*.txt
rm *.gz.csi
rm *.bed


# Cleaning folder
rm -r /scratch/TRACKPOSON-$out
