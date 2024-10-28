# The aim of to call insertions of a TE consensus on all samples and generate a bedfile of TIPs containing windows per individual

#!/bin/bash

#PBS -N trackposon     
#PBS -l nodes=1:ppn=6   
#PBS -q long            
#PBS -M abirami.soundiramourtty@univ-perp.fr  
#PBS -m e          

    
# Load tools
module load bowtie2/2.4.4
module load ncbi-blast+/2.12.0
module load samtools/1.13
module load bedtools2/2.30.0
module load mosdepth

source ~/.bashrc

# Define variables
FQ1=$file
FQ2=$(echo $file | sed -e "s/_1/_2/")
out=$(basename $file | awk -F "_" '{print $2}')
te=Cluster59.fa
DIR=/home/abirami/MASSANE/TRACKPOSON_Cluster59
DB=/home/abirami/MASSANE/MITEs_DB/clusters/mafft_Cluster59_seqs_trimal_gap_cons
blast_ref_database=/home/database/Beech_genomes/reference/Fagus_sylvatica_v3.fasta
perl_script=/home/abirami/MASSANE/find_insertion_point.pl
ref_genome_2kbpwindowsbed=/home/abirami/Fagus_sylvatica_v3_2kbwindows.bed
bam=/home/database/Beech_genomes/Massane/MAPPING_v3/BWA/$out-vs-Fagus_sylvatica_v3.sort.bam


# Work in local
mkdir /scratch/TRACKPOSON-$out
cd  /scratch/TRACKPOSON-$out

#copy data in local directecoty
cp $FQ1 .
cp $FQ2 .
cp $DB* .
cp $blast_ref_database* .
cp $ref_genome_2kbpwindowsbed .
fq1=$(basename $FQ1)
fq2=$(basename $FQ2)
db=$(basename $DB)
ref=$(basename $blast_ref_database)
wref=$(basename $ref_genome_2kbpwindowsbed)

# Mapping reads against TE reference
bowtie2 --time --end-to-end  -k 1 --very-fast -p 6 -x $DB  -1 $fq1 -2 $fq2  | samtools view -bS -@ 6 - > "$out"-vs-"$te".bam

# Keep only unmap reads with flag unmap/map and generate fasta 
samtools view "$out"-vs-"$te".bam | awk -F "\t" '{if ( ($1!~/^@/) && (($2==69) || ($2==133) || ($2==165) || ($2==181) || ($2==101) || ($2==117)) ) {print ">"$1"\n"$10}}' > $out-vs-$te.fa

# Blast fa against reference genome (Fagus v3) for identification insertion point
blastn -db $blast_ref_database -query $out-vs-$te.fa -out $out-vs-$te.fa.bl -num_threads 6 -evalue 1e-50

# Parse blast to find TE insertion point 
perl $perl_script $out-vs-$te.fa.bl $out-vs-$te

# Sort bed output
sort -k1,1 -k2,2n $out-vs-$te.bed > $out-vs-$te.sort.bed

# Coveragebed by 2kb windows
bedtools coverage -counts -nonamecheck -a $ref_genome_2kbpwindowsbed -b $out-vs-$te.sort.bed | awk -F "\t" '{if ($4>=2){print $0}}' > coveragebed_$out-vs-$te\_per2kb.bed

# Cleaning temporary files
cp coveragebed_$out-vs-$te\_per2kb.bed $DIR
rm $out-vs-$te.bam
rm $out-vs-$te.fa*
rm $out-vs-$te.bed
rm coveragebed_$out-vs-$te\_per2kb.bed 

# Cleaning
rm -r /scratch/TRACKPOSON-$out
