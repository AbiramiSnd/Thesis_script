# Author : Abirami SOUNDIRAMOURTTY
# PhD student, Laboratoire Génome et développement des plantes
# abirami.soundiramourtty@univ-perp.fr 

# The aim of this script is to filter raw ONT reads 
# NanoFilt is required.

# Select Pass reads

awk '$1 ~ /pass/ {print $0; getline; print $0; getline;print $0;getline;print $0}' sample.fastq > sample_pass_reads.fastq 

# Select reads with length >= 3kb

NanoFilt -l 3000  sample_pass_reads.fastq >  sample_pass_reads.3kb.fastq
