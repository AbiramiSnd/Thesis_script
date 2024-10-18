# Author : Abirami SOUNDIRAMOURTTY
# PhD student, Laboratoire Génome et développement des plantes (Perpignan, France)
# abirami.soundiramourtty@univ-perp.fr 

# Usage : The input used are ONT data of #354 individual and 
# Beech reference genome from VERZY project 
# The aim of this script is to detect structural variants (insertions and 
# deletions) using a mapping approach. Minimap2, bedtools and CigarSV 
# is required.    

# Mapping reads to a reference genome

minimap2 -t 8 --cs -cx map-ont Fagus_sylvatica_v3.fasta CYR_AAAA_ONT_1_PAG73847_pass_reads_3kb.fastq > CYR_AAAC_ONT_1_PAG73847_Fagus_sylvatica_v3_pass_reads_3kb.paf

# SV detection using CigarSV

python cigar_sv_fasta.py -i CYR_AAAC_ONT_1_PAG73847_Fagus_sylvatica_v3_pass_reads_3kb.paf  -r Fagus_sylvatica_v3.fasta  -l CYR_AAAA_ONT_1_PAG73847_pass_reads_3kb.fastq

# Filter 1 = select SV based on reads coverage, mapping quality and alignment length

# merge fragmented alignment of reads and calculate start and end of alignemnent 
cat CYR*.csv | awk -v OFS='\t' 'NR>1 {print $1"_"$6,$0}' | sort | groupBy -g 1 -c 4,5,3 -o min,max,distinct -full | awk -v OFS='\t' '{print $0,($4/$3)*100,(($3-$5)/$3)*100}' | cut -f 2-16,20,21 > sv_mutant_2-21.csv

# merge reads containing sv in 20 pb range as 1 SV 
awk -v OFS='\t' 'NR>1 {print $6,$8,$9,$1,$10,$11,$12,$16,$17}' sv_mutant_2-21.csv  | bedtools sort | bedtools merge -d 20 -c 4,5,6,4,7,8,9 -o count,distinct,mean,distinct,mean,mean,mean > sv_mutant_2-21.bed

# select sv based on start /end of alignment of reads and coverage of SV
awk -v OFS='\t' '$8>=40 && $6>=100 && (($4>=10 && $9<=15 && $10<=15) || ($4>=62)) { print $1,$2,$3,$6}'  sv_mutant_2-21.bed > sv_mutant_2-21_filtered.bed

# Filter 2 = Compare M/R to keep false negatives

# Determine shared SV between M/R

bedtools intersect -a sv_mutant_2-21_filtered.bed -b ../../Rév*/CIG*/sv_revertant_2-21_filtered.bed -wo  | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | awk '$9>max[$4,$5,$6]{max[$4,$5,$6]=$9; row[$4,$5,$6]=$0} END{for (i in row) print row[i]}' | bedtools sort |awk '$9>max[$1,$2,$3]{max[$1,$2,$3]=$9; row[$1,$2,$3]=$0} END{for (i in row) print row[i]}' | awk -v OFS="\t" '{print $1,$2,$3,$4}' | uniq > sh_sv_mutant_2-21_filtered.bed

# Intersection between initial cigarSV output and filtered SV to define false SV
bedtools intersect -a sv_mutant_2-21.bed -b sv_mutant_2-21_filtered.bed -v -wa | awk -v OFS="\t" '{print $1,$2,$3,$6}' > false_sv_mutant_2-21_filtered.bed

# Determine specific SV of M/R
bedtools intersect -a sv_mutant_2-21_filtered.bed -b sh_sv_mutant_2-21_filtered.bed -v -wa | awk -v OFS="\t" '{print $0}' > sp_sv_mutant_2-21_filtered.bed

# Considering high FN rate in detection, correction is required by comparing M/R
# Gather FN SV of revertant : Intersection between M filtered SV and false SV 
bedtools intersect -a sp_sv_mutant_2-21_filtered.bed -b false_sv_revertant_2-21_filtered.bed -wo  | awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' |  awk '$9>max[$1,$2,$3]{max[$1,$2,$3]=$9; row[$1,$2,$3]=$0} END{for (i in row) print row[i]}' | bedtools sort | awk '$9>max[$4,$5,$6]{max[$4,$5,$6]=$9; row[$4,$5,$6]=$0} END{for (i in row) print row[i]}' |awk -v OFS="\t" '{print $5,$6,$7,$8}' | uniq > sh2_sv_mutant_2-21_filtered.bed

bedtools intersect -a sv_mutant_2-21_filtered.bed -b sv_revertant_2-21_filtered.bed -wo  | gawk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | gawk '$9>max[$1,$2,$3]{max[$1,$2,$3]=$9; row[$1,$2,$3]=$0} END{for (i in row) print row[i]}' | bedtools sort | gawk '$9>max[$4,$5,$6]{max[$4,$5,$6]=$9; row[$4,$5,$6]=$0} END{for (i in row) print row[i]}' |  gawk -v OFS="\t" '{print $1,$2,$3,$4}' | uniq > sh_sv_mutant_2-21_filtered_tmp.bed

bedtools intersect -a sh_sv_mutant_2-21_filtered_tmp.bed -b sh_sv_revertant_2-21_filtered_tmp.bed -wo | gawk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | gawk '$9>max[$1,$2,$3]{max[$1,$2,$3]=$9; row[$1,$2,$3]=$0} END{for (i in row) print row[i]}' | bedtools sort | gawk '$9>max[$4,$5,$6]{max[$4,$5,$6]=$9; row[$4,$5,$6]=$0} END{for (i in row) print row[i]}' | gawk -v OFS="\t" '{print $1,$2,$3,$4}' | uniq > sh_sv_mutant_2-21_filtered.bed

# Corection of specific SVs
bedtools intersect -a sp_sv_mutant_2-21_filtered.bed -b sh2_sv_mutant_2-21_filtered.bed -v -wa | awk -v OFS="\t" '{print $0}' > sp2_sv_mutant_2-21_filtered.bed 
