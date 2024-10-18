# The aim of this script is to map mobilome reads to reference genome

# Mapping mobilome nanopore reads to reference genome in paf format
minimap2 -t 8 --cs -cx map-ont '/media/lgdp/TOSHIBA/Genome_ref_Fagus_v3/Fagus_sylvatica_v3.fasta' '/home/lgdp/Desktop/216MR/216MR.fastq' > 216MR_20_06.paf

# Filtering paf file into reduced format
awk -v OFS="\t" '{print $6,$8,$9,$1,$2,$3,$4,$5,$12}' 216MR_20_06.paf | bedtools sort > 216MR_20_06.paf.bed


