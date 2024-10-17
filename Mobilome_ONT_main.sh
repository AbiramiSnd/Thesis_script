#Mapping mobilome-seq reads in PAF format
minimap2 -t 8 --cs -cx map-ont '/media/lgdp/TOSHIBA/Genome_ref_Fagus_v3/Fagus_sylvatica_v3.fasta' '/home/lgdp/Desktop/194MR/194MR.fastq' > 194MR_20_06.paf

#select reads with 200 bp size; mapping quality of 10; and simplify mapping file
awk -v OFS="\t" '($5>=200 && $9>=10){print $1"_"$2"_"$3"_"$4"_"$5"_"$8,$6,$7,$7-$6,$8,$9}' 194MR_mobilome_minimap.paf.bed | bedtools sort | bedtools groupby -g 1 -c 1,4,2,3,6 -o count,mean,collapse,collapse,mean  | sed 's/_/\t/g' | awk -v OFS="\t" '{print $1"_"$2"_"$3""$4"_"$5""$6"_"$7"_"$8,$1"_"$2"_"$3,$4,$5,$6,$7,$8,$9,$11,$12,$13}' > 194MR_mobilome_minimap.paf_filt.bed

# python script to calculate repetition within reads, where if the read id is the same we count the number of time the reads maps to a locus
python ../repeats_mobilome0.py 194MR_mobilome_minimap.paf_filt.bed > 194MR_mobilome_minimap.paf_filt.bed.16_06.py_out

# Number of reads ic counted by incrementation, so we need to remove duplicated by selecting highest number of units
awk -v OFS="\t" '$6 >= max[$1,$2,$3,$4] {max[$1,$2,$3,$4]=$6;row[$1,$2,$3,$4]=$0} END { for (i in row) print row[i] }' BA2_mobilome_minimap_tmp.paf.py_out | bedtools sort > BA2_mobilome_minimap_tmp.paf.py_out2

# Select reads with Mapping quality >= 10; merge read coordinates mapping within 200bp overlap to calculate total coverage per locus with read alignments
awk -v OFS="\t" '($9>=10) {print $0}' 194MR_20_06.paf.bed | bedtools merge -d 200 -c 4,4 -o distinct,count_distinct | awk -v OFS="\t" '{print $1,$2,$3,$5}' | bedtools sort > all_194MR_mobilome.paf.merged200_rd.bed

# Select reads having atleast 3 repeat units
awk -v OFS="\t" '($6>2){print $0}' 194MR_mobilome_minimap_tmp.paf.py_out2 > 194MR_mobilome_minimap_tmp.paf.py_out3

# Intersection between position with total coverage and repeats containing locus to find regions hcovered by both and calculate a coverage table
bedtools intersect -a all_194MR_mobilome.paf.merged200_rd.bed -b 194MR_mobilome_minimap_tmp.paf.py_out3 -f 0.7 -loj  > all_194MR_mobilome.paf.merged200_rd_all_194MR_mobilome_mobilome_minimap_tmp.paf.py_out3_loj.bed

# Intersection with coverage table and TE annoation
bedtools intersect -a all_194MR_mobilome.paf.merged200_rd_all_194MR_mobilome_mobilome_minimap_tmp.paf.py_out3_loj.bed -b '/home/lgdp/Desktop/FsylCur4_refTEs_wclassif_wreliable.frags.no_microsat.bed' -wa -wb -f 0.2 | awk -v OFS="\t" '{print $0}' | sort |  uniq  | grep -v "Unclassified" > all_194MR_mobilome.paf.merged200_rd_all_194MR_mobilome_mobilome_minimap_tmp.paf.py_out3_loj-vs-TEs.bed
