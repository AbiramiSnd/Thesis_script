# The aim of this script is to map ONT mobilome reads to reference genome and simplify the PAF file, to retrieve mobilome long reads supporting RCA, 
# to count the number of reads mapping per locus and filter reads with repeats (having atleast 3 repeat units), to find per TE locus, the number of total read mapping and the number of reads supporting RCA

# Mapping mobilome nanopore reads to reference genome in paf format
minimap2 -t 8 --cs -cx map-ont '/media/lgdp/TOSHIBA/Genome_ref_Fagus_v3/Fagus_sylvatica_v3.fasta' '/home/lgdp/Desktop/194MR/194MR.fastq' > 194MR_20_06.paf

# Filtering paf file into reduced format
awk -v OFS="\t" '{print $6,$8,$9,$1,$2,$3,$4,$5,$12}' 194MR_20_06.paf | bedtools sort > 194MR_20_06.paf.bed

# Filter alignment for
#	- reads atleast 200 pb length
#	- MQ > 10
	
# bedtools groupby is to group reads and count the number of time the same read is mapping to a locus
# a unique ID is generated for each read mapping to a single locus by concatenating the columns about reads and mapping positions on reference genome : $1"_"$2"_"$3""$4"_"$5""$6"_"$7"_"$8,$1"_"$2"_"$3
# Format the file to have the below output
#
awk -v OFS="\t" '($5>=200 && $9>=10) {print $1"_"$2"_"$3"_"$4"_"$5"_"$8,$6,$7,$7-$6,$8,$9}' 194MR_20_06.paf.bed  | bedtools groupby -g 1 -c 1,4,2,3,6 -o count,mean,collapse,collapse,mean  | sed 's/_/\t/g' | awk -v OFS="\t" '{print $1"_"$2"_"$3""$4"_"$5""$6"_"$7"_"$8,$1"_"$2"_"$3,$4,$5,$6,$7,$8,$9,$11,$12,$13}' > 194MR_20_06.paf.filt.txt

# Custom python script is used to filter reads supporting RCA and is based on identifying reads mapping to the same locus multiple times
python ../repeats_mobilome0.py 194MR_20_06.paf.filt.txt > 194MR_mobilome_minimap.paf_filt.bed.20_06.py_out

# Removing duplicates from output and selecting reads mapping atleast 2 times to a same locus
awk -v OFS="\t" '$6 >= max[$1,$2,$3,$4] {max[$1,$2,$3,$4]=$6;row[$1,$2,$3,$4]=$0} END { for (i in row) print row[i] }' 194MR_mobilome_minimap.paf_filt.bed.20_06.py_out |  awk -v OFS="\t" '($6>=2) {print $0}' | bedtools sort | bedtools sort > 194MR_mobilome_minimap_tmp.paf.py_out2

# Selecting reads with MQ >=10
# Merging reads alignment positions overlapping over 200 bp and count the number of reads mapping
awk -v OFS="\t" '($9>=10) {print $0}' 194MR_20_06.paf.bed | bedtools merge -d 200 -c 4,4 -o distinct,count_distinct | awk -v OFS="\t" '{print $1,$2,$3,$5}' | bedtools sort > all_194MR_mobilome.paf.merged200_rd.bed

# Select reads having atleast 3 repeat units
awk -v OFS="\t" '($6>2){print $0}' 194MR_mobilome_minimap_tmp.paf.py_out2 > 194MR_mobilome_minimap_tmp.paf.py_out3

# Intersection
#	- file with locus from reference genome with total number of reads
#	- file with locus from reference genome with number of reads RCA supporting reads
bedtools intersect -a all_194MR_mobilome.paf.merged200_rd.bed -b 194MR_mobilome_minimap_tmp.paf.py_out3 -f 0.7 -loj  > all_194MR_mobilome.paf.merged200_rd_all_194MR_mobilome_mobilome_minimap_tmp.paf.py_out3_loj.bed

# Intersecting 
#	- file with total reads/RCA supporting reads per locus 
#	- filtered TE annotation
bedtools intersect -a all_194MR_mobilome.paf.merged200_rd_all_194MR_mobilome_mobilome_minimap_tmp.paf.py_out3_loj.bed -b '/home/lgdp/Desktop/FsylCur4_refTEs_wclassif_wreliable.frags.no_microsat.bed' -wa -wb -f 0.2 | awk -v OFS="\t" '{print $0}' | sort |  uniq  | grep -v "Unclassified" > all_194MR_mobilome.paf.merged200_rd_all_194MR_mobilome_mobilome_minimap_tmp.paf.py_out3_loj-vs-TEs.bed
