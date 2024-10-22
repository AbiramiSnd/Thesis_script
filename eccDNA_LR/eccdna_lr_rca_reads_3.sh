# The aim of this script is to retrieve mobilome long reads supporting RCA


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

