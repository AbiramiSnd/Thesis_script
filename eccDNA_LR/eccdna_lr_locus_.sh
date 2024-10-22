# The aim of this script is to count the number of reads mapping per locus and filter reads with repeats (having atleast 3 repeat units)

# Selecting reads with MQ >=10
# Merging reads alignment positions overlapping over 200 bp and count the number of reads mapping
awk -v OFS="\t" '($9>=10) {print $0}' 194MR_20_06.paf.bed | bedtools merge -d 200 -c 4,4 -o distinct,count_distinct | awk -v OFS="\t" '{print $1,$2,$3,$5}' | bedtools sort > all_194MR_mobilome.paf.merged200_rd.bed

# Select reads having atleast 3 repeat units
awk -v OFS="\t" '($6>2){print $0}' 194MR_mobilome_minimap_tmp.paf.py_out2 > 194MR_mobilome_minimap_tmp.paf.py_out3
