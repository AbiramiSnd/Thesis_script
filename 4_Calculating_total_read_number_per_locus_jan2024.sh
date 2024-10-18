# The aim of this script is to count the number of reads mapping per locus

# Selecting reads with MQ >=10
# Merging reads alignment positions overlapping over 200 bp and count the number of reads mapping
awk -v OFS="\t" '($9>=10) {print $0}' all_V25-001_mobilome.paf.bed | bedtools merge -d 200 -c 4,4 -o distinct,count_distinct | awk -v OFS="\t" '{print $1,$2,$3,$5}' | bedtools sort > all_V25-001_mobilome_minimap.paf.merged200_rd.bed
