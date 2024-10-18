# The aim of this script is to find per TE locus, the number of total read mapping and the number of reads supporting RCA

# Intersecting 
#	- file with locus from reference genome with total number of reads
#	- file with locus from reference genome with number of reads RCA supporting reads
bedtools intersect -a all_V25-001_mobilome_minimap.paf.merged200_rd.bed -b all_V25-001_mobilome_minimap_tmp.paf.py_out2 -f 0.75 -loj  > all_V25-001_mobilome_minimap.paf.merged200_rd_all_V25-001_mobilome_minimap_tmp.paf.py_out2_loj.bed


# Intersecting 
#	- file with total reads/RCA supporting reads per locus 
#	- filtered TE annotation

bedtools intersect -a all_V25-001_mobilome_minimap.paf.merged200_rd_all_V25-001_mobilome_minimap_tmp.paf.py_out2_loj.bed -b '/home/lgdp/Desktop/FsylCur4_refTEs_wclassif_wreliable.frags.no_microsat.bed' -wa -wb -f 0.25 | awk -v OFS="\t" '{print $0}' | sort |  uniq  | grep -v "Unclassified" > all_V25-001_mobilome_minimap.paf.merged200_rd_all_V25-001_mobilome_minimap_tmp.paf.py_out2_loj-vs-TEs.bed
