# The aim of this script is to find per TE locus, the number of total read mapping and the number of reads supporting RCA

# Intersection
#	- file with locus from reference genome with total number of reads
#	- file with locus from reference genome with number of reads RCA supporting reads
bedtools intersect -a all_194MR_mobilome.paf.merged200_rd.bed -b 194MR_mobilome_minimap_tmp.paf.py_out3 -f 0.7 -loj  > all_194MR_mobilome.paf.merged200_rd_all_194MR_mobilome_mobilome_minimap_tmp.paf.py_out3_loj.bed

# Intersecting 
#	- file with total reads/RCA supporting reads per locus 
#	- filtered TE annotation
bedtools intersect -a all_194MR_mobilome.paf.merged200_rd_all_194MR_mobilome_mobilome_minimap_tmp.paf.py_out3_loj.bed -b '/home/lgdp/Desktop/FsylCur4_refTEs_wclassif_wreliable.frags.no_microsat.bed' -wa -wb -f 0.2 | awk -v OFS="\t" '{print $0}' | sort |  uniq  | grep -v "Unclassified" > all_194MR_mobilome.paf.merged200_rd_all_194MR_mobilome_mobilome_minimap_tmp.paf.py_out3_loj-vs-TEs.bed

