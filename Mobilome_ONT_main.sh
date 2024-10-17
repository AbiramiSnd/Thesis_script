

awk -v OFS="\t" '($9>=10) {print $0}' 194MR_20_06.paf.bed | bedtools merge -d 200 -c 4,4 -o distinct,count_distinct | awk -v OFS="\t" '{print $1,$2,$3,$5}' | bedtools sort > all_194MR_mobilome.paf.merged200_rd.bed
awk -v OFS="\t" '($6>2){print $0}' 194MR_mobilome_minimap_tmp.paf.py_out2 > 194MR_mobilome_minimap_tmp.paf.py_out3

bedtools intersect -a all_194MR_mobilome.paf.merged200_rd.bed -b 194MR_mobilome_minimap_tmp.paf.py_out3 -f 0.7 -loj  > all_194MR_mobilome.paf.merged200_rd_all_194MR_mobilome_mobilome_minimap_tmp.paf.py_out3_loj.bed
bedtools intersect -a all_194MR_mobilome.paf.merged200_rd_all_194MR_mobilome_mobilome_minimap_tmp.paf.py_out3_loj.bed -b '/home/lgdp/Desktop/FsylCur4_refTEs_wclassif_wreliable.frags.no_microsat.bed' -wa -wb -f 0.2 | awk -v OFS="\t" '{print $0}' | sort |  uniq  | grep -v "Unclassified" > all_194MR_mobilome.paf.merged200_rd_all_194MR_mobilome_mobilome_minimap_tmp.paf.py_out3_loj-vs-TEs.bed
