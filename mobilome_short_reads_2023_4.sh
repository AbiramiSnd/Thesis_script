cd /media/lgdp/TOSHIBA/Verzy/Mobilome_Verzy/Resultats_eccFinder_Mapping/eccFinder_output_185M*
cat *num.bed | awk -v OFS="\t" '{a[$1"\t"$2"\t"$3]+=$4}END{for(i in a) print i,a[i]}' | sort | uniq | awk -v OFS="\t" '($4>=5){print $1,$2,$3}' | bedtools sort | uniq > ecc.sr.disc.split.num.bed
bedtools intersect -a ecc.sr.disc.split.num.bed -b align*/ecc.sr.sorted_185M.bed  -wa -wb | bedtools groupby -g 1,2,3 -c 7 -o count_distinct |  awk -v OFS="\t" '($4>=20 && ($3-$2)>=100){print $0}' > ecc.sr.disc.split.num.rd.bed
bedtools intersect -a ecc.sr.disc.split.num.rd.bed -b '/home/lgdp/Desktop/FsylCur4_refTEs_wclassif_wreliable.frags.no_microsat.bed' -wa -wb -f 0.2 -r | awk -v OFS="\t" ' ($5!=".") {print $0}' | sort |  uniq  | grep -v "Unclassified"   | awk -v OFS="\t" '{print $1,$2,$3,$4,$8}' | bedtools sort > ecc.sr.sorted_185M-TEs.bed

cd /media/lgdp/TOSHIBA/Verzy/Mobilome_Verzy/Resultats_eccFinder_Mapping/eccFinder_output_185R*
cat *num.bed | awk -v OFS="\t" '{a[$1"\t"$2"\t"$3]+=$4}END{for(i in a) print i,a[i]}' | sort  | uniq | awk -v OFS="\t" '($4>=5){print $1,$2,$3}' | bedtools sort | uniq > ecc.sr.disc.split.num.bed
bedtools intersect -a ecc.sr.disc.split.num.bed -b align*/ecc.sr.sorted_185R.bed  -wa -wb | bedtools groupby -g 1,2,3 -c 7 -o count_distinct |  awk -v OFS="\t" '($4>=20 && ($3-$2)>=100){print $0}' > ecc.sr.disc.split.num.rd.bed
bedtools intersect -a ecc.sr.disc.split.num.rd.bed -b '/home/lgdp/Desktop/FsylCur4_refTEs_wclassif_wreliable.frags.no_microsat.bed' -wa -wb -f 0.2 -r | awk -v OFS="\t" ' ($5!=".") {print $0}' | sort |  uniq  | grep -v "Unclassified"   | awk -v OFS="\t" '{print $1,$2,$3,$4,$8}' | bedtools sort > ecc.sr.sorted_185R-TEs.bed

