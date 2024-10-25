grep "Fsyl" matrice_final_LTR*.csv   |  sed -e 's/,/\t/g' | sed 's/NA/0/g' | sed -e 's/\t2/\t1/g' | sed 's/matrice_final_//g' | sed 's/.csv:/\t/g' > ../all_matrice_final_LTR.interval-TE.bed

cat used_matrix2/matrice_final_LTR*.csv | grep "^Fsyl" | sed 's/,/\t/g' | awk -v OFS="\t" '{print $1}' | sed 's/_/\t/g' | awk -v OFS="\t" '{print $1"_"$2"_"$3,$4,$5}' | sort | uniq | bedtools sort >  all_matrice_final_LTR.interval.pos.bed
awk -v OFS="\t" '{print $2,$1}' all_matrice_final_LTR.interval-TE.bed |  sed 's/_/\t/g' | awk -v OFS="\t" '{print $1"_"$2"_"$3,$4,$5,$6}' |  sort | uniq | bedtools sort > all_matrice_final_LTR.interval.pos_TE.bed


bedtools closest -a all_matrice_final_LTR.interval.pos.bed -b '/home/lgdp/Desktop/distribution_sv/Fagus_Sylvatica_gene_annot_fn2.bed' -wa -wb -D "a" | sort | uniq > all_matrice_final_LTR.interval.pos-gene_dist.bed


bedtools makewindows -g '/media/lgdp/TOSHIBA/Données_Fagus_sylvatica/Fagus_sylvatica_v3.30Scaff.length'  -w 100000 > Fagus_sylvatica_v3.30Scaff.length.100kb.bed

bedtools intersect -a Fagus_sylvatica_v3.30Scaff.length.100kb.bed -b all_matrice_final_LTR.interval.pos.bed -wa -wb | awk -v OFS="\t" '{count[$1"\t"$2"\t"$3]++} END {for (word in count) print word, count[word]}' > all_matrice_final_LTR.interval.pos_100kb_count_TIPs_without_TE.bed

bedtools intersect -a Fagus_sylvatica_v3.30Scaff.length.100kb.bed -b all_matrice_final_LTR.interval.pos_TE.bed -wa -wb | awk -v OFS="\t" '{count[$1"\t"$2"\t"$3]++} END {for (word in count) print word, count[word]}' > all_matrice_final_LTR.interval.pos_100kb_count_TIPs_with_TE.bed

bedtools intersect -a "/media/lgdp/Seagate Basic/P2/Abirami/BA2_mobilome/hotspot_TIPs_log2_filt6.bed"  -b '/home/lgdp/Desktop/GWAS_plink/all_matrice_final_LTR.interval.pos_TE.bed'  -wa -wb | bedtools groupby -g 1,2,3,4 -c 8,8 -o count_distinct,distinct > hotspot_TIPs_log2_filt6_TE_family.bed
bedtools intersect -a hotspot_TIPs_log2_filt6.bed -b '/home/lgdp/Desktop/distribution_sv/Fagus_Sylvatica_gene_annot_fn2.bed' -wa -wb > hotspot_TIPs_log2_filt6-gene_annot.bed

cat *.bed | bedtools sort | bedtools groupby -g 1,2,3 -c 4 -o mean  > all_trees_cov.100kb.regions.mean.bed


 sort  -k1,2n '/home/lgdp/Desktop/GWAS_plink/all_matrice_final_LTR.interval.pos_TE.bed' | bedtools groupby -g 1,2,3 -c 4,4 -o distinct,count   > all_matrice_final_LTR.interval.pos_TE.concat_TEfam.bed

bedtools intersect  -a '/home/lgdp/Desktop/GWAS_plink/Fagus_sylvatica_v3.30Scaff.length.100kb.bed' -b /home/lgdp/Desktop/GWAS_plink/all_matrice_final_LTR.interval.pos_TE.bed -wa -wb | bedtools groupby -g 1,2,3 -c 7 -o count_distinct > '/home/lgdp/Desktop/GWAS_plink/all_matrice_final_LTR.interval.pos_TE.concat_TEfam_distinct.100kb.bed'

awk -v OFS="\t" '($2-2000>=0) {print $1,$2-2000,$3+2000,$4} ($2-2000<0) {print $1,0,$3+2000,$4}' '/home/lgdp/Desktop/distribution_sv/Fagus_Sylvatica_gene_annot_fn2.bed' > Fagus_Sylvatica_gene_annot_2kb.bed

bedtools intersect -a Fagus_Sylvatica_gene_annot_2kb.bed -b '/home/lgdp/Desktop/GWAS_plink/all_matrice_final_LTR.interval.pos_TE.bed' -wa -wb | awk '{count[$1"\t"$2"\t"$3"\t"$4"\t"$5]++} END {for (word in count) print word, count[word]}' | sort -k4n > Fagus_Sylvatica_gene_annot_2kb_tips_count.bed
