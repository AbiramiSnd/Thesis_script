# The aim of this script is to concatenate all matrices and find position of TIPs and number of families involves.

# Select header of file with individual names for LTR TIPs
grep ""  LTR*/*final*/*csv | head -1  | sed 's/,/\t/g' | sed 's/LTR100\/final_cov2_LTR100\/matrice_final_LTR100.csv:/TE\tins/g' > all_LTRs_header.txt
# Reformat matrix for LTR TIPs
grep ""  LTR*/*final*/*csv | grep "Fsyl" | sed 's/,/\t/g' | sed 's/matrice_final_/\t/g' | sed 's/.csv:/\t/g'  | cut -f2- | sed 's/NA/0/g' | sed 's/\t2/\t1/g' > all_LTRs_ins_TE_ind.csv.bed
cat all_LTRs_header.txt all_LTRs_ins_TE_ind.csv.bed  > all_LTRs_ins_TE_ind2.csv.bed

# Select header of file with individual names for MITEs TIPs
grep ""  */*final*/*csv | head -1  | sed 's/,/\t/g' | sed 's/Cluster10\/final_cov2_Cluster10\/matrice_final_Cluster10.csv:/TE\tins/g' > all_MITEs_header.txt
# Reformat matrix for MITE TIPs
grep ""  */*final*/*csv | grep "Fsylvatica" | sed 's/,/\t/g' | sed 's/matrice_final_/\t/g' | sed 's/.csv:/\t/g'  | cut -f2- | sed 's/NA/0/g' | sed 's/\t2/\t1/g' > all_MITEs_ins_TE_ind.csv.bed
cat all_MITEs_header.txt all_MITEs_ins_TE_ind.csv.bed  > all_MITEs_ins_TE_ind2.csv.bed

# Output bedfile of position of TIP with LTR family involved
awk -v OFS="\t" '{print $2,$1}'  all_LTRs_ins_TE_ind2.csv.bed | sed 's/\(Fsylvatica_scaffold_[0-9]*\)_\([0-9]*\)_\([0-9]*\)/\1\t\2\t\3/g' | grep -v "TE" | bedtools sort > all_LTRs_ins_TE_ind2.csv_win.bed
# Output bedfile of position of TIP with the number of different number of families
awk -v OFS="\t" '{print $2,$1}'  all_LTRs_ins_TE_ind2.csv.bed | sed 's/\(Fsylvatica_scaffold_[0-9]*\)_\([0-9]*\)_\([0-9]*\)/\1\t\2\t\3/g' | grep -v "TE" | awk -v OFS="\t" '{count[$1"\t"$2"\t"$3]++} END {for (word in count) print word, count[word]}' | bedtools sort > all_LTRs_ins_TE_ind2.csv_count_fam.bed

# Output bedfile of position of TIP with MITE family involved
awk -v OFS="\t" '{print $2,$1}'  all_MITEs_ins_TE_ind2.csv.bed | sed 's/\(Fsylvatica_scaffold_[0-9]*\)_\([0-9]*\)_\([0-9]*\)/\1\t\2\t\3/g' | grep -v "TE" | bedtools sort > all_MITEs_ins_TE_ind2.csv_win.bed
# Output bedfile of position of TIP with the number of different number of families
awk -v OFS="\t" '{print $2,$1}'  all_MITEs_ins_TE_ind2.csv.bed | sed 's/\(Fsylvatica_scaffold_[0-9]*\)_\([0-9]*\)_\([0-9]*\)/\1\t\2\t\3/g' | grep -v "TE" | awk -v OFS="\t" '{count[$1"\t"$2"\t"$3]++} END {for (word in count) print word, count[word]}' | bedtools sort > all_MITEs_ins_TE_ind2.csv_count_fam.bed

# Remove LTR-TIPs windows opened by aberrant number of different families
bedtools intersect -a '/home/lgdp/Desktop/TRACKPOSON_test/TRACKPOSON_test2/LTR/all_LTRs_ins_TE_ind2.csv_win.bed'  -b all_LTRs_ins_TE_ind2.csv_win_grp_remove.bed  -v | awk '{print $1"_"$2"_"$3}' | sort | uniq  > all_LTRs_ins_TE_ind2.csv_win_removed.txt
grep -f all_LTRs_ins_TE_ind2.csv_win_removed.txt '/home/lgdp/Desktop/TRACKPOSON_test/TRACKPOSON_test2/LTR/all_LTRs_ins_TE_ind2.csv.bed' > /home/lgdp/Desktop/TRACKPOSON_test/TRACKPOSON_test2/LTR/all_LTRs_ins_TE_ind2.csv2.bed

# Remove MITE-TIPs windows opened by aberrant number of different families
bedtools intersect -a '/home/lgdp/Desktop/TRACKPOSON_test/TRACKPOSON_test3/MITE/all_MITEs_ins_TE_ind2.csv_win.bed'  -b all_MITEs_ins_TE_ind2.csv_win_grp_remove.bed  -v | awk '{print $1"_"$2"_"$3}' | sort | uniq > all_MITEs_ins_TE_ind2.csv_win_removed.txt
grep -f all_MITEs_ins_TE_ind2.csv_win_removed.txt '/home/lgdp/Desktop/TRACKPOSON_test/TRACKPOSON_test3/MITE/all_MITEs_ins_TE_ind2.csv.bed' > /home/lgdp/Desktop/TRACKPOSON_test/TRACKPOSON_test2/LTR/all_MITEs_ins_TE_ind2.csv2.bed
