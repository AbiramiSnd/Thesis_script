# The aim of this script is to loop analysis of TRACKPOSON outputs on all TEs at once
#!/usr/bin/bash

# Output folder
path=/home/lgdp/Desktop/TRACKPOSON_test/TRACKPOSON_test2/LTR

# Folder containing mosdepth coverage bedfile of with the total coverage of all windows 
all_tips=/home/lgdp/Desktop/TRACKPOSON_test/TRACKPOSON_test2/cov_TIPs

# Output of total coverage matrix
all_tips_file="all_tips_cov.csv"

for file in $all_tips/*regions.bed;do awk -F "\t" '{print $1"_"$2"_"$3}' $file; done  > $path/all_position_TIPs_cov_dup.names 
sort $path/all_position_TIPs_cov_dup.names | uniq > $path/all_position_TIPs_cov.names
rm $path/all_position_TIPs_cov_dup.names

# Calculate numerical matrix with total coverage with all possible windows
cd $path
R CMD BATCH "--args $all_tips $path/all_position_TIPs_cov.names $all_tips_file " $path/../all_cov_matrix.R

for i in $te
do
~/Desktop/TRACKPOSON_test/TRACKPOSON_test2/Analyse_pipeline.sh ${i/\//} $path $path/../all_position_TIPs_cov.names $path/$all_tips_file
done
