# The aim of this script is analyse TRACKPOSON bedfile outputs to generate a csv numerical matrix per TE family with individuals as columns and windows as rows.
#!/usr/bin/bash


#name of TE familly
csv_tip=$4
path=$2
te=$1
echo $te
out=$te
all_tips=$3
#only coveragebed analyzed
mkdir $path/$out/final_cov2_$te
cd $path/$out/final_cov2_$te

# Remove dead trees samples
rm $path/$out/*AAKOOSDC*  $path/$out/*AAKROSDC*   $path/$out/*AAOAOSDC*

# Store all insertions of the specific TE across samples
for file in $path/$out/*.fa_per*kb.bed;do awk -F "\t" '{print $1"_"$2"_"$3}' $file;done | sort -u > all_insertion_$te.names

# Number of TE insertion
wc -l all_insertion_$te.names

# Retrieve insertions where 6 reads are found across samples
awk -F "\t" '{print $1"_"$2"_"$3,$4}' $path/$out/*.fa_per*kb.bed | sort | awk -v OFS="\t" '{a[$1]+=$2;}END{for (i in a)print i,a[i]}' | awk -v OFS="\t" '($2>=6){print $1}' > all_position_cov5_$te.names 
wc -l all_position_cov5_$te.names

#reformat output
for file in $path/$out/*.fa_per*kb.bed ;do n=$(echo $file | sed -e "s/coveragebed_//"  | sed -e "s/-vs-[A-Z0-9_]*.fa_per*kb.bed//"); awk -F "\t" '{print $1"_"$2"_"$3}' $file > $n.txt;done

echo $csv_tip
R CMD BATCH "--args $te $path/$out /final_cov2_$te  $csv_tip" $path/../Analyse_pipeline.R
