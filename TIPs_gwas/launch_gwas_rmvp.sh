# The aim of this script is to run GWAS for all matrices using a loop
list=$(ls used_matrix2/*.csv)

for i in $list
do
./run_TE_gwas_rmvp.sh $i
done
