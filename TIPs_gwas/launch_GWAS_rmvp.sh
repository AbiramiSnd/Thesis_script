# The aim of this script is to run GWAS for all matrices using a loop
list=$(ls used_matrix2/*.csv)

for i in $list
do
./Run_TE_gwas_Rmvp.sh $i
done
