list=$(ls used_matrix2/*.csv)
#list=$(ls 500_LTR_2.csv)
for i in $list
do

./Run_TE_gwas_Rmvp.sh $i
#./Run_TE_gwas_Rmvp.sh $i
done
