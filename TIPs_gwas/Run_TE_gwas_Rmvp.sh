#!/bin/bash

te_tmp=$1
te=${te_tmp/used_matrix2\//}
path=/home/lgdp/Desktop/GWAS_plink
out=${te/.csv/}_MLM_ACP
#pheno=$path/table_beech_massane_phenotypes_2023_3_2
pheno=$path/nm_db_filtered_147_final_fevrier2024.txt
snp_vcf='/media/lgdp/P2/Abirami/vcf_SNP/Beech_147_complet_23062022.recode.only30Scaff.vcf'
#snp_vcf='/media/lgdp/TOSHIBA/Fagus_v4/MASSANE_BWA_Freebayes_167_allscaff_maf1_mm0.95.recode.vcf'
#vcf=$path/ltr2_matrix_reordered.vcf2
vcf=$path/$1

# reformat matrix
#echo $te

#R CMD BATCH "--args $path $snp_vcf " /home/lgdp/Desktop/GWAS_plink/rmvp_pop-struct.R 

mkdir $out  
cd $out
echo $out
R CMD BATCH "--args $path/$out $vcf $pheno $path/mvp.SNP.vcf $out $path" /home/lgdp/Desktop/GWAS_plink/rmvp_GWAS.R 

#rmdir /home/lgdp/Desktop/GWAS_plink/*MLM_ACP_MLM*
