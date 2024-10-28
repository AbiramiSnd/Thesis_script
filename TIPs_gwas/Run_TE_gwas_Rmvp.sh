# The aim of this script to launch GWAS on all TE numerical matrices in a folder using a loop
#!/bin/bash

te_tmp=$1
te=${te_tmp/used_matrix2\//}
path=/home/lgdp/Desktop/GWAS_plink
out=${te/.csv/}_MLM_ACP
pheno=$path/nm_db_filtered_147_final_fevrier2024.txt
snp_vcf='/media/lgdp/P2/Abirami/vcf_SNP/Beech_147_complet_23062022.recode.only30Scaff.vcf'
vcf=$path/$1

# Calculate population structure
R CMD BATCH "--args $path $snp_vcf " /home/lgdp/Desktop/GWAS_plink/rmvp_pop-struct.R 

mkdir $out  
cd $out

# Run gwas
R CMD BATCH "--args $path/$out $vcf $pheno $path/mvp.SNP.vcf $out $path" /home/lgdp/Desktop/GWAS_plink/rmvp_GWAS.R 
