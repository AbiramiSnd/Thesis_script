# The aim of this script is to remove redundancy of families by clustering TE consensus sequences

# Creating TE library database in blast format
makeblastdb -in FsylCur4_denovoLibTEs.fa  -title FsylCur4_denovoLibTEs.fa -dbtype nucl

# Blast of TE library vs TE library
blastn -query FsylCur4_denovoLibTEs.fa -db FsylCur4_denovoLibTEs.fa -out FsylCur4_denovoLibTEs.bl -outfmt 6 -evalue 1e-100

# Clustering TE library 
silix -i 0.9 -r 0.9 FsylCur4_denovoLibTEs.fa FsylCur4_denovoLibTEs.bl | awk -v OFS="\t" '{print $2,"Cluster_"$1}' > FsylCur4_denovoLibTEs_cluster.fa.clustr
