library(devtools)
library(rMVP)
library(dplyr)
library(stringr)
library(tidyverse)


args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])

### Generate numerical TIPs matrix
matrix=read.csv(file=args[2],header=T,row.names=1)


#matrix=read.table(file=args[2],sep="\t",header=T,row.names=1)
out_matrix=paste0(args[1],"/",args[5],".num")
ind=paste0(args[6],"/mvp.SNP.vcf.geno.ind")
paste0(args[6],"/mvp.SNP.vcf.geno.ind")
order_ind=read.table(ind,sep="\t")
#
#dim(order_ind)
#dim(matrix)

#matrix=matrix[ , order(order_ind[,1])]


#matrix=matrix[match(colnames(matrix),order_ind[["V1"]]), ]   
#matrix <- data[match(order_ind[["V1"]], colnames(matrix)), ]     

dim(matrix)


#
colnames_order=order_ind[["V1"]]

idx=match(colnames_order,colnames(matrix))
matrix_reordered=matrix[,idx]


matrix_reordered=matrix_reordered[grep("Fsylvatica_scaffold_[0-9]_|Fsylvatica_scaffold_1[0-9]_|Fsylvatica_scaffold_2[0-9]_|Fsylvatica_scaffold_30_", row.names(matrix_reordered)),]

write.table(matrix_reordered,out_matrix,quote=F,row.names = F,col.names = F,sep="\t")
dim(matrix_reordered)

### Generate map for TIPs matrix
#map=rownames(matrix_reordered)
#split_map=str_split_fixed(map, '_', 5)
#head(split_map)
#map=cbind(map,split_map)
#map=data.frame(map[,c(1,4,5)])  #TE per TE

#all_ins
#dim(matrix_reordered)
map=rownames(matrix_reordered)
head(map)
split_map=str_split_fixed(map, '_', 6)
head(split_map)
map=cbind(map,split_map)
map=data.frame(map[,c(1,4,5)])

#V4
#split_map=str_split_fixed(map, '_', 4)
#head(split_map)
#map=cbind(map,split_map)
#head(map)
#map=data.frame(map[,c(1,1,5)])


colnames(map)=c("TE","Chr","Pos")
map <- map %>%
  mutate(Chr =  gsub('_[0-9]*_[0-9]*$', '', Chr))
head(map)

write.table(map,"TE.map",quote=F,sep="\t",col.names = T,row.names = T)

### Read phenotype
phenotype=read.table(file =args[3] ,sep="\t",header=T,row.names = 1)

out_matrix
### Prepare mvp numerical matrix
MVP.Data(fileNum=out_matrix,
         filePhe=args[3],
         fileMap="TE.map",
         sep.num="\t",
         sep.map="\t", 
         sep.phe="\t",
         fileKin=FALSE,
         filePC=FALSE,
         auto_transpose = FALSE,
         #priority="memory",
         maxLine=100000,
         out="mvp.num"
)

### Read rmvp generated genotype matrix
genotype = attach.big.matrix("mvp.num.geno.desc")
dim(genotype)

### Read rmvp generated Hapmap
map=read.table(file="mvp.num.geno.map",header=T,sep="\t")
dim(map)
### Calculate kinship from mvp_geno_file
Kinship=MVP.Data.Kin(TRUE, mvp_prefix=args[4], out='mvp')
dim(Kinship)

### Calculate PCA from mvp_geno_file
MVP.Data.PC(TRUE, mvp_prefix=args[4], pcs.keep=5, out='mvp')
ACP=attach.big.matrix("mvp.pc.desc")
dim(ACP)

### Run GWAS
for(i in 2:ncol(phenotype)){
  GWAS=MVP(phe=phenotype[,c(1, i)], geno=genotype, map=map, K=Kinship,
           priority=speed, ncpus = 10, vc.method="EMMA",
           CV.FarmCPU=ACP,method = "FarmCPU", method.bin = "FaST-LMM",
           maxLoop = 5, permutation.threshold=TRUE, permutation.rep=100,
           file.output=T)
  GWAS=MVP(phe=phenotype[,c(1, i)], geno=genotype, map=map, K=Kinship,
           priority=speed, ncpus = 10, vc.method="EMMA",
           CV.MLM=ACP,method = "MLM", method.bin = "FaST-LMM",
           maxLoop = 5, permutation.threshold=TRUE, permutation.rep=100,
           file.output=T)
  gc()
} #method = "FarmCPU" CV.FarmCPU=ACP, CV.MLM=ACP


