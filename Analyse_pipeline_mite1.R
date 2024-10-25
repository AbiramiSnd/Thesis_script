#####
# Automatic analysis of TRACKPOSON output
#Create final matrix 
#####

library(pheatmap)
library(tidyverse)
`%nin%` <- negate(`%in%`)
#warning eliminations
options(warn=-1)

# Arguments retrieve
args <- commandArgs(trailingOnly = TRUE)

# TE family name
te<-as.character(args[1])
path<-as.character(args[2])
final_cov<-as.character(args[3])

# Read insertions files
inputFile<-try(system(paste0("ls ",args[2],"/cove*.bed"),intern=TRUE))
files <- unlist(strsplit(inputFile,split=" "))

M<-data.frame()
M=read.table(paste0(args[2],"/final_cov2_",args[1],"/all_position_cov5_",args[1],".names"),sep="\n",h=F)

for (i in 1:length(files)){
  res<-read.table(files[i],sep="\t",h=F)
  res$ID <- paste(res$V1,res$V2,res$V3,sep="_")
  tmp<-as.numeric(res$V4[match(M$V1, res$ID)]) 
  M<-cbind(M,tmp)
}
colnames(M)<-c("insertion",files)

paste0(args[2],"/final_cov2_",args[1],"/all_insertion_",args[1],".names")

pos<-paste0(args[2],"/final_cov2_",args[1],"/all_insertion_",args[1],".names")
POS<-read.table(pos,sep="\n",h=F)

M2<-M[which(M$insertion %in% POS$V1),]

colnames(M2)=gsub('coveragebed_',"", colnames(M2))
colnames(M2) <- gsub(paste0(args[2],"/"), '', colnames(M2), fixed=TRUE)
colnames(M2)=gsub('-vs-[A-Za-z0-9_.]*',"", colnames(M2))

# Read total local coverage table of all TIPs windows
N2<-read.csv(file=args[4])

# Compare 2 matrix of TIPs insertion coverage and total local coverage
M2<- data.frame(M2[,-1, drop = FALSE], row.names = M2[,1])
N2<- data.frame(N2[,-1, drop = FALSE], row.names = N2[,1])

# Extract from total coverage matrix, only required column
N2=subset(N2, select = -c(AAKOOSDC,AAKROSDC,AAOAOSDC))
N3=N2[rownames(N2) %in% rownames(M2), ]

M2_tmp=setdiff(colnames(N2), colnames(M2))
M2[M2_tmp] <- 0
M2=M2[, colnames(N2)]

# Read total local coverage, TIP insertion coverage matrices and generate a 3rd numerical matrix 
# to distinguish homozygous (2), heterozygous (1), no insertion (0), no data (NA)

M2[!is.na(M2) & (((M2/N3)*100) >= 70)]<-2
M2[!is.na(M2) & (((M2/N3)*100) <=  70)]<-1
M2[is.na(M2) & N3>=2]<-0
M2[is.na(M2) & N3<2]<-"NA"

write.csv(M2,file=paste0("matrice_final_",args[1],".csv"),row.names = T,quote=F,sep="\t",col.names=T)


