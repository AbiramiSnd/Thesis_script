#####
# Automatic analysis of TRACKPOSON output
#Create final matrix 
#####
library(pheatmap)
library(tidyverse)
#warning eliminations
options(warn=-1)

#arguments retrieve
args <- commandArgs(trailingOnly = TRUE)

#TE family name
files<-as.character(args[1])
args[0]
args[1]
args[2]



# TIPs coverage matrix
inputFile<-try(system(paste0("ls ",files,"/*all_sprite.regions.bed"),intern=TRUE)) #retro_lauren.regions.bed

#inputFile<-try(system(paste0("ls ",files,"/*retro_lauren.regions.bed"),intern=TRUE)) 
files <- unlist(strsplit(inputFile,split=" "))
length(files)
N<-data.frame()

N=read.table(paste0(args[2]),sep="\n",h=F)

head(files)
for (i in 1:length(files)){
  
  res<-read.table(files[i],sep="\t",h=F)
  res$ID <- paste(res$V1,res$V2,res$V3,sep="_")
  tmp<-as.numeric(res$V4[match(N$V1, res$ID)]) 
  N<-cbind(N,tmp)
  
}

colnames(N)<-c("insertion",files)
pos<-paste0(args[2])
POS<-read.table(pos,sep="\n",h=F)
N2<-N[which(N$insertion %in% POS$V1),]

colnames(N2) <- gsub(paste0("/home/lgdp/Desktop/TRACKPOSON_test/TRACKPOSON_test2/SPRITE/SPRITE1/cov_TIPs_2kb/"), '', colnames(N2), fixed=TRUE)
colnames(N2) <- gsub(paste0('.cov.all_sprite.regions.bed'), '', colnames(N2), fixed=TRUE)
#colnames(N2) <- gsub(paste0('.cov.retro_lauren.regions.bed'), '', colnames(N2), fixed=TRUE)
write.csv(N2,args[3],quote=F,row.names = F,col.names = T)
