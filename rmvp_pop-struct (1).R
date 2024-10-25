
library(devtools)
library(rMVP)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

setwd(args[1])

### read vcf and generate vcf for population structure
MVP.Data(fileVCF=args[2],
         fileKin=FALSE,
         filePC=FALSE,
         out="mvp.SNP.vcf")

