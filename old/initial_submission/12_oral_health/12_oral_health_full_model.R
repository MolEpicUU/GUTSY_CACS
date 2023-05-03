rm(list=ls())

#set the working director
set.seed(123)
input.path="./0_Data/"
output.path="./12_oral_health/"
dir.create(output.path, showWarnings = FALSE)
ncores=15

#load libraries
library(rio)
library(BiocParallel)

#import files
pheno.oral=import(paste(input.path,"./MOS.csv",sep=""))
res=import(paste(input.path,"./11_gut_oral_microbiome/oral_gut_microbiome.tsv",sep=""))
res$var.x=paste(res$hg3a.oral,res$maintax.oral,sep=".")

#load functions
source("./functions/ordinal_model_clust.se.R")

#preapare data
pheno.oral$Sex_k=as.factor(pheno.oral$Sex_k)
pheno.oral$Smoking_ctr=factor(pheno.oral$Smoking_ctr,level=c(1,4,2))
pheno.oral$Education_ctr=as.factor(pheno.oral$Education_ctr)
pheno.oral$Antibiotics=as.factor(pheno.oral$Antibiotics)

pheno.oral$eat=as.factor(pheno.oral$eat)
pheno.oral$tabacco=as.factor(pheno.oral$tabacco)
pheno.oral$brush=as.factor(pheno.oral$brush)

strep.id=res$var.x[as.numeric(res$q.value)<0.05]
pheno.oral$caries=pheno.oral$DiS32_DS32
outcome=c("caries","BoP","FS_32")

var=c("lopnrMOS","Age","Sex_k","Smoking_ctr","Education_ctr","shannon.oral","Plaque_score",
      "eat","tabacco","brush","Antibiotics","Family_ID")

pheno.oral=pheno.oral[,c(var,outcome,strep.id)]
pheno.oral=pheno.oral[complete.cases(pheno.oral),]

names(pheno.oral)=gsub(" ","_",names(pheno.oral))
strep.id=gsub(" ","_",strep.id)

pheno.oral[,strep.id]=apply(pheno.oral[,strep.id],2,log1p)
pheno.oral=pheno.oral[complete.cases(pheno.oral),]

var=c("Age","Sex_k",
         "Smoking_ctr",
         "Education_ctr",
         "eat","tabacco","brush",
         "shannon.oral","Plaque_score","Family_ID")

covari=paste(var,collapse="+")

pheno.oral$caries=as.factor(pheno.oral$caries)
pheno.oral$BoP=as.factor(pheno.oral$BoP)
pheno.oral$FS_32=as.factor(pheno.oral$FS_32)

system.time(res.outc2<-bplapply(outcome,function(y)
{
  
  system.time(xres<-bplapply(strep.id,function(x){
    ordinal.fun(y,x,covari,var,pheno.oral)
    
  }, BPPARAM = MulticoreParam(ncores)))
  
  xres <- do.call(rbind, xres)
  
  xres
  
}, BPPARAM = MulticoreParam(1)))

res.outc2 <- do.call(rbind, res.outc2)
res.outc2$q.value=p.adjust(res.outc2$p.value,"BH")

res.outc2=res.outc2[order(res.outc2$p.value),]

#Export data
write.table(res.outc2,file=paste(output.path,"/res_oral_health.tsv",sep=""),row.names=F,col.names=T,sep="\t")

