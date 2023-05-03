rm(list=ls())

#set the working director
set.seed(123)
input.path="/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/github/casc_microbiome/Demo/0_Data/"
output.path="/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/github/casc_microbiome/Demo/13_oral_health/"
dir.create(output.path, showWarnings = FALSE)
ncores=15

#load libraries
library(rio)
library(BiocParallel)

#import files
pheno.oral=import(paste(input.path,"/MOS.csv",sep=""))
res=import(paste(input.path,"../12_gut_oral_microbiome/oral_gut_microbiome.tsv",sep=""))
res$var.x=paste(res$hg3a.oral,res$maintax.oral,sep=".")

#load functions
source("/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/github/casc_microbiome/scripts/Demo/ordinal_model_clust.se.R")

#preapare data
pheno.oral$caries=pheno.oral$DiS32_DS32

pheno.oral$Sex_k=as.factor(pheno.oral$Sex_k)
pheno.oral$Smoking_ctr=factor(pheno.oral$Smoking_ctr,level=c(1,4,2))
pheno.oral$Education_ctr=as.factor(pheno.oral$Education_ctr)
pheno.oral$Antibiotics=as.factor(pheno.oral$Antibiotics)
pheno.oral$PPI_All=as.factor(pheno.oral$PPI_All)

pheno.oral$eat=as.factor(pheno.oral$eat)
pheno.oral$tabacco=as.factor(pheno.oral$tabacco)
pheno.oral$brush=as.factor(pheno.oral$brush)

pheno.oral$caries=as.factor(pheno.oral$caries)
pheno.oral$BoP=as.factor(pheno.oral$BoP)
pheno.oral$FS_32=as.factor(pheno.oral$FS_32)


strep.id=res$var.x[as.numeric(res$q.value)<0.05]
outcome=c("caries","BoP","FS_32")

var=c("lopnrMOS","Age","Sex_k","Smoking_ctr","Education_ctr","shannon.oral","Plaque_score",
      "eat","tabacco","brush","Antibiotics","BMI","PPI_All","Family_ID")

pheno.oral=pheno.oral[,c(var,outcome,strep.id)]
pheno.oral=pheno.oral[complete.cases(pheno.oral),]

names(pheno.oral)=gsub(" ","_",names(pheno.oral))
strep.id=gsub(" ","_",strep.id)

pheno.oral[,strep.id]=apply(pheno.oral[,strep.id],2,log1p)
pheno.oral=pheno.oral[complete.cases(pheno.oral),]

var=c("Age","Sex_k","Smoking_ctr","Education_ctr","shannon.oral","Plaque_score",
      "eat","tabacco","brush","Antibiotics","BMI","PPI_All","Family_ID")

covari=paste(var,collapse="+")

system.time(res.outc3<-bplapply(outcome,function(y)
{
  
  system.time(xres<-bplapply(strep.id,function(x){
    ordinal.fun(y,x,covari,var,pheno.oral)
    
  }, BPPARAM = MulticoreParam(ncores)))
  
  xres <- do.call(rbind, xres)
  
  xres
  
}, BPPARAM = MulticoreParam(1)))

res.outc3 <- do.call(rbind, res.outc3)
res.outc3$q.value=p.adjust(res.outc3$p.value,"BH")

res.outc3=res.outc3[order(res.outc3$p.value),]

#Export data
write.table(res.outc3,file=paste(output.path,"/res_oral_health_bmi_ppi_antib.tsv",sep=""),row.names=F,col.names=T,sep="\t")

