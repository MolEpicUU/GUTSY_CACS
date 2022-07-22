rm(list=ls())

#set the working director
set.seed(123)
input.path="/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/github/casc_microbiome/phenotypes/"
output.path="/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/github/casc_microbiome/results/9_metabolites/"
dir.create(output.path, showWarnings = FALSE)
ncores=10

#load libraries
library(rio)
library(ppcor)
library(BiocParallel)
library(ppcor)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape)
library(plyr)
library(ggpubr)
library(fBasics)

#import data
dades=import(paste(input.path,"pheno_2021_09_06_metab.tsv",sep=""))
res=import(paste(output.path,"/res_spearman_corr_metab_mod2.tsv",sep=""))

# load functions
source("/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/github/casc_microbiome/scripts/functions/linear_model.R")

#prepare data
pathway=unique(c("Sphingomyelins","Partially Characterized Molecules", "Androgenic Steroids", "Plasmalogen", "Secondary Bile Acid Metabolism","Acetylated Peptides","Primary Bile Acid Metabolism","Drug - Analgesics, Anesthetics"))
metab.hit=unique(res$var.y[res$SUB_PATHWAY%in%pathway])

dades$smokestatus=as.factor(dades$smokestatus)
dades$gender=as.factor(dades$gender)
dades$siteid=as.factor(dades$siteid)
dades$q134=as.factor(dades$q134)
dades$q005a=as.factor(dades$q005a)
dades$diab_treat=as.factor(dades$diab_treat)
dades$HBP_treat=as.factor(dades$HBP_treat)
dades$HC_treat=as.factor(dades$HC_treat)
dades$plate=as.factor(dades$plate)

strep=unique(grep("Streptococcus",res$var.x,value=T))

yi="casctot"
dades[,c(yi,strep)]=apply(dades[,c(yi,strep)],2,log1p)

dades=dades[,c(yi,strep,"agev1","gender","siteid","q005a","smokestatus",
                          "q134","diab_treat","HBP_treat","HC_treat",
                          "plate","log.fibrer","log.energi",
                          "shannon____mgs",metab.hit)]

dades <- dades[, colMeans(is.na(dades)) < 0.3]
metab.hit=grep("MET",names(dades),value=T)

var=c("agev1","gender","siteid","plate","q005a","smokestatus",
      "q134","diab_treat","HBP_treat","HC_treat",
      "log.fibrer","log.energi",
      "shannon____mgs",metab.hit)

covari=c("agev1","gender","siteid","siteid:plate","q005a","smokestatus",
         "q134","diab_treat","HBP_treat","HC_treat",
         "log.fibrer","log.energi",
         "shannon____mgs",metab.hit)

dades=dades[complete.cases(dades[,c(yi,strep,var,metab.hit)]),]

system.time(res4<-bplapply(strep,function(x){
  res4 <- fastlm.fun(x,yi,dades,covari)
  data.frame(res4)
}, BPPARAM = MulticoreParam(10)))
res4 <- do.call(rbind, res4)
res4=data.frame(res4,stringsAsFactors = F)
res4$q.value=p.adjust(res4$p.value,"BH")

res4=res4[order(res4$p.value),]


covari=c("agev1","gender","siteid","siteid:plate","q005a","smokestatus",
         "q134","diab_treat","HBP_treat","HC_treat",
         "log.fibrer","log.energi",
         "shannon____mgs")


system.time(res3<-bplapply(strep,function(x){
  res3 <- fastlm.fun(x,yi,dades,covari)
  data.frame(res3)
}, BPPARAM = MulticoreParam(10)))
res3 <- do.call(rbind, res3)
res3=data.frame(res3,stringsAsFactors = F)
res3$q.value=p.adjust(res3$p.value,"BH")

res3=res3[order(res3$p.value),]

#Export data
write.table(res4,file=paste(output.path,"/res_CACS_adj_metab.tsv",sep=""),col.names=T,row.names=F,sep="\t")

write.table(res3,file=paste(output.path,"/res_complete.cases_metab.tsv",sep=""),col.names=T,row.names=F,sep="\t")

