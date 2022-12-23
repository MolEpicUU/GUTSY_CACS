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
res_mod=import(paste(output.path,"../3_species_level/res_cacs_mgs_full_model.tsv",sep=""))
met_info=import(paste(input.path,"scapis_merged_annotations_batchnorm_clean.tsv",sep=""))

# load functions
source("/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/github/casc_microbiome/scripts/functions/spearman_correlation_androgenic_metab.R")

#prepare data
sel1=grep("HG3A.",res_mod[res_mod$q.value<0.05,"var.x"],value=T)
sel1=grep("Streptococcus",sel1,value=T)

dades$gender=as.factor(dades$gender)
dades$siteid=as.factor(dades$siteid)
dades$q005a=as.factor(dades$q005a)
dades$plate=as.factor(dades$plate)
dades$smokestatus=as.factor(dades$smokestatus)
dades$q134=as.factor(dades$q134)
dades$diab_treat=as.factor(dades$diab_treat)
dades$HBP_treat=as.factor(dades$HBP_treat)
dades$HC_treat=as.factor(dades$HC_treat)

names(met_info)[1]=c("var.y")

metab=met_info[which(met_info$SUB_PATHWAY%in%"Androgenic Steroids"),"var.y"]

res.andro.1=spearman.cor.fun(sel1,dades[dades$gender%in%1,],
                      covari1=c("agev1","siteid","siteid:plate","q005a","smokestatus",
                               "q134","diab_treat","HBP_treat","HC_treat",
                               "log.fibrer","log.energi",
                               "shannon____mgs"),
                      cores=ncores)

res.andro.1$sex="male"

res.andro.2=spearman.cor.fun(sel1,dades[dades$gender%in%2,],
                             covari1=c("agev1","siteid","siteid:plate","q005a","smokestatus",
                                       "q134","diab_treat","HBP_treat","HC_treat",
                                       "log.fibrer","log.energi",
                                       "shannon____mgs"),
                             cores=ncores)

res.andro.2$sex="female"

res.andro=data.frame(rbind(res.andro.1,res.andro.2),stringsAsFactors = F)
res.andro$q.value=p.adjust(res.andro$p.value,"BH")

res.andro.1=res.andro[res.andro$sex%in%"male",]
res.andro.2=res.andro[res.andro$sex%in%"female",]

res.andro_test=merge(res.andro.1,res.andro.2,by=c("var.x","var.y"))

z2=(1/2)*log((1+res.andro_test$estimate.y)/(1-res.andro_test$estimate.y))
z1=(1/2)*log((1+res.andro_test$estimate.x)/(1-res.andro_test$estimate.x))

chr.covari=c("siteid","q005a","smokestatus",
         "q134","diab_treat","HBP_treat","HC_treat","plate")
num.covari=c("agev1", "log.fibrer","log.energi",
             "shannon____mgs")

x=0
for (i in chr.covari){
 xx=length(levels(dades[,i]))-1
 x=x+xx
 }

res.andro_test$x2=(z2-z1)/(sqrt((1/sqrt(res.andro_test$n.y-x-3))^2+(1/sqrt(res.andro_test$n.x-x-3))^2))
res.andro_test$p.value_homo.test=2*(1-pnorm(abs(res.andro_test$x2)))
res.andro_test$q.value_homo.test=p.adjust(res.andro_test$p.value_homo.test,"BH")

res.andro_test=merge(res.andro_test,met_info,by="var.y")
res.andro_test=res.andro_test[,c(1:7,9:14,16:19,27,21,22)]
res.andro_test=res.andro_test[order(res.andro_test$p.value_homo.test),]

#Export data
write.table(res.andro_test,file=paste(output.path,"/res_sex_strata_androgenic_steroids.tsv",sep=""),col.names=T,row.names=F,sep="\t")

