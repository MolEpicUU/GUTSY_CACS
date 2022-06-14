rm(list=ls())

#set the working director
set.seed(123)
input.path="./"
output.path="4_taxonomic_enrichment/"
dir.create(output.path, showWarnings = FALSE)
ncores=5

#load libraries
library(rio)
library(BiocParallel)

dir.create(output.path, showWarnings = FALSE)
ncores=5

#import data
res_lm.mod1=import(paste(input.path,"/3_species_level/res_cacs_mgs_basic_model.tsv",sep=""))
dades.info=import(paste(input.path,"./0_Data/HG3.A.7_tax.xlsx",sep=""))

#load functions
source(paste(input.path,"./functions/fgsea.tax.R",sep=""))

#prepare data
res_lm.mod1$mgs=matrix(unlist(strsplit(res_lm.mod1$var.x,"____")),ncol=2,byrow=T)[,2]
names(dades.info)[1]="mgs"

#stratify for the direction of the estimate
res_lm.mod1.pos=res_lm.mod1[res_lm.mod1$estimate>0,]
res_lm.mod1.neg=res_lm.mod1[res_lm.mod1$estimate<0,]

#run the models
set.seed(123)
res_lm.mod1.pos=fgsea.tax.fun(res_lm.mod1.pos,dades.info,"genus","casctot")
res_lm.mod1.pos=res_lm.mod1.pos[order(res_lm.mod1.pos$pval),c(1:2,5:8,3,4)]
names(res_lm.mod1.pos)=c("Trait","Genus","log2err","ES","NES","size","p.value","q.value")

set.seed(123)
res_lm.mod1.neg=fgsea.tax.fun(res_lm.mod1.neg,dades.info,"genus","casctot")
res_lm.mod1.neg=res_lm.mod1.neg[order(res_lm.mod1.neg$pval),c(1:2,5:8,3,4)]
names(res_lm.mod1.neg)=c("Trait","Genus","log2err","ES","NES","size","p.value","q.value")

#merge together
res_lm.mod1.posneg=merge(res_lm.mod1.pos,res_lm.mod1.neg,by="Genus",all=T)
res_lm.mod1.posneg=res_lm.mod1.posneg[,!names(res_lm.mod1.posneg)%in%grep("Trait",names(res_lm.mod1.posneg),value=T)]

#clean the table and export
names(res_lm.mod1.posneg)=c("Genus","log2err.pos","ES.pos","NES.pos","size.pos","p.value.pos","q.value.pos","log2err.neg","ES.neg","NES.neg","size.neg","p.value.neg","q.value.neg")
res_lm.mod1.posneg[,c(2:4,8:10)]=round(res_lm.mod1.posneg[,c(2:4,8:10)],digit=5)
res_lm.mod1.posneg=res_lm.mod1.posneg[order(res_lm.mod1.posneg$p.value.pos),]

for(i in grep("value",names(res_lm.mod1.posneg),value=T))
{
  res_lm.mod1.posneg[,i]=ifelse(res_lm.mod1.posneg[,i]<0.001,formatC(res_lm.mod1.posneg[,i],format="e",digits=2),round(res_lm.mod1.posneg[,i],digit=3))
}

write.table(res_lm.mod1.posneg,file=paste(output.path,"res_enrichment_genus_pos_neg.tsv",sep=""),row.names=F,col.names=T,sep="\t")

