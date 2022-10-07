rm(list=ls())

#set the working director
set.seed(123)
input.path="./"
output.path="./10_inflammation_infection/"
dir.create(output.path, showWarnings = FALSE)
ncores=1

#load libraries
library(rio)
library(BiocParallel)

#import data
res_lm.mod2=import(paste(input.path,"/3_species_level/res_cacs_mgs_full_model.tsv",sep=""))
dades=import(paste(input.path,"./0_Data/pheno_2021_06_21.tsv",sep=""))

#load functions
source(paste(input.path,"/functions/linear_model.R",sep=""))

#preparate data
bact=res_lm.mod2[res_lm.mod2$q.value<0.05,"var.x"]

dades$gender=as.factor(dades$gender)
dades$smokestatus=as.factor(dades$smokestatus)
dades$siteid=as.factor(dades$siteid)
dades$q134=as.factor(dades$q134)
dades$q005a=as.factor(dades$q005a)
dades$diab_treat=as.factor(dades$diab_treat)
dades$HBP_treat=as.factor(dades$HBP_treat)
dades$HC_treat=as.factor(dades$HC_treat)
dades$plate=as.factor(dades$plate)

dades$log_CRP=log(dades$hscrp_res)
dades$log_neut=log(dades$neut_res)
dades$log_leukt=log(dades$lpk_res)

dades[,bact]=apply(dades[,bact],2,log1p)

yi=c("log_CRP","log_neut","log_leukt")

covari=c("agev1","gender","siteid","siteid:plate","q005a","smokestatus",
         "q134","diab_treat","HBP_treat","HC_treat","log.fibrer",
         "log.energi")

var=c("agev1","gender","siteid","q005a","smokestatus",
               "q134","diab_treat","HBP_treat","HC_treat","plate","log.fibrer",
               "log.energi", "bmi","ppi")

xdades=dades[,c(var,bact,yi)]

xdades=xdades[complete.cases(xdades[,var]),]

res<-bplapply(yi,function(i){
  
  xres<-bplapply(bact,function(x){
    xres <- fastlm.fun(x,i,xdades,covari,log=FALSE,rank.1=FALSE)
    data.frame(xres)
  }, BPPARAM = MulticoreParam(5))
  xres <- do.call(rbind, xres)
  
  return(xres)
}, BPPARAM = MulticoreParam(1))

res <- do.call(rbind, res)
res=data.frame(res,stringsAsFactors = F)
res$q.value=p.adjust(res$p.value,"BH")

bact2=unique(res$var.x[res$q.value<0.05])
covari=c("agev1","gender","siteid","siteid:plate","q005a","smokestatus",
         "q134","diab_treat","HBP_treat","HC_treat","log.fibrer","log.energi","bmi","ppi")

res2<-bplapply(yi,function(i){
  
  xres<-bplapply(bact2,function(x){
    xres <- fastlm.fun(x,i,xdades,covari,log=FALSE,rank.1=FALSE)
    data.frame(xres)
  }, BPPARAM = MulticoreParam(5))
  xres <- do.call(rbind, xres)
  
  return(xres)
}, BPPARAM = MulticoreParam(1))

res2 <- do.call(rbind, res2)
res2=data.frame(res2,stringsAsFactors = F)
res2$q.value=p.adjust(res2$p.value,"BH")

#export data
write.table(res,file=paste(output.path,"/res_inflammations_infections_full.model.tsv",sep=""),col.names=T,row.names=F,sep="\t")
write.table(res2,file=paste(output.path,"/res_inflammations_infections_full_bmi_ppi.tsv",sep=""),col.names=T,row.names=F,sep="\t")

