rm(list=ls())

#set the working director
set.seed(123)
input.path="./0_Data/"
output.path="./2_diversity/"
dir.create(output.path, showWarnings = FALSE)
ncores=5

#load libraries
library(rio)

#import data
dades=import(paste(input.path,"pheno_2021_06_21.tsv",sep=""))

#load functions
source("./functions/linear_model.R")

#prepare data
dades=dades[,c("casctot","agev1","gender","siteid",
               "q005a","plate","smokestatus",
               "q134","diab_treat","HBP_treat",
               "HC_treat","log.fibrer","log.energi",
               "shannon____mgs")]
               
yi="casctot"
x="shannon____mgs"

#log1p transform 
dades[,yi]=log1p(dades[,yi])

#covari basic model
covari=c("agev1","gender","siteid","siteid:plate","q005a")

dades$gender=as.factor(dades$gender)
dades$siteid=as.factor(dades$siteid)
dades$q005a=as.factor(dades$q005a)
dades$plate=as.factor(dades$plate)
dades$smokestatus=as.factor(dades$smokestatus)
dades$q134=as.factor(dades$q134)
dades$diab_treat=as.factor(dades$diab_treat)
dades$HBP_treat=as.factor(dades$HBP_treat)
dades$HC_treat=as.factor(dades$HC_treat)

#run the model
res.basic=fastlm.fun(x,yi,dades,covari=covari)

#covari full model
covari=c(covari,"smokestatus","q134","diab_treat",
         "HBP_treat","HC_treat","log.fibrer",
         "log.energi")

#run the model
res.full=fastlm.fun(x,yi,dades,covari=covari)

#export
write.table(res.basic,file=paste(output.path,"/res_alpha_diversity_basic_model.tsv",sep=""),col.names=T,row.names=F,sep="\t")
write.table(res.full,file=paste(output.path,"/res_alpha_diversity_full_model.tsv",sep=""),col.names=T,row.names=F,sep="\t")

