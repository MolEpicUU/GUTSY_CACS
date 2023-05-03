rm(list=ls())

#set the working director
set.seed(123)
input.path="./0_Data/"
output.path="./6_carotid/"
dir.create(output.path, showWarnings = FALSE)
ncores=15

#load libraries
library(rio)
library(BiocParallel)

#import data
dades=import(paste(input.path,"pheno_2021_06_21.tsv",sep=""))
import.data=import(paste(output.path,"../3_species_level/res_cacs_mgs_basic_model.tsv",sep=""))

#load functions
source("./functions/ordinal_model.R")

#prepare data
import.data=unique(import.data[import.data$q.value<0.05,"var.x"])
noms=grep("HG3A.",import.data,value=T)

dades=dades[,c("casctot","carotidplaque","agev1","gender","siteid","q005a","plate",noms,"shannon____mgs")]

var=c("agev1","gender","siteid","plate","q005a","shannon____mgs")

covari=c("agev1","gender","siteid","siteid:plate","q005a",
"shannon____mgs")

dades$gender=as.factor(dades$gender)
dades$siteid=as.factor(dades$siteid)
dades$q005a=as.factor(dades$q005a)
dades$plate=as.factor(dades$plate)
dades$carotidplaque=as.factor(dades$carotidplaque)

yi="carotidplaque"

dades[,noms]=apply(dades[,noms],2,log1p)

print("#### model ####")
print("ordinal model")
print(paste(yi,"~ microbiota+",paste(covari,collapse="+"),sep=""))

print("##### time #####")
system.time(res.ordinal<-bplapply(noms,function(x){
  res.ordinal<-ordinal.fun(yi,x,covari,var,dades)
  data.frame(res.ordinal)
}, BPPARAM = MulticoreParam(ncores)))
res <- do.call(rbind, res.ordinal)

res$q.value <- p.adjust(res$p.value, method = "BH")


write.table(res,file=paste(output.path,"/res_ordinal_carotidplaques_basic_model.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

print(Sys.time())

