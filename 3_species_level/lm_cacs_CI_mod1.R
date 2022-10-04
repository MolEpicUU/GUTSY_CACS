rm(list=ls())

#set the working director
set.seed(123)
input.path="./0_Data/"
output.path="./3_species_level/"
dir.create(output.path, showWarnings = FALSE)
ncores=5

#load libraries
library(rio)
library(BiocParallel)

#import data
dades=import(paste(input.path,"pheno_2021_06_21.tsv",sep=""))

#load functions
source("./functions/linear_model.R")

#prepare data
dades=dades[,c("casctot","agev1","gender","siteid","q005a","plate",grep("____mgs",names(dades),value=T),grep("HG3A.",names(dades),value=T))]

noms=grep("HG3A.",names(dades),value=T)


yi="casctot"
covari=c("agev1","gender","siteid","siteid:plate","q005a")

dades$gender=as.factor(dades$gender)
dades$siteid=as.factor(dades$siteid)
dades$q005a=as.factor(dades$q005a)
dades$plate=as.factor(dades$plate)

print("#### model ####")
print(paste(yi,"~ microbiota+",paste(covari,collapse="+"),sep=""))


dades[,c(yi,noms)]=apply(dades[,c(yi,noms)], 2, log1p)

print("##### time #####")
system.time(res<-bplapply(noms,function(x){
  res <- fastlm.fun(x,yi,dades,covari,log=FALSE,rank.1=FALSE)
  data.frame(res)
}, BPPARAM = MulticoreParam(ncores)))
res <- do.call(rbind, res)
res$q.value <- p.adjust(res$p.value, method = "BH")


#export
write.table(res,file=paste(output.path,"/res_cacs_mgs_basic_model.tsv",sep=""),col.names=T,row.names=F,sep="\t")


print(Sys.time())

