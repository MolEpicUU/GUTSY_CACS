rm(list=ls())

#set the working director
set.seed(123)
input.path="./0_Data/"
output.path="./5_sex_strat_species/"
dir.create(output.path, showWarnings = FALSE)
ncores=5

#load libraries
library(rio)
library(BiocParallel)

#import data
dades=import(paste(input.path,"pheno_2021_06_21.tsv",sep=""))
import.data=import(paste(output.path,"../3_species_level/res_cacs_mgs_basic_model.tsv",sep=""))


#load functions
source("./functions/linear_model.R")

#prepare data
import.data=unique(import.data[import.data$q.value<0.05,"var.x"])
noms=grep("HG3A.",names(dades),value=T)

dades=dades[,c("casctot","agev1","gender","siteid","q005a","plate","smokestatus","q134","diab_treat","HBP_treat","HC_treat","log.fibrer","log.energi",
               grep("____mgs",names(dades),value=T),
               noms)]

yi="casctot"
covari=c("agev1","siteid","siteid:plate","q005a","smokestatus","q134","diab_treat","HBP_treat","HC_treat","log.fibrer","log.energi","shannon____mgs")

dades$gender=as.factor(dades$gender)
dades$siteid=as.factor(dades$siteid)
dades$q005a=as.factor(dades$q005a)
dades$plate=as.factor(dades$plate)
dades$smokestatus=as.factor(dades$smokestatus)
dades$q134=as.factor(dades$q134)
dades$diab_treat=as.factor(dades$diab_treat)
dades$HBP_treat=as.factor(dades$HBP_treat)
dades$HC_treat=as.factor(dades$HC_treat)

print("#### model ####")
print(paste("lm model"))
print(paste(yi,"~ microbiota+",paste(covari,collapse="+"),sep=""))


dades[,c(yi,noms)]=apply(dades[,c(yi,noms)], 2, log1p)


dades1=dades[dades$gender%in%1,]

a1=dades1[,grep("HG3A.",names(dades1),value=T)]
a1=a1[,colSums(a1)!=0]
noms1=names(a1)
rm(a1)

dades2=dades[dades$gender%in%2,]

a1=dades2[,grep("HG3A.",names(dades2),value=T)]
noms2=names(a1)
rm(a1)

noms1.mgs=noms1
noms2.mgs=noms2

print(paste("microbiota1=",length(noms1)))
print(paste("microbiota2=",length(noms2)))

print(dim(dades1))
print(dim(dades2))


print("##### time #####")
system.time(res1<-bplapply(noms1,function(x){
  res1 <- fastlm.fun(x,yi,dades1,covari,log=FALSE,rank.1=FALSE)
  data.frame(res1)
}, BPPARAM = MulticoreParam(15)))
res1 <- do.call(rbind, res1)
res1$sex="male"


system.time(res2<-bplapply(noms2,function(x){
  res2 <- fastlm.fun(x,yi,dades2,covari,log=FALSE,rank.1=FALSE)
  data.frame(res2)
}, BPPARAM = MulticoreParam(15)))
res2 <- do.call(rbind, res2)
res2$sex="female"

res=data.frame(rbind(res1,res2),stringsAsFactors=F)
res$q.value <- p.adjust(res$p.value, method = "BH")

#export
write.table(res,file=paste(output.path,"/res_cacs_mgs_full_model_sex_strat.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

print(Sys.time())

