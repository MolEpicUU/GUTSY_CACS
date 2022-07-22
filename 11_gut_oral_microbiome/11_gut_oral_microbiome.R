rm(list=ls())

#set the working director
set.seed(123)
input.path="./0_Data/"
output.path="./11_gut_oral_microbiome/"
dir.create(output.path, showWarnings = FALSE)
ncores=1

#load libraries
library(rio)
library(BiocParallel)
library(fastDummies)
library(ppcor)
#import files
mos=import(paste(input.path,"/MOS.csv",sep=""))

#load function
source("./functions/spearman_clust.se.fun.R")

#preapare data

mos=mos[!mos$PPI_All%in%1,]

mos=mos[,c(grep("Streptoco",names(mos),value=T),"country_birth","Sex_k","Plate.x","Age","Family_ID")]
mos=mos[complete.cases(mos),]

covari=names(mos)[!names(mos)%in%c(grep("Strep",names(mos),value=T),"lopnrMOS")]
chr.covari=c(grep("country",covari,value=T),grep("Sex",covari,value=T),grep("Plate",covari,value=T))
names(mos)=gsub(" ",".",names(mos))

oral.strep=grep("Ho1B",names(mos),value=T)
gut.strep=grep("HG3A",names(mos),value=T)

mos$Sex_k=as.factor(mos$Sex_k)
mos$country_birth=as.factor(mos$country_birth)
mos$Plate.x=as.factor(mos$Plate.x)

mos[,!names(mos)%in%chr.covari]=apply(mos[,!names(mos)%in%chr.covari],2,rank)

res<-bplapply(oral.strep,function(x){
  x.s=matrix(unlist(strsplit(x," ")),ncol=2,byrow=T)[,2]
  x.s=substring(x.s,25,nchar(x.s))
  
  xx=sp.clust.se(mos,x,grep(x.s,gut.strep,value=T),covari,"Family_ID")
  
  data.frame(var.x=x,var.y=grep(x.s,gut.strep,value=T),
             coef=xx[1],
             p.value=xx[2],
             n=xx[3])
  
  
}, BPPARAM = MulticoreParam(ncores))

res <- do.call(rbind, res)
res$q.value=p.adjust(res$p.value,"BH")

# res$hg3a.oral=substring(res$var.x,1,9)
# res$maintax.oral=gsub("."," ",substring(res$var.x,11,nchar(res$var.x)),fixed=T)
# res$hg3a.gut=substring(res$var.y,1,9)
# res$maintax.gut=gsub("."," ",substring(res$var.y,11,nchar(res$var.x)),fixed=T)
# res=res[,c(10,9,8,7,5,3,4,6)]

write.table(res,file=paste(output.path,"/oral_gut_microbiome.tsv",sep=""),
            sep="\t",col.names=T,row.names=F)

