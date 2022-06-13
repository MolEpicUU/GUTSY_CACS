rm(list=ls())

#set the working director
set.seed(123)
input.path="./0_Data/"
output.path="./2_diversity/"
dir.create(output.path, showWarnings = FALSE)
ncores=16
permutations<-9999

#load libraries
library(rio)
library(vegan)
library(BiocParallel)

#import data
dades=import(paste(input.path,"pheno_2021_06_21.tsv",sep=""))
beta=import(paste(input.path,"bray_curtis.tsv",sep=""))

#prepare data
dades=dades[,!names(dades)%in%grep("HG3A",names(dades),value=T)]

rownames(beta)=colnames(beta)

beta=beta[rownames(beta)%in%dades$sample.id,names(beta)%in%dades$sample.id]
dades=dades[dades$sample.id%in%rownames(beta),]

dades=dades[ order(match(dades$sample.id, rownames(beta))), ]

print("identical1")
print(identical(dades$sample.id,rownames(beta)))

dades$casctot_cat=ifelse(dades$casctot==0,0,
                         ifelse(dades$casctot>0 & dades$casctot<=100,1,
                                ifelse(dades$casctot>100 & dades$casctot<=400,2,3
                                )))

dades$gender=as.factor(dades$gender)
dades$siteid=as.factor(dades$siteid)
dades$q005a=as.factor(dades$q005a)
dades$plate=as.factor(dades$plate)
dades$casctot_cat=as.factor(dades$casctot_cat)
dades$smokestatus=as.factor(dades$smokestatus)
dades$q134=as.factor(dades$q134)
dades$diab_treat=as.factor(dades$diab_treat)
dades$HBP_treat=as.factor(dades$HBP_treat)
dades$HC_treat=as.factor(dades$HC_treat)

var=c("sample.id","agev1","log.fibrer","log.energi","gender","siteid","q005a","plate","casctot_cat","diab_treat","HBP_treat","HC_treat","q134","smokestatus")
dades=dades[,var]
dades=dades[complete.cases(dades),]
beta=beta[dades$sample.id,dades$sample.id]

print("identical2")
print(identical(dades$sample.id,rownames(beta)))

print("adonis mod2")

#pairwise
dades.fun=dades[dades$casctot_cat==0| dades$casctot_cat==1,]
beta.fun<-beta[rownames(beta)%in%dades.fun$sample.id,names(beta)%in%dades.fun$sample.id]

print("identical3")
print(identical(dades$sample.id,rownames(beta)))

set.seed(12345) 
fit=adonis2(as.matrix(beta.fun)~casctot_cat+gender+agev1+siteid+siteid:plate+q005a+smokestatus+q134+diab_treat+HC_treat+HBP_treat+log.fibrer+log.energi,data=dades.fun,by="margin",parallel=ncores,permutations=permutations)
res=data.frame(var1="0",var2="1",R2=fit$R2[1],p.value=fit$"Pr(>F)"[1])

write.table(res,file=paste(output.path,"permanova_adonis2_dsbray_curtis_casctot_pairwise.mod2_0_1.tsv",sep=""),row.names=F,col.names=T,sep="\t")

print(Sys.time())



