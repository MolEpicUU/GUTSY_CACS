rm(list=ls())

#set the working director
set.seed(123)
input.path="./0_Data/"
output.path="./2_diversity/"
dir.create(output.path, showWarnings = FALSE)
ncores=10
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

dades$casctot_cat=ifelse(dades$casctot==0,0,
                         ifelse(dades$casctot>0 & dades$casctot<=100,1,
                                ifelse(dades$casctot>100 & dades$casctot<=400,2,3
                                )))

dades$gender=as.factor(dades$gender)
dades$siteid=as.factor(dades$siteid)
dades$q005a=as.factor(dades$q005a)
dades$plate=as.factor(dades$plate)
dades$casctot_cat=as.factor(dades$casctot_cat)


print("adonis mod1")

#pairwise
dades.fun=dades[dades$casctot_cat==1| dades$casctot_cat==2,]
beta.fun<-beta[rownames(beta)%in%dades.fun$sample.id,names(beta)%in%dades.fun$sample.id]

print("identical")
print(identical(rownames(beta),rownames(dades.fun$sample.id)))

set.seed(12345) 
fit=adonis2(as.matrix(beta.fun)~casctot_cat+gender+agev1+siteid+siteid:plate+q005a,data=dades.fun,by="margin",parallel=ncores,permutations=permutations)
res=data.frame(var1="1",var2="2",R2=fit$R2[1],p.value=fit$"Pr(>F)"[1])

write.table(res,file=paste(output.path,"permanova_adonis2_dsbray_curtis_casctot_pairwise.mod1_1_2.tsv",sep=""),row.names=F,col.names=T,sep="\t")

print(Sys.time())



