rm(list=ls())

#set the working director
set.seed(123)
input.path="./0_Data/"
output.path="./7_validation/"
dir.create(output.path, showWarnings = FALSE)
ncores=1


#load libraries
library(curatedMetagenomicData)
library(rio)
library(ggplot2)
library(ggpubr)
library(reshape)
library(mosaic)
library(FRGEpistasis)
library(vegan)
library(viridis)

#import data
pheno.v=import(paste(input.path,"validation_cohort.tsv",sep=""))
dades.v=data.frame(readRDS(paste(input.path,"ERP023788_chinese.MGS.comp.rds",sep="")),stringsAsFactors = F)
res=import(paste(output.path,"../3_species_level/res_cacs_mgs_basic_model.tsv",sep=""))

#prepare data
res$hg3a=matrix(unlist(strsplit(res$var.x,"____")),ncol=2,byrow=TRUE)[,2]
sel=grep("HG3A.",res[res$q.value<0.05,"hg3a"],value=T)

dades.v=dades.v[,names(dades.v)%in%sel]
dades.v$subject_id=rownames(dades.v)
data=merge(pheno.v,dades.v,by="subject_id")

data$disease=factor(data$disease, levels = c("healthy", "ACVD"))
data$gender=as.factor(data$gender)


noms2=grep("HG3A.",names(data),value=T)
exc=names(which(colSums(data[,grep("HG3A",names(data),value=T)])==0))

noms2=noms2[!noms2%in%exc]

#estimate shannon diversity index
data$shannon=diversity(data[,grep("HG3A.",names(data),value=T)],index="shannon")

data=data[!is.na(data$gender)==TRUE,]
data=data[!is.na(data$age)==TRUE,]

#run the models
xdata=data
res.val=NULL
for (i in 1:length(noms2))
{
  # rank-Based Inverse Normal Transformation
  data[,noms2[i]]=rankTransPheno(data[, noms2[i]], para_c=3/8)
  eval(parse(text=paste("fit=glm(disease~",noms2[i],"+age+gender+shannon,family='binomial',data=data)",sep="")))
  
  coef=summary(fit)$coefficient[noms2[i],]
  ci=confint(fit)[noms2[i],]
  
  a= try(if(table(data[data$disease%in%"ACVD",noms2[i]]!=0)["FALSE"]!=nrow(data[data$disease%in%"ACVD",]))
  {
    xres=data.frame(var.x=noms2[i],
                    prevalence.cases=ifelse(table(data[data$disease%in%"ACVD",noms2[i]]==0)<length(data[data$disease%in%"ACVD",noms2[i]]),
                                            table(data[data$disease%in%"ACVD",noms2[i]]!=0)["TRUE"]/sum(table(data[data$disease%in%"ACVD",noms2[i]]==0))*100,0),
                    prevalence.control=table(data[data$disease%in%"healthy",noms2[i]]!=0)["TRUE"]/sum(table(data[data$disease%in%"healthy",noms2[i]]==0))*100,
                    estimate=coef[1],
                    lower=ci[1],upper=ci[2],
                    se=coef[2],
                    n=length(fit$fitted.values),
                    p.value=coef[4],message=NA,estimate_cacs=res[res$hg3a%in%grep(noms2[i],res$hg3a,value=T),"estimate"])
  }else{
    xres=data.frame(var.x=noms2[i],
                    prevalence.cases=ifelse(table(data[data$disease%in%"ACVD",noms2[i]]==0)<length(data[data$disease%in%"ACVD",noms2[i]]),
                                            table(data[data$disease%in%"ACVD",noms2[i]]!=0)["TRUE"]/sum(table(data[data$disease%in%"ACVD",noms2[i]]==0))*100,0),
                    prevalence.control=table(data[data$disease%in%"healthy",noms2[i]]!=0)["TRUE"]/sum(table(data[data$disease%in%"healthy",noms2[i]]==0))*100,
                    estimate=NA,
                    lower=NA,upper=NA,
                    se=NA,
                    n=NA,
                    p.value=NA,message=NA,estimate_cacs=res[res$hg3a%in%grep(noms2[i],res$hg3a,value=T),"estimate"])
  })
  # }
  
  if(class(a)=="try-error") {
    xres=data.frame(var.x=noms2[i],
                    prevalence.cases=
                      ifelse(table(data[data$disease%in%"ACVD",noms2[i]]!=0)["TRUE"]==nrow(data[data$disease%in%"ACVD",]),100,
                             ifelse(table(data[data$disease%in%"ACVD",noms2[i]]==0)["TRUE"]<length(data[data$disease%in%"ACVD",noms2[i]],
                                                                                                   table(data[data$disease%in%"ACVD",noms2[i]]!=0)["TRUE"]/sum(table(data[data$disease%in%"ACVD",noms2[i]]==0))*100,0))),
                    prevalence.control=table(data[data$disease%in%"healthy",noms2[i]]!=0)["TRUE"]/sum(table(data[data$disease%in%"healthy",noms2[i]]==0))*100,
                    estimate=coef[1],
                    lower=ci[1],upper=ci[2],
                    se=coef[2],
                    n=length(fit$fitted.values),
                    p.value=coef[4],message=NA,estimate_cacs=res[res$hg3a%in%grep(noms2[i],res$hg3a,value=T),"estimate"])
  }
  
  res.val=rbind(res.val,xres)
  
  
}
res.val=unique(res.val)
res.val$q.value=p.adjust(res.val$p.value,"BH",n=length(na.omit(res.val$p.value)))
res.val$OR=exp(res.val$estimate)
res.val$lower.OR=exp(res.val$lower)
res.val$upper.OR=exp(res.val$upper)

#export results
write.table(res.val,file=paste(output.path,"/res_validation.tsv",sep="")
            ,row.names=F,col.names=T,sep="\t")
