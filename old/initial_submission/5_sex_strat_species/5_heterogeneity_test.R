rm(list=ls())

#set the working director
set.seed(123)
input.path="./"
output.path="./5_sex_strat_species/"
dir.create(output.path, showWarnings = FALSE)
ncores=5

#load libraries
library(rio)
library(BiocParallel)

#import data
dades=import(paste(input.path,"/0_Data/pheno_2021_06_21.tsv",sep=""))
res.sex1=import(paste(input.path,"5_sex_strat_species/res_cacs_mgs_basic_model_sex_strat.tsv",sep=""))
res.sex2=import(paste(input.path,"5_sex_strat_species/res_cacs_mgs_full_model_sex_strat.tsv",sep=""))
res_lm.mod1=import(paste(input.path,"3_species_level/res_cacs_mgs_basic_model.tsv",sep=""))
res_lm.mod2=import(paste(input.path,"3_species_level/res_cacs_mgs_full_model.tsv",sep=""))

#load functions
source(paste(input.path,"/functions/prevalence.fun2.R",sep=""))
source(paste(input.path,"/functions/beta.test.R",sep=""))

#prepare data
dades=dades[,c("casctot","agev1","gender","siteid",
               "q005a","plate",
               "smokestatus",
               "q134","diab_treat","HBP_treat","HC_treat",
               "log.fibrer","log.energi",
               grep("____mgs",names(dades),value=T),
               grep("HG3A.",names(dades),value=T))]

yi="casctot"
covari=c("agev1","siteid","q005a","plate",
         "shannon____mgs")

dades$gender=as.factor(dades$gender)
dades$siteid=as.factor(dades$siteid)
dades$q005a=as.factor(dades$q005a)
dades$plate=as.factor(dades$plate)
dades$smokestatus=as.factor(dades$smokestatus)
dades$q134=as.factor(dades$q134)
dades$diab_treat=as.factor(dades$diab_treat)
dades$HBP_treat=as.factor(dades$HBP_treat)
dades$HC_treat=as.factor(dades$HC_treat)

#select mgs
sel1=res_lm.mod1[which(res_lm.mod1$q.value<0.05),"var.x"]
sel2=res_lm.mod2[which(res_lm.mod2$q.value<0.05),"var.x"]


#Run the model1
system.time(res1<-bplapply(sel1,function(x){
  res1<- beta.test(dades,x,yi,covari,res.sex1)
  res1
}, BPPARAM = MulticoreParam(ncores)))
res1 <- do.call(rbind, res1)
res1=data.frame(res1,stringsAsFactors = F)
res1$q.value=p.adjust(res1$p.value,"BH")

covari=c("agev1","siteid","q005a","smokestatus",
         "q134","diab_treat","HBP_treat","HC_treat","plate","log.fibrer","log.energi","shannon____mgs")

#Run the model2
system.time(res2<-bplapply(sel2,function(x){
  res <- beta.test(dades,x,yi,covari,res.sex2)
  res
}, BPPARAM = MulticoreParam(10)))
res2 <- do.call(rbind, res2)
res2=data.frame(res2,stringsAsFactors = F)

res2$q.value=p.adjust(res2$p.value,"BH")

#export data
write.table(res1,file=paste(output.path,"heterogeneity_test_basic_model.tsv",sep=""),col.names=T,row.names=F,sep="\t")
write.table(res2,file=paste(output.path,"heterogeneity_test_full_model.tsv",sep=""),col.names=T,row.names=F,sep="\t")

