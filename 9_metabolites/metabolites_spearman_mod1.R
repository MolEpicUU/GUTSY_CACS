rm(list=ls())

#set the working director
set.seed(123)
input.path="/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/github/casc_microbiome/phenotypes/"
output.path="/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/github/casc_microbiome/results/9_metabolites/"
dir.create(output.path, showWarnings = FALSE)
ncores=10

#load libraries
library(rio)
library(BiocParallel)

#import data
dades=import(paste(input.path,"pheno_2021_09_06_metab.tsv",sep=""))
res_lm.mod1=import(paste(output.path,"../3_species_level/res_cacs_mgs_basic_model.tsv",sep=""))
met_info=import(paste(input.path,"scapis_merged_annotations_batchnorm_clean.tsv",sep=""))

#load functions
source("/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/github/casc_microbiome/scripts/functions/spearman_correlation.R")

#prepare data
sel1=grep("HG3A.",res_lm.mod1[res_lm.mod1$q.value<0.05,"var.x"],value=T)
sel1=grep("Streptococc",sel1,value=T)


dades$gender=as.factor(dades$gender)
dades$siteid=as.factor(dades$siteid)
dades$q005a=as.factor(dades$q005a)
dades$plate=as.factor(dades$plate)

met_info=met_info[,c("MET_ID","SUPER_PATHWAY","SUB_PATHWAY","CHEMICAL_NAME")]
names(met_info)[1]="var.y"


res1=spearman.cor.fun(sel1,dades,
		      var=c("agev1","gender","siteid","q005a","plate","shannon____mgs"),
                      covari1=c("agev1","gender","siteid","siteid:plate","q005a","shannon____mgs"),
                      cores=ncores)

res1=merge(res1,met_info,by="var.y")
#export data
write.table(res1,file=paste(output.path,"res_spearman_corr_metab_mod1.tsv",sep=""),col.names=T,row.names=F,sep="\t")

