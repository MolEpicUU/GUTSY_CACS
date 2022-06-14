rm(list=ls())

#set the working director
set.seed(123)
input.path="/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/github/casc_microbiome/results/"
output.path="/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/github/casc_microbiome/results/8_GMM_enrichment/"
dir.create(output.path, showWarnings = FALSE)
ncores=15

#load libraries
library(rio)
library(BiocParallel)

dir.create(output.path, showWarnings = FALSE)
ncores=5

#import data
res_lm.mod1=import(paste(input.path,"/3_species_level/res_cacs_mgs_basic_model.tsv",sep=""))
load(paste(input.path,"../phenotypes/MGS_HG3A.GMMs2MGS.RData",sep=""))
modules=import(paste(input.path,"../phenotypes/GMM_reference.csv",sep=""))
modules=modules[,1:4]

#load functions
source(paste(input.path,"../scripts/functions/fgsea.GMM.R",sep=""))

#prepare data
#pathway enrichment analyses
MGS_HG3A.GMMs2MGS=MGS_HG3A.GMMs2MGS[lapply(MGS_HG3A.GMMs2MGS,length)>0] ## you can use sapply,rapply
## Compute maximum length
max.length <- max(sapply(MGS_HG3A.GMMs2MGS, length))
## Add NA values to list elements
MGS_HG3A.GMMs2MGS <- lapply(MGS_HG3A.GMMs2MGS, function(v) { c(v, rep(NA, max.length-length(v)))})
## Rbind
MGS_HG3A.GMMs2MGS=data.frame(do.call(rbind, MGS_HG3A.GMMs2MGS),stringsAsFactors = F)

res_lm.mod1.pos=res_lm.mod1[res_lm.mod1$estimate>0,]
res_lm.mod1.neg=res_lm.mod1[res_lm.mod1$estimate<0,]

set.seed(1234)
res1.pos=fgsea.module.fun(res_lm.mod1.pos,"casctot")
res1.pos=res1.pos[order(res1.pos$pval),c(1,8,5,6,7,3,4,11:ncol(res1.pos))]
names(res1.pos)=c("Module","size","log2err","ES","NES","p.value","q.value","Name","HL1","HL2")

set.seed(1234)
res1.neg=fgsea.module.fun(res_lm.mod1.neg,"casctot")
res1.neg=res1.neg[order(res1.neg$pval),c(1,8,5,6,7,3,4,11:ncol(res1.neg))]
names(res1.neg)=c("Module","size","log2err","ES","NES","p.value","q.value","Name","HL1","HL2")


dades.info=import("/proj/sens2019512/SCAPIS_org/SCAPIS/final_release_CMv1/clean_files/HG3.A.7_tax.xlsx")
names(dades.info)[1]="mgs"

res_lm.mod1.pos$mgs=matrix(unlist(
  strsplit(res_lm.mod1.pos$var.x, "____")),ncol=2,byrow=T)[,2]

sig.module=res1.pos[res1.pos$q.value<0.05,"Module"]

gen=names(table(dades.info$genus))[!names(table(dades.info$genus))%in%"unclassified"]


system.time(res.tax<-bplapply(gen,function(i){
  xxx=fgsea.module.fun(res_lm.mod1.pos[res_lm.mod1.pos$mgs%in%dades.info$mgs[!dades.info$genus%in%i],],"casctot")
  xxx=xxx[xxx$pathway%in%sig.module,c(1,8,5,6,7,3,4,11:ncol(xxx))]
  xxx$ex.genus=i
  xxx=xxx[,c(ncol(xxx),1:(ncol(xxx)-1))]
  data.frame(xxx)
}, BPPARAM = MulticoreParam(ncores)))
res.tax <- do.call(rbind, res.tax)
res.tax=data.frame(res.tax,stringsAsFactors = F)

names(res.tax)=c("ex.genus","Module","size","log2err","ES","NES","p.value","q.value","Name","HL1","HL2")

res_lm.mod1.neg$mgs=matrix(unlist(
  strsplit(res_lm.mod1.neg$var.x, "____")),ncol=2,byrow=T)[,2]

sig.module=res1.neg[res1.neg$q.value<0.05,"Module"]

gen=names(table(dades.info$genus))[!names(table(dades.info$genus))%in%"unclassified"]


system.time(res.tax.neg<-bplapply(gen,function(i){
  xxx=fgsea.module.fun(res_lm.mod1.neg[res_lm.mod1.neg$mgs%in%dades.info$mgs[!dades.info$genus%in%i],],"casctot")
  xxx=xxx[xxx$pathway%in%sig.module,c(1,8,5,6,7,3,4,11:ncol(xxx))]
  xxx$ex.genus=i
  xxx=xxx[,c(ncol(xxx),1:(ncol(xxx)-1))]
  data.frame(xxx)
}, BPPARAM = MulticoreParam(ncores)))
res.tax.neg <- do.call(rbind, res.tax.neg)
res.tax.neg=data.frame(res.tax.neg,stringsAsFactors = F)

names(res.tax.neg)=c("ex.genus","Module","size","log2err","ES","NES","p.value","q.value","Name","HL1","HL2")


#export data
write.table(res1.pos,file=paste(output.path,"/enrichment_GMM_pos.tsv",sep=""),
            row.names=F,col.names=T,sep="\t")
write.table(res1.neg,file=paste(output.path,"/enrichment_GMM_neg.tsv",sep=""),
            row.names=F,col.names=T,sep="\t")
write.table(res.tax,file=paste(output.path,"/enrichment_GMM_leave_one_out_pos.tsv",sep=""),
            row.names=F,col.names=T,sep="\t")
write.table(res.tax.neg,file=paste(output.path,"/enrichment_GMM_leave_one_out_neg.tsv",sep=""),
            row.names=F,col.names=T,sep="\t")





