rm(list=ls())
#set the working director

output.folder="./0_Data"
dir.create(output.folder, showWarnings = FALSE)

set.seed(123)

#load libraries
library(vegan)
library(scales)
library(FRGEpistasis)
library(xlsx)

#simulate CACS distribution
mu=60
var=44200
n=100

k=(mu^2)/(var - mu)

#simulate species distribution
y <- rnbinom(n = n, mu = mu, size =k)

mu=0.18
var=1.4
k=(mu^2)/(var - mu)

x=as.data.frame(matrix(abs(jitter(rnbinom(n*50,mu=mu,size=k),1)), ncol = 50))

rho=seq(0.01,0.1,by=0.01)

z=NULL
for(i in 1:10)
{
  z=cbind(z,abs(rho[i]*abs(jitter(y,1)) + sqrt(1-(rho[i]^2))*x[,i]))
}

x[,1:10]=z
names(x)[1:9]=paste("species____HG3A.0",1:9,sep="")
names(x)[2:7]=paste("Streptococcus.",names(x)[2:7],sep="")
names(x)[10:ncol(x)]=paste("species____HG3A.",10:ncol(x),sep="")

x$sample.id=paste("sample.gut.",1:nrow(x),sep="")
x=x[,c(ncol(x),1:(ncol(x)-1))]

mgs=x

hscrp_res=rankTransPheno(abs(0.9*log1p(abs(jitter(y,1))) + sqrt(1-(0.9^2))*rnorm(n,2.3,4.1)),para_c=3/8)
neut_res=rankTransPheno(abs(0.5*log1p(abs(jitter(y,1))) + sqrt(1-(0.5^2))*rnorm(n,3.1,1.1)),para_c=3/8)
lpk_res=rankTransPheno(abs(0.5*log1p(abs(jitter(y,1))) + sqrt(1-(0.5^2))*rnorm(n,5.6,2.3)),para_c=3/8)

#simulate the other traits
pheno=data.frame(sample.id=paste("sample.gut.",1:nrow(x),sep=""),
                 casctot=y,
                 agev1=rnorm(n,57.4,4.3),
                 gender = sample(c("1", "2"), n, replace = T),
                 q005a = sample(c("country1", "country2", "country3"), n, replace = T),
                 siteid = sample(c("site1", "site2"), n, replace = T),
                 smokestatus=sample(c("never","former","current"),n,replace=T),
                 q134=sample(c("sedentary","moderate","regular-regular","regular"),n,replace = T),
                 diab_treat=sample(c("0","1"),n,replace = T),
                 HBP_treat=sample(c("0","1"),n,replace = T),
                 HC_treat=sample(c("0","1"),n,replace = T),
                 plate = sample(c("plate1", "plate2", "plate3"), n, replace = T),
                 log.fibrer=rnorm(n,2.9,0.5),
                 log.energi=rnorm(n,7.4,0.4),
                 carotidplaque=sample(c("none","unilateral","bilateral"),n,replace=T),
                 
                 hscrp_res=hscrp_res+abs(round(min(hscrp_res)))+1,
                 neut_res=neut_res+abs(round(min(neut_res)))+1,
                 lpk_res=lpk_res+abs(round(min(lpk_res)))+1,
                 bmi=rnorm(n,27.1,4.5),
                 ppi=sample(c("0","1"),n,replace=T),
                 stringsAsFactors = F)

#estimate shannon diversity index
pheno$shannon____mgs=diversity(x[,grep("species.",names(x),value=T)],index="shannon")

#merge mgs
pheno=merge(pheno,mgs,by="sample.id")
mgs=mgs[sample(n,n*40/100),]
rownames(mgs)=mgs$sample.id
mgs=data.frame(t(mgs[,-1]),stringsAsFactors = F)
mgs$mgs=rownames(mgs)

#estimate bray-curtis matrix

BC=data.frame(as.matrix(
  vegdist(pheno[,grep("species",names(pheno),
                               value=T)],method="bray")),
  stringsAsFactors=F)

rownames(BC)=names(BC)=pheno$sample.id

#simulate metabolomic data
met=pheno[,grep("species",names(pheno),value=T)]+rnorm(n*30)
names(met)=paste("metabolite..MET.",1:ncol(met),sep="")

met <- apply(met, 2, function(x) {
  
  sample <- sample(1:2, 1)
  scale <- list(rescale(x, c(0, 5)), 
                rescale(x, c(5, 0)))
  scale[[sample]]
  
})
met=data.frame(met,stringsAsFactors = F)
met$sample.id=paste("sample.gut.",1:nrow(met),sep="")

met=merge(pheno,met,by="sample.id")

set.seed(888)
met=met[sample(1:nrow(met),n-100),]

#simulate genus groups
HG3A=matrix(unlist(strsplit(rownames(mgs),"____")),ncol=2,byrow=T)[,2]
MainTax=matrix(unlist(strsplit(rownames(mgs),"____")),ncol=2,byrow=T)[,1]
genus=c(rep("genus1",10),
  rep("genus2",15),
  rep("genus3",25))
                 
dades.info=data.frame(cbind(HG3A,MainTax,genus),stringsAsFactors=F)

#simulate functional modules
 module1=c(dades.info$HG3A[1],dades.info$HG3A[50],
 sample(dades.info$HG3A[which(dades.info$HG3A%in%c(dades.info$HG3A[1],dades.info$HG3A[50],"HG3A.07")==F)],8))
 module2=c(sample(dades.info$HG3A[2:10],6),sample(dades.info$HG3A[which(dades.info$HG3A%in%c(module1,"HG3A.07")==F)],3),"HG3A.14")
module3=sample(dades.info$HG3A[which(dades.info$HG3A%in%c(module1,module2,"HG3A.07")==F)],20)
module4=dades.info$HG3A[which(dades.info$HG3A%in%c(module1,module2,module3)==F)]
 
MGS_HG3A.GMMs2MGS <- list(module1 = module1,
               module2 = module2,
               module3 = module3,
               module4=module4)
               
modu=data.frame(Module=c("module1","module2","module3","module4"),
Name=c("name1","name2","name3","name4"),
HL1=c("HL1_1","HL1_2","HL1_3","HL1_4"),
HL2=c("HL2_1","HL2_2","HL2_3","HL2_4"),stringsAsFactors=F)


#simulate metabolomic subpathways

subpathway1=sample(names(met)[2:11],5)
subpathway2=names(met)[which(names(met)[1:11]%in%c("sample.id",subpathway1)==F)]
subpathway3= sample(names(met)[12:50],20)
subpathway4= names(met)[which(names(met)%in%c(subpathway1,subpathway2,subpathway3,"sample.id")==F)]

subpathway <- list(subpathway1 = subpathway1,
               subpathway2 = subpathway2,
               subpathway3 = subpathway3,
               subpathway4=subpathway4)

#simulate mgs on the validation data

x=as.data.frame(matrix(abs(jitter(rnbinom(n*50,mu=mu,size=k))), ncol = 50))
mgs.v=x

for(i in 1:4)
{
  mgs.v[,i]=abs(rho[i]*abs(jitter(pheno$casctot,1)) + sqrt(1-(rho[i]^2))*x[,i])
}

mgs.v=mgs.v[,c(1:5,10:50)]
names(mgs.v)=c("HG3A.01","HG3A.02","HG3A.03","HG3A.04","HG3A.06","HG3A.07","HG3A.08","HG3A.09",paste("HG3A.",10:47,sep=""))
rownames(mgs.v)=mgs.v$sample.id=paste("sample.v.",1:nrow(mgs.v),sep="")
mgs.v=mgs.v[,c(ncol(mgs.v),1:(ncol(mgs.v)-1))]

#simulate the other traits
pheno.v=data.frame(subject_id=paste("sample.v.",1:nrow(mgs.v),sep=""),
                 age=rnorm(n,61,10),
                 gender = sample(c("2", "1"), n, replace = T),
                 disease=ifelse(pheno$casctot==0,"healthy","ACVD"),
                 shannon____mgs=diversity(mgs.v[,grep("HG3A.",names(mgs.v),value=T)],index="shannon"),
                 stringsAsFactors = F)


mgs.g=data.frame(apply(mgs[,grep("species",names(mgs),value=T)],2,function(x)  jitter(x)),stringsAsFactors = F)

#simulate the oral mgs
rho=seq(0.003,0.0075,by=0.001)
x=as.data.frame(matrix(abs(jitter(rnbinom(n*50,mu=jitter(mu),size=jitter(k)),1)), ncol = 50))
mgs.o=x
for(i in 1:4)
{
  mgs.o[,i]=abs(rho[i]*abs(jitter(y,1)) + sqrt(1-(rho[i]^2))*x[,i])
}
names(mgs.o)=c("species____HG3A.1","species____HG3A.2","species____HG3A.3","species____HG3A.4",paste("species____HG3A.",6:47,sep=""))
mgs.g=mgs.o$sample.id=paste("sample.o.",1:nrow(mgs.o),sep="")
mgs.o=mgs.o[,c(ncol(mgs.o),1:(ncol(mgs.o)-1))]

#simulate the other traits
pheno.o=data.frame(sample.id=paste("sample.o.",1:nrow(mgs.o),sep=""),
                   agev1=rnorm(n,53,8),
                   country_birth = sample(c("country1", "country2"), n, replace = T),
                   gender = sample(c("1", "2"), n, replace = T),
                   shannon____mgs.oral=diversity(mgs.o[,grep("species",names(mgs.o),value=T)],index="shannon"),
                   shannon____mgs.gut=diversity(mgs[,grep("species",names(mgs),value=T)],index="shannon"),
                   plate = sample(c("plate1", "plate2", "plate3"), n, replace = T),
                   antibiotics=sample(c("0","1"),n,replace = T),
                   education=sample(c("primary","secundary","university"),n,replace=T),
                   eat=sample(c("0","1"),n,replace = T),
                   tabacco=sample(c("0","1"),n,replace = T),
                   brush=sample(c("0","1"),n,replace = T),
                   plaque_score=rbeta(n,1,2)*1.5,
                   stringsAsFactors = F)

#simulate the oral phenotypes
#caries
mu=2.5
var=32.5
k=(mu^2)/(var-mu)
x=rnbinom(1:n,mu=mu,size=k)
rho=0.001
pheno.o$caries=abs(rho*abs(jitter(y,1)) + sqrt(1-(rho^2))*x)


#FS
mu=31.1
var=461.4
k=(mu^2)/(var-mu)
rho=0.003
x=rnbinom(1:n,mu=mu,size=k)
pheno.o$FS=abs(rho*abs(jitter(y,1)) + sqrt(1-(rho^2))*x)

#BoP
rho=0.8
x=rbeta(n,2,5)*100
pheno.o$BoP=rankTransPheno(abs(rho*abs(jitter(y,1)) + sqrt(1-(rho^2))*rnorm(n,24.6,17.2)),para_c=3/8)

#export data

write.table(BC,file=paste(output.folder,"/bray_curtis.tsv",sep=""),row.names=F,col.names=T,sep="\t")
write.table(pheno,file=paste(output.folder,"/pheno_2021_06_21.tsv",sep=""),row.names=F,col.names=T,sep="\t")
write.table(pheno.v,file=paste(output.folder,"/validation_cohort.tsv",sep=""),row.names=F,col.names=T,sep="\t")
write.table(pheno.o,file=paste(output.folder,"/pheno_oral.tsv",sep=""),row.names=F,col.names=T,sep="\t")
write.table(mgs,file=paste(output.folder,"/upugut_batch1.mgsComp_tax.txt",sep=""),row.names=F,col.names=T,sep="\t")
saveRDS(mgs.v,file=paste(output.folder,"/ERP023788_chinese.MGS.comp.rds",sep=""))
write.table(mgs.g,file=paste(output.folder,"/mgs_gut.tsv",sep=""),row.names=F,col.names=T,sep="\t")
write.table(mgs.g,file=paste(output.folder,"/mgs_oral.tsv",sep=""),row.names=F,col.names=T,sep="\t")
write.table(met,file=paste(output.folder,"/metabolites.tsv",sep=""),row.names=F,col.names=T,sep="\t")
write.xlsx(dades.info,file=paste(output.folder,"/HG3.A.7_tax.xlsx",sep=""),row.names=F,col.names=T)
save(MGS_HG3A.GMMs2MGS,file=paste(output.folder,"/MGS_HG3A.GMMs2MGS.RData",sep=""))
write.table(modu,file=paste(output.folder,"/GMM_reference.csv",sep=""),row.names=F,col.names=T,sep=",")
save(subpathway,file=paste(output.folder,"/metabolic_subpathways.RData",sep=""))

