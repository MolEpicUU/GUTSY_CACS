rm(list=ls())

print(Sys.time())

#set the working director
set.seed(123)
input.path="./0_Data/"
output.path="./1_simulations/"
dir.create(output.path, showWarnings = FALSE)
ncores=5

##load libraries
library(compareGroups)
library(pscl)
require(MASS)
library(ggplot2)

library(VennDiagram)


library(rio)
library(car)
library(caret)
library(robust)
library(ordinal)
library(censReg)
library(BiocParallel)
library(reshape2)
library(ggplot2)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(sandwich)
library(lmtest)
library(pscl)

# load data
dades=import("/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/final_results_cacs_microbiota/phenotype/manuscriptV1/random_data_covar_1batch.tsv")


# regression functions

two.stage.log.lm.fun<-function(x,y,dades,covari=NULL){
  xdades=dades[,c(x,y,covari)]
  xdades[,paste(y,"_cat",sep="")]=ifelse(xdades[,y]==0,0,1)
  xdades2=xdades[xdades[,y]!=0,]
  covari=paste("+",paste(covari,collapse="+"),sep="")
  if(covari=="+"){
    covari=NULL
  }
  
  tryCatch({
    log.fit <- eval(parse(text=paste("glm(",y,"_cat~",x,covari,",data=xdades,family='binomial')",sep='')))
    lm.fit=eval(parse(text=paste("lm(",y,"~",x,covari,",data=xdades2)")))
    data.frame(var.x=x,var.y=y,
               estimate.zero = summary(log.fit)$coeff[x,"Estimate"],
               se.zero = summary(log.fit)$coeff[x,"Std. Error"],
               p.value.zero = summary(log.fit)$coeff[x,"Pr(>|z|)"],
               n.zero=length(log.fit$fitted.values),
               aic.zero=AIC(log.fit),
               estimate.count = summary(lm.fit)$coeff[x,"Estimate"],
               se.count = summary(lm.fit)$coeff[x,"Std. Error"],
               p.value.count = summary(lm.fit)$coeff[x,"Pr(>|t|)"],
               n.count = length(lm.fit$fitted.values),
               aic.count=AIC(lm.fit))
  }, error = function(e) {
    
    data.frame(var.x=x,var.y=y,estimate.zero = NA,se.zero = NA,
               p.value.zero = NA, n.zero=NA,aic.zero=NA,estimate.count=NA,se.count=NA,
               p.value.count=NA,n.count = NA,aic.count=NA)
    
    
  })
}


log.op="log1p"

noms=grep("X.",names(dades),value=T)
yi="casctot"
dades[,c(yi,noms)]=apply(dades[,c(yi,noms)],2,log1p)


print("#### model ####")
print(paste("two.stage.log.lm",log.op))

print("##### time #####")
system.time(res<-bplapply(noms,function(x){
  res <- two.stage.log.lm.fun(x,yi,dades,covari=NULL)
  data.frame(res)
}, BPPARAM = MulticoreParam(ncores)))
res <- do.call(rbind, res)


write.table(res,file=paste(output.path,"res_two.stage.log.lm.fun_",log.op,"_1batch.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

print(Sys.time())

