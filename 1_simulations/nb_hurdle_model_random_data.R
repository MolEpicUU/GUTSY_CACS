rm(list=ls())

print(Sys.time())

#set the working director
set.seed(123)
input.path="./0_Data/"
output.path="./1_simulations/"
dir.create(output.path, showWarnings = FALSE)
ncores=5
set.seed(123)

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
dades=import(paste(input.path,"/random_data_covar_1batch.tsv",sep=""))

# regression functions

hurdle.nb.fun<-function(x,y,dades,covari=NULL){
  xdades=dades[,c(x,y,covari)]
  covari=paste("+",paste(covari,collapse="+"),sep="")
  if(covari=="+"){
    covari=NULL
  }
  
  tryCatch({
    hurdle.nb <- eval(parse(text=paste("hurdle(",y,"~",x,covari,",data=xdades, dist = 'negbin',zero.dist = 'binomial')")))
    
    data.frame(var.x=x,var.y=y,
               estimate.zero = summary(hurdle.nb)$coeff$zero[x,"Estimate"],
               se.zero = summary(hurdle.nb)$coeff$zero[x,"Std. Error"],
               p.value.zero = summary(hurdle.nb)$coeff$zero[x,"Pr(>|z|)"],
               estimate.count = summary(hurdle.nb)$coeff$count[x,"Estimate"],
               se.count = summary(hurdle.nb)$coeff$count[x,"Std. Error"],
               p.value.count = summary(hurdle.nb)$coeff$count[x,"Pr(>|z|)"],
               n = length(hurdle.nb$fitted.values),aic=AIC(hurdle.nb))
    
  }, error = function(e) {
    
    data.frame(var.x=x,var.y=y,estimate.zero = NA,se.zero = NA,
               p.value.zero = NA, estimate.count=NA,se.count=NA,
               p.value.count=NA,n = NA,aic=NA)
    
    
  })
}


log.op="log1p"

#noms=c(grep("MGS",names(dades),value=T),grep("^K[0-9]",names(dades),value=T))
noms=grep("X.",names(dades),value=T)

yi="casctot"

dades[,noms]=apply(dades[,noms],2,log1p)


print("#### model ####")
print(paste("hurdle.nb.fun",log.op))


print("##### time #####")
system.time(res<-bplapply(noms,function(x){
  res <- hurdle.nb.fun(x,yi,dades,covari=NULL)
  data.frame(res)
}, BPPARAM = MulticoreParam(ncores)))
res <- do.call(rbind, res)


write.table(res,file=paste(output.path,"res_hurdle.nb.fun_",log.op,"_1batch.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

print(Sys.time())


