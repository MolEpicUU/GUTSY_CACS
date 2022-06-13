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
library(caret)
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

# load data
dades=import(paste(input.path,"/random_data_covar_1batch.tsv",sep=""))

# regression functions

nb.fun<-function(x,y,dades,covari=NULL,log=TRUE){
  xdades=dades[,c(x,y,covari)]
  covari=paste("+",paste(covari,collapse="+"),sep="")
  if(covari=="+"){
    covari=NULL
  }
  if(log==TRUE){
    xdades[,x]<-log(xdades[,x]+1)
  }
  
  tryCatch({
    mod.nb <- eval(parse(text=paste("glm.nb(",y,"~",x,covari,",data=xdades)")))
    fit=coeftest(mod.nb)
    data.frame(var.x=x,var.y=y,estimate = fit[x,"Estimate"],
               se = fit[x,"Std. Error"],
               p.value = fit[x,"Pr(>|z|)"],
               n = length(mod.nb$fitted.values),aic=AIC(mod.nb))
    
  }, error = function(e) {
    
    
    data.frame(var.x=x,var.y=y,estimate = NA,
               se = NA,
               p.value = NA, n = NA,aic=NA)
    
    
  })
}


noms2=grep("X.",names(dades),value=T)

dades[,noms2]=apply(dades[,noms2],2,log1p)


log.op="log1p"

yi="casctot"

print("#### model ####")
print(paste("nb.fun",log.op))

print("##### time #####")
system.time(res<-bplapply(noms2,function(x){
  res <- nb.fun(x,yi,dades,covari=NULL,log=FALSE)
  data.frame(res)
}, BPPARAM = MulticoreParam(ncores)))
res <- do.call(rbind, res)

#export data
write.table(res,file=paste(output.path,"/nb_",log.op,"_1batch.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

print(Sys.time())


