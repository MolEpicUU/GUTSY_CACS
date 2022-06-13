rm(list=ls())

#set the working director
set.seed(123)
input.path="./0_Data/"
output.path="./1_simulations/"
dir.create(output.path, showWarnings = FALSE)
ncores=5

#load libraries
require(MASS)
library(boot)
library(rio)
library(car)
library(caret)
library(robust)
library(BiocParallel)
library(sandwich)
library(lmtest)

# load data
dades=import(paste(input.path,"/random_data_covar_1batch.tsv",sep=""))

# regression functions
fastrob.fun <- function(x,y,dades,covari=NULL,log=T) {
  xdades=dades[,c(x,y,covari)]
  covari=paste("+",paste(covari,collapse="+"),sep="")
  if(covari=="+"){
    covari=NULL
  }
  if(log==T){
    xdades[,x]<-log(xdades[,x]+1)
    xdades[,y]<-log(xdades[,y]+1)
  }
  tryCatch({
    fit <- eval(parse(text=paste("lm(",y,"~",x,covari,",data=xdades)")))
    coef <- coeftest(fit, df=Inf,vcovHC(fit,type="HC0"))
    # conf <- confint(coef) # ci
    coef <- as.matrix(coef)
    coef <- coef[grep(x,rownames(coef)), ]
    data.frame(var.x=x,var.y=y, estimate = coef["Estimate"],
               se = coef["Std. Error"],
               p.value = coef["Pr(>|z|)"],
               n = length(fit$fitted.values))
  }, error = function(e) {
    
    data.frame(var.x = NA, var.y = NA, estimate = NA, se = NA, p.value = NA, n =
                 NA)
    
  })
  
}

log.op="log1p"

#variables to assess
noms=grep("X.",names(dades),value=T)

#outcome
yi="casctot"

#remove individuals with NAs in the outcome
dades=dades[which(is.na(dades[,yi])==F),]

#log1p transform 
dades[,c(yi,noms)]=apply(dades[,c(yi,noms)], 2, log1p)

print("#### model ####")
print(paste("fastrob.fun",log.op))


print("##### time #####")
system.time(res<-bplapply(noms,function(x){
  res <- fastrob.fun(x,yi,dades,covari=NULL,log=FALSE)
  data.frame(res)
}, BPPARAM = MulticoreParam(5)))
res <- do.call(rbind, res)



write.table(res,file=paste(output.path,"/res_fastrobust.fun_CASC_",log.op,"_1batch.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

print(Sys.time())

