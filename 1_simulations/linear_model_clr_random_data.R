rm(list=ls())

print(Sys.time())

#set the working director
set.seed(123)
input.path="./0_Data/"
output.path="./1_simulations/"
dir.create(output.path, showWarnings = FALSE)
ncores=5

##load libraries
library(rio)
library(BiocParallel)
library(microbiome)

#load data
dades=rio::import(paste(input.path,"/random_data_clr_covar_1batch.tsv",sep=""))

# regression functions
fastlm.fun <- function(x,y,dades,covari=NULL,log=FALSE,rank.1=FALSE,ties.meth="keep") {
  tryCatch({
    xdades=dades[,c(x,y,covari)]
    covari=paste("+",paste(covari,collapse="+"),sep="")
    if(covari=="+"){
      covari=NULL
    }
    if(log==TRUE){
      xdades[,x]<-log(xdades[,x]+1)
      xdades[,y]<-log(xdades[,y]+1)
    }
    if(rank.1==TRUE){
      xdades[,x] <- rank(xdades[,x], ties.meth)
      xdades[,y] <- rank(xdades[,y], ties.meth)
      
    }
    
    fit <- eval(parse(text=paste("lm(",y,"~",x,covari,", data = xdades)")))
    coef<-summary(fit)$coefficients[x,]
    lower= coef["Estimate"] - qt(0.975, df = fit$df) * coef["Std. Error"]
    upper= coef["Estimate"] + qt(0.975, df = fit$df) * coef["Std. Error"]
    data.frame(var.x=x,var.y=y, estimate = coef["Estimate"],
               lower=lower,
               upper=upper,
               se = coef["Std. Error"],
               p.value = coef["Pr(>|t|)"],
               n = length(fit$fitted.values),
               aic=AIC(fit))
    
  }, error = function(e) {
    
    data.frame(var.x=x,var.y=y, estimate = NA,
               lower=NA,upper=NA,se = NA,
               p.value = NA,
               n = NA,aic=NA)
  })
}

log.op="clr"
rank.op=""

#variables to assess
noms=grep("X.",names(dades),value=T)

#outcome
yi="casctot"


print("#### model ####")
print(paste("lm model",log.op,rank.op))

dades$casctot=log1p(dades$casctot)

#run the models
print("##### time #####")
system.time(res.lm<-bplapply(noms,function(x){
  res.lm <- fastlm.fun(x,yi,dades,covari=NULL,log=FALSE,rank.1=FALSE)
  data.frame(res.lm)
}, BPPARAM = MulticoreParam(ncores)))
res.lm <- do.call(rbind, res.lm)

write.table(res.lm,file=paste(output.path,"/res_fastlm.fun_CASC_",log.op,rank.op,"_1batch.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

print(Sys.time())

