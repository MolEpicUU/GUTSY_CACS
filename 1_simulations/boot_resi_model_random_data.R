rm(list=ls())

#set the working director
set.seed(123)
input.path="./0_Data/"
output.path="./1_simulations/"
dir.create(output.path, showWarnings = FALSE)
ncores=10

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

#load data
dades=import(paste(input.path,"/random_data_covar_1batch.tsv",sep=""))


# regression functions
fastboot.fun <- function(x,y,dades,covari=NULL,log=TRUE,R=1000,meth="residual"){

  covari=paste("+",paste(covari,collapse="+"),sep="")
  if(covari=="+"){
    covari=NULL
  }
  
  tryCatch({
    fit <- eval(parse(text=paste("lm(",y,"~ ",x,covari,", data = dades)",sep=""))
    )
    boot <- Boot(fit, R = R,method=meth)
    # # conf <- confint(boot)
    boot <- as.data.frame(summary(boot))
    coef <- boot[x, ]
    
    
    coef$p.value <- 2 * pnorm(-abs(coef[,"bootMed"] / coef[, "bootSE"]))
    n <-  length(fit$fitted.values)
    data.frame(var.x=x,var.y=y,estimate = coef[,"bootMed"],
               se = coef[, "bootSE"],
               p.value = coef$p.value,
               n = length(fit$fitted.values))
  }, error = function(e) {
    
    data.frame(var.x = x, var.y = y, estimate = NA, se = NA, p.value = NA, n = NA
    )
    
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
print(paste("boots residual",log.op))


print("##### time #####")
system.time(res<-bplapply(noms,function(x){
  res <- fastboot.fun(x,yi,dades,covari=NULL,log=FALSE,R=1000,meth="residual")
  data.frame(res)
}, BPPARAM = MulticoreParam(15)))
res <- do.call(rbind, res)
res$q.value <- p.adjust(res$p.value, method = "BH")


write.table(res, file=paste(output.path,"/res_fastboot_resi.fun_CASC_",log.op,"_1batch.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

print(Sys.time())

