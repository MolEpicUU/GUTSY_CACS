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
library(boot)

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
fastboot.fun <- function(x,y,dades,covari=NULL,log=TRUE,R=1000) {
  covari=paste("+",paste(covari,collapse="+"),sep="")
  if(covari=="+"){
    covari=NULL
  }
  
  tryCatch({
    coef_se.boot.fun<-function(dades,index){
      xdades=dades[index,c(y,x)]
      model2<-eval(parse(text=paste("glm.nb(",y,"~.,data=xdades,maxit=100)",sep="")))
      c(coef(model2)[x], summary(model2)$coefficient[x,"Std. Error"])
    }
    
    boot <- boot(dades,coef_se.boot.fun, R = R)
    # # conf <- confint(boot)
    boot <- as.data.frame(summary(boot))
    rownames(boot)=c("Estimate","SE")
    #    coef$p.value <- 2 * pt(-abs(coef[,"bootMed"] / coef[, "bootSE"]), fit$df)
    boot$p.value <- 2 * pnorm(-abs(boot["Estimate","bootMed"] / boot["SE", "bootMed"]))
    n <-  nrow(dades)
    data.frame(var.x=x,var.y=y,
               estimate = boot["Estimate","bootMed"],
               se = boot["SE", "bootMed"],
               p.value = boot$p.value[1],
               n = nrow(dades))
  }, error = function(e) {
    
    data.frame(var.x = x, var.y = y, estimate = NA, se = NA, p.value = NA, n = NA)
    
  })
}

log.op="log1p"
noms=grep("X.",names(dades),value=T)
yi="casctot"

dades[,noms]=apply(dades[,noms], 2, log1p)

print("#### model ####")
print(paste("fastboot.fun.nb",log.op))

print("##### time #####")
system.time(res<-bplapply(noms,function(x){
  res <- fastboot.fun(x,yi,dades,covari=NULL,log=FALSE,R=1000)
  data.frame(res)
}, BPPARAM = MulticoreParam(ncores)))
res <- do.call(rbind, res)


#export data
write.table(res,file=paste(output.path,"/res_nb.boot.",log.op,"_1batch.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

print(Sys.time())


