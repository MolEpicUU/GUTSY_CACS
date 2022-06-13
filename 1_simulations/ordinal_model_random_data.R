rm(list=ls())

print(Sys.time())

#set the working director
set.seed(123)
input.path="/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/github/casc_microbiome/Demo/0_Data/"
output.path="/proj/nobackup/sens2019512/wharf/ssayols/ssayols-sens2019512/github/casc_microbiome/Demo/1_simulations/"
dir.create(output.path, showWarnings = FALSE)
ncores=5


##load libraries
library(rio)
library(BiocParallel)
library(ordinal)

# load data
dades=import(paste(input.path,"/random_data_covar_1batch.tsv",sep=""))

# regression functions
ordinal.fun <- function(y, x,data) {
  
  tryCatch({
    fit <- eval(parse(text=paste("clm(as.factor(",y,") ~ ",x,", data = data)",sep="")))
    coef <- summary(fit)$coefficients
    # ci <- confint(fit)
    if (is.finite(coef[x, 2])) {
      
      data.frame(var.x=x,var.y=y,estimate = coef[x, 1],
                 # lower=ci["x",1],
                 # upper=ci["x",2],
                 se=coef[x,2],p.value = coef[x, 4],n=length(fit$fitted.values),aic=AIC(fit))
      
    } else {
      
      data.frame(var.x=x,var.y=y,estimate = NA, 
                 # lower=NA,upper=NA,
                 se=NA,p.value = NA,n=NA,aic=NA)
    }
    
  }, error = function(e) {
    
    data.frame(var.x=x,var.y=y,estimate = NA, 
               # lower=NA,upper=NA,
               se=NA,p.value = NA,n=NA,aic=NA)
    
  })
  
}

log.op="log1p"
rank.op=""

#variables to assess
noms=grep("X.",names(dades),value=T)

#outcome
yi="casctot"

print("#### model ####")
print(paste("lm model",log.op,rank.op))

#remove individuals with NAs in the outcome
dades=dades[which(is.na(dades[,yi])==F),]

#log1p transform 
dades[,noms]=apply(dades[,noms], 2, log1p)

#run the models
print("##### time #####")
system.time(res.ordinal<-bplapply(noms,function(x){
  res.ordinal<-ordinal.fun(yi,x,dades)
  data.frame(res.ordinal)
}, BPPARAM = MulticoreParam(ncores)))
res.ordinal <- do.call(rbind, res.ordinal)

write.table(res.ordinal,file=paste(output.path,"/res_ordinal.fun_CASC_",log.op,rank.op,"_1batch.tsv",sep=""),
            col.names = T,row.names = F,sep="\t")

print(Sys.time())

