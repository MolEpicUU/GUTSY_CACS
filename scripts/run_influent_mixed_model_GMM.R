library(rio)
library(BiocParallel)
library(lmerTest)
library(lmtest)
library(nnet)

args <- commandArgs(trailingOnly = TRUE)
num <- as.numeric(args[1])


cores <- 15
set.seed(1)


data <- import("processed/data.tsv")
main_sig <- import("results/GMM_main_significant.tsv")

 # main model

 cacstot <- data$cacstot

mgsid=main_sig$id[num]
 GMM <- data[, mgsid]
covari=c("age", "sex", "country", "site_plate", "smoke", "pa", "carb","protein","fiber", "sbp", "dbp", "chol", "hdl", "ldl", "tg", "diab", "bmi", "chol_med", "bp_med", "diab_med")
 cov <- data[, covari]
family<-data$family

print(num)
print(mgsid)

 # linear regression function
infl.fun <- function(y, x, z,re,covari,mgsid) {

  tryCatch({
        covari=paste("+",paste(covari,collapse="+"),sep="")
   y <- log1p(y)
     x <- scale(log1p(x))
     data <- data.frame(y = y, x = x, z,re)
     data <- data[which(complete.cases(data)), ]

fit <- eval(parse(text=paste("lmer(y~x",covari,"+(1|re), data = data)")))
print(mgsid)
inf_gm1 <- influence(fit, ncores = getOption("mc.cores",cores))
beta=dfbetas(inf_gm1)
nb=which.is.max(abs(beta[,2]))
fit2 <- eval(parse(text=paste("lmer(y~x",covari,"+(1|re), data = data[-nb,])")))
    coef <- summary(fit)$coefficients
    lower <- coef[, 1] + qt(0.025, coef[, 3]) * coef[, 2]
    upper <- coef[, 1] + qt(0.975, coef[, 3]) * coef[, 2]

    coef2 <- summary(fit2)$coefficients
    lower2 <- coef2[, 1] + qt(0.025, coef2[, 3]) * coef2[, 2]
    upper2 <- coef2[, 1] + qt(0.975, coef2[, 3]) * coef2[, 2]

    data.frame(mgsid=mgsid,estimate = coef[2, 1], lower = lower[2], upper = upper[2], se = coef[2, 2], p.value = coef[2, 5], n = nrow(data), estimate.infl=coef2[2, 1], lower.infl = lower2[2], upper.infl = upper2[2], se.infl = coef2[2, 2], p.value.infl = coef2[2, 5], n.infl = nrow(data[-nb,]),rm.val=nb,converg.infl=table(inf_gm1$converged),message = NA)

  }, warning = function(w) {

    data.frame(mgsid=mgsid,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA,estimate.infl = NA, lower.infl = NA, upper.infl = NA, se.infl = NA, p.value.infl = NA, n.infl = NA,converg.infl=NA,rm.val=NA,message = paste("Warning:", w$message))

  }, error = function(e) {

    data.frame(mgsid=mgsid,estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, estimate.infl = NA, lower.infl = NA, upper.infl = NA, se.infl = NA, p.value.infl = NA, n.infl = NA,converg.infl=NA,rm.val=NA,message = paste("Error:", e$message))

  })

}

infl=infl.fun(cacstot,GMM,cov,family,covari,mgsid)

export(infl, paste("results/GMM_influential_",num,".tsv",sep=""))

