
# load libraries

library(rio)
library(lmerTest)

# import data

data <- import("processed/data.tsv")

# linear mixed regression function

lmer.fun <- function(y, x, z, re) {
  
  tryCatch({
    
    y <- log1p(y)
    data <- data.frame(y = y, x = x, z, re)
    data <- data[which(complete.cases(data)), ]
    fit <- lmer(y ~ . - re + (1 | re), data)
    coef <- summary(fit)$coefficients
    lower <- coef[, 1] + qt(0.025, coef[, 3]) * coef[, 2]
    upper <- coef[, 1] + qt(0.975, coef[, 3]) * coef[, 2]
    data.frame(estimate = coef[2, 1], lower = lower[2], upper = upper[2], se = coef[2, 2], p.value = coef[2, 5], n = nrow(data), message = NA)
    
  }, warning = function(w) {
    
    data.frame(estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Warning:", w$message))
    
  }, error = function(e) {
    
    data.frame(estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Error:", e$message))
    
  })
  
}

# basic model

cacstot <- data$cacstot
alpha <- data[, c("shannon", "simpson", "invsimpson", "chao")]
xcov<-c("age", "sex", "country", "site_plate")
cov <- data[, xcov]
family <- data$family

basic <- lapply(colnames(alpha), function(x) lmer.fun(cacstot, alpha[, x], cov,family))
basic <- do.call(rbind, basic)
res <- data.frame(alpha =  colnames(alpha),model="basic", basic)

# main model

xcov2 <-c( "smoke", "pa", "carb", "protein", "fiber", "sbp", "dbp", "chol", "hdl", "ldl", "tg", "diab", "bmi", "chol_med", "bp_med", "diab_med")


for (i in xcov2)
{
  cov2=c(xcov,i)
  
  cov=data[,cov2]
  
  xxx <- lapply(colnames(alpha), function(x) lmer.fun(cacstot, alpha[, x], cov,family))
  xxx <- do.call(rbind, xxx)
  res <- rbind(res,data.frame(alpha = colnames(alpha),model=paste0("basic+",i),xxx))
}
 
# make and export tables
  
res <- data.frame(res,stringsAsFactors=F)
export(res, "results/alpha_stepwise.tsv")
  
sessionInfo()
  
