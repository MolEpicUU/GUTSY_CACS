
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
cov <- data[, c("age", "sex", "country", "site_plate")]
family <- data$family

basic <- lapply(colnames(alpha), function(x) lmer.fun(cacstot, alpha[, x], cov,family))
basic <- do.call(rbind, basic)
colnames(basic) <- paste0("basic_", colnames(basic))

# main model

cov <- data[, c( "smoke", "pa", "carb", "protein", "fiber", "sbp", "dbp", "chol", "hdl", "ldl", "tg", "diab", "bmi", "chol_med", "bp_med", "diab_med")]


main <- lapply(colnames(alpha), function(x) lmer.fun(cacstot, alpha[, x], cov,family))
main <- do.call(rbind, main)
colnames(main) <- paste0("main_", colnames(main))

# make and export tables

res <- data.frame(alpha = colnames(alpha), basic, main)
export(res, "results/alpha.tsv")

sessionInfo()

