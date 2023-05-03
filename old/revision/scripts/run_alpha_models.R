# load libraries

library(rio)

# import data

data <- import("processed/data.tsv")

# linear regression function

lm.fun <- function(y, x, z) {
  
  tryCatch({
    
    y <- log1p(y)
    data <- data.frame(y = y, x = x, z)
    data <- data[which(complete.cases(data)), ]
    fit <- lm(y ~ ., data)
    coef <- summary(fit)$coefficients
    ci <- confint(fit)
    data.frame(estimate = coef[2, 1], lower = ci[2, 1], upper = ci[2, 2], se = coef[2, 2], p.value = coef[2, 4], n = nrow(data), message = NA)
    
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

basic <- lapply(colnames(alpha), function(x) lm.fun(cacstot, alpha[, x], cov))
basic <- do.call(rbind, basic)
colnames(basic) <- paste0("basic_", colnames(basic))

# main model

cov <- data[, c("age", "sex", "country", "site_plate", "smoke", "pa", "carb", "protein", "fiber", "sbp", "dbp", "chol", "hdl", "ldl", "tg", "diab", "bmi", "chol_med", "bp_med", "diab_med")]

main <- lapply(colnames(alpha), function(x) lm.fun(cacstot, alpha[, x], cov))
main <- do.call(rbind, main)
colnames(main) <- paste0("main_", colnames(main))

# make and export tables

res <- data.frame(alpha = colnames(alpha), basic, main)
export(res, "results/alpha.tsv")

sessionInfo()
