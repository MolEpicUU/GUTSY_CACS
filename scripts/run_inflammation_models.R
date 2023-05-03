# load libraries and set options

library(rio)
library(BiocParallel)
library(lmerTest)

cores <- 16
set.seed(1)

# import data

data <- import("processed/data.tsv")
info <- import("raw/HG3.A.7_tax.xlsx")
colnames(info)[1:3] <- c("id", "name", "level")
main_sig <- import("results/main_significant.tsv")
id <- main_sig$id
id <- id[which(id %in% c("HG3A.1967","HG3A.0270","HG3A.1800")==F)]
name <- info$name[match(id, info$id)]
level <- info$level[match(id, info$id)]

# main model

mgs <- data[, match(id, colnames(data))]
cov <- data[, c("age", "sex", "country", "site_plate", "smoke", "pa", "carb", "protein", "fiber", "sbp", "dbp", "chol", "hdl", "ldl", "tg", "diab", "bmi", "chol_med", "bp_med", "diab_med")]
family<-data$family

# linear regression function

lmer.fun <- function(y, x, z,re) {
  
  tryCatch({
    
    x <- scale(log1p(x))
    data <- data.frame(y = y, x = x, z,re)
    data <- data[which(complete.cases(data)), ]
    fit <- lmer(y ~ .-re+(1|re), data)
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

# crp

crp <- data$crp
crp <- bplapply(colnames(mgs), function(x) lmer.fun(crp, mgs[, x], cov,family), BPPARAM = MulticoreParam(cores))
crp <- do.call(rbind, crp)
print(paste("main model warnings and errors:", sum(!is.na(crp$message))))
crp$q.value <- p.adjust(crp$p.value, method = "BH", n = sum(!is.na(crp$p.value)))
crp <- crp[, c("n", "estimate", "se", "p.value", "q.value")]
colnames(crp) <- paste0("crp_", colnames(crp))


# neut

neut <- data$neut
neut <- bplapply(colnames(mgs), function(x) lmer.fun(neut, mgs[, x], cov,family), BPPARAM = MulticoreParam(cores))
neut <- do.call(rbind, neut)
print(paste("main model warnings and errors:", sum(!is.na(neut$message))))
neut$q.value <- p.adjust(neut$p.value, method = "BH", n = sum(!is.na(neut$p.value)))
neut <- neut[, c("n", "estimate", "se", "p.value", "q.value")]
colnames(neut) <- paste0("neut_", colnames(neut))

# leuk

leuk <- data$leuk
leuk <- bplapply(colnames(mgs), function(x) lmer.fun(leuk, mgs[, x], cov,family), BPPARAM = MulticoreParam(cores))
leuk <- do.call(rbind, leuk)
print(paste("main model warnings and errors:", sum(!is.na(leuk$message))))
leuk$q.value <- p.adjust(leuk$p.value, method = "BH", n = sum(!is.na(leuk$p.value)))
leuk <- leuk[, c("n", "estimate", "se", "p.value", "q.value")]
colnames(leuk) <- paste0("leuk_", colnames(leuk))

res <- data.frame(name, id, level, crp, neut, leuk)
export(res, "results/inflammation.tsv")

sessionInfo()
