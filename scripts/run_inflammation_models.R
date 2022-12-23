# load libraries and set options

library(rio)
library(BiocParallel)

cores <- 16
set.seed(1)

# import data

data <- import("processed/data.tsv")
info <- import("raw/HG3.A.7_tax.xlsx")
main_sig <- import("results/main_significant.tsv")

colnames(info)[1:3] <- c("id", "name", "level")
id <- main_sig$id
id <- id[which(!id == "HG3A.1967")]
name <- info$name[match(id, info$id)]
level <- info$level[match(id, info$id)]

# main model

mgs <- data[, match(id, colnames(data))]
cov <- data[, c("age", "sex", "country", "site_plate", "smoke", "pa", "carb", "protein", "fiber", "sbp", "dbp", "chol", "hdl", "ldl", "tg", "diab", "bmi", "chol_med", "bp_med", "diab_med")]

# linear regression function

lm.fun <- function(y, x, z) {
  
  tryCatch({
    
    x <- scale(log1p(x))
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

# crp

crp <- data$crp
crp <- bplapply(colnames(mgs), function(x) lm.fun(crp, mgs[, x], cov), BPPARAM = MulticoreParam(cores))
crp <- do.call(rbind, crp)
print(paste("main model warnings and errors:", sum(!is.na(crp$message))))
crp$q.value <- p.adjust(crp$p.value, method = "BH", n = sum(!is.na(crp$p.value)))
crp <- crp[, c("n", "estimate", "se", "p.value", "q.value")]
colnames(crp) <- paste0("crp_", colnames(crp))


# neut

neut <- data$neut
neut <- bplapply(colnames(mgs), function(x) lm.fun(neut, mgs[, x], cov), BPPARAM = MulticoreParam(cores))
neut <- do.call(rbind, neut)
print(paste("main model warnings and errors:", sum(!is.na(neut$message))))
neut$q.value <- p.adjust(neut$p.value, method = "BH", n = sum(!is.na(neut$p.value)))
neut <- neut[, c("n", "estimate", "se", "p.value", "q.value")]
colnames(neut) <- paste0("neut_", colnames(neut))

# leuk

leuk <- data$leuk
leuk <- bplapply(colnames(mgs), function(x) lm.fun(leuk, mgs[, x], cov), BPPARAM = MulticoreParam(cores))
leuk <- do.call(rbind, leuk)
print(paste("main model warnings and errors:", sum(!is.na(leuk$message))))
leuk$q.value <- p.adjust(leuk$p.value, method = "BH", n = sum(!is.na(leuk$p.value)))
leuk <- leuk[, c("n", "estimate", "se", "p.value", "q.value")]
colnames(leuk) <- paste0("leuk_", colnames(leuk))

res <- data.frame(name, id, level, crp, neut, leuk)
export(res, "results/inflammation.tsv")

sessionInfo()