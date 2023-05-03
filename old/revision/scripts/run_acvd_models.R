# load libraries and set options

library(rio)
library(BiocParallel)
library(MASS)

cores <- 16
set.seed(1)

# import and clean data

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

# logistic regression function

logist.fun <- function(y, x, z) {
  
  tryCatch({
    
    x <- scale(log1p(x))
    data <- data.frame(y = y, x = x, z)
    data <- data[which(complete.cases(data)), ]
    fit <- glm(y ~ ., data, family = "binomial")
    coef <- summary(fit)$coefficients
    lower <- coef[, 1] + qnorm(0.025) * coef[, 2]
    upper <- coef[, 1] + qnorm(0.975) * coef[, 2]
    data.frame(or = exp(coef[2, 1]), lower = exp(lower[2]), upper = exp(upper[2]), se = coef[2, 2], p.value = coef[2, 4], n = nrow(data), message = NA)
    
  }, warning = function(w) {
    
    data.frame(or = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Warning:", w$message))
    
  }, error = function(e) {
    
    data.frame(or = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Error:", e$message))
    
  })
  
}

# athero

athero <- as.numeric(factor(data$athero, levels = c("No", "Yes"))) - 1
athero <- bplapply(colnames(mgs), function(x) logist.fun(athero, mgs[, x], cov), BPPARAM = MulticoreParam(cores))
athero <- do.call(rbind, athero)
athero$q.value <- p.adjust(athero$p.value, method = "BH", n = sum(!is.na(athero$p.value)))
athero <- athero[, c("n", "or", "lower", "upper", "p.value", "q.value", "message")]
colnames(athero) <- paste0("athero_", colnames(athero))

# stenosis

stenosis <- as.numeric(factor(data$stenosis, levels = c("No", "Yes"))) - 1
stenosis <- bplapply(colnames(mgs), function(x) logist.fun(stenosis, mgs[, x], cov), BPPARAM = MulticoreParam(cores))
stenosis <- do.call(rbind, stenosis)
stenosis$q.value <- p.adjust(stenosis$p.value, method = "BH", n = sum(!is.na(stenosis$p.value)))
stenosis <- stenosis[, c("n", "or", "lower", "upper", "p.value", "q.value", "message")]
colnames(stenosis) <- paste0("stenosis_", colnames(stenosis))

# ordinal regression function

ordinal.fun <- function(y, x, z) {
  
  tryCatch({
    
    x <- scale(log1p(x))
    data <- data.frame(y = y, x = x, z)
    data <- data[which(complete.cases(data)), ]
    fit <- polr(ordered(y) ~ ., data = data, Hess = T)
    coef <- summary(fit)$coefficients
    lower <- coef[, 1] + qnorm(0.025) * coef[, 2]
    upper <- coef[, 1] + qnorm(0.975) * coef[, 2]
    p <- 2 * pnorm(-abs(coef[, 3]))
    data.frame(or = exp(coef[1, 1]), lower = exp(lower[1]), upper = exp(upper[1]), se = coef[1, 2], p.value = p[1], n = nrow(data), message = NA)
    
  }, warning = function(w) {
    
    data.frame(or = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Warning:", w$message))
    
  }, error = function(e) {
    
    data.frame(or = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Error:", e$message))
    
  })
  
}

# duke

duke <- as.numeric(gsub("DUKE_", "", data$duke))
duke[which(duke == 5)] <- 4
duke[which(duke == 6)] <- 4
duke <- bplapply(colnames(mgs), function(x) ordinal.fun(duke, mgs[, x], cov), BPPARAM = MulticoreParam(cores))
duke <- do.call(rbind, duke)
duke$q.value <- p.adjust(duke$p.value, method = "BH", n = sum(!is.na(duke$p.value)))
duke <- duke[, c("n", "or", "lower", "upper", "p.value", "q.value", "message")]
colnames(duke) <- paste0("duke_", colnames(duke))

# sis

sis <- as.numeric(gsub("SIS_", "", data$sis))
sis <- bplapply(colnames(mgs), function(x) ordinal.fun(sis, mgs[, x], cov), BPPARAM = MulticoreParam(cores))
sis <- do.call(rbind, sis)
sis$q.value <- p.adjust(sis$p.value, method = "BH", n = sum(!is.na(sis$p.value)))
sis <- sis[, c("n", "or", "lower", "upper", "p.value", "q.value", "message")]
colnames(sis) <- paste0("sis_", colnames(sis))

# plaque

plaque <- as.numeric(factor(data$plaque, levels = c("None", "Unilateral", "Bilateral")))
plaque <- bplapply(colnames(mgs), function(x) ordinal.fun(plaque, mgs[, x], cov), BPPARAM = MulticoreParam(cores))
plaque <- do.call(rbind, plaque)
plaque$q.value <- p.adjust(plaque$p.value, method = "BH", n = sum(!is.na(plaque$p.value)))
plaque <- plaque[, c("n", "or", "lower", "upper", "p.value", "q.value", "message")]
colnames(plaque) <- paste0("plaque_", colnames(plaque))

res <- data.frame(name, id, level, athero, stenosis, duke, sis, plaque)
export(res, "results/acvd.tsv")

sessionInfo()