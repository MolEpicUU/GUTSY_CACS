# load libraries and set options

library(rio)
library(BiocParallel)

cores <- 16

# import and clean data

data <- import("processed/data.tsv")
info <- import("raw/HG3.A.7_tax.xlsx")
colnames(info)[1:3] <- c("id", "name", "level")

# linear regression function

lm.fun <- function(y, x, z) {
  
  tryCatch({
    
    y <- log1p(y)
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

# basic model

cacstot <- data$cacstot
mgs <- data[, which(grepl("HG3A", colnames(data)))]
cov <- data[, c("age", "sex", "country", "site_plate")]

basic <- bplapply(colnames(mgs), function(x) lm.fun(cacstot, mgs[, x], cov), BPPARAM = MulticoreParam(cores))
basic <- do.call(rbind, basic)
basic$id <- colnames(mgs)
basic$name <- info$name[match(basic$id, info$id)]
basic$level <- info$level[match(basic$id, info$id)]
basic$prev <- apply(mgs[which(complete.cases(cov)), ], 2, function(x) sum(x > 0) / length(x))
basic$q.value <- p.adjust(basic$p.value, method = "BH", n = sum(!is.na(basic$p.value)))

# make and export tables

basic <- merge(basic, info[, c("id", "subspecies", "species", "genus", "family", "order", "class", "phylum", "superkingdom")], by = "id")
basic <- basic[, c("name", "id", "level", "prev", "n", "estimate", "lower", "upper", "se", "p.value", "q.value", "subspecies", "species", "genus", "family", "order", "class", "phylum", "superkingdom")]
basic <- basic[order(basic$p.value), ]
export(basic, "results/basic.tsv")

basic_sig <- basic[which(basic$q.value < 0.05), ]
export(basic_sig, "results/basic_significant.tsv")

# main model

mgs <- data[, which(colnames(data) %in% basic_sig$id)]
cov <- data[, c("age", "sex", "country", "site_plate", "smoke", "pa", "carb", "protein", "fiber", "sbp", "dbp", "chol", "hdl", "ldl", "tg", "diab", "bmi", "chol_med", "bp_med", "diab_med")]

main <- bplapply(colnames(mgs), function(x) lm.fun(cacstot, mgs[, x], cov), BPPARAM = MulticoreParam(cores))
main <- do.call(rbind, main)
main$id <- colnames(mgs)
main$name <- info$name[match(main$id, info$id)]
main$level <- info$level[match(main$id, info$id)]
main$prev <- apply(mgs[which(complete.cases(cov)), ], 2, function(x) sum(x > 0) / length(x))
main$q.value <- p.adjust(main$p.value, method = "BH", n = sum(!is.na(main$p.value)))
main <- main[, c("name", "id", "level", "prev", "n", "estimate", "lower", "upper", "se", "p.value", "q.value", "message")]

# influential observations

infl.fun <- function(y, x, z) {
  
  tryCatch({
    
    y <- log1p(y)
    x <- scale(log1p(x))
    data <- data.frame(y = y, x = x, z)
    data <- data[which(complete.cases(data)), ]
    fit <- lm(y ~ ., data)
    inf <- influence.measures(fit)
    dfbeta <- ifelse(sum(inf$is.inf[, 2]) > 0, paste(inf$infmat[, 2][which(inf$is.inf[, 2])], collapse = ","), NA)
    fit <- lm(y ~ ., data[which(!inf$is.inf[, 2]), ])
    coef <- summary(fit)$coefficients
    ci <- confint(fit)
    data.frame(dfbeta = dfbeta, estimate = coef[2, 1], lower = ci[2, 1], upper = ci[2, 2], se = coef[2, 2], p.value = coef[2, 4], n = length(fit$residuals), message = NA)
    
  }, warning = function(w) {
    
    data.frame(dfbeta = dfbeta, estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Warning:", w$message))
    
  }, error = function(e) {
    
    data.frame(dfbeta = dfbeta, estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Error:", e$message))
    
  })
  
}

infl <- bplapply(colnames(mgs), function(x) infl.fun(cacstot, mgs[, x], cov), BPPARAM = MulticoreParam(cores))
infl <- do.call(rbind, infl)
infl$q.value <- p.adjust(infl$p.value, method = "BH", n = sum(!is.na(infl$p.value)))
infl <- infl[, c("n", "dfbeta", "estimate", "lower", "upper", "se", "p.value", "q.value", "message")]
colnames(infl) <- paste0("influential_", colnames(infl))

# sex interaction

interaction.fun <- function(y, x, z) {
  
  tryCatch({
    
    y <- log1p(y)
    x <- scale(log1p(x))
    data <- data.frame(y = y, x = x, z)
    data <- data[which(complete.cases(data)), ]
    fit <- lm(paste0("y ~ ", paste("sex *", colnames(z), collapse = " + "), " - sex * sex + sex * x"), data)
    coef <- summary(fit)$coefficients
    ci <- confint(fit)
    data.frame(estimate = coef[nrow(coef), 1], lower = ci[nrow(ci), 1], upper = ci[nrow(ci), 2], se = coef[nrow(coef), 2], p.value = coef[nrow(coef), 4], n = nrow(data), message = NA)
    
  }, warning = function(w) {
    
    data.frame(estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Warning:", w$message))
    
  }, error = function(e) {
    
    data.frame(estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Error:", e$message))
    
  })
  
}

interaction <- bplapply(colnames(mgs), function(x) interaction.fun(cacstot, mgs[, x], cov), BPPARAM = MulticoreParam(cores))
interaction <- do.call(rbind, interaction)
interaction$q.value <- p.adjust(interaction$p.value, method = "BH", n = sum(!is.na(interaction$p.value)))
interaction <- interaction[, c("estimate", "p.value", "q.value", "message")]
colnames(interaction) <- paste0("interaction_", colnames(interaction))

# make and export tables

main <- data.frame(main, infl, interaction)
main <- main[order(main$p.value), ]
export(main, "results/main.tsv")

main_sig <- main[which(main$q.value < 0.05), ]
export(main_sig, "results/main_significant.tsv")

sessionInfo()
