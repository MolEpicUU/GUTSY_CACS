# load libraries and set options

library(rio)
library(BiocParallel)
library(lmerTest)
library(ggpubr)

cores <- 16

# import and clean data

data <- import("processed/data.tsv")
info <- import("raw/HG3.A.7_tax.xlsx")
colnames(info)[1:3] <- c("id", "name", "level")

# linear regression function

lmer.fun <- function(y, x, z, re) {
 
  tryCatch({
   
    y <- log1p(y)
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

# basic model

cacstot <- data$cacstot
mgs <- data[, which(grepl("HG3A", colnames(data)))]
cov <- data[, c("age", "sex", "country", "site_plate")]
family<-data$family

basic <- bplapply(colnames(mgs), function(x) lmer.fun(cacstot, mgs[, x], cov,family), BPPARAM = MulticoreParam(cores))
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

main <- bplapply(colnames(mgs), function(x) lmer.fun(cacstot, mgs[, x], cov,family), BPPARAM = MulticoreParam(cores))
main <- do.call(rbind, main)
main$id <- colnames(mgs)
main$name <- info$name[match(main$id, info$id)]
main$level <- info$level[match(main$id, info$id)]
main$prev <- apply(mgs[which(complete.cases(cov)), ], 2, function(x) sum(x > 0) / length(x))
main$q.value <- p.adjust(main$p.value, method = "BH", n = sum(!is.na(main$p.value)))
main <- main[, c("name", "id", "level", "prev", "n", "estimate", "lower", "upper", "se", "p.value", "q.value", "message")]

# sex interaction

interaction.fun <- function(y, x, z,re) {
 
  tryCatch({
   
    y <- log1p(y)
    x <- scale(log1p(x))
    data <- data.frame(y = y, x = x, z,re)
    data <- data[which(complete.cases(data)), ]
    fit <- lmer(paste0("y ~ ", paste("sex *", colnames(z), collapse = " + "), " - sex * sex +(1|re)+ sex * x"), data)
    coef <- summary(fit)$coefficients
    lower <- coef[, 1] + qt(0.025, coef[, 3]) * coef[, 2]
    upper <- coef[, 1] + qt(0.975, coef[, 3]) * coef[, 2]

    data.frame(estimate = coef[nrow(coef), 1], lower = lower[nrow(coef)], upper = upper[nrow(coef)], se = coef[nrow(coef), 2], p.value = coef[nrow(coef), 5], n = nrow(data), message = NA)
   
  }, warning = function(w) {
   
    data.frame(estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Warning:", w$message))
   
  }, error = function(e) {
   
    data.frame(estimate = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Error:", e$message))
   
  })
 
}

interaction <- bplapply(colnames(mgs), function(x) interaction.fun(cacstot, mgs[, x], cov,family), BPPARAM = MulticoreParam(cores))
interaction <- do.call(rbind, interaction)
interaction$q.value <- p.adjust(interaction$p.value, method = "BH", n = sum(!is.na(interaction$p.value)))
interaction <- interaction[, c("estimate", "p.value", "q.value", "message")]
colnames(interaction) <- paste0("interaction_", colnames(interaction))

# make and export tables

main <- data.frame(main, interaction)
main <- main[order(main$p.value), ]
export(main, "results/main.tsv")

main_sig <- main[which(main$q.value < 0.05), ]
export(main_sig, "results/main_significant.tsv")

sessionInfo()

