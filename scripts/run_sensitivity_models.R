# load libraries and set options

library(rio)
library(BiocParallel)
library(lmerTest)
library(lmtest)
library(ordinal)

cores <- 16
set.seed(1)

# import and clean data

data <- import("processed/data.tsv")
info <- import("raw/HG3.A.7_tax.xlsx")
main_sig <- import("results/main_significant.tsv")

colnames(info)[1:3] <- c("id", "name", "level")
id <- main_sig$id

id <- id[which(main_sig$id %in% c("HG3A.1967","HG3A.1800","HG3A.0270")==F) ]
name <- info$name[match(id, info$id)]
level <- info$level[match(id, info$id)]

# main model

cacstot <- data$cacstot
mgs <- data[, match(id, colnames(data))]
cov <- data[, c("age", "sex", "country", "site_plate", "smoke", "pa", "carb", "protein", "fiber", "sbp", "dbp", "chol", "hdl", "ldl", "tg", "diab", "bmi", "chol_med", "bp_med", "diab_med")]
family<-data$family

# linear regression function

lmer.fun <- function(y, x, z,re) {
  
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

# ppi

ppi <- bplapply(colnames(mgs), function(x) lmer.fun(cacstot[which(data$ppi == "No")], mgs[which(data$ppi == "No"), x], cov[which(data$ppi == "No"), ],family[which(data$ppi=="No")]), BPPARAM = MulticoreParam(cores))
ppi <- do.call(rbind, ppi)
ppi$q.value <- p.adjust(ppi$p.value, method = "BH", n = sum(!is.na(ppi$p.value)))
ppi <- ppi[, c("n", "estimate", "se", "p.value", "q.value")]
colnames(ppi) <- paste0("ppi_", colnames(ppi)) 

# crohn

crohn <- bplapply(colnames(mgs), function(x) lmer.fun(cacstot[which(data$crohn == "No")], mgs[which(data$crohn == "No"), x], cov[which(data$crohn == "No"), ],family[which(data$crohn=="No")]), BPPARAM = MulticoreParam(cores))
crohn <- do.call(rbind, crohn)
crohn$q.value <- p.adjust(crohn$p.value, method = "BH", n = sum(!is.na(crohn$p.value)))
crohn <- crohn[, c("n", "estimate", "se", "p.value", "q.value")]
colnames(crohn) <- paste0("crohn_", colnames(crohn)) 

# antibiotics

antib <- ifelse(data$narrow == "Yes" | data$broad == "Yes", "Yes", "No")
antib <- bplapply(colnames(mgs), function(x) lmer.fun(cacstot[which(antib == "No")], mgs[which(antib == "No"), x], cov[which(antib == "No"), ],family[which(antib=="No")]), BPPARAM = MulticoreParam(cores))
antib <- do.call(rbind, antib)
antib$q.value <- p.adjust(antib$p.value, method = "BH", n = sum(!is.na(antib$p.value)))
antib <- antib[, c("n", "estimate", "se", "p.value", "q.value")]
colnames(antib) <- paste0("antib_", colnames(antib)) 

# shannon

shannon <- bplapply(colnames(mgs), function(x) lmer.fun(cacstot, mgs[, x], data.frame(cov, data$shannon),family), BPPARAM = MulticoreParam(cores))
shannon <- do.call(rbind, shannon)
shannon$q.value <- p.adjust(shannon$p.value, method = "BH", n = sum(!is.na(shannon$p.value)))
shannon <- shannon[, c("n", "estimate", "se", "p.value", "q.value")]
colnames(shannon) <- paste0("shannon_", colnames(shannon)) 

# make and export table

res <- data.frame(name, id, level, ppi, crohn, antib, shannon)
export(res, "results/sensitivity.tsv")

sessionInfo()
