# load libraries and set options

library(rio)
library(BiocParallel)
library(lmerTest)
library(ordinal)

cores <- 16
set.seed(1)

# import and clean data

data <- import("processed/mods.tsv")
gut_info <- import("raw/HG3.A.7_tax.xlsx")
oral_info <- import("processed/mods_tax.tsv")
main_sig <- import("results/main_significant.tsv")

id <- main_sig$id
id <- id[which(!id == "HG3A.1967")]
colnames(gut_info)[1:3] <- c("id", "name", "level")
colnames(oral_info)[1:3] <- c("id", "name", "level")

# calculate overlap between gut and oral

overlap <- lapply(id, function(x) {
  
  if (any(oral_info$species %in% gut_info$species[which(gut_info$id == x)])) {
    
    species <- gut_info$species[which(gut_info$id == x)]
    oral <- oral_info$id[which(oral_info$species == species)]
    res <- data.frame(gut = x, species = species, oral = oral)
    res[which(!res$species == "unclassified"), ]
    
  }
  
})
overlap <- do.call(rbind, overlap)

# spearman mixed correlation function

spearman.fun <- function(y, x, z, re) {
  
  tryCatch({
    
    data <- data.frame(y = y, x = x, z, re = re)
    data <- data[which(complete.cases(data)), ]
    re2 <- data$re
    data <- model.matrix(~., data[, -ncol(data)])[, -1]
    data <- as.data.frame(apply(data, 2, function(x) rank(x, na = "keep")))
    data$re <- re2
    res1 <- residuals(lmer(y ~ . - x - re + (1 | re), data), type = "response")
    res2 <- residuals(lmer(x ~ . - y - re + (1 | re), data), type = "response")
    fit <- lm(res1 ~ res2)
    estimate <- sqrt(summary(fit)$r.squared) * sign(summary(fit)$coefficients[2, 1])
    p <- summary(lmer(y ~ . - re + (1 | re), data))$coefficients[2, 5]
    data.frame(estimate = estimate, p.value = p, n = nrow(data), message = NA)
    
  }, warning = function(w) {
    
    data.frame(estimate = NA, p.value = NA, n = NA, message = paste("Warning:", w$message))
    
  }, error = function(e) {
    
    data.frame(estimate = NA, p.value = NA, n = NA, message = paste("Error:", e$message))
    
  })
  
}

# gut-oral

cov <- data[, c("age", "sex", "country", "plate")]
gut_oral <- bplapply(1:nrow(overlap), function(i) spearman.fun(data[, overlap$gut[i]], data[, overlap$oral[i]], cov, data$family), BPPARAM = MulticoreParam(cores))
gut_oral <- do.call(rbind, gut_oral)
gut_oral$q.value <- p.adjust(gut_oral$p.value, method = "BH", n = sum(!is.na(gut_oral$p.value)))
gut_oral$gut_id <- overlap$gut
gut_oral$gut_name <- gut_info$name[match(gut_oral$gut_id, gut_info$id)]
gut_oral$gut_level <- gut_info$level[match(gut_oral$gut_id, gut_info$id)]
gut_oral$oral_id <- overlap$oral
gut_oral$oral_name <- oral_info$name[match(gut_oral$oral_id, oral_info$id)]
gut_oral$oral_level <- oral_info$level[match(gut_oral$oral_id, oral_info$id)]
gut_oral <- gut_oral[, c("gut_name", "gut_id", "gut_level", "oral_name", "oral_id", "oral_level", "n", "estimate", "p.value", "q.value", "message")]
export(gut_oral, "results/gut_oral.tsv")

# dental hygiene

id <- gut_oral$oral_id[which(gut_oral$q.value < 0.05)]
cov <- data[, c("age", "sex", "smoke", "plaque", "education", "lasthour")]
mgs <- data[, match(id, colnames(data))]

# ordinal mixed regression function

ordinal.fun <- function(y, x, z, re) {
  
  tryCatch({
    
    x <- scale(log1p(x))
    data <- data.frame(y = y, x = x, z, re)
    data <- data[which(complete.cases(data)), ]
    fit <- clmm(ordered(y) ~ . - re + (1 | re), data = data)
    coef <- summary(fit)$coefficients
    lower <- coef["x", 1] + qnorm(0.025) * coef["x", 2]
    upper <- coef["x", 1] + qnorm(0.975) * coef["x", 2]
    data.frame(or = exp(coef["x", 1]), lower = exp(lower), upper = exp(upper), se = coef["x", 2], p.value = coef["x", 4], n = nrow(data), message = NA)
    
  }, warning = function(w) {
    
    data.frame(or = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Warning:", w$message))
    
  }, error = function(e) {
    
    data.frame(or = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Error:", e$message))
    
  })
  
}

# caries

caries <- bplapply(colnames(mgs), function(x) ordinal.fun(data$caries, mgs[, x], cov, data$family), BPPARAM = MulticoreParam(cores))
caries <- do.call(rbind, caries)
caries$q.value <- p.adjust(caries$p.value, method = "BH", n = sum(!is.na(caries$p.value)))
caries <- caries[, c("n", "or", "lower", "upper", "p.value", "q.value", "message")]
colnames(caries) <- paste0("caries_", colnames(caries))

# fillings

fillings <- cut(rank(data$fillings, "keep"), 10, labels = F)
fillings <- bplapply(colnames(mgs), function(x) ordinal.fun(fillings, mgs[, x], cov, data$family), BPPARAM = MulticoreParam(cores))
fillings <- do.call(rbind, fillings)
fillings$q.value <- p.adjust(fillings$p.value, method = "BH", n = sum(!is.na(fillings$p.value)))
fillings <- fillings[, c("n", "or", "lower", "upper", "p.value", "q.value", "message")]
colnames(fillings) <- paste0("fillings_", colnames(fillings))

# linear mixed model function

lmer.fun <- function(y, x, z, re) {
  
  tryCatch({
    
    x <- scale(log1p(x))
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

# gingivitis

gingivitis <- bplapply(colnames(mgs), function(x) lmer.fun(data$gingivitis, mgs[, x], cov, data$family), BPPARAM = MulticoreParam(cores))
gingivitis <- do.call(rbind, gingivitis)
gingivitis$q.value <- p.adjust(gingivitis$p.value, method = "BH", n = sum(!is.na(gingivitis$p.value)))
gingivitis <- gingivitis[, c("n", "estimate", "lower", "upper", "p.value", "q.value", "message")]
colnames(gingivitis) <- paste0("gingivitis_", colnames(gingivitis))

# make and export table

res <- data.frame(name = gut_oral$oral_name[which(gut_oral$q.value < 0.05)], id = id, level = gut_oral$oral_level[which(gut_oral$q.value < 0.05)], caries, fillings, gingivitis)
export(res, "results/oral_health.tsv")

sessionInfo()
