# load libraries and set options

library(rio)
library(BiocParallel)
library(logistf)

cores <- 16
set.seed(1)

# import and clean data

mgs <- import("raw/ERP023788_chinese.MGS.comp.rds")
pheno <- import("raw/validation_cohort.tsv")
info <- import("raw/HG3.A.7_tax.xlsx")
main_sig <- import("results/main_significant.tsv")

mgs <- mgs[, which(colSums(mgs > 0) > 10)]
colnames(info)[1:3] <- c("id", "name", "level")
id <- main_sig$id
id <- id[which(!id == "HG3A.1967")]
id <- id[which(id %in% colnames(mgs))]
name <- info$name[match(id, info$id)]
level <- info$level[match(id, info$id)]
pheno <- pheno[which(complete.cases(pheno[, c("disease", "age", "gender")])), ]
pheno <- pheno[match(rownames(mgs), pheno$subject_id), ]
main_sig <- main_sig[match(id, main_sig$id), ]

# case-control

logistf.fun <- function(y, x, z) {
  
  tryCatch({
    
    x <- scale(log1p(x))
    data <- data.frame(y = y, x = x, z)
    data <- data[which(complete.cases(data)), ]
    fit <- logistf(y ~ ., data, control= logistf.control(maxit = 10000))
    coef <- fit$coefficients
    se <- sqrt(diag(vcov(fit)))
    lower <- fit$ci.lower
    upper <- fit$ci.upper
    p <- fit$prob
    data.frame(or = exp(coef[2]), lower = exp(lower[2]), upper = exp(upper[2]), se = se[2], p.value = p[2], n = nrow(data), message = NA)
    
  }, warning = function(w) {
    
    data.frame(or = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Warning:", w$message))
    
  }, error = function(e) {
    
    data.frame(or = NA, lower = NA, upper = NA, se = NA, p.value = NA, n = NA, message = paste("Error:", e$message))
    
  })
  
}

cvd <- ifelse(pheno$disease == "ACVD", 1, 0)
mgs <- mgs[, match(id, colnames(mgs))]
cov <- pheno[, c("age", "gender")]

casecontrol <- bplapply(colnames(mgs), function(x) logistf.fun(cvd, mgs[, x], cov), BPPARAM = MulticoreParam(cores))
casecontrol <- do.call(rbind, casecontrol)
casecontrol$id <- id
casecontrol$name <- name
casecontrol$level <- level
casecontrol$prev.cases <- apply(mgs[which(cvd == 1), ], 2, function(x) sum(x > 0) / length(x))
casecontrol$prev.controls <- apply(mgs[which(cvd == 0), ], 2, function(x) sum(x > 0) / length(x))
casecontrol$q.value <- p.adjust(casecontrol$p.value, method = "BH", n = sum(!is.na(casecontrol$p.value)))
casecontrol$direction <- NA
casecontrol$direction[which(main_sig$estimate >= 0 & casecontrol$or >= 1)] <- "++"
casecontrol$direction[which(main_sig$estimate >= 0 & casecontrol$or < 1)] <- "+-"
casecontrol$direction[which(main_sig$estimate < 0 & casecontrol$or >= 1)] <- "-+"
casecontrol$direction[which(main_sig$estimate < 0 & casecontrol$or < 1)] <- "--"
casecontrol <- casecontrol[, c("name", "id", "level", "prev.cases", "prev.controls", "or", "lower", "upper", "p.value", "q.value", "direction", "n", "message")]

# export table

export(casecontrol, "results/casecontrol.tsv")

sessionInfo()