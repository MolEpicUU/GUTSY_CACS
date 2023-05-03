# load libraries and set options

library(rio)
library(GUniFrac)

set.seed(1)

# import data

key <- import("raw/id_conversion.txt")
data <- import("processed/data.tsv")
tax_pheno <- import("raw/pheno_MGS_shannon_bray_curtis_MGP_4839_upp_4980_malmo.tsv")
beta <- import("raw/bray_curtis_dsMGS_4838_upp_4980_malmo.tsv")

data=data[!data$family1==0,]

key$pheno_id <- key$export_id
key$pheno_id[which(grepl("2-", key$subject_id))] <- key$subject_id[which(grepl("2-", key$subject_id))]
colnames(tax_pheno)[2] <- "pheno_id"
tax_pheno$scapis_id <- key$subject_id[match(tax_pheno$pheno_id, key$pheno_id)]
colnames(beta) <- tax_pheno$scapis_id[match(colnames(beta), tax_pheno$sample.id)]
rownames(beta) <- colnames(beta)
id <- data$scapis_id[which(data$scapis_id %in% rownames(beta))]
beta <- beta[match(id, rownames(beta)), match(id, colnames(beta))]
data <- data[match(id, data$scapis_id), ]

print("identical beta and data")
identical(colnames(beta),data$scapis_id)

# dmanova function

dmanova.fun <- function(beta, x, z) {
  
  tryCatch({
    
    x <- log1p(x)
    data <- data.frame(z, x)
    beta <- beta[which(complete.cases(data)), which(complete.cases(data))]
    data <- data[which(complete.cases(data)), ]
    res <- dmanova(as.dist(beta) ~ ., data = data)$aov.tab
   print(res)
    data.frame(r.squared = res[1, 5], p.value = res[1, 6], n = nrow(data), message = NA)
    
  }, warning = function(w) {
    
    data.frame(r.squared = NA, p.value = NA, n = NA, message = paste("Warning:", w$message))
    
  }, error = function(e) {
    
    data.frame(r.squared = NA, p.value = NA, n = NA, message = paste("Error:", e$message))
    
  })
  
}

# basic model

cacstot <- data$cacstot
cov <- data[, c("age", "sex", "country", "site_plate")]
basic <- dmanova.fun(beta, cacstot, cov)
colnames(basic) <- paste0("basic_", colnames(basic))

# main model

cov <- data[, c("age", "sex", "country", "site_plate", "smoke", "pa", "carb", "protein", "fiber", "sbp", "dbp", "chol", "hdl", "ldl", "tg", "diab", "bmi", "chol_med", "bp_med", "diab_med")]
main <- dmanova.fun(beta, cacstot, cov)
colnames(main) <- paste0("main_", colnames(main))

# export data

res <- cbind(basic, main)
export(res, "results/beta.tsv")

sessionInfo()
