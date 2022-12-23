# load libraries and set options

library(rio)
library(fgsea)

set.seed(1)

# import and clean data

data <- import("processed/data.tsv")
info <- import("raw/HG3.A.7_tax.xlsx")
gmm <- import("raw/MGS_HG3A.GMMs2MGS.RData")
gmm_info <- import("raw/GMM_reference.csv")
basic <- import("results/basic.tsv")

colnames(info)[1:3] <- c("id", "name", "level")
colnames(gmm_info) <- c("id", "name", "hl1", "hl2")

# enrichment function

gsea.fun <- function(res, pathways) {
  
  stats <- -res$p.value
  names(stats) <- res$id
  
  pathways <- lapply(pathways, function(pathway) {
    
    pathway[which(pathway %in% names(stats))]
    
  })
  
  gsea1 <- tryCatch({
    
    gsea1 <- as.data.frame(fgsea(pathways, rank(stats[which(res$estimate >= 0)], na = "keep"), scoreType = "pos", eps = 0))
    data.frame(pathway = gsea1$pathway, estimate = gsea1$NES, p.value = gsea1$pval, size = gsea1$size, leading = sapply(gsea1$leadingEdge, function(x) paste(x, collapse = ";")), message = NA)
    
  }, warning = function(w) {
    
    data.frame(pathway = names(pathways), estimate = NA, p.value = NA, size = NA, leading = NA, message = w$message)
    
  }, error = function(e) {
    
    data.frame(pathway = names(pathways), estimate = NA, p.value = NA, size = NA, leading = NA, message = e$message)
    
  })
  
  gsea2 <- tryCatch({
    
    gsea2 <- as.data.frame(fgsea(pathways, rank(stats[which(res$estimate < 0)], na = "keep"), scoreType = "pos", eps = 0))
    data.frame(pathway = gsea2$pathway, estimate = gsea2$NES, p.value = gsea2$pval, size = gsea2$size, leading = sapply(gsea2$leadingEdge, function(x) paste(x, collapse = ";")), message = NA)
    
  }, warning = function(w) {
    
    data.frame(pathway = names(pathways), estimate = NA, p.value = NA, size = NA, leading = NA, message = w$message)
    
  }, error = function(e) {
    
    data.frame(pathway = names(pathways), estimate = NA, p.value = NA, size = NA, leading = NA, message = e$message)
    
  })
  
  q <- p.adjust(c(gsea1$p.value, gsea2$p.value), method = "BH", n = sum(!is.na(c(gsea1$p.value, gsea2$p.value))))
  gsea1$q.value <- q[1:nrow(gsea1)]
  gsea2$q.value <- q[-(1:nrow(gsea1))]
  id <- unique(c(gsea1$pathway, gsea2$pathway)[order(q)])
  gsea1 <- gsea1[match(id, gsea1$pathway), c("size", "estimate", "p.value", "q.value", "message")]
  gsea2 <- gsea2[match(id, gsea2$pathway), c("size", "estimate", "p.value", "q.value", "message")]
  colnames(gsea1) <- paste0("positive_", colnames(gsea1))
  colnames(gsea2) <- paste0("negative_", colnames(gsea2))
  data.frame(id, gsea1, gsea2)
  
}

# gmm

gmm.res <- gsea.fun(basic, gmm)
gmm.res <- data.frame(gmm.res, gmm_info[match(gmm.res$id, gmm_info$id), c("name", "hl1", "hl2")])
export(gmm.res, "results/gmm.tsv")

# genera

genera <- lapply(unique(info$genus), function(genus) {
  
  info$id[which(info$genus == genus)]
  
})
names(genera) <- unique(info$genus)
genera <- genera[which(!names(genera) == "unclassified")]
genera.res <- gsea.fun(basic, genera)
export(genera.res, "results/genera.tsv")

# leave-one-out

loo <- lapply(names(genera), function(genus) {
  
  basic <- basic[which(!basic$id %in% genera[[genus]]), ]
  gmm.res <- gsea.fun(basic, gmm)
  data.frame(genus = genus, gmm.res, gmm_info[match(gmm.res$id, gmm_info$id), c("name", "hl1", "hl2")])
  
})
loo <- do.call(rbind, loo)
loo <- loo[which(loo$id %in% gmm.res$id[which(gmm.res$positive_q.value < 0.05| gmm.res$negative_q.value < 0.05)]), ]
export(loo, "results/loo.tsv")
