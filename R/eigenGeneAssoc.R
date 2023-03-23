# =======================================
# = eigengene association with lifespan =
# =======================================

args <- commandArgs(TRUE) # args <- c("../qg/female.sig.gene.exp.RData", "../figureData/female.wgcna.RData", "../eqtl/dgrp.common.fam", female.$trait.pheno)

# 1. extract data
# ============================================================

load(args[1])
blup.old <- t(blup.18c.data[match(intersect(sig.gene.name.18c, sig.gene.name.25c), sig.gene.name), ])
colnames(blup.old) <- intersect(sig.gene.name.18c, sig.gene.name.25c)
blup.young <- t(blup.25c.data[match(intersect(sig.gene.name.18c, sig.gene.name.25c), sig.gene.name), ])
colnames(blup.young) <- intersect(sig.gene.name.18c, sig.gene.name.25c)
geno.line.order <- read.table(args[3], as.is = TRUE, header = F)[, 1]
rownames(blup.young) <- geno.line.order
rownames(blup.old) <- geno.line.order

# remove incomplete data upfront
# ============================================================

na.idx <- sort(unique(c(which(is.na(rowSums(blup.old))), which(is.na(rowSums(blup.young))))))
blup.old <- blup.old[-na.idx, ]
blup.young <- blup.young[-na.idx, ]

# retain only those with variance > 0.01
# ============================================================

sig.var <- which(apply(blup.old, 2, var) >= 0.01 & apply(blup.young, 2, var) >= 0.01)
blup.old <- blup.old[, sig.var]
blup.young <- blup.young[, sig.var]

# get modules
load(args[2])

pheno <- read.table(args[4], header = T, as.is = T)
rownames(pheno) <- pheno[, 1]

common.line <- intersect(rownames(pheno), rownames(blup.young))
pheno <- pheno[common.line, 3]
blup.young <- blup.young[common.line, ]
blup.old <- blup.old[common.line, ]

# pca on each modules
# ============================================================

unique.module <- setdiff(sort(unique(blup.25c.tree)), 0)
module.assoc.young <- matrix(NA, ncol = 4, nrow = length(unique.module))
module.assoc.old <- module.assoc.young
module.assoc.diff <- module.assoc.young

for (i in 1:length(unique.module)) {
  
  module.exp.young <- blup.young[, which(blup.25c.tree == unique.module[i])]
  pca.module.young <- predict(prcomp(module.exp.young, scale = T))[, 1]
  module.assoc.young[i, ] <- c(unique.module[i], sum(blup.25c.tree == unique.module[i]), unlist(cor.test(pca.module.young, pheno))[c(4, 3)])
  
  module.exp.old <- blup.old[, which(blup.25c.tree == unique.module[i])]
  pca.module.old <- predict(prcomp(module.exp.old, scale = T))[, 1]
  module.assoc.old[i, ] <- c(unique.module[i], sum(blup.25c.tree == unique.module[i]), unlist(cor.test(pca.module.old, pheno))[c(4, 3)])
  
  module.exp.diff <- blup.old[, which(blup.25c.tree == unique.module[i])] - blup.young[, which(blup.25c.tree == unique.module[i])]
  pca.module.diff <- predict(prcomp(module.exp.diff, scale = T))[, 1]
  module.assoc.diff[i, ] <- c(unique.module[i], sum(blup.25c.tree == unique.module[i]), unlist(cor.test(pca.module.diff, pheno))[c(4, 3)])
  
}

save(module.assoc.young, module.assoc.old, module.assoc.diff, file = args[5])
