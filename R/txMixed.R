# ==================
# = tx mixed model =
# ==================

# module load GCC/11.2.0 OpenMPI/4.1.1 R/4.1.2
args = commandArgs(TRUE) # args <- c("../qg/female.sig.gene.exp.RData", "female.25c.mean.pheno", "../eqtl/dgrp.common.fam")
library("sommer")

# read data
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

# phenotype
pheno <- read.table(args[2], header = T, as.is = T)
rownames(pheno) <- pheno[, 1]

common.line <- intersect(rownames(pheno), rownames(blup.young))
pheno <- pheno[common.line, 3]
blup.young <- blup.young[common.line, ]
blup.old <- blup.old[common.line, ]

# remove low variance
# ============================================================

sig.var <- which(apply(blup.old, 2, var) >= 0.01 & apply(blup.young, 2, var) >= 0.01)
blup.old <- blup.old[, sig.var]
blup.young <- blup.young[, sig.var]

# compute K matrix
# ============================================================

pheno.data <- data.frame(y = pheno, id <- factor(1:length(pheno)))

blup.young.scaled <- scale(blup.young)
k.young <- blup.young.scaled %*% t(blup.young.scaled)
k.young <- k.young/mean(diag(k.young))
rownames(k.young) <- 1:length(pheno)
colnames(k.young) <- 1:length(pheno)
young.mix <- mmer(y ~ 1, random =~ vsr(id, Gu = k.young), rcov =~ units, data = pheno.data)
young.mix.h2 <- vpredict(young.mix, h1 ~ V1/(V1+V2))

blup.old.scaled <- scale(blup.old)
k.old <- blup.old.scaled %*% t(blup.old.scaled)
k.old <- k.old/mean(diag(k.old))
rownames(k.old) <- 1:length(pheno)
colnames(k.old) <- 1:length(pheno)
old.mix <- mmer(y ~ 1, random =~ vsr(id, Gu = k.old), rcov =~ units, data = pheno.data)
old.mix.h2 <- vpredict(old.mix, h1 ~ V1/(V1+V2))

blup.diff.scaled <- scale(blup.old - blup.young)
k.diff <- blup.diff.scaled %*% t(blup.diff.scaled)
k.diff <- k.diff/mean(diag(k.diff))
rownames(k.diff) <- 1:length(pheno)
colnames(k.diff) <- 1:length(pheno)
diff.mix <- mmer(y ~ 1, random =~ vsr(id, Gu = k.diff), rcov =~ units, data = pheno.data)
diff.mix.h2 <- vpredict(diff.mix, h1 ~ V1/(V1+V2))

save(young.mix.h2, old.mix.h2, diff.mix.h2, file = args[4])
