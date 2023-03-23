# ===================
# = tx qtt analysis =
# ===================

args = commandArgs(TRUE) # args <- c("../qg/male.sig.gene.exp.RData", "male.25c.mean.pheno", "../eqtl/dgrp.common.fam")

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

# qtt cor
qtt.cor.young <- numeric(ncol(blup.young))
qtt.cor.old <- qtt.cor.young
qtt.cor.pval.young <- qtt.cor.young
qtt.cor.pval.old <- qtt.cor.young
qtt.cor.diff <- qtt.cor.young
qtt.cor.pval.diff <- qtt.cor.young

for (i in 1:length(qtt.cor.young)) {
	
	qtt.cor.young[i] <- cor(blup.young[, i], pheno)
	qtt.cor.pval.young[i] <- cor.test(blup.young[, i], pheno)$p.value
	qtt.cor.old[i] <- cor(blup.old[, i], pheno)
	qtt.cor.pval.old[i] <- cor.test(blup.old[, i], pheno)$p.value
	qtt.cor.diff[i] <- cor(blup.old[, i] - blup.young[, i], pheno)
	qtt.cor.pval.diff[i] <- cor.test(blup.old[, i] - blup.young[, i], pheno)$p.value
	
}

gene.name <- colnames(blup.young)

save(gene.name, qtt.cor.young, qtt.cor.old, qtt.cor.pval.young, qtt.cor.pval.old, qtt.cor.diff, qtt.cor.pval.diff, file = args[4])
