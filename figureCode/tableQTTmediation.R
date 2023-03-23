# ==============================
# = mediation analysis results =
# ==============================

args <- commandArgs(TRUE) # args <- c("../figureData/female.eqtl.gwas.trio.mediation.RData", "../figureData/male.eqtl.gwas.trio.mediation.RData", "../figureData/gene.info")

library("openxlsx")

# read go table
# ============================================================

gene.info <- read.table(args[3], header = FALSE, as.is = TRUE)
rownames(gene.info) <- gene.info[, 1]

# load data
load(args[1])
female.res <- cbind(trio.res[, 1], gene.info[trio.res[, 2], ], trio.res[, 3:ncol(trio.res)])
load(args[2])
male.res <- cbind(trio.res[, 1], gene.info[trio.res[, 2], ], trio.res[, 3:ncol(trio.res)])

# col names
col.names <- c("snp", "gene", "symbol", "type", "location", "trait", "other-snps", "other-snps-r2", "discovery", "snp-trait-eff", "snp-trait-pval", "snp-young-exp-eff", "snp-young-exp-pval", "snp-old-exp-eff", "snp-old-exp-pval", "snp-diff-exp-eff", "snp-diff-exp-pval", "young-trait-cor", "young-trait-cor-pval", "old-trait-cor", "old-trait-cor-pval", "diff-trait-cor", "diff-trait-cor-pval", paste(rep(c("young", "old", "diff"), each = 16), rep(c("acme-est", "acme-low", "acme-high", "acme-pval", "ade-est", "ade-low", "ade-high", "ade-pval", "prop-est", "prop-low", "prop-high", "prop-pval", "total-est", "total-low", "total-high", "total-pval"), 3), sep = "-"))

# make table
# ============================================================

colnames(female.res) <- col.names
colnames(male.res) <- col.names

# write.table
# ============================================================

write.xlsx(list('female' = female.res, 'male' = male.res), file = args[4], colNames = TRUE, rowNames = FALSE)

