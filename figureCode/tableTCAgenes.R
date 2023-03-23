# ===================
# = TCA cycle genes =
# ===================

args <- commandArgs(TRUE) # args <- c("../figureData/female.wgcna.RData", "../figureData/male.wgcna.RData", "../figureData/gene_go.table", "../figureData/fbgn.kegg", "../figureData/gene.info")

go.id <- "GO:0006099"
kegg.id <- "dme00020"

library("openxlsx")

# read go table
# ============================================================

annot <- read.table(args[3], header = TRUE, as.is = TRUE, sep = " ") # annotation table
gene.info <- read.table(args[5], header = FALSE, as.is = TRUE)
rownames(gene.info) <- gene.info[, 1]

# read kegg table

kegg <- read.table(args[4], header = FALSE, as.is = TRUE)

# find genes in go and kegg
# ============================================================

bp.go.genes <- annot[grepl(go.id, annot[, 3]), 1]
kegg.genes <- kegg[kegg[, 2] == kegg.id, 1]

# make table
# ============================================================

both.genes <- sort(unique(c(bp.go.genes, kegg.genes)))
final.table <- cbind(both.genes, ifelse(both.genes %in% bp.go.genes, "x", ""), ifelse(both.genes %in% kegg.genes, "x", ""), gene.info[both.genes, 2])[, c(1, 4, 2, 3)]

# write.table(final.table, file = "../figureData/tca.genes.txt", col.names = FALSE, sep = " ", row.names = F, quote = F)

# write.table
# ============================================================

colnames(final.table) <- c("FBgn", "symbol", "GO:0006099", "dme00020")

write.xlsx(list('module' = final.table), file = args[6], colNames = TRUE, rowNames = FALSE)

