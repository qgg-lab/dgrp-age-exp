# ====================
# = network analysis =
# ====================

args <- commandArgs(TRUE) # args <- c("../figureData/female.wgcna.RData", "../figureData/gene_go.table", "../figureData/fbgn.kegg", "../figureData/gene.info", "../figureData/kegg.id")

library("openxlsx")
library("GO.db")

# read go table
# ============================================================

annot <- read.table(args[2], header = TRUE, as.is = TRUE, sep = " ") # annotation table
gene.info <- read.table(args[4], header = FALSE, as.is = TRUE)
rownames(gene.info) <- gene.info[, 1]
kegg <- read.table(args[3], header = FALSE, as.is = TRUE)
kegg.id <- read.table(args[5], header = FALSE, as.is = TRUE, sep = "\t")     
rownames(kegg.id) <- kegg.id[, 1]

# read female wgcna
# ============================================================

load(args[1])
gene.id <- colnames(blup.25c.cor)

bp.annot <- annot[annot[, 1] %in% gene.id & annot[, 3] != "", c(1, 3)]
bp.gene.id <- bp.annot[, 1]
bp.count <- table(unlist(strsplit(bp.annot[, 2], split = ",")))
unique.bp <- names(bp.count[bp.count > 20])

kegg <- kegg[kegg[, 1] %in% gene.id, ]
kegg.gene.id <- unique(kegg[, 1])
kegg.count <- table(kegg[, 2])
unique.kegg <- names(kegg.count[kegg.count > 20])

unique.bp <- setdiff(unique.bp, "GO:0006099")
unique.kegg <- setdiff(unique.kegg, "dme00020")

# loop through network
# ============================================================

net.diff <- matrix(ncol = 7, nrow = length(unique.bp) + length(unique.kegg) + 1)
rownames(net.diff) <- c(unique.bp, unique.kegg, "GO:0006099 + dme00020")

for (i in 1:length(unique.bp)) {
	
	this.go <- unique.bp[i]
	this.go.genes <- bp.annot[grepl(this.go, bp.annot[, 2]), 1]
	this.young.cor <- blup.25c.cor[this.go.genes, this.go.genes]
	this.old.cor <- blup.18c.cor[this.go.genes, this.go.genes]
  
  
	
	net.diff[i, ] <- c(this.go, Term(as.character(this.go)), length(this.go.genes), sum(abs(this.young.cor[lower.tri(this.young.cor, diag = F)]) > 0.5),
												sum(abs(this.old.cor[lower.tri(this.old.cor, diag = F)]) > 0.5),
                        mean(abs(this.young.cor[lower.tri(this.young.cor, diag = F)])),
                        mean(abs(this.old.cor[lower.tri(this.old.cor, diag = F)])))
	
}

for (i in 1:length(unique.kegg)) {
  
  this.kegg <- unique.kegg[i]
  this.kegg.genes <- kegg[kegg[, 2] == this.kegg, 1]
	this.young.cor <- blup.25c.cor[this.kegg.genes, this.kegg.genes]
	this.old.cor <- blup.18c.cor[this.kegg.genes, this.kegg.genes]
	net.diff[length(unique.bp) + i, ] <- c(this.kegg, kegg.id[this.kegg, 2], length(this.kegg.genes), sum(abs(this.young.cor[lower.tri(this.young.cor, diag = F)]) > 0.5),
												sum(abs(this.old.cor[lower.tri(this.old.cor, diag = F)]) > 0.5),
                        mean(abs(this.young.cor[lower.tri(this.young.cor, diag = F)])),
                        mean(abs(this.old.cor[lower.tri(this.old.cor, diag = F)])))
  
}

this.tca <- "GO:0006099 + dme00020"
this.tca.genes <- unique(c(kegg[kegg[, 2] == "dme00020", 1], bp.annot[grepl("GO:0006099", bp.annot[, 2]), 1]))
this.young.cor <- blup.25c.cor[this.tca.genes, this.tca.genes]
this.old.cor <- blup.18c.cor[this.tca.genes, this.tca.genes]
net.diff[length(unique.bp) + length(unique.kegg) + 1, ] <- c(this.tca, "TCA cycle", length(this.tca.genes), sum(abs(this.young.cor[lower.tri(this.young.cor, diag = F)]) > 0.5),
											sum(abs(this.old.cor[lower.tri(this.old.cor, diag = F)]) > 0.5),
                      mean(abs(this.young.cor[lower.tri(this.young.cor, diag = F)])),
                      mean(abs(this.old.cor[lower.tri(this.old.cor, diag = F)])))


net.diff <- net.diff[(as.numeric(net.diff[, 4]) > 10 | as.numeric(net.diff[, 5]) > 10) & !is.na(net.diff[, 2]), ]
write.xlsx(list('final' = net.diff), file = args[6], colNames = TRUE, rowNames = FALSE)

