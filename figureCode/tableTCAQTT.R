# ================
# = qtt table (TCA) =
# ================

args <- commandArgs(TRUE) # args <- c("../figureData/female.qtt.RData", "../figureData/male.qtt.RData", "../figureData/gene_go.table", "../figureData/fbgn.kegg", "../figureData/gene.info")

go.id <- "GO:0006099"
kegg.id <- "dme00020"

library("openxlsx")

# read data
# ============================================================

annot <- read.table(args[3], header = TRUE, as.is = TRUE, sep = " ") # annotation table
gene.info <- read.table(args[5], header = FALSE, as.is = TRUE)
rownames(gene.info) <- gene.info[, 1]

# read kegg table

kegg <- read.table(args[4], header = FALSE, as.is = TRUE)

bp.go.genes <- annot[grepl(go.id, annot[, 3]), 1]
kegg.genes <- kegg[kegg[, 2] == kegg.id, 1]

tca.genes <- unique(c(bp.go.genes, kegg.genes))

# female
# ============================================================

load(args[1])
names(qtt.cor.young) <- gene.name
names(qtt.cor.old) <- gene.name
names(qtt.cor.diff) <- gene.name

tca.cor <- cbind(qtt.cor.young[tca.genes], qtt.cor.old[tca.genes], qtt.cor.diff[tca.genes])

load(args[2])
names(qtt.cor.young) <- gene.name
names(qtt.cor.old) <- gene.name
names(qtt.cor.diff) <- gene.name

tca.cor <- cbind(tca.cor, qtt.cor.young[tca.genes], qtt.cor.old[tca.genes], qtt.cor.diff[tca.genes])

not.all.na <- which(apply(tca.cor, 1, function(x){sum(is.na(x))}) < 6)

final.table <- cbind(tca.genes[not.all.na], gene.info[tca.genes[not.all.na], 2], tca.cor[not.all.na, ])

write.xlsx(list('module' = final.table), file = args[6], colNames = TRUE, rowNames = FALSE)
