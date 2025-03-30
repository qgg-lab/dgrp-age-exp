# =========================================
# = testis specificity index and pathways =
# =========================================


args <- commandArgs(TRUE) # args = c("../figureData/gene_go.table", "../figureData/fbgn.kegg", "../figureData/testis.specificity.index.2021")
library("RColorBrewer")

# prepare file to plot
# ============================================================

file.width = 150
cairo_pdf(file = args[6], width = file.width/25.4, height = file.width/25.4*0.35, family = args[7])
layout(matrix(c(1, 2), ncol = 2, byrow = T), widths = c(1, 0.5))

# load data for pathways
# ============================================================

annot <- read.table(args[1], header = TRUE, as.is = TRUE, sep = " ") # 
kegg <- read.table(args[2], header = FALSE, as.is = TRUE)

# data for testis
# ============================================================

testis <- read.table(args[3], header = FALSE, as.is = TRUE)
rownames(testis) <- testis[, 6]

# GO ids:
# ============================================================

# female GO:0006120 -> mitochondrial electron transport, NADH to ubiquinone
# female dme00020 -> tca cycle
# male dme00190 -> Oxidative phosphorylation
# male dme01040 ->	Biosynthesis of unsaturated fatty acids

bp.go.genes <- annot[grepl(go.id, annot[, 3]), 1]
kegg.genes <- kegg[kegg[, 2] == kegg.id, 1]

testis.index <- rbind(rbind(rbind(rbind(cbind("all", testis[, 9]), cbind("TCA", testis[kegg[kegg[, 2] == "dme00020", 1], 9])), cbind("Mitochondrial", testis[annot[grepl("GO:0006120", annot[, 3]), 1], 9])), cbind("Oxi. Phosph.", testis[kegg[kegg[, 2] == "dme00190", 1], 9])), cbind("Fatty acids", testis[kegg[kegg[, 2] == "dme01040", 1], 9]))

boxplot(as.numeric(testis.index[, 2]) ~ testis.index[, 1], ylab = "Testis Specificity Index", xlab = "")




