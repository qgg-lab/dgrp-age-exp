# ===============
# = TCA network =
# ===============

args <- commandArgs(TRUE) # args <- c("../figureData/female.wgcna.RData", "../figureData/gene_go.table", "../figureData/fbgn.kegg", "../figureData/gene.info")

library("igraph")
library("RColorBrewer")
library("fields") #to use designer.colors
library("diagram")

go.id <- "GO:0006099"
kegg.id <- "dme00020"

module.col <- c("grey80", brewer.pal(8, "Dark2"), brewer.pal(9, "Set1"))

# function to re-order and get mean connectivity across module
# order based on module connectivity then connectivity within module
# ============================================================

orderCor <- function(cor.mat, module) {
	
	module <- as.character(module)
	unique.module <- unique(module)
	module.mean.corr <- numeric(length(unique.module))
	names(module.mean.corr) <- unique.module
	gene.conn <- numeric(nrow(cor.mat))
	names(gene.conn) <- rownames(cor.mat)
	
	for (i in unique.module) {
		
    this.module.idx <- which(module == i)
    if (length(this.module.idx) > 1) {
		  this.cor <- cor.mat[which(module == i), which(module == i)]
      module.mean.corr[i] <- mean(abs(this.cor[lower.tri(this.cor, diag = FALSE)]))
		  gene.conn[rownames(this.cor)] <- rowMeans(this.cor)
    } else {
      module.mean.corr[i] <- 0
      gene.conn[this.module.idx] <- 0
    }
		
	}
	
	gene.order <- order(module.mean.corr[module], gene.conn)
	cor.mat <- cor.mat[gene.order, gene.order]
	module <- module[gene.order]
	
  # # get mean cross module correlation
  # for (i in 1:(length(unique.module) - 1)) {
  #   for (j in (i+1):length(unique.module)) {
  #
  #     module1 <- which(module == unique.module[i])
  #     module2 <- which(module == unique.module[j])
  #     mean.cross.cor <- mean(cor.mat[module1, module2])
  #     cor.mat[module1, module2] <- mean.cross.cor
  #     cor.mat[module2, module1] <- mean.cross.cor
  #
  #   }
  # }
  	
	return(list(cor = cor.mat, mod = module, gene.name = colnames(cor.mat)))
	
}

# prepare file
# ============================================================

file.width = 89
cairo_pdf(file = args[5], width = file.width/25.4, height = file.width/25.4, family = args[6])
par(las = 1, tcl = -0.2, mar = c(0.5, 0.5, 0.5, 0.5), ps = 7, lwd = 0.5)
layout(mat = matrix(c(1, 2, 3, 4, 5, 6), ncol = 3, byrow = TRUE), widths = c(1, 0.25, 1))

# read go table
# ============================================================

annot <- read.table(args[2], header = TRUE, as.is = TRUE, sep = " ") # annotation table
gene.info <- read.table(args[4], header = FALSE, as.is = TRUE)
rownames(gene.info) <- gene.info[, 1]

# read kegg table

kegg <- read.table(args[3], header = FALSE, as.is = TRUE)

# read female wgcna
# ============================================================

load(args[1])
gene.id <- colnames(blup.25c.cor)

mf.annot <- annot[annot[, 1] %in% gene.id & annot[, 2] != "", c(1, 2)]
bp.annot <- annot[annot[, 1] %in% gene.id & annot[, 3] != "", c(1, 3)]
cc.annot <- annot[annot[, 1] %in% gene.id & annot[, 4] != "", c(1, 4)]

mf.gene.id <- mf.annot[, 1]
bp.gene.id <- bp.annot[, 1]
cc.gene.id <- cc.annot[, 1]

unique.mf <- unique(unlist(strsplit(mf.annot[, 2], split = ",")))
unique.bp <- unique(unlist(strsplit(bp.annot[, 2], split = ",")))
unique.cc <- unique(unlist(strsplit(cc.annot[, 2], split = ",")))

# genes
# ============================================================

bp.go.gene <- bp.annot[grepl(go.id, bp.annot[, 2]), 1]
kegg.gene <- intersect(kegg[kegg[, 2] == kegg.id, 1], gene.id)
all.tca.gene <- unique(c(bp.go.gene, kegg.gene))

young.cor <- blup.25c.cor[all.tca.gene, all.tca.gene]
old.cor <- blup.18c.cor[all.tca.gene, all.tca.gene]

# scsbetaA
cbind(young.cor["FBgn0037643", ], old.cor["FBgn0037643", ], gene.info[colnames(young.cor), 2])
# cg5214 FBgn0037891
cbind(young.cor["FBgn0037891", ], old.cor["FBgn0037891", ], gene.info[colnames(young.cor), 2])
# kdn FBgn0261955
cbind(young.cor["FBgn0261955", ], old.cor["FBgn0261955", ], gene.info[colnames(young.cor), 2])

# get network
# ============================================================

young.edge <- young.cor
young.gene.id <- colnames(young.edge)
colnames(young.edge) <- gene.info[colnames(young.edge), 2]
rownames(young.edge) <- gene.info[rownames(young.edge), 2]
diag(young.edge) <- 0
young.edge[abs(young.edge) < 0.5] <- 0
sum(young.edge["ScsbetaA", ] > 0)
young.graph <- graph.adjacency(young.edge, mode = "undirected", weighted = TRUE)

old.edge <- old.cor
old.gene.id <- colnames(old.edge)
colnames(old.edge) <- gene.info[colnames(old.edge), 2]
rownames(old.edge) <- gene.info[rownames(old.edge), 2]
diag(old.edge) <- 0
old.edge[abs(old.edge) < 0.5] <- 0
old.graph <- graph.adjacency(old.edge, mode = "undirected", weighted = TRUE)

old.edge <- old.cor
diag(old.edge) <- 0
old.edge[abs(old.edge) < 0.5] <- 0

set.seed(1)
old.layout <- layout_with_fr(old.graph)


# correlation matrix
# ============================================================
names(blup.18c.tree) <- names(blup.25c.tree)
young.order.cor <- orderCor(young.cor, blup.25c.tree[all.tca.gene])
young.order.module <- young.order.cor$mod
names(young.order.module) <- rownames(young.order.cor$cor)

old.order <- rownames(young.order.cor$cor)
new.order <- rownames(young.order.cor$cor)
first.gene.idx <- which(old.order == "FBgn0052026")
second.gene.idx <- which(old.order == "FBgn0261955")
new.idx <- 1:length(new.order)

new.order[first.gene.idx] <- old.order[second.gene.idx]
new.order[second.gene.idx] <- old.order[first.gene.idx]

young.order.cor$cor <- young.order.cor$cor[new.order, new.order]
young.order.module <- young.order.module[new.order]

old.order.cor <- old.cor[new.order, new.order]
old.order.module <- blup.18c.tree[new.order]

cat(young.order.module, "\n")
cat(old.order.module, "\n")

par(las = 1, tcl = -0.2, mar = c(1, 1, 1, 1), ps = 7, lwd = 0.5)

plot(c(0, length(all.tca.gene)), c(0, length(all.tca.gene)), type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1)
#image(1:length(all.tca.gene), 1:length(all.tca.gene), young.order.cor$cor, col = colorRampPalette(c(brewer.pal(9, "Blues")[9], "white", brewer.pal(9, "Reds")[9]))(101), useRaster = F, add = TRUE)
image(1:length(all.tca.gene), 1:length(all.tca.gene), young.order.cor$cor, col = rev(designer.colors(n = 101, col = brewer.pal(9, "RdBu"))), useRaster = F, add = TRUE)
# unique.module <- unique(young.order.cor[[2]])
# module.start <- sapply(split(1:length(young.order.cor[[2]]), factor(young.order.cor[[2]], levels = unique.module)), min)
# module.end <- sapply(split(1:length(young.order.cor[[2]]), factor(young.order.cor[[2]], levels = unique.module)), max)
# for (i in 1:length(module.start)) {
#
# 	rect(module.start[i] - 0.5, module.start[i] - 0.5, module.end[i] + 0.5, module.end[i] + 0.5, lwd = 1)
#
# }


module1 <- which(young.order.module == as.character(old.order.module) & young.order.module == "1")
rect(min(module1) - 0.5, min(module1) - 0.5, max(module1) + 0.5, max(module1) + 0.5, lwd = 1)



text(1:length(all.tca.gene) + 1.5, rep(0, length(all.tca.gene)), gene.info[rownames(young.order.cor$cor), 2], font = 3, srt = 45, pos = 2, xpd = TRUE, col = ifelse(gene.info[rownames(young.order.cor$cor), 2] %in% c("kdn", "ScsbetaA"), "red", "black"), cex = 0.75)
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)
text(0, length(all.tca.gene) + 2, "Young", pos = 4, cex = 7/par("ps")/par("cex"), xpd = TRUE)


par(las = 1, tcl = -0.2, mar = c(0, 0, 0, 0), ps = 7, lwd = 0.5)
plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "'")
#arrows(0, 0.5, 1, 0.5, length = 0.08, lwd = 2, col = "grey50")
straightarrow(c(0, 0.5), c(0.92, 0.5), arr.pos = 1, lwd = 1, arr.length = 0.15, arr.width = 0.1)

# scale bar
cor20 <- seq(-1, 1, length = 19)
n.col <- length(cor20)
par(las = 1, tcl = -0.2, mar = c(1, 9, 0, 9), ps = 7, lwd = 0.5)
image((1:length(cor20))/20, c(0.72, 0.74), matrix(rep(cor20, 2), ncol = 2, byrow = F), col = rev(designer.colors(n = 20, col = brewer.pal(9, "RdBu"))), add = TRUE, useRaster = F)
rect(1/40, 0.71, 19/20 + 1/40, 0.75, lwd = 0.5)
text(0.5, 0.85, "Correlation", pos = 1, cex = 6/par("ps")/par("cex"), xpd = TRUE)
segments(1/40, 0.71, 1/40, 0.70, lwd = 0.5)
text(1/40, 0.67, "-1", cex = 6/par("ps")/par("cex"))
segments(19/20 + 1/40, 0.71, 19/20 + 1/40, 0.70, lwd = 0.5)
text(19/20 + 1/40, 0.67, "1", cex = 6/par("ps")/par("cex"))
segments(10/20, 0.71, 10/20, 0.70, lwd = 0.5)
text(10/20, 0.67, "0", cex = 6/par("ps")/par("cex"))




# old, must preserve gene order
par(las = 1, tcl = -0.2, mar = c(1, 1, 1, 1), ps = 7, lwd = 0.5)
plot(c(0, length(all.tca.gene)), c(0, length(all.tca.gene)), type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1)
image(1:length(all.tca.gene), 1:length(all.tca.gene), old.order.cor, col = rev(designer.colors(n = 101, col = brewer.pal(9, "RdBu"))), useRaster = F, add = TRUE)
# unique.module <- unique(old.order.module)
# module.start <- sapply(split(1:length(old.order.module), factor(old.order.module, levels = unique.module)), min)
# module.end <- sapply(split(1:length(old.order.module), factor(old.order.module, levels = unique.module)), max)
#
# for (i in 1:length(module.start)) {
#
#   rect(module.start[i] - 0.5, module.start[i] - 0.5, module.end[i] + 0.5, module.end[i] + 0.5, lwd = 1)
#
# }

text(1:length(all.tca.gene) + 1.5, rep(0, length(all.tca.gene)), gene.info[rownames(young.order.cor$cor), 2], font = 3, srt = 45, pos = 2, xpd = TRUE, col = ifelse(gene.info[rownames(young.order.cor$cor), 2] %in% c("kdn", "ScsbetaA"), "red", "black"), cex = 0.75)
text(0, length(all.tca.gene) + 2, "Aged", pos = 4, cex = 7/par("ps")/par("cex"), xpd = TRUE)
rect(min(module1) - 0.5, min(module1) - 0.5, max(module1) + 0.5, max(module1) + 0.5, lwd = 1)


# network graph
# ============================================================

# young

par(las = 1, tcl = -0.2, mar = c(0.5, 0.5, 0.5, 0.5), ps = 7, lwd = 0.5)
plot.igraph(young.graph, layout = old.layout, vertex.color = ifelse(blup.25c.tree[old.gene.id] == 1 & blup.18c.tree[old.gene.id] == 1, brewer.pal(12, "Set3")[5], brewer.pal(12, "Set3")[6]), vertex.size = 20, vertex.label.color = "black", vertex.frame.color = NA, edge.width = 0.5, vertex.label.font = 3, vertex.label.cex = 2.2/par("ps")/par("cex"), vertex.label.family = args[6])
text(min(par('usr')[1:2]) + 0.2*(par('usr')[2] - par('usr')[1]), max(par('usr')[3:4]), "Young", pos = 1, adj = 2, cex = 7/par("ps")/par("cex"))
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# arrow
par(las = 1, tcl = -0.2, mar = c(0, 0, 0, 0), ps = 7, lwd = 0.5)
plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "'")
#arrows(0, 0.5, 1, 0.5, length = 0.08, lwd = 2, col = "grey50")
straightarrow(c(0, 0.5), c(0.92, 0.5), arr.pos = 1, lwd = 1, arr.length = 0.15, arr.width = 0.1)

par(las = 1, tcl = -0.2, mar = c(0.5, 0.5, 0.5, 0.5), ps = 7, lwd = 0.5)
plot.igraph(old.graph, layout = old.layout, vertex.color = ifelse(blup.25c.tree[old.gene.id] == 1 & blup.18c.tree[old.gene.id] == 1, brewer.pal(12, "Set3")[5], brewer.pal(12, "Set3")[6]), vertex.size = 20, vertex.label.color = "black", vertex.frame.color = NA, edge.width = 0.5, vertex.label.font = 3, vertex.label.cex = 2.2/par("ps")/par("cex"), vertex.label.family = args[6])

text(min(par('usr')[1:2]) + 0.2*(par('usr')[2] - par('usr')[1]), max(par('usr')[3:4]), "Aged", pos = 1, adj = 2, cex = 7/par("ps")/par("cex"))

dev.off()
