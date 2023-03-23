# =======================================
# = make figure TCA cycle GSEA in QTT   =
# =======================================

args <- commandArgs(TRUE) # args = c("../figureData/female.young.qtt.gsea.example.RData", "../figureData/female.old.qtt.gsea.example.RData", "../figureData/female.qtt.RData", "../figureData/gene_go.table", "../figureData/fbgn.kegg")
library("RColorBrewer")

# prepare file to plot
# ============================================================

file.width = 150
cairo_pdf(file = args[6], width = file.width/25.4, height = file.width/25.4*0.35, family = args[7])
layout(matrix(c(1, 2), ncol = 2, byrow = T), widths = c(1, 0.5))

# load data for young
# ============================================================

load(args[1])
go.young <- gsea.example[[1]]

load(args[2])
go.old <- gsea.example[[1]]

par(las = 1, tcl = -0.2, mai = c(0.20, 0.20, 0.05, 0.1)*file.width/25.4/2/1.5, ps = 7, lwd = 0.5, xpd = T)

plot(c(0, length(go.young$p)), c(-0.2, 0.8), type = "n", axes = FALSE, xlab = "", ylab = "")
points((1:length(go.young$p)), cumsum(go.young$p), type = "s", lwd = 1, col = brewer.pal(9, "Blues")[9])
points((1:length(go.old$p)), cumsum(go.old$p), type = "s", lwd = 1, col = brewer.pal(9, "Reds")[9])

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), c(1, 4000, 7041), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 0.8, 0.2), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
segments(0, -0.12, length(go.young$p), -0.12, col = brewer.pal(9, "Blues")[9])
segments(0, -0.15, length(go.old$p), -0.15, col = brewer.pal(9, "Reds")[9])
segments(go.young$S, -0.12, go.young$S, -0.10, col = brewer.pal(9, "Blues")[9])
segments(go.old$S, -0.15, go.old$S, -0.17, col = brewer.pal(9, "Reds")[9])
text(2000, 0.6, "Young", col = brewer.pal(9, "Blues")[9], pos = 3, cex = 7/par("ps")/par("cex"))
text(2000, 0.34, "Old", col = brewer.pal(9, "Reds")[9], pos = 3, cex = 7/par("ps")/par("cex"))
text(6000, 0.6, "GO: TCA cycle", cex = 7/par("ps")/par("cex"))
title(ylab = "Enrichment score", mgp = c(1, 0.3, 0), cex.lab = 8/par("ps")/par("cex"))
title(xlab = expression(paste("" + "", "" %<-% "", "Rank of correlation with lifespan", "" %->% "", "" - "")), mgp = c(0.5, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(2000, 0.15, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], family = "Arial")

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# cumulative plots
# ============================================================

cumfrac <- function(x, from, to, int, genes) {
	
	x.cut <- cut(x, breaks = seq(from, to, int))
	frac <- c(0, cumsum(table(x.cut))/length(x))
	x.pos <- seq(from, to, int)
  
  int.genes <- intersect(genes, names(x))
  mark <- matrix(ncol = 2, nrow = length(int.genes))
  
  for (i in 1:length(int.genes)) {
    
    mark[i, 1] <- x[int.genes[i]]
    mark[i, 2] <- sum(x <= x[int.genes[i]])/length(x)
    
  }
  
  
	return(list(x = x.pos, frac = frac, mark = mark))
	
}

go.id <- "GO:0006099"
kegg.id <- "dme00020"

annot <- read.table(args[4], header = TRUE, as.is = TRUE, sep = " ") # annotation table

# read kegg table

kegg <- read.table(args[5], header = FALSE, as.is = TRUE)

# find genes in go and kegg
# ============================================================

bp.go.genes <- annot[grepl(go.id, annot[, 3]), 1]
kegg.genes <- kegg[kegg[, 2] == kegg.id, 1]

# both genes
# ============================================================

both.genes <- sort(unique(c(bp.go.genes, kegg.genes)))

# cumulative
# ============================================================

load(args[3])

names(qtt.cor.young) <- gene.name
names(qtt.cor.old) <- gene.name
qtt.young.frac <- cumfrac(qtt.cor.young, -0.36, 0.36, 0.01, both.genes)
qtt.old.frac <- cumfrac(qtt.cor.old, -0.36, 0.36, 0.01, both.genes)


plot(c(-0.36, 0.36), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")

abline(h = 0.5, lwd = 1, col = "grey80", xpd = F)
abline(v = 0, lwd = 1, col = "grey80", xpd = F)

lines(qtt.young.frac$x, qtt.young.frac$frac, lwd = 1, col = brewer.pal(9, "Set1")[2])
lines(qtt.old.frac$x, qtt.old.frac$frac, lwd = 1, col = brewer.pal(9, "Set1")[1])

points(qtt.young.frac$mark, pch = 3, col = brewer.pal(9, "Set1")[2], cex = 0.3)
points(qtt.old.frac$mark, pch = 3, col = brewer.pal(9, "Set1")[1], cex = 0.3)

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), seq(-0.3, 0.3, 0.2), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = "Correlation with lifespan", mgp = c(0.5, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Cumulative fraction", mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(0.2, 0.3, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], family = "Arial")
legend("topleft", lwd = 1, col = brewer.pal(9, "Set1")[c(2, 1)],
       legend = c("Young", "Old"), bty = "n", seg.len = 0.8,
       x.intersp = 0.5, y.intersp = 0.5, cex = 6/par("ps")/par("cex"))

text(grconvertX(file.width/25.4/3*2 + 0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

dev.off()


