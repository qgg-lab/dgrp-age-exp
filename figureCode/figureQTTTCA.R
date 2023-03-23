# ================
# = qtt analysis =
# ================

args <- commandArgs(TRUE) # args <- c("../figureData/female.phototaxis.decline.qtt.RData", "../figureData/female.fecundity.decline.qtt.RData", "../figureData/female.lifespan.qtt.RData", "../figureData/male.phototaxis.decline.qtt.RData", "../figureData/male.speed.decline.qtt.RData", "../figureData/male.endurance.decline.qtt.RData", "../figureData/male.lifespan.qtt.RData", "../figureData/tca.genes.txt", "../figureData/female.lifespan.qtt.young.gsea.example.RData", "../figureData/female.lifespan.qtt.old.gsea.example.RData")

library("RColorBrewer")
library("fields") #to use designer.colors

# prepare file
# ============================================================

file.width = 89
cairo_pdf(file = args[11], width = file.width/25.4, height = file.width/25.4*1.6, family = args[12])
layout(mat = matrix(c(1, 2), ncol = 1, byrow = TRUE), heights = c(2.5, 1))

# get TCA genes
tca.genes <- read.table(args[8], header = FALSE, sep = " ", as.is = TRUE)

# get all TCA genes present
load(args[1])
female.gene.name <- gene.name
load(args[4])
male.gene.name <- gene.name
all.gene.name <- unique(c(female.gene.name, male.gene.name))

tca.genes <- tca.genes[tca.genes[, 1] %in% all.gene.name, ]
tca.genes <- tca.genes[order(tca.genes[, 2]), ]

# positive association
# ============================================================

# young and aged expression with female lifespan

load(args[3])
tca.genes <- cbind(tca.genes, qtt.cor.young[match(tca.genes[, 1], gene.name)])
tca.genes <- cbind(tca.genes, qtt.cor.old[match(tca.genes[, 1], gene.name)])

# decline in speed in males when taking difference
load(args[5])
tca.genes <- cbind(tca.genes, qtt.cor.diff[match(tca.genes[, 1], gene.name)])

# negative association
# ============================================================

# decline in phototaxis in young males
load(args[4])
tca.genes <- cbind(tca.genes, qtt.cor.young[match(tca.genes[, 1], gene.name)])

# decline in climbing speed in young males
load(args[5])
tca.genes <- cbind(tca.genes, qtt.cor.young[match(tca.genes[, 1], gene.name)])

colnames(tca.genes)[5:9] <- c("female-young-exp-lifespan", "female-aged-exp-lifespan", "male-diff-exp-speed", "male-young-exp-phototaxis", "male-young-exp-speed")

# plot heatmap
par(las = 1, tcl = -0.2, mai = c(0.02, 0.1, 0.25, 0.15)*file.width/25.4, ps = 7, lwd = 0.5)
plot(c(0, 9.5), c(0, nrow(tca.genes)), type = "n", axes = FALSE, xlab = "", ylab = "")
text(rep(0.5, nrow(tca.genes)), nrow(tca.genes):1, tca.genes[, 2], font = 3, pos = 2, cex = 6/par("ps")/par("cex"), xpd = TRUE)
text(rep(1.2, nrow(tca.genes)), nrow(tca.genes):1 + 0.1, tca.genes[, 3], pos = 2, cex = 6/par("ps")/par("cex"))
text(rep(2, nrow(tca.genes)), nrow(tca.genes):1 + 0.1, tca.genes[, 4], pos = 2, cex = 6/par("ps")/par("cex"))

text(0.8, nrow(tca.genes) + 1.5, "GO", cex = 6/par("ps")/par("cex"), xpd = TRUE)
text(1.6, nrow(tca.genes) + 1.5, "KEGG", cex = 6/par("ps")/par("cex"), xpd = TRUE)

text(0.5, nrow(tca.genes) + 4, "TCA cycle genes", xpd = TRUE, cex = 6/par("ps")/par("cex"), font = 2)

# plot positive associations
pos.tca.cor <- t(tca.genes[, 5:7])
pos.tca.cor[is.na(pos.tca.cor)] <- 2
image(3:5, 1:nrow(tca.genes), pos.tca.cor, col = c(rev(designer.colors(n = 100, col = brewer.pal(9, "RdBu"))), "grey80"), useRaster = F, add = TRUE, breaks = c(seq(-0.32, 0.32, 0.0064), 10))

# plot negative associations
neg.tca.cor <- t(tca.genes[, 8:9])
neg.tca.cor[is.na(neg.tca.cor)] <- 2
image(7:8, 1:nrow(tca.genes), neg.tca.cor, col = c(rev(designer.colors(n = 100, col = brewer.pal(9, "RdBu"))), "grey80"), useRaster = F, add = TRUE, breaks = c(seq(-0.32, 0.32, 0.0064), 10))

# text
par(lheight = 0.9)

text(3, 45, "\u2640", xpd = TRUE, srt = 45, pos = 3, cex = 6/par("ps")/par("cex"), family = "Arial")
text(4.5, 51.2, "Young expression vs lifespan", xpd = TRUE, srt = 45, pos = 3, cex = 6/par("ps")/par("cex"))
text(4, 45, "\u2640", xpd = TRUE, srt = 45, pos = 3, cex = 6/par("ps")/par("cex"), family = "Arial")
text(5.45, 51.1, "Aged expression vs lifespan", xpd = TRUE, srt = 45, pos = 3, cex = 6/par("ps")/par("cex"))
text(5, 45, "\u2642", xpd = TRUE, srt = 45, pos = 3, cex = 6/par("ps")/par("cex"), family = "Arial")
text(6.84, 50.45, "Aged - Young expression vs   \n climbing speed decline", xpd = TRUE, srt = 45, pos = 3, cex = 6/par("ps")/par("cex"))

text(7, 45, "\u2642", xpd = TRUE, srt = 45, pos = 3, cex = 6/par("ps")/par("cex"), family = "Arial")
text(8.5, 49.8, "Young expression vs       \n phototaxis decline", xpd = TRUE, srt = 45, pos = 3, cex = 6/par("ps")/par("cex"))

text(8, 45, "\u2642", xpd = TRUE, srt = 45, pos = 3, cex = 6/par("ps")/par("cex"), family = "Arial")
text(9.75, 50.6, "Young expression vs                \n climbing speed decline", xpd = TRUE, srt = 45, pos = 3, cex = 6/par("ps")/par("cex"))


# plot correlation bar
# par(las = 1, tcl = -0.2, mai = c(0.3, 0.05, 0.25, 0.02)*file.width/25.4, ps = 7, lwd = 0.5)

cor100 <- seq(-0.32, 0.32, by = 0.0064)
n.col <- length(cor100)
#plot(c(9, 10), c(0, length(cor100) + 1), type = "n", axes = FALSE, xlab = "", ylab = "")
image(c(9.1, 9.3), seq(10, 30, by = 0.2), matrix(rep(cor100, 2), nrow = 2, byrow = T), col = rev(designer.colors(n = 101, col = brewer.pal(9, "RdBu"))), add = TRUE, useRaster = T, xpd = T)
rect(9, 9.9, 9.4, 30.1, lwd = 0.5)
segments(9.4, c(9.9, 20, 30.1), 9.5, c(9.9, 20, 30.1))
text(9.3, c(9.9, 20, 30.1), c(-0.32, 0, 0.32), cex = 6/par("ps")/par("cex"), pos = 4, xpd = T)
text(9.5, 30.2, "Correlation", cex = 6/par("ps")/par("cex"), pos = 3, xpd = T)

image(c(9.1, 9.3), seq(5, 6, by = 1), matrix(rep(2, 4), nrow = 2, byrow = T), col = c(rev(designer.colors(n = 100, col = brewer.pal(9, "RdBu"))), "grey80"), useRaster = T, add = TRUE, breaks = c(seq(-0.32, 0.32, 0.0064), 10), xpd = TRUE)
rect(9, 4.5, 9.4, 6.5, lwd = 0.5)
segments(9.4, 5.5, 9.5, 5.5)
text(9.3, 5.5, "missing", cex = 6/par("ps")/par("cex"), pos = 4, xpd = T)

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# gsea plot
# ============================================================

load(args[9])
go.young <- gsea.example[[1]]

load(args[10])
go.old <- gsea.example[[1]]

par(las = 1, tcl = -0.2, mai = c(0.20, 0.25, 0.05, 0.1)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = T)

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
text(2000, 0.34, "Aged", col = brewer.pal(9, "Reds")[9], pos = 3, cex = 7/par("ps")/par("cex"))
text(6000, 0.6, "GO: TCA cycle", cex = 7/par("ps")/par("cex"))
title(ylab = "Enrichment score", mgp = c(1, 0.3, 0), cex.lab = 8/par("ps")/par("cex"))
title(xlab = expression(paste("" + "", "" %<-% "", "Rank of correlation with lifespan", "" %->% "", "" - "")), mgp = c(0.5, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(4800, 0.6, "\u2640", cex = 15/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], family = "Arial")

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

dev.off()
