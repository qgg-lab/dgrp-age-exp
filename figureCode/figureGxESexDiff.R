# =======================================
# = make figure to summarize QG results =
# =======================================

args <- commandArgs(TRUE) # args = c("../figureData/female.adj.qgGxE.RData", "../figureData/male.adj.qgGxE.RData", "Myriad Pro")
library("RColorBrewer")

# prepare file to plot
# ============================================================

file.width = 89
cairo_pdf(file = args[3], width = file.width/25.4/2, height = file.width/25.4/2, family = args[4])
par(las = 1, tcl = -0.2, mai = c(0.2, 0.2, 0.1, 0.02)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = T)

# load data
# ============================================================

load(args[1]); female.gxe <- gxe; rownames(female.gxe) <- gene.name
load(args[2]); male.gxe <- gxe; rownames(male.gxe) <- gene.name

# calculate gxa and fdr
# ============================================================

female.gxa <- female.gxe[, 9]/(female.gxe[, 9] + female.gxe[, 10])
female.gxa.fdr <- p.adjust(female.gxe[, 13], method = "BH")
female.g.fdr <- p.adjust(female.gxe[, 12], method = "BH")
female.gxa <- female.gxa[female.gxa.fdr <= 0.05 | female.g.fdr <= 0.05]

male.gxa <- male.gxe[, 9]/(male.gxe[, 9] + male.gxe[, 10])
male.gxa.fdr <- p.adjust(male.gxe[, 13], method = "BH")
male.g.fdr <- p.adjust(male.gxe[, 12], method = "BH")
male.gxa <- male.gxa[male.gxa.fdr <= 0.05 | male.g.fdr <= 0.05]

common.gene <- intersect(names(female.gxa), names(male.gxa))

# scatter plot
# ============================================================

plot(c(0, 1), c(0, 1), axes = FALSE, type = "n", xlab = "", ylab = "")
abline(a = 0, b = 1, xpd = F)
points(female.gxa[common.gene], male.gxa[common.gene], col = brewer.pal(9, "Purples")[9], cex = 0.5)
title(xlab = expression(paste(italic("GxA"), " in females")), mgp = c(0.6, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic("GxA"), " in males")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
eff.cor <- formatC(cor.test(female.gxa[common.gene], male.gxa[common.gene], method = "spearman")$estimate, format = "f", digits = 2)
text(0.1, 1.1, parse(text = paste("paste(italic(\"r\"), \" = \", ", eff.cor, ")")), cex = 7/par("ps")/par("cex"))

dev.off()


