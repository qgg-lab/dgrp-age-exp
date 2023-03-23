# =======================================
# = make figure to summarize QG results =
# =======================================

args <- commandArgs(TRUE) # args = c("../figureData/female.adj.qgAgeEff.RData", "../figureData/male.adj.qgAgeEff.RData", "../figureData/gene.exp.mean.RData", "../figureData/female.age.effect.gsea.example.RData", "../figureData/male.age.effect.gsea.example.RData", "../figure/figureAgingSumamry.pdf", "Myriad Pro")
library("RColorBrewer")
library("diagram")

# prepare file to plot
# ============================================================

file.width = 89
cairo_pdf(file = args[6], width = file.width/25.4, height = file.width/25.4*1.5, family = args[7])
layout(matrix(c(1, 3, 2, 4, 5, 6, 5, 7), byrow = T, ncol = 2, nrow = 4), widths = c(0.5, 0.5), heights = c(0.5, 0.5, 0.45*0.5, 0.55*0.5))

# load data for female
# ============================================================

load(args[1]); load(args[3]);

# count number of genes

cat("female expressed NTRs: ", sum(grepl("XLOC", gene.name)), "\n", sep = "")
cat("female number of genes with significant age effect and FC > 2: ", sum(p.adjust(age.effect[, 6], method = "BH") <= 0.05 & abs(age.effect[, 4]) >= 1), " out of ", nrow(age.effect), "\n", sep = "")
cat("of these, number of NTRs: ", sum(p.adjust(age.effect[, 6], method = "BH") <= 0.05 & abs(age.effect[, 4]) >= 1 & grepl("XLOC", gene.name)), "\n", sep = "")

cat("female number of genes with significant age effect: ", sum(p.adjust(age.effect[, 6], method = "BH") <= 0.05), " out of ", nrow(age.effect), "\n", sep = "")
cat("of these, number of NTRs: ", sum(p.adjust(age.effect[, 6], method = "BH") <= 0.05 & grepl("XLOC", gene.name)), "\n", sep = "")

# plot histogram of effect
# cap age effect
age.effect.cap <- age.effect[, 4]
age.effect.cap[age.effect.cap < -5] <- -5.01
age.effect.cap[age.effect.cap > 5] <- 5.1

par(las = 1, tcl = -0.2, mai = c(0.15, 0.2, 0.1, 0.1)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = F)

hist(age.effect.cap, ylim = c(0, 2500), col = "white", axes = FALSE, xlab = "", ylab = "", main = "", breaks = seq(-5.1, 5.1, 0.1), lwd = 0.1)
hist(age.effect.cap[p.adjust(age.effect[, 6], method = "BH") <= 0.05], breaks = seq(-5.1, 5.1, 0.1), lwd = 0.1, col = brewer.pal(9, "Reds")[9], add = TRUE)
axis(side = 1, mgp = c(0.8, 0.05, 0), at = c(-5, -2, 0, 2, 5), label = c("< -5", -2, 0, 2, "> 5"), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(2, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = expression(paste(Delta, "log", ""[2], "TPM (aged - young)")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Number of genes", mgp = c(1.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
legend(-6, 2800, pch = 22, pt.bg = c(brewer.pal(9, "Reds")[9]),
       legend = c("FDR = 0.05"), bty = "n",
       x.intersp = 0.5, y.intersp = 0.6, cex = 7/par("ps")/par("cex"))
box(bty = "l", lwd = 0.5)
text(2, 1500, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 4, family = "Arial")
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# plot female old versus young

plot(c(-5, 15), c(-5, 15), type = "n", axes = FALSE)
abline(a = 0, b = 1)
points(female.young.mean[p.adjust(age.effect[, 6], method = "BH") > 0.05 | abs(age.effect[, 4]) < 1], female.old.mean[p.adjust(age.effect[, 6], method = "BH") > 0.05 | abs(age.effect[, 4]) < 1], col = "grey60", cex = 0.5)
points(female.young.mean[p.adjust(age.effect[, 6], method = "BH") <= 0.05 & abs(age.effect[, 4]) >= 1], female.old.mean[p.adjust(age.effect[, 6], method = "BH") <= 0.05 & abs(age.effect[, 4]) >= 1], col = brewer.pal(9, "Reds")[9], cex = 0.5)

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = seq(-5, 15, 5), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(-5, 15, 5), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = expression(paste("Expression in young flies (log", ""[2], "TPM)")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste("Expression in aged flies (log", ""[2], "TPM)")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(12, -2, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 3, family = "Arial")

legend(-7, 17, pch = 1, col = brewer.pal(9, "Reds")[9], x.intersp = 0.5, legend = expression(paste("FDR = 0.05\nFold change" >= "2")), bty = "n", cex = 7/par("ps")/par("cex"), adj = c(0, 1))

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# load data for male
# ============================================================

load(args[2])

# count number of genes
cat("male expressed NTRs: ", sum(grepl("XLOC", gene.name)), "\n", sep = "")
cat("male number of genes with significant age effect and FC > 2: ", sum(p.adjust(age.effect[, 6], method = "BH") <= 0.05 & abs(age.effect[, 4]) >= 1), " out of ", nrow(age.effect), "\n", sep = "")
cat("of these, number of NTRs: ", sum(p.adjust(age.effect[, 6], method = "BH") <= 0.05 & abs(age.effect[, 4]) >= 1 & grepl("XLOC", gene.name)), "\n", sep = "")

cat("male number of genes with significant age effect: ", sum(p.adjust(age.effect[, 6], method = "BH") <= 0.05), " out of ", nrow(age.effect), "\n", sep = "")
cat("of these, number of NTRs: ", sum(p.adjust(age.effect[, 6], method = "BH") <= 0.05 & grepl("XLOC", gene.name)), "\n", sep = "")

# plot histogram of effect
# cap age effect
age.effect.cap <- age.effect[, 4]
age.effect.cap[age.effect.cap < -5] <- -5.01
age.effect.cap[age.effect.cap > 5] <- 5.1

par(las = 1, tcl = -0.2, mai = c(0.15, 0.2, 0.1, 0.1)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = F)

hist(age.effect.cap, ylim = c(0, 2500), col = "white", axes = FALSE, xlab = "", ylab = "", main = "", breaks = seq(-5.1, 5.1, 0.1), lwd = 0.1)
hist(age.effect.cap[p.adjust(age.effect[, 6], method = "BH") <= 0.05], breaks = seq(-5.1, 5.1, 0.1), lwd = 0.1, col = brewer.pal(9, "Blues")[9], add = TRUE)
axis(side = 1, mgp = c(0.8, 0.05, 0), at = c(-5, -2, 0, 2, 5), label = c("< -5", -2, 0, 2, "> 5"), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(2, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = expression(paste(Delta, "log", ""[2], "TPM (aged - young)")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Number of genes", mgp = c(1.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
legend(-6, 2800, pch = 22, pt.bg = c(brewer.pal(9, "Blues")[9]),
       legend = c("FDR = 0.05"), bty = "n",
       x.intersp = 0.5, y.intersp = 0.6, cex = 7/par("ps")/par("cex"))
box(bty = "l", lwd = 0.5)
text(2, 1500, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 4, family = "Arial")
text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)



# plot male old versus young
# par(las = 1, tcl = -0.2, mai = c(0.15, 0.16, 0.1, 0.1)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = F)

plot(c(-5, 15), c(-5, 15), type = "n", axes = FALSE)
abline(a = 0, b = 1)
points(male.young.mean[p.adjust(age.effect[, 6], method = "BH") > 0.05 | abs(age.effect[, 4]) < 1], male.old.mean[p.adjust(age.effect[, 6], method = "BH") > 0.05 | abs(age.effect[, 4]) < 1], col = "grey60", cex = 0.5)
points(male.young.mean[p.adjust(age.effect[, 6], method = "BH") <= 0.05 & abs(age.effect[, 4]) >= 1], male.old.mean[p.adjust(age.effect[, 6], method = "BH") <= 0.05 & abs(age.effect[, 4]) >= 1], col = brewer.pal(9, "Blues")[9], cex = 0.5)

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), at = seq(-5, 15, 5), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(-5, 15, 5), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = expression(paste("Expression in young flies (log", ""[2], "TPM)")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste("Expression in aged flies (log", ""[2], "TPM)")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(12, -2, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 3, family = "Arial")

legend(-7, 17, pch = 1, col = brewer.pal(9, "Blues")[9], x.intersp = 0.5, legend = expression(paste("FDR = 0.05\nFold change" >= "2")), bty = "n", cex = 7/par("ps")/par("cex"), adj = c(0, 1))

text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# female versus male common aging effect
# ============================================================

load(args[1]); female.eff <- age.effect[, 4]; names(female.eff) <- gene.name
load(args[2]); male.eff <- age.effect[, 4]; names(male.eff) <- gene.name

# find genes in TCA
load(args[4])

common.gene <- intersect(names(female.eff), names(male.eff))
#par(las = 1, tcl = -0.2, mai = c(0.15, 0.16, 0.1, 0.1)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = F)

tca.gene <- intersect(common.gene, unique(c(gsea.example[[1]]$gene.set, gsea.example[[3]]$gene.set)))

plot(c(-7, 4), c(-7, 4), type = "n", axes = FALSE)
abline(a = 0, b = 1)
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
points(female.eff[common.gene], male.eff[common.gene], col = brewer.pal(9, "Set1")[9], cex = 0.5)
points(female.eff[tca.gene], male.eff[tca.gene], col = brewer.pal(9, "Purples")[9], cex = 0.5)
title(xlab = expression(paste(Delta, "log", ""[2], "TPM (aged - young) in females")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(Delta, "log", ""[2], "TPM (aged - young) in males")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
eff.cor <- formatC(cor.test(female.eff[common.gene], male.eff[common.gene], method = "spearman")$estimate, format = "f", digits = 2)
text(-5, 4, parse(text = paste("paste(italic(\"r\"), \" = \", ", eff.cor, ")")), cex = 7/par("ps")/par("cex"))
legend("topleft", pch = 1, legend = c("TCA cycle genes"), col = brewer.pal(9, "Purples")[9], bty = "n", cex = 7/par("ps")/par("cex"), x.intersp = 0.5, pt.cex = 1)

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("e")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# gsea example, there are two of them BP> GO:0006120, mitochondrial electron transport, NADH to ubiquinone
# and BP> GO:0007606 sensory perception of chemical stimulus
# ============================================================

# GSEA in females
load(args[4])
go1.female <- gsea.example[[1]]
go2.female <- gsea.example[[2]]
load(args[5])
go1.male <- gsea.example[[1]]
go2.male <- gsea.example[[2]]
par(las = 1, tcl = -0.2, mai = c(0.05, 0.16, 0.01, 0.1)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = T)

plot(c(0, length(go1.male$p)), c(-1.2, 0.4), type = "n", axes = FALSE, xlab = "", ylab = "")
points((1:length(go1.female$p)), cumsum(go1.female$p), type = "s", lwd = 1, col = brewer.pal(9, "Reds")[9])
points(1:length(go1.male$p), cumsum(go1.male$p), type = "s", lwd = 1, col = brewer.pal(9, "Blues")[9])
axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), seq(0, 8000, 4000), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
segments(0, -0.95, length(go1.female$p), -0.95, col = brewer.pal(9, "Reds")[9])
segments(0, -1, length(go1.male$p), -1, col = brewer.pal(9, "Blues")[9])
segments(go1.female$S, -0.95, go1.female$S, -0.90, col = brewer.pal(9, "Reds")[9])
segments(go1.male$S, -1, go1.male$S, -1.05, col = brewer.pal(9, "Blues")[9])
text(length(go1.male$p), 0, "\u2642", cex = 8/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 3, family = "Arial")
text(length(go1.female$p), 0, "\u2640", cex = 8/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 3, family = "Arial")
text(3000, 0.2, "GO: TCA cycle", cex = 7/par("ps")/par("cex"))

text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("f")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)
text(-2100, -1.1, "Score", srt = 90, cex = 7/par("ps")/par("cex"))



par(las = 1, tcl = -0.2, mai = c(0.15, 0.16, 0.04, 0.1)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = T)

plot(c(0, length(go2.male$p)), c(-0.2, 0.6), type = "n", axes = FALSE, xlab = "", ylab = "")
points((1:length(go2.female$p)), cumsum(go2.female$p), type = "s", lwd = 1, col = brewer.pal(9, "Reds")[9])
points(1:length(go2.male$p), cumsum(go2.male$p), type = "s", lwd = 1, col = brewer.pal(9, "Blues")[9])
axis(side = 1, lwd = 0.5, mgp = c(0.8, 0.1, 0), seq(0, 1200, 300), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), seq(0, 0.6, 0.3), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
segments(0, -0.1, length(go2.female$p), -0.1, col = brewer.pal(9, "Reds")[9])
segments(0, -0.13, length(go2.male$p), -0.13, col = brewer.pal(9, "Blues")[9])
segments(go2.female$S, -0.1, go2.female$S, -0.08, col = brewer.pal(9, "Reds")[9])
segments(go2.male$S, -0.13, go2.male$S, -0.15, col = brewer.pal(9, "Blues")[9])
text(length(go2.male$p), 0, "\u2642", cex = 8/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 3, family = "Arial")
text(length(go2.female$p), 0, "\u2640", cex = 8/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 3, family = "Arial")
text(800, 0.5, "KEGG: Longevity\nregulating pathway", cex = 7/par("ps")/par("cex"))

text(-2100, 0, "Enrichment", srt = 90, cex = 7/par("ps")/par("cex"))

title(xlab = expression(paste("" + "", "" %<-% "", Delta, "log", ""[2], "TPM rank", "" %->% "", "" - "")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))


dev.off()


