# =======================================
# = make figure to summarize QG results =
# =======================================

args <- commandArgs(TRUE) # args = c("../figureData/female.adj.qgGxE.RData", "../figureData/male.adj.qgGxE.RData", "../figureData/female.adj.qgVarHet.RData", "../figureData/female.adj.qgSingleTemp.RData", "../figureData/male.adj.qgVarHet.RData", "../figureData/male.adj.qgSingleTemp.RData", "../figure/figureGxE.pdf", "Myriad Pro")
library("RColorBrewer")

# prepare file to plot
# ============================================================

file.width = 89
cairo_pdf(file = args[7], width = file.width/25.4, height = file.width/25.4 * 0.9, family = args[8])
par(las = 1, tcl = -0.2, mai = c(0.15, 0.25, 0.1, 0.02)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = F, mfrow = c(2, 2))

# load data
# ============================================================

load(args[1]); female.gxe <- gxe;
load(args[2]); male.gxe <- gxe;

# calculate gxa and fdr
# ============================================================

female.gxa <- female.gxe[, 9]/(female.gxe[, 9] + female.gxe[, 10])
female.gxa.fdr <- p.adjust(female.gxe[, 13], method = "BH")
female.g.fdr <- p.adjust(female.gxe[, 12], method = "BH")

male.gxa <- male.gxe[, 9]/(male.gxe[, 9] + male.gxe[, 10])
male.gxa.fdr <- p.adjust(male.gxe[, 13], method = "BH")
male.g.fdr <- p.adjust(male.gxe[, 12], method = "BH")

# count number of genes
cat("female number of genes with significant gxa effect: ", sum(p.adjust(female.gxe[, 13], method = "BH") <= 0.05), " out of ", nrow(female.gxe), "\n", sep = "")
cat("male number of genes with significant gxa effect: ", sum(p.adjust(male.gxe[, 13], method = "BH") <= 0.05), " out of ", nrow(male.gxe), "\n", sep = "")

cat("female number of genes with significant g effect: ", sum(p.adjust(female.gxe[, 12], method = "BH") <= 0.05), " out of ", nrow(female.gxe), "\n", sep = "")
cat("male number of genes with significant g effect: ", sum(p.adjust(male.gxe[, 12], method = "BH") <= 0.05), " out of ", nrow(male.gxe), "\n", sep = "")

cat("female number of genes with significant g or gxe effect: ", sum(p.adjust(female.gxe[, 12], method = "BH") <= 0.05 | p.adjust(female.gxe[, 13], method = "BH") <= 0.05), " out of ", nrow(female.gxe), "\n", sep = "")
cat("male number of genes with significant g or gxe effect: ", sum(p.adjust(male.gxe[, 12], method = "BH") <= 0.05 | p.adjust(male.gxe[, 13], method = "BH") <= 0.0), " out of ", nrow(male.gxe), "\n", sep = "")

# female histogram
# ============================================================

female.hist.plot <- hist(female.gxa[female.gxa.fdr <= 0.05 | female.g.fdr <= 0.05], breaks = seq(0, 1, 0.05), grey = "grey80", axes = FALSE, xlab = "", ylab = "", main = "", ylim = c(0, 1400))
hist(female.gxa[female.gxa.fdr <= 0.05], breaks = seq(0, 1, 0.05), col = brewer.pal(9, "Reds")[9], add = TRUE)
axis(side = 1, at = seq(0, 1, 0.2), mgp = c(0.8, -0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = expression(italic("GxA")), mgp = c(0.7, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Number of genes", mgp = c(1.6, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
legend(-0.1, 1600, pch = 22, pt.bg = c(brewer.pal(9, "Reds")[9]),
       legend = c(expression(paste("FDR (", italic("GxA"), " > 0) = 0.05"))), bty = "n",
       x.intersp = 0.5, y.intersp = 0.6, cex = 7/par("ps")/par("cex"))
box(bty = "l", lwd = 0.5)
axis(side = 2, mgp = c(2, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
text(0.7, 400, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 4, family = "Arial")
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# male histogram
# ============================================================

male.hist.plot <- hist(male.gxa[male.gxa.fdr <= 0.05 | male.g.fdr <= 0.05], breaks = seq(0, 1, 0.05), grey = "grey80", axes = FALSE, xlab = "", ylab = "", main = "", ylim = c(0, 1400))
hist(male.gxa[male.gxa.fdr <= 0.05], breaks = seq(0, 1, 0.05), col = brewer.pal(9, "Blues")[9], add = TRUE)
axis(side = 1, at = seq(0, 1, 0.2), mgp = c(0.8, -0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = expression(italic("GxA")), mgp = c(0.7, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Number of genes", mgp = c(1.6, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
legend(-0.1, 1600, pch = 22, pt.bg = c(brewer.pal(9, "Blues")[9]),
       legend = c(expression(paste("FDR (", italic("GxA"), " > 0) = 0.05"))), bty = "n",
       x.intersp = 0.5, y.intersp = 0.6, cex = 7/par("ps")/par("cex"))
box(bty = "l", lwd = 0.5)
axis(side = 2, mgp = c(2, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
text(0.7, 400, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 4, family = "Arial")
text(grconvertX(0.05 + file.width/2/25.4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


par(mai = c(0.16, 0.2, 0.1, 0.02)*file.width/25.4/2)

# female variance het
# ============================================================

load(args[3])
rownames(var.het) <- gene.name; female.var.het <- var.het
load(args[4])
rownames(single.temp) <- gene.name; female.single.temp <- single.temp

female.fdr.single.g <- p.adjust(as.numeric(female.var.het[, 12]), method = "BH")
female.fdr.single.e <- p.adjust(as.numeric(female.var.het[, 11]), method = "BH")

female.fdr.young.g <- p.adjust(as.numeric(female.single.temp[, 12]), method = "BH")
female.fdr.old.g <- p.adjust(as.numeric(female.single.temp[, 9]), method = "BH")

sig.female.g <- female.fdr.single.g < 0.05 & ( ((as.numeric(female.var.het[, 7])/as.numeric(female.var.het[, 8]) > 2) & (female.fdr.young.g < 0.05)) | ((as.numeric(female.var.het[, 8])/as.numeric(female.var.het[, 7]) > 2) & (female.fdr.old.g < 0.05)) )

plot(c(0, 3.2), c(0, 3.2), axes = FALSE, type = "n", xlab = "", ylab = "")

points(sqrt(as.numeric(female.var.het[!sig.female.g, 7])), sqrt(as.numeric(female.var.het[!sig.female.g, 8])), col = rgb(col2rgb("grey50")[1, 1], col2rgb("grey50")[2, 1], col2rgb("grey50")[3, 1], 40, max = 255), pch = 16, cex = 0.5)

segments(-0.5, -0.5, 2.5, 2.5, lty = 2, col = "grey20")

points(sqrt(as.numeric(female.var.het[sig.female.g, 7])), sqrt(as.numeric(female.var.het[sig.female.g, 8])), col = rgb(col2rgb(brewer.pal(9, "Reds")[9])[1, 1], col2rgb(brewer.pal(9, "Reds")[9])[2, 1], col2rgb(brewer.pal(9, "Reds")[9])[3, 1], 120, max = 255), pch = 16, cex = 0.5)

lines(c(0, 3), c(0, 3)*sqrt(2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Reds")[9])
lines(c(0, 3), c(0, 3)*sqrt(1/2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Reds")[9])

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")

title(xlab = expression(paste(italic(sigma[G]), " in young flies")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(sigma[G]), " in aged flies")), mgp = c(0.6, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(2.8, 2.4, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 3, family = "Arial")
text(2, 0.2, expression(paste(italic(sigma[G]^2), " (young/aged)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])
text(1, 2.9, expression(paste(italic(sigma[G]^2), " (aged/young)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

cat("In females, there are", sum(sig.female.g, na.rm = T), "genes with variance het.\n")
cat("among which,", sum(sig.female.g & (as.numeric(female.var.het[, 7])/as.numeric(female.var.het[, 8]) > 2), na.rm = T), "higher var in young.\n")
cat("among which,", sum(sig.female.g & (as.numeric(female.var.het[, 7])/as.numeric(female.var.het[, 8]) < 1/2), na.rm = T), "higher var in old.\n")

# male
load(args[5])
rownames(var.het) <- gene.name; male.var.het <- var.het
load(args[6])
rownames(single.temp) <- gene.name; male.single.temp <- single.temp

male.fdr.single.g <- p.adjust(as.numeric(male.var.het[, 12]), method = "BH")
male.fdr.single.e <- p.adjust(as.numeric(male.var.het[, 11]), method = "BH")

male.fdr.young.g <- p.adjust(as.numeric(male.single.temp[, 12]), method = "BH")
male.fdr.old.g <- p.adjust(as.numeric(male.single.temp[, 9]), method = "BH")

sig.male.g <- male.fdr.single.g < 0.05 & ( ((as.numeric(male.var.het[, 7])/as.numeric(male.var.het[, 8]) > 2) & (male.fdr.young.g < 0.05)) | ((as.numeric(male.var.het[, 8])/as.numeric(male.var.het[, 7]) > 2) & (male.fdr.old.g < 0.05)) )

plot(c(0, 3.2), c(0, 3.2), axes = FALSE, type = "n", xlab = "", ylab = "")

points(sqrt(as.numeric(male.var.het[!sig.male.g, 7])), sqrt(as.numeric(male.var.het[!sig.male.g, 8])), col = rgb(col2rgb("grey50")[1, 1], col2rgb("grey50")[2, 1], col2rgb("grey50")[3, 1], 40, max = 255), pch = 16, cex = 0.5)

segments(-0.5, -0.5, 2.5, 2.5, lty = 2, col = "grey20")

points(sqrt(as.numeric(male.var.het[sig.male.g, 7])), sqrt(as.numeric(male.var.het[sig.male.g, 8])), col = rgb(col2rgb(brewer.pal(9, "Blues")[9])[1, 1], col2rgb(brewer.pal(9, "Blues")[9])[2, 1], col2rgb(brewer.pal(9, "Blues")[9])[3, 1], 120, max = 255), pch = 16, cex = 0.5)

lines(c(0, 3), c(0, 3)*sqrt(2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Blues")[9])
lines(c(0, 3), c(0, 3)*sqrt(1/2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Blues")[9])

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 3, 1), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")

title(xlab = expression(paste(italic(sigma[G]), " in young flies")), mgp = c(0.8, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(sigma[G]), " in aged flies")), mgp = c(0.6, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(2.8, 2.4, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 3, family = "Arial")
text(2, 0.2, expression(paste(italic(sigma[G]^2), " (young/aged)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])
text(1, 2.9, expression(paste(italic(sigma[G]^2), " (aged/young)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])

text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

cat("In males, there are", sum(sig.male.g, na.rm = T), "genes with variance het.\n")
cat("among which,", sum(sig.male.g & (as.numeric(male.var.het[, 7])/as.numeric(male.var.het[, 8]) > 2), na.rm = T), "higher var in young.\n")
cat("among which,", sum(sig.male.g & (as.numeric(male.var.het[, 7])/as.numeric(male.var.het[, 8]) < 1/2), na.rm = T), "higher var in old.\n")


dev.off()


