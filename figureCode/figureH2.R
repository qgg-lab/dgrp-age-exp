# =======================================
# = make figure to summarize QG results =
# =======================================

args <- commandArgs(TRUE) # args = c("../figureData/female.adj.qgSingleTemp.RData", "../figureData/male.adj.qgSingleTemp.RData", "../figure/figureH2.pdf", "Myriad Pro")
library("RColorBrewer")

# prepare file to plot
# ============================================================

file.width = 89
cairo_pdf(file = args[3], width = file.width/25.4 * 1.5, height = file.width/25.4 * 0.9, family = args[4])
par(las = 1, tcl = -0.2, mai = c(0.15, 0.25, 0.1, 0.02)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = F, mfrow = c(2, 3))

# load data
# ============================================================

load(args[1]); female.var <- single.temp;
load(args[2]); male.var <- single.temp;

# calculate H2 and fdr
# ============================================================

female.H2.old <- female.var[, 7]/(female.var[, 7] + female.var[, 8])
female.H2.old.fdr <- p.adjust(female.var[, 9], method = "BH")

female.H2.young <- female.var[, 10]/(female.var[, 10] + female.var[, 11])
female.H2.young.fdr <- p.adjust(female.var[, 12], method = "BH")

# female young histogram
# ============================================================

female.hist.plot <- hist(female.H2.young, breaks = seq(0, 1, 0.05), col = "grey80", axes = FALSE, xlab = "", ylab = "", main = "", ylim = c(0, 2000))
hist(female.H2.young[female.H2.young.fdr <= 0.05], breaks = seq(0, 1, 0.05), col = brewer.pal(9, "Reds")[9], add = TRUE)
axis(side = 1, at = seq(0, 1, 0.2), mgp = c(0.8, -0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = expression(paste(italic(H^2), " in young flies")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Number of genes", mgp = c(1.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
legend(-0.1, 2250, pch = 22, pt.bg = c(brewer.pal(9, "Reds")[9]),
       legend = c(expression(paste("FDR (", italic(H^2), " > 0)"))), bty = "n",
       x.intersp = 0.5, y.intersp = 0.6, cex = 7/par("ps")/par("cex"))
box(bty = "l", lwd = 0.5)
axis(side = 2, mgp = c(2, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
text(0.7, 1000, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 4, family = "Arial")
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# female old histogram
# ============================================================

female.hist.plot <- hist(female.H2.old, breaks = seq(0, 1, 0.05), col = "grey80", axes = FALSE, xlab = "", ylab = "", main = "", ylim = c(0, 2000))
hist(female.H2.old[female.H2.old.fdr <= 0.05], breaks = seq(0, 1, 0.05), col = brewer.pal(9, "Reds")[9], add = TRUE)
axis(side = 1, at = seq(0, 1, 0.2), mgp = c(0.8, -0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = expression(paste(italic(H^2), " in aged flies")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Number of genes", mgp = c(1.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
legend(-0.1, 2250, pch = 22, pt.bg = c(brewer.pal(9, "Reds")[9]),
       legend = c(expression(paste("FDR (", italic(H^2), " > 0)"))), bty = "n",
       x.intersp = 0.5, y.intersp = 0.6, cex = 7/par("ps")/par("cex"))
box(bty = "l", lwd = 0.5)
axis(side = 2, mgp = c(2, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
text(0.7, 1000, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 4, family = "Arial")
text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# scatter plot
# ============================================================

plot(c(0, 1), c(0, 1), axes = FALSE, type = "n", xlab = "", ylab = "")
abline(a = 0, b = 1, xpd = F)
points(female.H2.young, female.H2.old, col = brewer.pal(9, "Reds")[9], cex = 0.5)
title(xlab = expression(paste(italic(H^2), " in young flies")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(H^2), " in aged flies")), mgp = c(1.3, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
text(grconvertX(0.05 + file.width/25.4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)
text(0.7, 0.15, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 4, family = "Arial")

# calculate H2 and fdr in males
# ============================================================

male.H2.old <- male.var[, 7]/(male.var[, 7] + male.var[, 8])
male.H2.old.fdr <- p.adjust(male.var[, 9], method = "BH")

male.H2.young <- male.var[, 10]/(male.var[, 10] + male.var[, 11])
male.H2.young.fdr <- p.adjust(male.var[, 12], method = "BH")

# male young histogram
# ============================================================

male.hist.plot <- hist(male.H2.young, breaks = seq(0, 1, 0.05), col = "grey80", axes = FALSE, xlab = "", ylab = "", main = "", ylim = c(0, 2000))
hist(male.H2.young[male.H2.young.fdr <= 0.05], breaks = seq(0, 1, 0.05), col = brewer.pal(9, "Blues")[9], add = TRUE)
axis(side = 1, at = seq(0, 1, 0.2), mgp = c(0.8, -0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = expression(paste(italic(H^2), " in young flies")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Number of genes", mgp = c(1.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
legend(-0.1, 2250, pch = 22, pt.bg = c(brewer.pal(9, "Blues")[9]),
       legend = c(expression(paste("FDR (", italic(H^2), " > 0)"))), bty = "n",
       x.intersp = 0.5, y.intersp = 0.6, cex = 7/par("ps")/par("cex"))
box(bty = "l", lwd = 0.5)
axis(side = 2, mgp = c(2, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
text(0.7, 1000, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 4, family = "Arial")
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# male old histogram
# ============================================================

male.hist.plot <- hist(male.H2.old, breaks = seq(0, 1, 0.05), col = "grey80", axes = FALSE, xlab = "", ylab = "", main = "", ylim = c(0, 2000))
hist(male.H2.old[male.H2.old.fdr <= 0.05], breaks = seq(0, 1, 0.05), col = brewer.pal(9, "Blues")[9], add = TRUE)
axis(side = 1, at = seq(0, 1, 0.2), mgp = c(0.8, -0.1, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(xlab = expression(paste(italic(H^2), " in aged flies")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Number of genes", mgp = c(1.8, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
legend(-0.1, 2250, pch = 22, pt.bg = c(brewer.pal(9, "Blues")[9]),
       legend = c(expression(paste("FDR (", italic(H^2), " > 0)"))), bty = "n",
       x.intersp = 0.5, y.intersp = 0.6, cex = 7/par("ps")/par("cex"))
box(bty = "l", lwd = 0.5)
axis(side = 2, mgp = c(2, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
text(0.7, 1000, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 4, family = "Arial")
text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("e")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# scatter plot
# ============================================================

plot(c(0, 1), c(0, 1), axes = FALSE, type = "n", xlab = "", ylab = "")
abline(a = 0, b = 1, xpd = F)
points(male.H2.young, male.H2.old, col = brewer.pal(9, "Blues")[9], cex = 0.5)
title(xlab = expression(paste(italic(H^2), " in young flies")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(H^2), " in aged flies")), mgp = c(1.3, 0, 0), cex.lab = 7/par("ps")/par("cex"))
axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")
text(grconvertX(0.05 + file.width/25.4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("f")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)
text(0.7, 0.15, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 4, family = "Arial")

dev.off()


