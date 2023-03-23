# =======================================
# = make figure to summarize QG results =
# =======================================

args <- commandArgs(TRUE) # args = c("../figureData/female.adj.qgVarHet.RData", "../figureData/male.adj.qgVarHet.RData", "Myriad Pro")
library("RColorBrewer")

# prepare file to plot
# ============================================================

file.width = 89
cairo_pdf(file = args[3], width = file.width/25.4, height = file.width/25.4 * 0.45, family = args[4])
par(las = 1, tcl = -0.2, mai = c(0.16, 0.2, 0.1, 0.02)*file.width/25.4/2, ps = 7, lwd = 0.5, xpd = F, mfrow = c(1, 2))

# female variance het
# ============================================================

load(args[1])
rownames(var.het) <- gene.name; female.var.het <- var.het

female.fdr.single.e <- p.adjust(as.numeric(female.var.het[, 11]), method = "BH")

sig.female.e <- female.fdr.single.e < 0.05 & ( (as.numeric(female.var.het[, 9])/as.numeric(female.var.het[, 10]) > 2)  | (as.numeric(female.var.het[, 10])/as.numeric(female.var.het[, 9]) > 2) )

plot(c(0, 6), c(0, 6), axes = FALSE, type = "n", xlab = "", ylab = "")

points(sqrt(as.numeric(female.var.het[!sig.female.e, 10])), sqrt(as.numeric(female.var.het[!sig.female.e, 9])), col = rgb(col2rgb("grey50")[1, 1], col2rgb("grey50")[2, 1], col2rgb("grey50")[3, 1], 40, max = 255), pch = 16, cex = 0.5)

segments(-0.5, -0.5, 6, 6, lty = 2, col = "grey20")

points(sqrt(as.numeric(female.var.het[sig.female.e, 10])), sqrt(as.numeric(female.var.het[sig.female.e, 9])), col = rgb(col2rgb(brewer.pal(9, "Reds")[9])[1, 1], col2rgb(brewer.pal(9, "Reds")[9])[2, 1], col2rgb(brewer.pal(9, "Reds")[9])[3, 1], 120, max = 255), pch = 16, cex = 0.5)

lines(c(0, 6), c(0, 6)*sqrt(2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Reds")[9])
lines(c(0, 6), c(0, 6)*sqrt(1/2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Reds")[9])

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = seq(0, 6, 2), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 6, 2), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")

title(xlab = expression(paste(italic(sigma[e]), " in young flies")), mgp = c(0.5, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(sigma[e]), " in aged flies")), mgp = c(0.6, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(5, 4, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 3, family = "Arial")
text(4, 0.2, expression(paste(italic(sigma[e]^2), " (young/aged)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])
text(1.8, 5.5, expression(paste(italic(sigma[e]^2), " (aged/young)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

cat("In females, there are", sum(sig.female.e, na.rm = T), "genes with variance het.\n")
cat("among which,", sum(sig.female.e & (as.numeric(female.var.het[, 10])/as.numeric(female.var.het[, 9]) > 2), na.rm = T), "higher var in young.\n")
cat("among which,", sum(sig.female.e & (as.numeric(female.var.het[, 10])/as.numeric(female.var.het[, 9]) < 1/2), na.rm = T), "higher var in old.\n")

# male
load(args[2])
rownames(var.het) <- gene.name; male.var.het <- var.het

male.fdr.single.e <- p.adjust(as.numeric(male.var.het[, 11]), method = "BH")

sig.male.e <- male.fdr.single.e < 0.05 & ( (as.numeric(male.var.het[, 9])/as.numeric(male.var.het[, 10]) > 2)  | (as.numeric(male.var.het[, 10])/as.numeric(male.var.het[, 9]) > 2) )

plot(c(0, 6), c(0, 6), axes = FALSE, type = "n", xlab = "", ylab = "")

points(sqrt(as.numeric(male.var.het[!sig.male.e, 10])), sqrt(as.numeric(male.var.het[!sig.male.e, 9])), col = rgb(col2rgb("grey50")[1, 1], col2rgb("grey50")[2, 1], col2rgb("grey50")[3, 1], 40, max = 255), pch = 16, cex = 0.5)

segments(-0.5, -0.5, 6, 6, lty = 2, col = "grey20")

points(sqrt(as.numeric(male.var.het[sig.male.e, 10])), sqrt(as.numeric(male.var.het[sig.male.e, 9])), col = rgb(col2rgb(brewer.pal(9, "Blues")[9])[1, 1], col2rgb(brewer.pal(9, "Blues")[9])[2, 1], col2rgb(brewer.pal(9, "Blues")[9])[3, 1], 120, max = 255), pch = 16, cex = 0.5)

lines(c(0, 6), c(0, 6)*sqrt(2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Blues")[9])
lines(c(0, 6), c(0, 6)*sqrt(1/2), lwd = 0.5, lty = 2, col = brewer.pal(9, "Blues")[9])

axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), at = seq(0, 6, 2), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), at = seq(0, 6, 2), cex.axis = 7/par("ps")/par("cex"))
box(bty = "l")

title(xlab = expression(paste(italic(sigma[e]), " in young flies")), mgp = c(0.5, 0, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic(sigma[e]), " in aged flies")), mgp = c(0.6, 0, 0), cex.lab = 7/par("ps")/par("cex"))
text(5, 4, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 3, family = "Arial")
text(4, 0.2, expression(paste(italic(sigma[e]^2), " (young/aged)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])
text(1.8, 5.5, expression(paste(italic(sigma[e]^2), " (aged/young)", " > 2")), cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

text(grconvertX(0.05 + file.width/2/25.4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

cat("In males, there are", sum(sig.male.e, na.rm = T), "genes with variance het.\n")
cat("among which,", sum(sig.male.e & (as.numeric(male.var.het[, 10])/as.numeric(male.var.het[, 9]) > 2), na.rm = T), "higher var in young.\n")
cat("among which,", sum(sig.male.e & (as.numeric(male.var.het[, 10])/as.numeric(male.var.het[, 9]) < 1/2), na.rm = T), "higher var in old.\n")


dev.off()


