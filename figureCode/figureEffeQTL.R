# ============================
# = figure for response eQTL =
# ============================

args <- commandArgs(TRUE) # args <- c("../figureData/female.eqtl.eff.RData", "../figureData/male.eqtl.eff.RData")
library("RColorBrewer")
library("plotrix")

# read data
# ============================================================

load(args[1])
female.eqtl.pair <- eqtl.pair
load(args[2])
male.eqtl.pair <- eqtl.pair

# set up file
# ============================================================

file.width = 100
cairo_pdf(file = args[3], width = file.width/25.4, height = file.width/25.4*0.5, family = args[4])
par(las = 1, tcl = -0.2, mfrow = c(1, 2), mar = c(2, 2, 0.5, 0.3), ps = 7, lwd = 0.5)


# ============================================================
female.col <- rep(brewer.pal(9, "Purples")[9], nrow(female.eqtl.pair))
female.col[female.eqtl.pair[, 3] == "18c"] <- brewer.pal(9, "Blues")[9]
female.col[female.eqtl.pair[, 3] == "25c"] <- brewer.pal(9, "Reds")[9]

cat("range of effects (standardized):", range(as.numeric(female.eqtl.pair[, 4])/sqrt(as.numeric(female.eqtl.pair[, 6]))), range(as.numeric(female.eqtl.pair[, 7])/sqrt(as.numeric(female.eqtl.pair[, 9]))), "\n")

cat("there are a total of", nrow(female.eqtl.pair), "eQTL-gene pairs.\n")
print(table(female.eqtl.pair[,3]))

plot(as.numeric(female.eqtl.pair[, 7])/sqrt(as.numeric(female.eqtl.pair[, 9])), as.numeric(female.eqtl.pair[, 4])/sqrt(as.numeric(female.eqtl.pair[, 6])), col = female.col, cex = 0.2, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), xlab = "", ylab = "", axes = FALSE, lwd = 0.2)

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(-1, 1, 0.5), cex.axis = 7/par("ps")/par("cex"))
axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = -0.5, cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, at = seq(-1, 1, 0.5), cex.axis = 7/par("ps")/par("cex"))

title(xlab = "Standardized effect in young flies", mgp = c(0.7, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Standardized effect in aged flies", mgp = c(1.15, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")

legend("topleft", col = c(brewer.pal(9, "Reds")[9], brewer.pal(9, "Blues")[9], brewer.pal(9, "Purples")[9]), pch = 1, legend = c("Young", "Aged", "Both"), y.intersp = 0.6, x.intersp = 0.5, pt.lwd = 1, cex = 6/par("ps")/par("cex"), bty = "n")

text(1, -0.5, "\u2640", cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 1, family = "Arial")
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# male
# ============================================================

male.col <- rep(brewer.pal(9, "Purples")[9], nrow(male.eqtl.pair))
male.col[male.eqtl.pair[, 3] == "18c"] <- brewer.pal(9, "Blues")[9]
male.col[male.eqtl.pair[, 3] == "25c"] <- brewer.pal(9, "Reds")[9]

cat("range of effects (standardized):", range(as.numeric(male.eqtl.pair[, 4])/sqrt(as.numeric(male.eqtl.pair[, 6]))), range(as.numeric(male.eqtl.pair[, 7])/sqrt(as.numeric(male.eqtl.pair[, 9]))), "\n")

cat("there are a total of", nrow(male.eqtl.pair), "eQTL-gene pairs.\n")
print(table(male.eqtl.pair[,3]))

plot(as.numeric(male.eqtl.pair[, 7])/sqrt(as.numeric(male.eqtl.pair[, 9])), as.numeric(male.eqtl.pair[, 4])/sqrt(as.numeric(male.eqtl.pair[, 6])), col = male.col, cex = 0.2, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), xlab = "", ylab = "", axes = FALSE, lwd = 0.2)

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(-1, 1, 0.5), cex.axis = 7/par("ps")/par("cex"))
axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = -0.5, cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, at = seq(-1, 1, 0.5), cex.axis = 7/par("ps")/par("cex"))

title(xlab = "Standardized effect in young flies", mgp = c(0.7, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Standardized effect in aged flies", mgp = c(1.15, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")

legend("topleft", col = c(brewer.pal(9, "Reds")[9], brewer.pal(9, "Blues")[9], brewer.pal(9, "Purples")[9]), pch = 1, legend = c("Young", "Aged", "Both"), y.intersp = 0.6, x.intersp = 0.5, pt.lwd = 1, cex = 6/par("ps")/par("cex"), bty = "n")

text(1, -0.5, "\u2642", cex = 16/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 1, family = "Arial")
text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


dev.off()
