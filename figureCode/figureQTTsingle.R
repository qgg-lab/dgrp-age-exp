# ================
# = qtt analysis =
# ================

args <- commandArgs(TRUE) # args <- c("../figureData/female.phototaxis.decline.qtt.RData", "../figureData/female.fecundity.decline.qtt.RData", "../figureData/female.lifespan.qtt.RData", "../figureData/male.phototaxis.decline.qtt.RData", "../figureData/male.speed.decline.qtt.RData", "../figureData/male.endurance.decline.qtt.RData", "../figureData/male.lifespan.qtt.RData")

library("RColorBrewer")


# prepare file
# ============================================================

file.width = 89
cairo_pdf(file = args[8], width = file.width/25.4, height = file.width/25.4, family = args[9])
par(las = 1, tcl = -0.2, mai = c(0.15, 0.1, 0.04, 0.01)*file.width/25.4, ps = 7, lwd = 0.5)
layout(mat = matrix(c(1, 2), ncol = 1, byrow = TRUE))

# female traits
# ============================================================

load(args[1])
phototaxis.box <- boxplot(cbind(qtt.cor.young, qtt.cor.old, qtt.cor.diff), plot = F, range = 0)
load(args[2])
fecundity.box <- boxplot(cbind(qtt.cor.young, qtt.cor.old, qtt.cor.diff), plot = F, range = 0)
load(args[3])
lifespan.box <- boxplot(cbind(qtt.cor.young, qtt.cor.old, qtt.cor.diff), plot = F, range = 0)

box.col <- c(brewer.pal(9, "Blues")[9], brewer.pal(9, "Reds")[9], brewer.pal(9, "Purples")[9])

plot(c(0, 16), c(-0.5, 0.5), type = "n", axes = FALSE, xlab = "", ylab = "")
bxp(phototaxis.box, show.names = FALSE, at = 1:3, add = TRUE, axes = FALSE, boxfill = box.col, medlwd = 0.5, boxcol = box.col, whiskcol = box.col, staplecol = box.col, medcol = "white")

bxp(fecundity.box, show.names = FALSE, at = 5:7, add = TRUE, axes = FALSE, boxfill = box.col, medlwd = 0.5, boxcol = box.col, whiskcol = box.col, staplecol = box.col, medcol = "white")

bxp(lifespan.box, show.names = FALSE, at = 9:11, add = TRUE, axes = FALSE, boxfill = box.col, medlwd = 0.5, boxcol = box.col, whiskcol = box.col, staplecol = box.col, medcol = "white")

axis(side = 1, at = c(1:3, 5:7, 9:11), labels = FALSE, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, -0.1, 0), lwd = 0.5)
axis(side = 2, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0.3, 0), lwd = 0.5)
box(bty = "l")

text(1:3 + 0.52, -0.61, c("Young", "Aged", "Aged - Young"), srt = 45, xpd = TRUE, pos = 2, col = box.col)
text(5:7 + 0.52, -0.61, c("Young", "Aged", "Aged - Young"), srt = 45, xpd = TRUE, pos = 2, col = box.col)
text(9:11 + 0.52, -0.61, c("Young", "Aged", "Aged - Young"), srt = 45, xpd = TRUE, pos = 2, col = box.col)

text(2, 0.5, "Phototaxis decline", cex = 6/par("ps")/par("cex"))
text(6, 0.5, "Fecundity decline", cex = 6/par("ps")/par("cex"))
text(10, 0.5, "Lifespan", cex = 6/par("ps")/par("cex"))

title(ylab = "Correlation", cex.lab = 7/par("ps")/par("cex"), mgp = c(1, 0, 0))
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)
# males
# ============================================================

load(args[4])
phototaxis.box <- boxplot(cbind(qtt.cor.young, qtt.cor.old, qtt.cor.diff), plot = F, range = 0)
load(args[5])
speed.box <- boxplot(cbind(qtt.cor.young, qtt.cor.old, qtt.cor.diff), plot = F, range = 0)
load(args[6])
endurance.box <- boxplot(cbind(qtt.cor.young, qtt.cor.old, qtt.cor.diff), plot = F, range = 0)
load(args[7])
lifespan.box <- boxplot(cbind(qtt.cor.young, qtt.cor.old, qtt.cor.diff), plot = F, range = 0)

box.col <- c(brewer.pal(9, "Blues")[9], brewer.pal(9, "Reds")[9], brewer.pal(9, "Purples")[9])

plot(c(0, 16), c(-0.5, 0.5), type = "n", axes = FALSE, xlab = "", ylab = "")
bxp(phototaxis.box, show.names = FALSE, at = 1:3, add = TRUE, axes = FALSE, boxfill = box.col, medlwd = 0.5, boxcol = box.col, whiskcol = box.col, staplecol = box.col, medcol = "white")

bxp(speed.box, show.names = FALSE, at = 5:7, add = TRUE, axes = FALSE, boxfill = box.col, medlwd = 0.5, boxcol = box.col, whiskcol = box.col, staplecol = box.col, medcol = "white")

bxp(endurance.box, show.names = FALSE, at = 9:11, add = TRUE, axes = FALSE, boxfill = box.col, medlwd = 0.5, boxcol = box.col, whiskcol = box.col, staplecol = box.col, medcol = "white")

bxp(lifespan.box, show.names = FALSE, at = 13:15, add = TRUE, axes = FALSE, boxfill = box.col, medlwd = 0.5, boxcol = box.col, whiskcol = box.col, staplecol = box.col, medcol = "white")

axis(side = 1, at = c(1:3, 5:7, 9:11, 13:15), labels = FALSE, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, -0.1, 0), lwd = 0.5)
axis(side = 2, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0.3, 0), lwd = 0.5)
box(bty = "l")

text(1:3 + 0.52, -0.61, c("Young", "Aged", "Aged - Young"), srt = 45, xpd = TRUE, pos = 2, col = box.col)
text(5:7 + 0.52, -0.61, c("Young", "Aged", "Aged - Young"), srt = 45, xpd = TRUE, pos = 2, col = box.col)
text(9:11 + 0.52, -0.61, c("Young", "Aged", "Aged - Young"), srt = 45, xpd = TRUE, pos = 2, col = box.col)
text(13:15 + 0.52, -0.61, c("Young", "Aged", "Aged - Young"), srt = 45, xpd = TRUE, pos = 2, col = box.col)

text(2, 0.5, "Phototaxis decline", cex = 6/par("ps")/par("cex"))
text(6, 0.5, "Climbing speed\ndecline", cex = 6/par("ps")/par("cex"), xpd = T)
text(10, 0.5, "Climbing endurance\ndecline", cex = 6/par("ps")/par("cex"), xpd = T)
text(14, 0.5, "Lifespan", cex = 6/par("ps")/par("cex"))

title(ylab = "Correlation", cex.lab = 7/par("ps")/par("cex"), mgp = c(1, 0, 0))
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# # function to calculate cumulative fraction
# # ============================================================
#
# cumfrac <- function(x, from, to, int) {
#
# 	x.cut <- cut(x, breaks = seq(from, to, int))
# 	frac <- c(0, cumsum(table(x.cut))/length(x))
# 	x.pos <- seq(from, to, int)
# 	return(list(x = x.pos, frac = frac))
#
# }
#
# # female
# # ============================================================
#
# load(args[1])
#
#
# #
# #
# # qtt.young.hist <- hist(qtt.cor.young, breaks = seq(-0.36, 0.36, 0.02), plot = FALSE)$counts
# # qtt.old.hist <- hist(qtt.cor.old, breaks = seq(-0.36, 0.36, 0.02), plot = FALSE)$counts
# # qtt.diff.hist <- hist(qtt.cor.diff, breaks = seq(-0.36, 0.36, 0.02), plot = FALSE)$counts
# #
# # barplot(rbind(qtt.young.hist, qtt.old.hist, qtt.diff.hist), xlim = c(0, 144), ylim = c(-100, 1300), beside = TRUE, col = c(brewer.pal(9, "Set1")[c(2, 1, 9)]), border = NA, axes = FALSE, xlab = "", ylab = "")
# # axis(side = 1, at = c(12, 32, 52, 72, 92, 112, 132) + 0.5, labels = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3), cex.axis = 6/par("ps")/par("cex"), mgp = c(2, -0.1, 0), lwd = 0.5)
# # axis(side = 2, cex.axis = 6/par("ps")/par("cex"), mgp = c(2, 0.3, 0), lwd = 0.5)
# # box(bty = "l")
# #
# # title(xlab = "Life span (days)", cex.lab = 7/par("ps")/par("cex"), mgp = c(0.5, 0, 0))
# # title(ylab = "Number of flies", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.2, 0, 0))
#
# # cumulative fraction
# qtt.young.frac <- cumfrac(qtt.cor.young, -0.36, 0.36, 0.01)
# qtt.old.frac <- cumfrac(qtt.cor.old, -0.36, 0.36, 0.01)
# qtt.diff.frac <- cumfrac(qtt.cor.diff, -0.36, 0.36, 0.01)
#
# range(qtt.cor.young); range(qtt.cor.old); range(qtt.cor.diff)
#
# plot(c(-0.36, 0.36), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
# lines(qtt.young.frac$x, qtt.young.frac$frac, lwd = 1, col = brewer.pal(9, "Set1")[2])
# lines(qtt.old.frac$x, qtt.old.frac$frac, lwd = 1, col = brewer.pal(9, "Set1")[1])
# lines(qtt.diff.frac$x, qtt.diff.frac$frac, lwd = 1, col = brewer.pal(9, "Set1")[9])
#
# axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), seq(-0.3, 0.3, 0.2), cex.axis = 7/par("ps")/par("cex"))
# axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
# box(bty = "l")
# title(xlab = "Correlation with lifespan", mgp = c(0.5, 0, 0), cex.lab = 7/par("ps")/par("cex"))
# title(ylab = "Cumulative fraction", mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
# text(0.2, 0.3, "\u2640", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], family = "Arial")
# legend("topleft", lwd = 1, col = brewer.pal(9, "Set1")[c(2, 1, 9)],
#        legend = c("Young", "Old", "Old - Young"), bty = "n", seg.len = 0.8,
#        x.intersp = 0.5, y.intersp = 0.5, cex = 6/par("ps")/par("cex"))
#
#
# text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)
#
#
# # males
# load(args[2])
# qtt.young.frac <- cumfrac(qtt.cor.young, -0.36, 0.36, 0.01)
# qtt.old.frac <- cumfrac(qtt.cor.old, -0.36, 0.36, 0.01)
# qtt.diff.frac <- cumfrac(qtt.cor.diff, -0.36, 0.36, 0.01)
# range(qtt.cor.young); range(qtt.cor.old); range(qtt.cor.diff)
#
# plot(c(-0.36, 0.36), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
# lines(qtt.young.frac$x, qtt.young.frac$frac, lwd = 1, col = brewer.pal(9, "Set1")[2])
# lines(qtt.old.frac$x, qtt.old.frac$frac, lwd = 1, col = brewer.pal(9, "Set1")[1])
# lines(qtt.diff.frac$x, qtt.diff.frac$frac, lwd = 1, col = brewer.pal(9, "Set1")[9])
#
# axis(side = 1, lwd = 0.5, mgp = c(0.8, -0.1, 0), seq(-0.3, 0.3, 0.2), cex.axis = 7/par("ps")/par("cex"))
# axis(side = 2, lwd = 0.5, mgp = c(0.8, 0.3, 0), cex.axis = 7/par("ps")/par("cex"))
# box(bty = "l")
# title(xlab = "Correlation with lifespan", mgp = c(0.5, 0, 0), cex.lab = 7/par("ps")/par("cex"))
# title(ylab = "Cumulative fraction", mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
# text(0.2, 0.3, "\u2642", cex = 20/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], family = "Arial")
# legend("topleft", lwd = 1, col = brewer.pal(9, "Set1")[c(2, 1, 9)],
#        legend = c("Young", "Old", "Old - Young"), bty = "n", seg.len = 0.8,
#        x.intersp = 0.5, y.intersp = 0.5, cex = 6/par("ps")/par("cex"))
#
#
# text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

dev.off()
