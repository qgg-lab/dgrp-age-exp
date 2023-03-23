# =================================
# = PCA and batch correction, etc =
# =================================

args <- commandArgs(TRUE) # args <- c("female.pca.RData", "female.sva.adjust.RData", "male.pca.RData", "male.sva.adjust.RData", "figurePCA.pdf", "Myriad Pro")

library("RColorBrewer")


# prepare file to plot
# ============================================================

file.width = 135
cairo_pdf(file = args[5], width = file.width/25.4, height = file.width * 1.2/25.4, family = args[6])
layout(mat = matrix(c(1, 2, 3, 4, 4, 4, 5, 6, 7, 8, 8, 8), ncol = 3, nrow = 4, byrow = TRUE))
par(las = 1, tcl = -0.2, mar = c(2, 2.5, 1, 0.9), ps = 7, lwd = 0.5)
fc.col <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))

# 1. female PCA versus age
# ============================================================

load(args[1])

plot(c(-160, 160), c(-160, 160), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)

points(all.pca.top2, col = ifelse(age == "3wk", brewer.pal(9, "Blues")[9], brewer.pal(9, "Reds")[9]), pch = 1, cex = 0.5)

legend("topleft", bty = "n", pch = 1, col = brewer.pal(9, "Blues")[9], legend = expression(paste("3wk")), cex = 8/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5, text.col = brewer.pal(9, "Blues")[9])
legend("bottomright", bty = "n", pch = 1, col = brewer.pal(9, "Reds")[9], legend = expression(paste("baseline")), cex = 8/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5, text.col = brewer.pal(9, "Reds")[9])

axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5, at = seq(-150, 150, 50)[seq(1, 7, 2)])
axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5, at = seq(-150, 150, 50)[seq(2, 6, 2)])
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
title(ylab = paste("PC2 (", formatC((all.pca.eigen/sum(all.pca.eigen))[2]*100, digits = 2, format = "f"), "% variance)", sep = ""), cex.lab = 7/par("ps")/par("cex"), mgp = c(1.5, 0, 0))
title(xlab = paste("PC1 (", formatC((all.pca.eigen/sum(all.pca.eigen))[1]*100, digits = 2, format = "f"), "% variance)", sep = ""), cex.lab = 7/par("ps")/par("cex"), mgp = c(0.8, 0, 0))
title(main = "Females (both ages)", cex.main = 7/par("ps")/par("cex"))

box(bty = "l")

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# 2. 3wk PCA versus fc
# ============================================================

plot(c(-160, 160), c(-160, 160), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)

points(pca.3wk.top2, col = fc.col[as.numeric(factor(flow.cell[,2][age == "3wk"]))], pch = 1, cex = 0.5)

axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5, at = seq(-150, 150, 50)[seq(1, 7, 2)])
axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5, at = seq(-150, 150, 50)[seq(2, 6, 2)])
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
title(ylab = paste("PC2 (", formatC((pca.3wk.eigen/sum(pca.3wk.eigen))[2]*100, digits = 2, format = "f"), "% variance)", sep = ""), cex.lab = 7/par("ps")/par("cex"), mgp = c(1.5, 0, 0))
title(xlab = paste("PC1 (", formatC((pca.3wk.eigen/sum(pca.3wk.eigen))[1]*100, digits = 2, format = "f"), "% variance)", sep = ""), cex.lab = 7/par("ps")/par("cex"), mgp = c(0.8, 0, 0))
title(main = expression(bold(paste("Females (3wk)"))), cex.main = 7/par("ps")/par("cex"))
col.levels <- sort(as.numeric(unique(factor(flow.cell[,2][age == "3wk"]))))
points(seq(-40, by = 13, length = length(col.levels)), rep(150, length(col.levels)), col = fc.col[col.levels], pch = 1, cex = 0.6, xpd = TRUE)
text(-180, 150, "Flow cell", pos = 4, cex = 7/par("ps")/par("cex"))

box(bty = "l")

text(grconvertX(0.05 + file.width/25.4/3, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# 3. baseline PCA versus fc
# ============================================================

plot(c(-160, 160), c(-160, 160), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)

points(pca.baseline.top2, col = fc.col[as.numeric(factor(flow.cell[,2][age == "baseline"]))], pch = 1, cex = 0.5)

axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5, at = seq(-150, 150, 50)[seq(1, 7, 2)])
axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5, at = seq(-150, 150, 50)[seq(2, 6, 2)])
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
title(ylab = paste("PC2 (", formatC((pca.baseline.eigen/sum(pca.baseline.eigen))[2]*100, digits = 2, format = "f"), "% variance)", sep = ""), cex.lab = 7/par("ps")/par("cex"), mgp = c(1.5, 0, 0))
title(xlab = paste("PC1 (", formatC((pca.baseline.eigen/sum(pca.baseline.eigen))[1]*100, digits = 2, format = "f"), "% variance)", sep = ""), cex.lab = 7/par("ps")/par("cex"), mgp = c(0.8, 0, 0))
title(main = expression(bold(paste("Females (baseline)"))), cex.main = 7/par("ps")/par("cex"))
col.levels <- sort(as.numeric(unique(factor(flow.cell[,2][age == "baseline"]))))
points(seq(-40, by = 13, length = length(col.levels)), rep(150, length(col.levels)), col = fc.col[col.levels], pch = 1, cex = 0.6, xpd = TRUE)
text(-180, 150, "Flow cell", pos = 4, cex = 7/par("ps")/par("cex"))

box(bty = "l")

text(grconvertX(0.05 + file.width/25.4/3*2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# 4. boxplot for sva versus fc
# ============================================================
load(args[2])

plot(c(0, 19), c(-0.15, 0.15), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)

bp.3wk <- boxplot(sva.adj$sv[age == "3wk", 1] ~ flow.cell[,2][age == "3wk"], plot = FALSE, range = 0)
col.levels <- sort(as.numeric(unique(factor(flow.cell[,2][age == "3wk"]))))

bxp(bp.3wk, width = rep(0.5, length(col.levels)), add = TRUE, axes = FALSE, at = 1:length(bp.3wk$names), boxfill = fc.col[col.levels], medlwd = 0.5, boxcol = fc.col[col.levels], whiskcol = fc.col[col.levels], staplecol = fc.col[col.levels], whisklty = 1)

bp.baseline <- boxplot(sva.adj$sv[age == "baseline", 1] ~ flow.cell[,2][age == "baseline"], plot = FALSE, range = 0)
col.levels <- sort(as.numeric(unique(factor(flow.cell[,2][age == "baseline"]))))

bxp(bp.baseline, width = rep(0.5, length(col.levels)), add = TRUE, axes = FALSE, at = (length(bp.3wk$names)+1):(length(bp.3wk$names)+length(bp.baseline$names)), boxfill = fc.col[col.levels], medlwd = 0.5, boxcol = fc.col[col.levels], whiskcol = fc.col[col.levels], staplecol = fc.col[col.levels], whisklty = 1)


axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, -0.3))
title(ylab = "1st surrogate variable", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.5, 0, 0))

segments(1, -0.13, length(bp.3wk$names), -0.13, lwd = 1, col = brewer.pal(9, "Blues")[9])
text((1+length(bp.3wk$names))/2, -0.15, expression(paste("3wk flow cell")), xpd = TRUE, cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])


segments(length(bp.3wk$names) + 1, -0.13, length(bp.3wk$names) + length(bp.baseline$names), -0.13, lwd = 1, col = brewer.pal(9, "Reds")[9])
text(length(bp.baseline$names)/2+length(bp.3wk$names) + 1, -0.15, expression(paste("baseline flow cell")), xpd = TRUE, cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])


text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)
title(main = "Females", cex.main = 7/par("ps")/par("cex"), line = -1)


# 5. female PCA versus age
# ============================================================

load(args[3])

plot(c(-160, 160), c(-160, 160), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)

points(all.pca.top2, col = ifelse(age == "3wk", brewer.pal(9, "Blues")[9], brewer.pal(9, "Reds")[9]), pch = 1, cex = 0.5)

legend("topleft", bty = "n", pch = 1, col = brewer.pal(9, "Blues")[9], legend = expression(paste("3wk")), cex = 8/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5, text.col = brewer.pal(9, "Blues")[9])
legend("bottomright", bty = "n", pch = 1, col = brewer.pal(9, "Reds")[9], legend = expression(paste("baseline")), cex = 8/par("ps")/par("cex"), x.intersp = 0.5, y.intersp = 0.5, text.col = brewer.pal(9, "Reds")[9])

axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5, at = seq(-150, 150, 50)[seq(1, 7, 2)])
axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5, at = seq(-150, 150, 50)[seq(2, 6, 2)])
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
title(ylab = paste("PC2 (", formatC((all.pca.eigen/sum(all.pca.eigen))[2]*100, digits = 2, format = "f"), "% variance)", sep = ""), cex.lab = 7/par("ps")/par("cex"), mgp = c(1.5, 0, 0))
title(xlab = paste("PC1 (", formatC((all.pca.eigen/sum(all.pca.eigen))[1]*100, digits = 2, format = "f"), "% variance)", sep = ""), cex.lab = 7/par("ps")/par("cex"), mgp = c(0.8, 0, 0))
title(main = "Males (both ages)", cex.main = 7/par("ps")/par("cex"))

box(bty = "l")

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("e")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# 6. 3wk PCA versus fc
# ============================================================

plot(c(-160, 160), c(-160, 160), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)

points(pca.3wk.top2, col = fc.col[as.numeric(factor(flow.cell[,2][age == "3wk"]))], pch = 1, cex = 0.5)

axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5, at = seq(-150, 150, 50)[seq(1, 7, 2)])
axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5, at = seq(-150, 150, 50)[seq(2, 6, 2)])
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
title(ylab = paste("PC2 (", formatC((pca.3wk.eigen/sum(pca.3wk.eigen))[2]*100, digits = 2, format = "f"), "% variance)", sep = ""), cex.lab = 7/par("ps")/par("cex"), mgp = c(1.5, 0, 0))
title(xlab = paste("PC1 (", formatC((pca.3wk.eigen/sum(pca.3wk.eigen))[1]*100, digits = 2, format = "f"), "% variance)", sep = ""), cex.lab = 7/par("ps")/par("cex"), mgp = c(0.8, 0, 0))
title(main = expression(bold(paste("Males (3wk)"))), cex.main = 7/par("ps")/par("cex"))
col.levels <- sort(as.numeric(unique(factor(flow.cell[,2][age == "3wk"]))))
points(seq(-40, by = 13, length = length(col.levels)), rep(150, length(col.levels)), col = fc.col[col.levels], pch = 1, cex = 0.6, xpd = TRUE)
text(-180, 150, "Flow cell", pos = 4, cex = 7/par("ps")/par("cex"))

box(bty = "l")

text(grconvertX(0.05 + file.width/25.4/3, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("f")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# 7. baseline PCA versus fc
# ============================================================

plot(c(-160, 160), c(-160, 160), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)

points(pca.baseline.top2, col = fc.col[as.numeric(factor(flow.cell[,2][age == "baseline"]))], pch = 1, cex = 0.5)

axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5, at = seq(-150, 150, 50)[seq(1, 7, 2)])
axis(side = 1, mgp = c(1, 0, 0), cex.axis = 7/par("ps")/par("cex"), lwd = 0.5, at = seq(-150, 150, 50)[seq(2, 6, 2)])
axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, 0))
title(ylab = paste("PC2 (", formatC((pca.baseline.eigen/sum(pca.baseline.eigen))[2]*100, digits = 2, format = "f"), "% variance)", sep = ""), cex.lab = 7/par("ps")/par("cex"), mgp = c(1.5, 0, 0))
title(xlab = paste("PC1 (", formatC((pca.baseline.eigen/sum(pca.baseline.eigen))[1]*100, digits = 2, format = "f"), "% variance)", sep = ""), cex.lab = 7/par("ps")/par("cex"), mgp = c(0.8, 0, 0))
title(main = expression(bold(paste("Males (baseline)"))), cex.main = 7/par("ps")/par("cex"))
col.levels <- sort(as.numeric(unique(factor(flow.cell[,2][age == "baseline"]))))
points(seq(-40, by = 13, length = length(col.levels)), rep(150, length(col.levels)), col = fc.col[col.levels], pch = 1, cex = 0.6, xpd = TRUE)
text(-180, 150, "Flow cell", pos = 4, cex = 7/par("ps")/par("cex"))

box(bty = "l")

text(grconvertX(0.05 + file.width/25.4/3*2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("g")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# 8. boxplot for sva versus fc
# ============================================================
load(args[4])

plot(c(0, 19), c(-0.15, 0.15), type = "n", xlab = "", ylab = "", main = "", axes = FALSE)

bp.3wk <- boxplot(sva.adj$sv[age == "3wk", 1] ~ flow.cell[,2][age == "3wk"], plot = FALSE, range = 0)
col.levels <- sort(as.numeric(unique(factor(flow.cell[,2][age == "3wk"]))))

bxp(bp.3wk, width = rep(0.5, length(col.levels)), add = TRUE, axes = FALSE, at = 1:length(bp.3wk$names), boxfill = fc.col[col.levels], medlwd = 0.5, boxcol = fc.col[col.levels], whiskcol = fc.col[col.levels], staplecol = fc.col[col.levels], whisklty = 1)

bp.baseline <- boxplot(sva.adj$sv[age == "baseline", 1] ~ flow.cell[,2][age == "baseline"], plot = FALSE, range = 0)
col.levels <- sort(as.numeric(unique(factor(flow.cell[,2][age == "baseline"]))))

bxp(bp.baseline, width = rep(0.5, length(col.levels)), add = TRUE, axes = FALSE, at = (length(bp.3wk$names)+1):(length(bp.3wk$names)+length(bp.baseline$names)), boxfill = fc.col[col.levels], medlwd = 0.5, boxcol = fc.col[col.levels], whiskcol = fc.col[col.levels], staplecol = fc.col[col.levels], whisklty = 1)


axis(side = 2, lwd = 0.5, cex.axis = 7/par("ps")/par("cex"), mgp = c(2.5, 0.3, -0.3))
title(ylab = "1st surrogate variable", cex.lab = 7/par("ps")/par("cex"), mgp = c(1.5, 0, 0))

segments(1, -0.13, length(bp.3wk$names), -0.13, lwd = 1, col = brewer.pal(9, "Blues")[9])
text((1+length(bp.3wk$names))/2, -0.15, expression(paste("3wk flow cell")), xpd = TRUE, cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])


segments(length(bp.3wk$names) + 1, -0.13, length(bp.3wk$names) + length(bp.baseline$names), -0.13, lwd = 1, col = brewer.pal(9, "Reds")[9])
text(length(bp.baseline$names)/2+length(bp.3wk$names) + 1, -0.15, expression(paste("baseline flow cell")), xpd = TRUE, cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])


text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("h")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)
title(main = "Males", cex.main = 7/par("ps")/par("cex"), line = -1)

dev.off()
