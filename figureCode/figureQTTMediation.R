# =================================
# = figure for mediation analysis =
# =================================

args = commandArgs(TRUE) # args <- c("../figureData/male.eqtl.gwas.trio.mediation.example.RData")

# 3R_12642719, 0 = TCTCCACCAAGC, 2 = T
# 3R_12642725, 0 = C, 2 = * (r2 = 1 with 3R_12642719)
# 2R_4440030, 0 = C, 2 = T
# 2R_22212424, 0 = G, 2 = A
# FBgn0260798: Gprk1 2R:-:4322295-4471060
# FBgn0025186: ari-2 2R:+:22212835-22217264
# FBgn0026207: mbo 3R:-:12676487-12679052
# load data

load(args[1])
library("RColorBrewer")
library("grid")
library("diagram")

# prepare file
# ============================================================

file.width = 89
cairo_pdf(file = args[2], width = file.width/25.4, height = file.width*1.25/25.4, family = args[3])
par(las = 1, tcl = -0.2, mai = c(0.15, 0.2, 0.1, 0.02)*file.width/25.4/2, ps = 7, lwd = 0.5)
layout(mat = matrix(c(1, 1, 2, 3, 4, 5), ncol = 2, byrow = TRUE), height = c(0.5, 1, 1))
set.seed(1)

# plot gene models
# ============================================================

par(mai = c(0.1, 0.1, 0.01, 0.01)*file.width/25.4/2)

plot(c(0, 9.2), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
points(1, 0.5, pch = 21, bg = brewer.pal(9, "Set1")[1], cex = 2)
rect(4, 0.35, 5, 0.65, col = brewer.pal(9, "Set1")[3])
text(4.5, 0.5, expression(italic("mbo")), cex = 7/par("ps")/par("cex"), col = "white")
rect(7.5, 0.35, 9.5, 0.65, col = brewer.pal(9, "Set1")[8])
par(lheight = 0.7)
text(8.5, 0.5, "Climbing speed\ndecline", cex = 7/par("ps")/par("cex"), col = "white")
#arrows(1.12, 0.55, 4, 0.55, length = 0.04)
# straightarrow(c(1.12, 0.55), c(3.85, 0.55), arr.pos = 1, lwd = 0.5, arr.length = 0.15, arr.width = 0.1)
segments(1.12, 0.55, 4, 0.55)


text(0.1, 0.5, expression(paste(italic("3R"), ":12642725")), cex = 6/par("ps")/par("cex"), xpd = TRUE)

text(2.5, 0.7, "eQTL-gene association (b)", cex = 6/par("ps")/par("cex"))
# segments(1.12, 0.45, 4, 0.45, lty = 2)
# arrows(3.99, 0.45, 4, 0.45, length = 0.04)
straightarrow(c(1.12, 0.45), c(3.85, 0.45), arr.pos = 1, lty = 1, lwd = 0.5, arr.length = 0.15, arr.width = 0.1)

straightarrow(c(5, 0.45), c(7.35, 0.45), arr.pos = 1, lty = 1, lwd = 0.5, arr.length = 0.15, arr.width = 0.1)

#straightarrow(c(5, 0.55), c(7.35, 0.55), arr.pos = 1, lty = 1, lwd = 0.5, arr.length = 0.15, arr.width = 0.1)
segments(5, 0.55, 7.5, 0.55)

curvedarrow(from = c(1, 0.3), to = c(8.5, 0.3), curve = 0.04, arr.pos = 1, lwd = 0.5, arr.length = 0, arr.width = 0)
text(4.5, -0.1, "QTL-trait association (c)", cex = 6/par("ps")/par("cex"), xpd = TRUE)

text(6.2, 0.7, "Gene-trait association (d)", cex = 6/par("ps")/par("cex"), xpd = TRUE)

text(4.5, 0.28, "QTL-trait association mediated by expression (e)", cex = 6/par("ps")/par("cex"), xpd = TRUE)
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)


# association for eQTL
par(las = 1, tcl = -0.2, mai = c(0.2, 0.22, 0.12, 0.02)*file.width/25.4/2, ps = 7, lwd = 0.5)
plot(c(0, 2), c(3.8, 4.5), type = "n", axes = F, xlab = "", ylab = "")

points(select.geno["3R_12642725", ]/2 + 0.5 + runif(length(select.geno["3R_12642725", ]), -0.25, 0.25), select.young.exp[, "FBgn0026207"], col = brewer.pal(9, "Set1")[3])

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(0.5, 1.5), label = c("ref", "del"), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(ylab = expression(paste(italic("mbo"), " expression in young males")), mgp = c(1.3, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = expression(paste(italic("3R"), ":12642725 genotype")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))

segments(c(0.4, 1.4), sapply(split(select.young.exp[, "FBgn0026207"], select.geno["3R_12642725", ]), mean, na.rm = T), c(0.6, 1.6), sapply(split(select.young.exp[, "FBgn0026207"], select.geno["3R_12642725", ]), mean, na.rm = T), lwd = 2)


eqtl.p <- formatC(summary(lm(select.young.exp[, "FBgn0026207"] ~ select.geno["3R_12642725", ]))$coefficients[2, 4], format = "e", digits = 2)
eqtl.eff <- formatC(summary(lm(select.young.exp[, "FBgn0026207"] ~ select.geno["3R_12642725", ]))$coefficients[2, 1]*2, format = "f", digits = 2)


text(1, 4.43, parse(text = paste("paste(beta, \" = \",", eqtl.eff, ",  \" (\", italic(P), \" = \", ", eqtl.p, ", \")\")")), pos = 3, cex = 6/par("ps")/par("cex"))

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# association for QTL

plot(c(0, 2), c(-0.2, 2.8), type = "n", axes = F, xlab = "", ylab = "")

points(select.geno["3R_12642725", ]/2 + 0.5 + runif(length(select.geno["3R_12642725", ]), -0.25, 0.25), select.pheno[, "male.speed.decline"], col = brewer.pal(9, "Set1")[8])

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(0.5, 1.5), label = c("ref", "del"), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(ylab = "Climbing speed decline in males", mgp = c(1.3, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")
title(xlab = expression(paste(italic("3R"), ":12642725 genotype")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))

segments(c(0.4, 1.4), sapply(split(select.pheno[, "male.speed.decline"], select.geno["3R_12642725", ]), mean, na.rm = T), c(0.6, 1.6), sapply(split(select.pheno[, "male.speed.decline"], select.geno["3R_12642725", ]), mean, na.rm = T), lwd = 2)

qtl.p <- formatC(summary(lm(select.pheno[, "male.speed.decline"] ~ select.geno["3R_12642725", ]))$coefficients[2, 4], format = "e", digits = 2)
qtl.eff <- formatC(summary(lm(select.pheno[, "male.speed.decline"] ~ select.geno["3R_12642725", ]))$coefficients[2, 1]*2, format = "f", digits = 2)
cat(qtl.eff, "\n")
summary(lm(select.pheno[, "male.speed.decline"] ~ select.geno["3R_12642725", ]))
# text(1, 2.5, parse(text = paste("paste(italic(P), \" = \", ", qtl.p, ")")), pos = 3, cex = 6/par("ps")/par("cex"))

text(1, 2.52, parse(text = paste("paste(beta, \" = \",", qtl.eff, ",  \" (\", italic(P), \" = \", ", qtl.p, ", \")\")")), pos = 3, cex = 6/par("ps")/par("cex"))


text(grconvertX(0.05 + file.width/2/25.4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# scatter plot, expression versus trait

plot(c(3.8, 4.5), c(-0.2, 2.8), type = "n", axes = F, xlab = "", ylab = "")
points(select.young.exp[, "FBgn0026207"], select.pheno[, "male.speed.decline"], col = "grey50")
axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(3.8, 4.4, 0.2), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
title(ylab = "Climbing speed decline in males", mgp = c(1.3, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(xlab = expression(paste(italic("mbo"), " expression in young males")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")

cor.est <- formatC(cor.test(select.young.exp[, "FBgn0026207"], select.pheno[, "male.speed.decline"])$estimate, format = "f", digits = 2)
cor.p <- formatC(cor.test(select.young.exp[, "FBgn0026207"], select.pheno[, "male.speed.decline"])$p.value, format = "e", digits = 2)
#text(4.15, 2.55, parse(text = paste("paste(italic(P), \" = \", ", cor.p, ")")), pos = 3, cex = 6/par("ps")/par("cex"))
cat("\n\np values and r are hard coded *************\n\n \n")
text(4.15, 2.55, expression(paste(italic(r), " = 0.32 (", italic(P), " = 2.81e-04)")), pos = 3, cex = 6/par("ps")/par("cex"))
text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# after adjustment

pheno.exp.fit <- lm(select.pheno[, "male.speed.decline"] ~ select.young.exp[, "FBgn0026207"])
fit.line <- rownames(pheno.exp.fit$model)
adj.pheno <- mean(pheno.exp.fit$model[, 1]) + residuals(pheno.exp.fit)

plot(c(0, 2), c(-0.2, 2.8), type = "n", axes = F, xlab = "", ylab = "")
points(select.geno["3R_12642725", match(fit.line, common.lines)]/2 + 0.5 + runif(length(select.geno["3R_12642725", match(fit.line, common.lines)]), -0.25, 0.25), adj.pheno, col = brewer.pal(9, "Set1")[8])

segments(c(0.4, 1.4), sapply(split(adj.pheno, select.geno["3R_12642725", match(fit.line, common.lines)]), mean, na.rm = T), c(0.6, 1.6), sapply(split(adj.pheno, select.geno["3R_12642725", match(fit.line, common.lines)]), mean, na.rm = T), lwd = 2)

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(0.5, 1.5), label = c("ref", "del"), cex.axis = 7/par("ps")/par("cex"))
axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))
#title(ylab = expression(paste("Climbing speed decline in males\n(adjusted for ", italic("mbo"), " expression)")), mgp = c(1.3, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = "Climbing speed decline in males", mgp = c(2.1, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
title(ylab = expression(paste("(adjusted for ", italic("mbo"), " expression)")), mgp = c(1.3, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))

box(bty = "l")

cat("\n\np values and r are hard coded *************\n\n \n")

adj.fit <- lm(adj.pheno ~ select.geno["3R_12642725", match(fit.line, common.lines)])
qtl.p <- formatC(summary(adj.fit)$coefficients[2, 4], format = "e", digits = 2)
qtl.eff <- formatC(summary(adj.fit)$coefficients[2, 1] * 2, format = "f", digits = 2)

#
#text(1, 2.52, parse(text = paste("paste(beta, \" = \",", qtl.eff, ",  \" (\", italic(P), \" = \", ", qtl.p, ", \")\")")), pos = 3, cex = 6/par("ps")/par("cex"))
text(1, 2.52, expression(paste(beta, " = 0.28 (", italic(P), " = 6.31e-03)")), pos = 3, cex = 6/par("ps")/par("cex"))

title(xlab = expression(paste(italic("3R"), ":12642725 genotype")), mgp = c(1, 0, 0), cex.lab = 7/par("ps")/par("cex"))

text(grconvertX(0.05 + file.width/2/25.4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("e")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)



dev.off()
