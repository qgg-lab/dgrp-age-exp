# =================================
# = figure for mediation analysis =
# =================================

args = commandArgs(TRUE) # args <- c("../figureData/female.example.mediation.RData")

# 3R_4680760, 0 = TGTACACTTTTCTTTCTAGATACAAAC, 2 = T
# FBgn0263346 FBgn0026620
# 3R3R:+:4659579-4712193  3R:-:4733654-4749056
# load data
load(args[1])
library("RColorBrewer")

# prepare file
# ============================================================

file.width = 89
cairo_pdf(file = args[2], width = file.width/25.4, height = file.width*1.25/25.4, family = args[3])
par(las = 1, tcl = -0.2, mai = c(0.15, 0.2, 0.1, 0.02)*file.width/25.4/2, ps = 7, lwd = 0.5)
layout(mat = matrix(c(1, 1, 2, 3, 4, 5), ncol = 2, byrow = TRUE), height = c(0.5, 1, 1))
set.seed(1)

# plot gene models
# ============================================================

par(mai = c(0.1, 0.01, 0.01, 0.01)*file.width/25.4/2)
plot(c(4650000, 4760000), c(0.2, 1.2), type = "n", axes = FALSE, xlab = "", ylab = "")
segments(4650000, 0.5, 4760000, 0.5, lwd = 1, col = "grey30")
rect(4659579, 0.4, 4712193, 0.6, col = brewer.pal(9, "Set1")[1])
rect(4733654, 0.4, 4749056, 0.6, col = brewer.pal(9, "Set1")[2])
text((4659579 + 4712193)/2, 0.5, expression(italic("smash")), col = "white", cex = 7/par("ps")/par("cex"))
text((4733654 + 4749056)/2, 0.5, expression(italic("tacc")), col = "white", cex = 7/par("ps")/par("cex"))
points(4680760, 0.65, pch = 6, cex = 0.8)
segments(4650000, 0.9, 4660000, 0.9, lwd = 1, col = "grey30")
text(4660000, 0.9, "10kb", cex = 7/par("ps")/par("cex"), pos = 4)
text(4650000, 0.2, expression(paste("Chromosome ", italic("3R"))), cex = 7/par("ps")/par("cex"), pos = 4)
text(4680760, 0.72, "27 bp deletion", cex = 7/par("ps")/par("cex"), pos = 4)

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# association in young and old emails with cis gene
par(las = 1, tcl = -0.2, mai = c(0.15, 0.2, 0.1, 0.02)*file.width/25.4/2, ps = 7, lwd = 0.5)

plot(c(0, 5), c(3.5, 6.5), type = "n", axes = F, xlab = "", ylab = "")

points(snp.geno/2 + 0.5 + runif(length(snp.geno), -0.25, 0.25), cis.gene.exp1, col = brewer.pal(9, "Reds")[6])

cis1.p <- formatC(summary(lm(cis.gene.exp1 ~ snp.geno))$coefficients[2, 4], format = "e", digits = 2)
text(1, 6.3, parse(text = paste("paste(italic(P), \" = \", ", cis1.p, ")")), pos = 3, cex = 6/par("ps")/par("cex"))


points(snp.geno/2 + 3.5 + runif(length(snp.geno), -0.25, 0.25), cis.gene.exp2, col = brewer.pal(9, "Reds")[9])

cis2.p <- formatC(summary(lm(cis.gene.exp2 ~ snp.geno))$coefficients[2, 4], format = "e", digits = 2)
text(4, 6.3, parse(text = paste("paste(italic(P), \" = \", ", cis2.p, ")")), pos = 3, cex = 6/par("ps")/par("cex"))


axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(0.5, 1.5), label = c("ref", "del"), cex.axis = 7/par("ps")/par("cex"))
axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(3.5, 4.5), label = c("ref", "del"), cex.axis = 7/par("ps")/par("cex"))
text(1, 3.5, "Young", pos = 3, cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[6])
text(4, 3.5, "Aged", pos = 3, cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = expression(paste(italic("smash"), " expression")), mgp = c(1.3, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# trans

plot(c(0, 5), c(6.2, 7), type = "n", axes = F, xlab = "", ylab = "")

points(snp.geno/2 + 0.5 + runif(length(snp.geno), -0.25, 0.25), trans.gene.exp1, col = brewer.pal(9, "Reds")[6])
trans1.p <- formatC(summary(lm(trans.gene.exp1 ~ snp.geno))$coefficients[2, 4], format = "e", digits = 2)
text(1, 6.94, parse(text = paste("paste(italic(P), \" = \", ", trans1.p, ")")), pos = 3, cex = 6/par("ps")/par("cex"))

points(snp.geno/2 + 3.5 + runif(length(snp.geno), -0.25, 0.25), trans.gene.exp2, col = brewer.pal(9, "Reds")[9])

#trans2.p <- formatC(summary(lm(trans.gene.exp2 ~ snp.geno))$coefficients[2, 4], format = "e", digits = 2)
#text(4, 6.94, parse(text = paste("paste(italic(P), \" = \", ", trans2.p, ")")), pos = 3, cex = 6/par("ps")/par("cex"))
#hard coded
text(4, 6.94, expression(paste(italic(P), " = 6.35e-03")), pos = 3, cex = 6/par("ps")/par("cex"))

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(0.5, 1.5), label = c("ref", "del"), cex.axis = 7/par("ps")/par("cex"))
axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(3.5, 4.5), label = c("ref", "del"), cex.axis = 7/par("ps")/par("cex"))
text(1, 6.2, "Young", pos = 3, cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[6])
text(4, 6.2, "Aged", pos = 3, cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = expression(paste(italic("tacc"), " expression")), mgp = c(1.3, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")
text(grconvertX(0.05 + file.width/2/25.4, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

## conditiona analysis
# ============================================================


## smash | tacc
plot(c(0, 5), c(3.5, 6.5), type = "n", axes = F, xlab = "", ylab = "")

data1 <- na.omit(data.frame(cis = cis.gene.exp1, trans = trans.gene.exp1, snp = snp.geno))
cond.trans.fit1 <- lm(cis ~ trans, data1)
cond.cis.exp1 <- mean(data1$cis) + residuals(cond.trans.fit1)
with(data1, points(snp/2 + 0.5 + runif(length(snp), -0.25, 0.25), cond.cis.exp1, col = brewer.pal(9, "Reds")[6]))

cis1.p <- formatC(summary(lm(cond.cis.exp1 ~ data1$snp))$coefficients[2, 4], format = "e", digits = 2)
text(1, 6.3, parse(text = paste("paste(italic(P), \" = \", ", cis1.p, ")")), pos = 3, cex = 6/par("ps")/par("cex"))

data2 <- na.omit(data.frame(cis = cis.gene.exp2, trans = trans.gene.exp2, snp = snp.geno))
cond.trans.fit2 <- lm(cis ~ trans, data2)
cond.cis.exp2 <- mean(data2$cis) + residuals(cond.trans.fit2)
with(data2, points(snp/2 + 3.5 + runif(length(snp), -0.25, 0.25), cond.cis.exp2, col = brewer.pal(9, "Reds")[9]))

cis2.p <- formatC(summary(lm(cond.cis.exp2 ~ data2$snp))$coefficients[2, 4], format = "e", digits = 2)
text(4, 6.3, parse(text = paste("paste(italic(P), \" = \", ", cis2.p, ")")), pos = 3, cex = 6/par("ps")/par("cex"))

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(0.5, 1.5), label = c("ref", "del"), cex.axis = 7/par("ps")/par("cex"))
axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(3.5, 4.5), label = c("ref", "del"), cex.axis = 7/par("ps")/par("cex"))
text(1, 3.5, "Young", pos = 3, cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[6])
text(4, 3.5, "Aged", pos = 3, cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = expression(paste(italic("smash"), "|", italic("tacc"), " expression")), mgp = c(1.3, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("d")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)



plot(c(0, 5), c(6.2, 7), type = "n", axes = F, xlab = "", ylab = "")

data1 <- na.omit(data.frame(cis = cis.gene.exp1, trans = trans.gene.exp1, snp = snp.geno))
cond.cis.fit1 <- lm(trans ~ cis, data1)
cond.trans.exp1 <- mean(data1$trans) + residuals(cond.cis.fit1)
with(data1, points(snp/2 + 0.5 + runif(length(snp), -0.25, 0.25), cond.trans.exp1, col = brewer.pal(9, "Reds")[6]))

trans1.p <- formatC(summary(lm(cond.trans.exp1 ~ data1$snp))$coefficients[2, 4], format = "f", digits = 2)
text(1, 6.94, parse(text = paste("paste(italic(P), \" = \", ", trans1.p, ")")), pos = 3, cex = 6/par("ps")/par("cex"))

data2 <- na.omit(data.frame(cis = cis.gene.exp2, trans = trans.gene.exp2, snp = snp.geno))
cond.cis.fit2 <- lm(trans ~ cis, data2)
cond.trans.exp2 <- mean(data2$trans) + residuals(cond.cis.fit2)
with(data2, points(snp/2 + 3.5 + runif(length(snp), -0.25, 0.25), cond.trans.exp2, col = brewer.pal(9, "Reds")[9]))

trans2.p <- formatC(summary(lm(cond.trans.exp2 ~ data2$snp))$coefficients[2, 4], format = "f", digits = 2)
text(4, 6.94, parse(text = paste("paste(italic(P), \" = \", ", trans2.p, ")")), pos = 3, cex = 6/par("ps")/par("cex"))

axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(0.5, 1.5), label = c("ref", "del"), cex.axis = 7/par("ps")/par("cex"))
axis(side = 1, lwd = 0.5, mgp = c(0.8, 0, 0), at = seq(3.5, 4.5), label = c("ref", "del"), cex.axis = 7/par("ps")/par("cex"))
text(1, 6.2, "Young", pos = 3, cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[6])
text(4, 6.2, "Aged", pos = 3, cex = 7/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

axis(side = 2, mgp = c(0.8, 0.3, 0), lwd = 0.5, cex.axis = 7/par("ps")/par("cex"))

title(ylab = expression(paste(italic("tacc"), "|", italic("smash"), " expression")), mgp = c(1.3, 0.3, 0), cex.lab = 7/par("ps")/par("cex"))
box(bty = "l")





text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("e")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)



dev.off()
