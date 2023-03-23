# ============================
# = figure for response eQTL =
# ============================

args <- commandArgs(TRUE) # args <- c("../figureData/female.18c.fdr05.eqtl.table.txt", "../figureData/female.25c.fdr05.eqtl.table.txt", "../figureData/male.18c.fdr05.eqtl.table.txt", "../figureData/male.25c.fdr05.eqtl.table.txt", "../figureData/gene.info", "../figureData/female.18c.genvar.id", "../figureData/female.25c.genvar.id", "../figureData/male.18c.genvar.id", "../figureData/male.25c.genvar.id", "../figureData/female.adj.qgGxE.RData", "../figureData/male.adj.qgGxE.RData", "../figureData/female.tf.summary.RData", "../figureData/male.tf.summary.RData")
library("RColorBrewer")
library("plotrix")

# read data
# ============================================================

female.18c <- read.table(args[1], header = FALSE, as.is = TRUE)
female.25c <- read.table(args[2], header = FALSE, as.is = TRUE)
male.18c <- read.table(args[3], header = FALSE, as.is = TRUE)
male.25c <- read.table(args[4], header = FALSE, as.is = TRUE)
rownames(female.18c) <- female.18c[, 1]
rownames(female.25c) <- female.25c[, 1]
rownames(male.18c) <- male.18c[, 1]
rownames(male.25c) <- male.25c[, 1]


cat("eGenes in young females: ", nrow(female.25c), "\n")
cat("eGenes in old females: ", nrow(female.18c), "\n")
cat("eGenes in both females: ", length(intersect(female.25c[, 1], female.18c[, 1])), "\n")

cat("eGenes in young males: ", nrow(male.25c), "\n")
cat("eGenes in old males: ", nrow(male.18c), "\n")
cat("eGenes in both males: ", length(intersect(male.25c[, 1], male.18c[, 1])), "\n")



gene.info <- read.table(args[5], header = FALSE, as.is = TRUE, quote = "\"", row.names = 1)

female.18c.genvar <- scan(args[6], what = "", quiet = TRUE)
female.25c.genvar <- scan(args[7], what = "", quiet = TRUE)
male.18c.genvar <- scan(args[8], what = "", quiet = TRUE)
male.25c.genvar <- scan(args[9], what = "", quiet = TRUE)
load(args[10])
female.gxe <- gxe
load(args[11])
male.gxe <- gxe


load(args[12])

female.tf.fold.enrich <- tf.fold.enrich[tf.list != "mE1_TFBS_HSA" & grepl("TFBS", tf.list), ]
female.tf.fold.enrich.p <- tf.fold.enrich.p[tf.list != "mE1_TFBS_HSA" & grepl("TFBS", tf.list), ]

load(args[13])
male.tf.fold.enrich <- tf.fold.enrich[tf.list != "mE1_TFBS_HSA" & grepl("TFBS", tf.list), ]
male.tf.fold.enrich.p <- tf.fold.enrich.p[tf.list != "mE1_TFBS_HSA" & grepl("TFBS", tf.list), ]

tf.list <- tf.list[tf.list != "mE1_TFBS_HSA" & grepl("TFBS", tf.list)]

range(female.tf.fold.enrich)
range(male.tf.fold.enrich)

cat("females share eqtl:", sum(apply(cbind(female.18c[intersect(female.25c[, 1], female.18c[, 1]), c(2, 4)], female.25c[intersect(female.25c[, 1], female.18c[, 1]), c(2, 4)]), 1, function(x) { return( length(intersect(unlist(strsplit(x[1], split = ",")), unlist(strsplit(x[4], split = ",")))) + length(intersect(unlist(strsplit(x[2], split = ",")), unlist(strsplit(x[3], split = ",")))) ) }) > 0), "\n")

cat("males share eqtl:", sum(apply(cbind(male.18c[intersect(male.25c[, 1], male.18c[, 1]), c(2, 4)], male.25c[intersect(male.25c[, 1], male.18c[, 1]), c(2, 4)]), 1, function(x) { return( length(intersect(unlist(strsplit(x[1], split = ",")), unlist(strsplit(x[4], split = ",")))) + length(intersect(unlist(strsplit(x[2], split = ",")), unlist(strsplit(x[3], split = ",")))) ) }) > 0), "\n")

# set up file
# ============================================================

file.width = 89
cairo_pdf(file = args[14], width = file.width/25.4, height = file.width/25.4, family = args[15])
par(las = 1, tcl = -0.2, mar = c(0.5, 0.5, 0.5, 0.5), ps = 7, lwd = 0.5)
layout(mat = matrix(c(1, 2, 3, 3), ncol = 2, byrow = TRUE), height = c(1, 1.2))


# make venn diagram for females
# ============================================================

plot(c(0, 1), c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
draw.circle(0.3, 0.5, 0.3, col = paste(brewer.pal(9, "Reds")[6], "AA", sep = ""))
draw.circle(0.6, 0.5, 0.3, col = paste(brewer.pal(9, "Reds")[9], "AA", sep = ""))
text(0.05, 0.8, "\u2640", cex = 12/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9], pos = 3, family = "Arial")
text(-0.02, 0.05, "Young eGenes", cex = 7/par("ps")/par("cex"), pos = 4)
text(0.95, 0.93, "Aged eGenes", cex = 7/par("ps")/par("cex"), pos = 2)

# young specific
text(0, 0.5, formatC(length(setdiff(female.25c[, 1], female.18c[, 1])), big.mark = ","), col = "white", cex = 7/par("ps")/par("cex"), pos = 4)
text(0.45, 0.5, formatC(length(intersect(female.25c[, 1], female.18c[, 1])), big.mark = ","), col = "white", cex = 7/par("ps")/par("cex"))
text(0.9, 0.5, formatC(length(setdiff(female.18c[, 1], female.25c[, 1])), big.mark = ","), col = "white", cex = 7/par("ps")/par("cex"), pos = 2)

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("a")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# make venn diagram for males
# ============================================================
plot(c(0, 1), c(0, 1), axes = F, type = "n", xlab = "", ylab = "")
draw.circle(0.3, 0.5, 0.3, col = paste(brewer.pal(9, "Blues")[6], "AA", sep = ""))
draw.circle(0.6, 0.5, 0.3, col = paste(brewer.pal(9, "Blues")[9], "AA", sep = ""))
text(0.05, 0.8, "\u2642", cex = 12/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9], pos = 3, family = "Arial")
text(-0.02, 0.05, "Young eGenes", cex = 7/par("ps")/par("cex"), pos = 4)
text(0.95, 0.93, "Aged eGenes", cex = 7/par("ps")/par("cex"), pos = 2)

# young specific
text(0, 0.5, formatC(length(setdiff(male.25c[, 1], male.18c[, 1])), big.mark = ","), col = "white", cex = 7/par("ps")/par("cex"), pos = 4)
text(0.45, 0.5, formatC(length(intersect(male.25c[, 1], male.18c[, 1])), big.mark = ","), col = "white", cex = 7/par("ps")/par("cex"))
text(0.9, 0.5, formatC(length(setdiff(male.18c[, 1], male.25c[, 1])), big.mark = ","), col = "white", cex = 7/par("ps")/par("cex"), pos = 2)

text(grconvertX(0.05 + file.width/25.4/2, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("b")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

# TF enrichment
# ============================================================

par(mar = c(0, 0.2, 0.1, 0))

tf.name <- gsub(".*_TFBS_", "", tf.list)

plot(c(-10, length(tf.name) + 1), c(0, 14), type = "n", axes = FALSE, xlab = "", ylab = "")

image(1:length(tf.name), 1:8, cbind(male.tf.fold.enrich[, 3], female.tf.fold.enrich[, 3], NA, male.tf.fold.enrich[, 2], female.tf.fold.enrich[, 2], NA, male.tf.fold.enrich[, 1], female.tf.fold.enrich[, 1]), add = T, breaks = c(seq(0, 1, 0.2), 2:4), col = c(brewer.pal(9, "Blues")[c(8, 6, 4, 2, 1)], brewer.pal(9, "Reds")[c(1, 4, 8)]))

segments(0.5, seq(0.5, 8.5, 1), length(tf.name) + 0.5, seq(0.5, 8.5, 1))
segments(seq(0.5, length(tf.name) + 0.5, 1), 0.5, seq(0.5, length(tf.name) + 0.5, 1), 2.5)
segments(seq(0.5, length(tf.name) + 0.5, 1), 3.5, seq(0.5, length(tf.name) + 0.5, 1), 5.5)
segments(seq(0.5, length(tf.name) + 0.5, 1), 6.5, seq(0.5, length(tf.name) + 0.5, 1), 8.5)

for (i in 1:length(tf.name)) {
	
	this.tf <- tf.name[i]
	text(i - 0.75, 8.8, parse(text = paste("italic(\"", this.tf, "\")", sep = "")), srt = 60, pos = 4, cex = 7/par("ps")/par("cex"))
	#text(i - 0.6, 9, parse(text = paste("expression(italic(", this.tf, "))", sep = "")), pos = 4)
	
}

p.adj <- matrix(p.adjust(c(female.tf.fold.enrich.p, male.tf.fold.enrich.p), method = "BH"), ncol = 6, byrow = FALSE)
female.tf.fold.enrich.p <- p.adj[, 1:3]
male.tf.fold.enrich.p <- p.adj[, 4:6]


text(which(male.tf.fold.enrich.p[, 3] < 0.05), 1, "**", cex = 7/par("ps")/par("cex"))
text(which(male.tf.fold.enrich.p[, 3] >= 0.05 & male.tf.fold.enrich.p[, 3] < 0.1), 1, "*", cex = 7/par("ps")/par("cex"))
text(which(female.tf.fold.enrich.p[, 3] < 0.05), 2, "**", cex = 7/par("ps")/par("cex"))
text(which(female.tf.fold.enrich.p[, 3] >= 0.05 & female.tf.fold.enrich.p[, 3] < 0.1), 2, "*", cex = 7/par("ps")/par("cex"))

text(which(male.tf.fold.enrich.p[, 2] < 0.05), 4, "**", cex = 7/par("ps")/par("cex"))
text(which(male.tf.fold.enrich.p[, 2] >= 0.05 & male.tf.fold.enrich.p[, 3] < 0.1), 4, "*", cex = 7/par("ps")/par("cex"))
text(which(female.tf.fold.enrich.p[, 2] < 0.05), 5, "**", cex = 7/par("ps")/par("cex"))
text(which(female.tf.fold.enrich.p[, 2] >= 0.05 & female.tf.fold.enrich.p[, 3] < 0.1), 5, "*", cex = 7/par("ps")/par("cex"))

text(which(male.tf.fold.enrich.p[, 1] < 0.05), 7, "**", cex = 7/par("ps")/par("cex"))
text(which(male.tf.fold.enrich.p[, 1] >= 0.05 & male.tf.fold.enrich.p[, 3] < 0.1), 7, "*", cex = 7/par("ps")/par("cex"))
text(which(female.tf.fold.enrich.p[, 1] < 0.05), 8, "**", cex = 7/par("ps")/par("cex"))
text(which(female.tf.fold.enrich.p[, 1] >= 0.05 & female.tf.fold.enrich.p[, 3] < 0.1), 8, "*", cex = 7/par("ps")/par("cex"))

text(-11.7, 7.5, "Young only vs random", cex = 7/par("ps")/par("cex"), pos = 4)
text(-11.7, 4.5, "Common vs young only", cex = 7/par("ps")/par("cex"), pos = 4)
text(-11.7, 1.5, "Aged only vs young only", cex = 7/par("ps")/par("cex"), pos = 4)


text(0, 1, "\u2642", family = "Arial", cex = 5/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])
text(0, 2, "\u2640", family = "Arial", cex = 5/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

text(0, 4, "\u2642", family = "Arial", cex = 5/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])
text(0, 5, "\u2640", family = "Arial", cex = 5/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

text(0, 7, "\u2642", family = "Arial", cex = 5/par("ps")/par("cex"), col = brewer.pal(9, "Blues")[9])
text(0, 8, "\u2640", family = "Arial", cex = 5/par("ps")/par("cex"), col = brewer.pal(9, "Reds")[9])

text(grconvertX(0.05, from = "inches", to = "user"), grconvertY(1, from = "nfc", to = "user"), expression(bold("c")), cex = 9/par("ps")/par("cex"), pos = 1, xpd = TRUE)

image(-10:-3, c(12, 13), matrix(c(seq(0.1, 0.9, 0.2), seq(1.5, 3.5, 1), rep(-0.1, 8)), ncol = 2), add = T, breaks = c(seq(-0.2, 1, 0.2), 2:4), col = c("white", brewer.pal(9, "Blues")[c(8, 6, 4, 2, 1)], brewer.pal(9, "Reds")[c(1, 4, 8)]))
segments(-10.5, 11.5, -2.5, 11.5)
segments(-10.5, 12.5, -2.5, 12.5)
segments(seq(-10.5, -2.5, 1), 11.5, seq(-10.5, -2.5, 1), 12.5)
text(seq(-10.5, -2.5, 1) + 1, 11, c(seq(0, 1, 0.2), 2:4), srt = 60, pos = 2)
text(-6.5, 13.2, "Fold enrichment", cex = 7/par("ps")/par("cex"))

text(7, 13.5, "* FDR = 0.10", cex = 7/par("ps")/par("cex"), pos = 2)
text(7, 12.25, "** FDR = 0.05", cex = 7/par("ps")/par("cex"), pos = 2)


dev.off()
