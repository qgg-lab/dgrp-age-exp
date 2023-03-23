# ======================================
# = table for eqtl and model selection =
# ======================================

args <- commandArgs(TRUE) # args <- c("../figureData/female.18c.fdr05.eqtl.table.txt", "../figureData/female.25c.fdr05.eqtl.table.txt", "../figureData/male.18c.fdr05.eqtl.table.txt", "../figureData/male.25c.fdr05.eqtl.table.txt", "../figureData/gene.info", "../figureData/female.18c.genvar.id", "../figureData/female.25c.genvar.id", "../figureData/male.18c.genvar.id", "../figureData/male.25c.genvar.id", "../figureData/female.adj.qgGxE.RData", "../figureData/male.adj.qgGxE.RData")
library(openxlsx)

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

gene.info <- read.table(args[5], header = FALSE, as.is = TRUE, quote = "\"", row.names = 1)

female.18c.genvar <- scan(args[6], what = "", quiet = TRUE)
female.25c.genvar <- scan(args[7], what = "", quiet = TRUE)
male.18c.genvar <- scan(args[8], what = "", quiet = TRUE)
male.25c.genvar <- scan(args[9], what = "", quiet = TRUE)
load(args[10])
female.gxe <- gxe
load(args[11])
male.gxe <- gxe


# combine data
# ============================================================

final.table <- cbind(c(rep("Female old", nrow(female.18c)), rep("Female young", nrow(female.25c)),
											 rep("Male old", nrow(male.18c)), rep("Male young", nrow(male.25c))),
											 rbind(female.18c[, -3], female.25c[, -3], male.18c[, -3], male.25c[, -3]))
final.table <- cbind(final.table[, 1:2], gene.info[final.table[, 2], ], final.table[, -(1:2)])	

# output table
# ============================================================

cat("number of selected eqtls:\n")
cat(table(sapply(strsplit(final.table[, 6], split = ","), length)), "\n")

cat("number of eqtls mapped among shared gen var genes in females:\n")
cat("18C:", length(intersect(intersect(female.18c.genvar, female.25c.genvar), female.18c[, 1])), "\n")
cat("25C:", length(intersect(intersect(female.18c.genvar, female.25c.genvar), female.25c[, 1])), "\n")
cat("bothC:", length(intersect(intersect(intersect(female.18c.genvar, female.25c.genvar), female.25c[, 1]), female.18c[, 1])), "\n")
cat("share eqtl:", sum(apply(cbind(female.18c[intersect(intersect(intersect(female.18c.genvar, female.25c.genvar), female.25c[, 1]), female.18c[, 1]), c(2, 4)], female.25c[intersect(intersect(intersect(female.18c.genvar, female.25c.genvar), female.25c[, 1]), female.18c[, 1]), c(2, 4)]), 1, function(x) { return( length(intersect(unlist(strsplit(x[1], split = ",")), unlist(strsplit(x[4], split = ",")))) + length(intersect(unlist(strsplit(x[2], split = ",")), unlist(strsplit(x[3], split = ",")))) ) }) > 0), "\n")

female.gene.share.eqtl <- intersect(intersect(intersect(female.18c.genvar, female.25c.genvar), female.25c[, 1]), female.18c[, 1])[apply(cbind(female.18c[intersect(intersect(intersect(female.18c.genvar, female.25c.genvar), female.25c[, 1]), female.18c[, 1]), c(2, 4)], female.25c[intersect(intersect(intersect(female.18c.genvar, female.25c.genvar), female.25c[, 1]), female.18c[, 1]), c(2, 4)]), 1, function(x) { return( length(intersect(unlist(strsplit(x[1], split = ",")), unlist(strsplit(x[4], split = ",")))) + length(intersect(unlist(strsplit(x[2], split = ",")), unlist(strsplit(x[3], split = ",")))) ) }) > 0]

female.int.fdr <- p.adjust(as.numeric(female.gxe[, 10]), method = "BH")
cat("share eqtl and GxE sig:", length(intersect(gene.name[female.int.fdr < 0.05], female.gene.share.eqtl)), "\n")



cat("number of eqtls mapped among shared gen var genes in males:\n")
cat("18C:", length(intersect(intersect(male.18c.genvar, male.25c.genvar), male.18c[, 1])), "\n")
cat("25C:", length(intersect(intersect(male.18c.genvar, male.25c.genvar), male.25c[, 1])), "\n")
cat("bothC:", length(intersect(intersect(intersect(male.18c.genvar, male.25c.genvar), male.25c[, 1]), male.18c[, 1])), "\n")
cat("share eqtl:", sum(apply(cbind(male.18c[intersect(intersect(intersect(male.18c.genvar, male.25c.genvar), male.25c[, 1]), male.18c[, 1]), c(2, 4)], male.25c[intersect(intersect(intersect(male.18c.genvar, male.25c.genvar), male.25c[, 1]), male.18c[, 1]), c(2, 4)]), 1, function(x) { return( length(intersect(unlist(strsplit(x[1], split = ",")), unlist(strsplit(x[4], split = ",")))) + length(intersect(unlist(strsplit(x[2], split = ",")), unlist(strsplit(x[3], split = ",")))) ) }) > 0), "\n")

male.gene.share.eqtl <- intersect(intersect(intersect(male.18c.genvar, male.25c.genvar), male.25c[, 1]), male.18c[, 1])[apply(cbind(male.18c[intersect(intersect(intersect(male.18c.genvar, male.25c.genvar), male.25c[, 1]), male.18c[, 1]), c(2, 4)], male.25c[intersect(intersect(intersect(male.18c.genvar, male.25c.genvar), male.25c[, 1]), male.18c[, 1]), c(2, 4)]), 1, function(x) { return( length(intersect(unlist(strsplit(x[1], split = ",")), unlist(strsplit(x[4], split = ",")))) + length(intersect(unlist(strsplit(x[2], split = ",")), unlist(strsplit(x[3], split = ",")))) ) }) > 0]

male.int.fdr <- p.adjust(as.numeric(male.gxe[, 10]), method = "BH")
cat("share eqtl and GxE sig:", length(intersect(gene.name[male.int.fdr < 0.05], male.gene.share.eqtl)), "\n")





colnames(final.table) <- c("condition", "FBgn", "symbol", "type", "pos", "selected.eqtl", "eqtl", "pval", "eff")

final.table[nchar(final.table[, 7]) > 32767, 7] <- "too many to display"
final.table[nchar(final.table[, 8]) > 32767, 8] <- "too many to display"
final.table[nchar(final.table[, 9]) > 32767, 9] <- "too many to display"

                   

# write excel file
# ============================================================

write.xlsx(list('table' = final.table), file = args[12], colNames = TRUE, rowNames = FALSE)
