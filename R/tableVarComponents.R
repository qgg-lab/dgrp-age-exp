# ============================================================
# = make a table for variance components as well as GxE etc. =
# ============================================================

args <- commandArgs(TRUE) # args <- c("qg/female.adj.qgSingleTemp.RData", "qg/male.adj.qgSingleTemp.RData", "qg/female.adj.qgVarHet.RData", "qg/male.adj.qgVarHet.RData", "qg/female.adj.qgGxE.RData", "qg/male.adj.qgGxE.RData", "joinExp/gene.info", "report/Table_VarComp.female.xlsx", "report/Table_VarComp.male.xlsx")
library(openxlsx)

# read data
# ============================================================
load(args[1])
female.single.temp <- as.data.frame(single.temp, stringsAsFactors = FALSE)
load(args[3])
female.var.het <- as.data.frame(var.het, stringsAsFactors = FALSE)
load(args[5])
female.gxe <- as.data.frame(gxe, stringsAsFactors = FALSE)

gene.info <- read.table(args[7], header = FALSE, as.is = TRUE, quote = "\"", row.names = 1)

# combine data
# ============================================================

female.table <- cbind(female.single.temp$wolba.var.line.25c, female.single.temp$wolba.var.e.25c, female.single.temp$wolba.p.line.25c, p.adjust(female.single.temp$wolba.p.line.25c, method = "BH"),
                      female.single.temp$wolba.var.line.18c, female.single.temp$wolba.var.e.18c, female.single.temp$wolba.p.line.18c, p.adjust(female.single.temp$wolba.p.line.18c, method = "BH"),
                      as.numeric(female.var.het$wolba.p.single.g), p.adjust(as.numeric(female.var.het$wolba.p.single.g), method = "BH"),
                      as.numeric(female.var.het$wolba.p.single.e), p.adjust(as.numeric(female.var.het$wolba.p.single.e), method = "BH"),
                      female.gxe$wolba.var.g, female.gxe$wolba.var.gxe, female.gxe$wolba.p.g, p.adjust(female.gxe$wolba.p.g, method = "BH"),
                      female.gxe$wolba.p.gxe, p.adjust(female.gxe$wolba.p.gxe, method = "BH"))

# output numbers
# ============================================================
female.int.fdr <- p.adjust(as.numeric(female.gxe[, 10]), method = "BH")
female.var.het.fdr <- p.adjust(as.numeric(var.het[, 12]), method = "BH")

cat("there are ", sum(female.int.fdr < 0.05 & female.var.het.fdr < 0.05, na.rm = TRUE), " genes with int FDR < 0.05 & var.het FDR < 0.05 in females.\n")

final.table.female <- cbind(gene.name, gene.info[gene.name, ], female.table)

colnames(final.table.female) <- c("FBgn", "symbol", "type", "pos",
                           "female.var.line.25c", "female.var.e.25c", "female.p.line.25c", "female.fdr.line.25c", 
                           "female.var.line.18c", "female.var.e.18c", "female.p.line.18c", "female.fdr.line.18c",
                           "female.p.single.g", "female.fdr.single.g", "female.p.single.e", "female.fdr.single.e",
                           "female.var.line", "female.var.int", "female.p.g", "female.fdr.g", "female.p.gxe", "female.fdr.gxe")
# write excel file
# ============================================================

write.xlsx(final.table.female, file = args[8], col.names = TRUE, row.names = FALSE)


# For males
# ============================================================
load(args[2])
male.single.temp <- as.data.frame(single.temp, stringsAsFactors = FALSE)
load(args[4])
male.var.het <- as.data.frame(var.het, stringsAsFactors = FALSE)
load(args[6])
male.gxe <- as.data.frame(gxe, stringsAsFactors = FALSE)

male.table <- cbind(male.single.temp$wolba.var.line.25c, male.single.temp$wolba.var.e.25c, male.single.temp$wolba.p.line.25c, p.adjust(male.single.temp$wolba.p.line.25c, method = "BH"),
                    male.single.temp$wolba.var.line.18c, male.single.temp$wolba.var.e.18c, male.single.temp$wolba.p.line.18c, p.adjust(male.single.temp$wolba.p.line.18c, method = "BH"),
                    as.numeric(male.var.het$wolba.p.single.g), p.adjust(as.numeric(male.var.het$wolba.p.single.g), method = "BH"),
                    as.numeric(male.var.het$wolba.p.single.e), p.adjust(as.numeric(male.var.het$wolba.p.single.e), method = "BH"),
                    male.gxe$wolba.var.g, male.gxe$wolba.var.gxe, male.gxe$wolba.p.g, p.adjust(male.gxe$wolba.p.g, method = "BH"),
                    male.gxe$wolba.p.gxe, p.adjust(male.gxe$wolba.p.gxe, method = "BH"))

male.int.fdr <- p.adjust(as.numeric(male.gxe[, 10]), method = "BH")
male.var.het.fdr <- p.adjust(as.numeric(var.het[, 12]), method = "BH")

cat("there are ", sum(male.int.fdr < 0.05 & male.var.het.fdr < 0.05, na.rm = TRUE), " genes with int FDR < 0.05 & var.het FDR < 0.05 in males.\n")

final.table.male <- cbind(gene.name, gene.info[gene.name, ], male.table)

colnames(final.table.male) <- c("FBgn", "symbol", "type", "pos",
                           "male.var.line.25c", "male.var.e.25c", "male.p.line.25c", "male.fdr.line.25c", 
                           "male.var.line.18c", "male.var.e.18c", "male.p.line.18c", "male.fdr.line.18c",
                           "male.p.single.g", "male.fdr.single.g", "male.p.single.e", "male.fdr.single.e",
                           "male.var.line", "male.var.int", "male.p.g", "male.fdr.g", "male.p.gxe", "male.fdr.gxe")

# write excel file
# ============================================================

write.xlsx(final.table.male, file = args[9], col.names = TRUE, row.names = FALSE)
