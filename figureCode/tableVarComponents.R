# ============================================================
# = make a table for variance components as well as GxE etc. =
# ============================================================

args <- commandArgs(TRUE) # args <- c("../figureData/female.adj.qgSingleTemp.RData", "../figureData/male.adj.qgSingleTemp.RData", "../figureData/female.adj.qgVarHet.RData", "../figureData/male.adj.qgVarHet.RData", "../figureData/female.adj.qgGxE.RData", "../figureData/male.adj.qgGxE.RData", "../figureData/gene.info", "../figureData/female.adj.qgAgeEff.RData", "../figureData/male.adj.qgAgeEff.RData", "../figure/Table_VarComp.xlsx")
library(openxlsx)

# read data
# ============================================================
load(args[1])
female.single.temp <- as.data.frame(single.temp, stringsAsFactors = FALSE)
female.gene.name <- gene.name
load(args[2])
male.single.temp <- as.data.frame(single.temp, stringsAsFactors = FALSE)
male.gene.name <- gene.name

load(args[3])
female.var.het <- as.data.frame(var.het, stringsAsFactors = FALSE)
load(args[4])
male.var.het <- as.data.frame(var.het, stringsAsFactors = FALSE)

load(args[5])
female.gxe <- as.data.frame(gxe, stringsAsFactors = FALSE)
load(args[6])
male.gxe <- as.data.frame(gxe, stringsAsFactors = FALSE)

gene.info <- read.table(args[7], header = FALSE, as.is = TRUE, quote = "\"", row.names = 1)

load(args[8])
female.age.effect <- age.effect
load(args[9])
male.age.effect <- age.effect



# combine data
# ============================================================

female.table <- cbind(female.age.effect[, 4:6], p.adjust(female.age.effect[, 6], method = "BH"),
											female.single.temp$wolba.var.line.25c, female.single.temp$wolba.var.e.25c, female.single.temp$wolba.var.line.25c/(female.single.temp$wolba.var.line.25c + female.single.temp$wolba.var.e.25c), female.single.temp$wolba.p.line.25c, p.adjust(female.single.temp$wolba.p.line.25c, method = "BH"),
                      female.single.temp$wolba.var.line.18c, female.single.temp$wolba.var.e.18c, female.single.temp$wolba.var.line.18c/(female.single.temp$wolba.var.line.18c + female.single.temp$wolba.var.e.18c), female.single.temp$wolba.p.line.18c, p.adjust(female.single.temp$wolba.p.line.18c, method = "BH"),
                      as.numeric(female.var.het$wolba.p.single.g), p.adjust(as.numeric(female.var.het$wolba.p.single.g), method = "BH"),
                      as.numeric(female.var.het$wolba.p.single.e), p.adjust(as.numeric(female.var.het$wolba.p.single.e), method = "BH"),
                      female.gxe$wolba.var.g, female.gxe$wolba.var.gxe, female.gxe$wolba.var.gxe/(female.gxe$wolba.var.gxe + female.gxe$wolba.var.g), female.gxe$wolba.p.g, p.adjust(female.gxe$wolba.p.g, method = "BH"),
                      female.gxe$wolba.p.gxe, p.adjust(female.gxe$wolba.p.gxe, method = "BH"))

cat("there are ", sum(female.table[, 8] < 0.05, na.rm = TRUE), " genes with int FDR < 0.05 at 25c (young) in females.\n")
cat("there are ", sum(female.table[, 12] < 0.05, na.rm = TRUE), " genes with int FDR < 0.05 at 18c (old) in females.\n")

# add proportion due to rank order change
# ============================================================

female.sigma.young <- sqrt(as.numeric(female.table[, 5]))
female.sigma.old <- sqrt(as.numeric(female.table[, 10]))
female.sigma2.g <- as.numeric(female.table[, 19])
female.sigma2.ga <- as.numeric(female.table[, 20])

female.young.sig <- as.numeric(female.table[, 9] <= 0.05)
female.old.sig <- as.numeric(female.table[, 14] <= 0.05)

female.change.rank <- female.sigma.young * female.sigma.old * (1 - female.sigma2.g/(female.sigma.young * female.sigma.old))
female.change.rank <- ifelse(female.change.rank < 0, 0, female.change.rank)
female.change.rank[!female.young.sig | !female.old.sig] <- NA

female.change.var <- (female.sigma.young - female.sigma.old)^2/2
female.change.var[!female.young.sig | !female.old.sig] <- NA

female.rank.percent <- female.change.rank/(female.change.rank + female.change.var)


male.table <- cbind(male.age.effect[, 4:6], p.adjust(male.age.effect[, 6], method = "BH"),
											male.single.temp$wolba.var.line.25c, male.single.temp$wolba.var.e.25c, male.single.temp$wolba.var.line.25c/(male.single.temp$wolba.var.line.25c + male.single.temp$wolba.var.e.25c), male.single.temp$wolba.p.line.25c, p.adjust(male.single.temp$wolba.p.line.25c, method = "BH"),
                      male.single.temp$wolba.var.line.18c, male.single.temp$wolba.var.e.18c, male.single.temp$wolba.var.line.18c/(male.single.temp$wolba.var.line.18c + male.single.temp$wolba.var.e.18c), male.single.temp$wolba.p.line.18c, p.adjust(male.single.temp$wolba.p.line.18c, method = "BH"),
                      as.numeric(male.var.het$wolba.p.single.g), p.adjust(as.numeric(male.var.het$wolba.p.single.g), method = "BH"),
                      as.numeric(male.var.het$wolba.p.single.e), p.adjust(as.numeric(male.var.het$wolba.p.single.e), method = "BH"),
                      male.gxe$wolba.var.g, male.gxe$wolba.var.gxe, male.gxe$wolba.var.gxe/(male.gxe$wolba.var.gxe + male.gxe$wolba.var.g), male.gxe$wolba.p.g, p.adjust(male.gxe$wolba.p.g, method = "BH"),
                      male.gxe$wolba.p.gxe, p.adjust(male.gxe$wolba.p.gxe, method = "BH"))


cat("there are ", sum(male.table[, 8] < 0.05, na.rm = TRUE), " genes with int FDR < 0.05 at 25c (young) in males.\n")
cat("there are ", sum(male.table[, 12] < 0.05, na.rm = TRUE), " genes with int FDR < 0.05 at 18c (old) in males.\n")

# add proportion due to rank order change
# ============================================================

male.sigma.young <- sqrt(as.numeric(male.table[, 5]))
male.sigma.old <- sqrt(as.numeric(male.table[, 10]))
male.sigma2.g <- as.numeric(male.table[, 19])
male.sigma2.ga <- as.numeric(male.table[, 20])

male.young.sig <- as.numeric(male.table[, 9] <= 0.05)
male.old.sig <- as.numeric(male.table[, 14] <= 0.05)

male.change.rank <- male.sigma.young * male.sigma.old * (1 - male.sigma2.g/(male.sigma.young * male.sigma.old))
male.change.rank <- ifelse(male.change.rank < 0, 0, male.change.rank)
male.change.rank[!male.young.sig | !male.old.sig] <- NA

male.change.var <- (male.sigma.young - male.sigma.old)^2/2
male.change.var[!male.young.sig | !male.old.sig] <- NA

male.rank.percent <- male.change.rank/(male.change.rank + male.change.var)

# output numbers
# ============================================================

female.int.fdr <- p.adjust(as.numeric(female.gxe$wolba.p.gxe), method = "BH")
female.var.het.fdr <- p.adjust(as.numeric(female.var.het$wolba.p.single.g), method = "BH")

cat("there are ", sum(female.int.fdr < 0.05 & female.var.het.fdr < 0.05, na.rm = TRUE), " genes with int FDR < 0.05 & var.het FDR < 0.05 in females.\n")

male.int.fdr <- p.adjust(as.numeric(male.gxe$wolba.p.gxe), method = "BH")
male.var.het.fdr <- p.adjust(as.numeric(male.var.het$wolba.p.single.g), method = "BH")

cat("there are ", sum(male.int.fdr < 0.05 & male.var.het.fdr < 0.05, na.rm = TRUE), " genes with int FDR < 0.05 & var.het FDR < 0.05 in males.\n")

final.female.table <- cbind(female.gene.name, gene.info[female.gene.name, ], female.table, female.change.rank, female.change.var, female.rank.percent)
final.male.table <- cbind(male.gene.name, gene.info[male.gene.name, ], male.table, male.change.rank, male.change.var, male.rank.percent)



colnames(final.female.table) <- c("FBgn", "symbol", "type", "pos", "female.age.eff", "female.age.eff.sd", "female.age.eff.p", "female.age.eff.fdr",
                           "female.var.line.25c", "female.var.e.25c", "female.25c.h2", "female.p.line.25c", "female.fdr.line.25c", 
                           "female.var.line.18c", "female.var.e.18c", "female.18c.h2", "female.p.line.18c", "female.fdr.line.18c",
                           "female.p.single.g", "female.fdr.single.g", "female.p.single.e", "female.fdr.single.e",
                           "female.var.line", "female.var.int", "female.gxe", "female.p.g", "female.fdr.g", "female.p.gxe", "female.fdr.gxe", "female.change.rank", "female.change.var", "female.rank.percent")
colnames(final.male.table) <- c("FBgn", "symbol", "type", "pos", "male.age.eff", "male.age.eff.sd", "male.age.eff.p", "male.age.eff.fdr",
                           "male.var.line.25c", "male.var.e.25c", "male.25c.h2", "male.p.line.25c", "male.fdr.line.25c", 
                           "male.var.line.18c", "male.var.e.18c", "male.18c.h2", "male.p.line.18c", "male.fdr.line.18c",
                           "male.p.single.g", "male.fdr.single.g", "male.p.single.e", "male.fdr.single.e",
                           "male.var.line", "male.var.int", "male.gxe", "male.p.g", "male.fdr.g", "male.p.gxe", "male.fdr.gxe", "male.change.rank", "male.change.var", "male.rank.percent")
# write excel file
# ============================================================

write.xlsx(list('Female' = final.female.table, 'Male' = final.male.table), file = args[10], colNames = TRUE, rowNames = FALSE)

