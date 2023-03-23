# ================
# = qtt analysis =
# ================

args <- commandArgs(TRUE) # args <- c("../figureData/female.phototaxis.decline.qtt.RData", "../figureData/female.fecundity.decline.qtt.RData", "../figureData/female.lifespan.qtt.RData", "../figureData/male.phototaxis.decline.qtt.RData", "../figureData/male.speed.decline.qtt.RData", "../figureData/male.endurance.decline.qtt.RData", "../figureData/male.lifespan.qtt.RData", "../figureData/gene.info")

library("openxlsx")

gene.info <- read.table(args[8], header = FALSE, as.is = TRUE, quote = "\"", row.names = 1)

# female traits
# ============================================================

load(args[1])

final.female.table <- cbind(qtt.cor.young, qtt.cor.pval.young, qtt.cor.old, qtt.cor.pval.old, qtt.cor.diff, qtt.cor.pval.diff)

cat("female phototaxis decline:\n")
cat("young", unlist(gene.info[gene.name[which.max(abs(qtt.cor.young))], ]), qtt.cor.young[which.max(abs(qtt.cor.young))], qtt.cor.pval.young[which.max(abs(qtt.cor.young))], "\n")
cat("old", unlist(gene.info[gene.name[which.max(abs(qtt.cor.old))], ]), qtt.cor.old[which.max(abs(qtt.cor.old))], qtt.cor.pval.old[which.max(abs(qtt.cor.old))], "\n")
cat("diff", unlist(gene.info[gene.name[which.max(abs(qtt.cor.diff))], ]), qtt.cor.diff[which.max(abs(qtt.cor.diff))], qtt.cor.pval.diff[which.max(abs(qtt.cor.diff))], "\n")
cat("\n")


load(args[2])

final.female.table <- cbind(final.female.table, qtt.cor.young, qtt.cor.pval.young, qtt.cor.old, qtt.cor.pval.old, qtt.cor.diff, qtt.cor.pval.diff)

cat("female fecundity decline:\n")
cat("young", unlist(gene.info[gene.name[which.max(abs(qtt.cor.young))], ]), qtt.cor.young[which.max(abs(qtt.cor.young))], qtt.cor.pval.young[which.max(abs(qtt.cor.young))], "\n")
cat("old", unlist(gene.info[gene.name[which.max(abs(qtt.cor.old))], ]), qtt.cor.old[which.max(abs(qtt.cor.old))], qtt.cor.pval.old[which.max(abs(qtt.cor.old))], "\n")
cat("diff", unlist(gene.info[gene.name[which.max(abs(qtt.cor.diff))], ]), qtt.cor.diff[which.max(abs(qtt.cor.diff))], qtt.cor.pval.diff[which.max(abs(qtt.cor.diff))], "\n")
cat("\n")


load(args[3])

final.female.table <- cbind(final.female.table, qtt.cor.young, qtt.cor.pval.young, qtt.cor.old, qtt.cor.pval.old, qtt.cor.diff, qtt.cor.pval.diff)


cat("female lifespan:\n")
cat("young", unlist(gene.info[gene.name[which.max(abs(qtt.cor.young))], ]), qtt.cor.young[which.max(abs(qtt.cor.young))], qtt.cor.pval.young[which.max(abs(qtt.cor.young))], "\n")
cat("old", unlist(gene.info[gene.name[which.max(abs(qtt.cor.old))], ]), qtt.cor.old[which.max(abs(qtt.cor.old))], qtt.cor.pval.old[which.max(abs(qtt.cor.old))], "\n")
cat("diff", unlist(gene.info[gene.name[which.max(abs(qtt.cor.diff))], ]), qtt.cor.diff[which.max(abs(qtt.cor.diff))], qtt.cor.pval.diff[which.max(abs(qtt.cor.diff))], "\n")
cat("\n")


final.female.table <- cbind(gene.name, gene.info[gene.name, ], final.female.table)

colnames(final.female.table) <- c("FBgn", "symbol", "type", "pos", paste(rep(c("phototaxis", "fecundity", "lifespan"), each = 6), rep(c("cor(young)", "pval(young)", "cor(aged)", "pval(aged)", "cor(diff)", "pval(diff)"), 3), sep = "-"))

# male traits
# ============================================================

load(args[4])

final.male.table <- cbind(qtt.cor.young, qtt.cor.pval.young, qtt.cor.old, qtt.cor.pval.old, qtt.cor.diff, qtt.cor.pval.diff)

cat("male phototaxis:\n")
cat("young", unlist(gene.info[gene.name[which.max(abs(qtt.cor.young))], ]), qtt.cor.young[which.max(abs(qtt.cor.young))], qtt.cor.pval.young[which.max(abs(qtt.cor.young))], "\n")
cat("old", unlist(gene.info[gene.name[which.max(abs(qtt.cor.old))], ]), qtt.cor.old[which.max(abs(qtt.cor.old))], qtt.cor.pval.old[which.max(abs(qtt.cor.old))], "\n")
cat("diff", unlist(gene.info[gene.name[which.max(abs(qtt.cor.diff))], ]), qtt.cor.diff[which.max(abs(qtt.cor.diff))], qtt.cor.pval.diff[which.max(abs(qtt.cor.diff))], "\n")
cat("\n")


load(args[5])

final.male.table <- cbind(final.male.table, qtt.cor.young, qtt.cor.pval.young, qtt.cor.old, qtt.cor.pval.old, qtt.cor.diff, qtt.cor.pval.diff)

cat("male speed:\n")
cat("young", unlist(gene.info[gene.name[which.max(abs(qtt.cor.young))], ]), qtt.cor.young[which.max(abs(qtt.cor.young))], qtt.cor.pval.young[which.max(abs(qtt.cor.young))], "\n")
cat("old", unlist(gene.info[gene.name[which.max(abs(qtt.cor.old))], ]), qtt.cor.old[which.max(abs(qtt.cor.old))], qtt.cor.pval.old[which.max(abs(qtt.cor.old))], "\n")
cat("diff", unlist(gene.info[gene.name[which.max(abs(qtt.cor.diff))], ]), qtt.cor.diff[which.max(abs(qtt.cor.diff))], qtt.cor.pval.diff[which.max(abs(qtt.cor.diff))], "\n")
cat("\n")

load(args[6])

final.male.table <- cbind(final.male.table, qtt.cor.young, qtt.cor.pval.young, qtt.cor.old, qtt.cor.pval.old, qtt.cor.diff, qtt.cor.pval.diff)

cat("male endurance:\n")
cat("young", unlist(gene.info[gene.name[which.max(abs(qtt.cor.young))], ]), qtt.cor.young[which.max(abs(qtt.cor.young))], qtt.cor.pval.young[which.max(abs(qtt.cor.young))], "\n")
cat("old", unlist(gene.info[gene.name[which.max(abs(qtt.cor.old))], ]), qtt.cor.old[which.max(abs(qtt.cor.old))], qtt.cor.pval.old[which.max(abs(qtt.cor.old))], "\n")
cat("diff", unlist(gene.info[gene.name[which.max(abs(qtt.cor.diff))], ]), qtt.cor.diff[which.max(abs(qtt.cor.diff))], qtt.cor.pval.diff[which.max(abs(qtt.cor.diff))], "\n")
cat("\n")

load(args[7])

final.male.table <- cbind(final.male.table, qtt.cor.young, qtt.cor.pval.young, qtt.cor.old, qtt.cor.pval.old, qtt.cor.diff, qtt.cor.pval.diff)

cat("male lifespan:\n")
cat("young", unlist(gene.info[gene.name[which.max(abs(qtt.cor.young))], ]), qtt.cor.young[which.max(abs(qtt.cor.young))], qtt.cor.pval.young[which.max(abs(qtt.cor.young))], "\n")
cat("old", unlist(gene.info[gene.name[which.max(abs(qtt.cor.old))], ]), qtt.cor.old[which.max(abs(qtt.cor.old))], qtt.cor.pval.old[which.max(abs(qtt.cor.old))], "\n")
cat("diff", unlist(gene.info[gene.name[which.max(abs(qtt.cor.diff))], ]), qtt.cor.diff[which.max(abs(qtt.cor.diff))], qtt.cor.pval.diff[which.max(abs(qtt.cor.diff))], "\n")
cat("\n")

final.male.table <- cbind(gene.name, gene.info[gene.name, ], final.male.table)

colnames(final.male.table) <- c("FBgn", "symbol", "type", "pos", paste(rep(c("phototaxis", "speed", "endurance", "lifespan"), each = 6), rep(c("cor(young)", "pval(young)", "cor(aged)", "pval(aged)", "cor(diff)", "pval(diff)"), 4), sep = "-"))

write.xlsx(list('Female' = final.female.table, 'Male' = final.male.table), file = args[9], colNames = TRUE, rowNames = FALSE)

