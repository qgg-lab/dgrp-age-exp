# ================================
# = table for mediation analysis =
# ================================

args = commandArgs(TRUE) # args <- c("../figureData/female.young.trio.mediation.RData", "../figureData/female.old.trio.mediation.RData", "../figureData/male.young.trio.mediation.RData", "../figureData/male.old.trio.mediation.RData", "../figureData/gene.info")
library("openxlsx")

# female young
load(args[1])
female.young.trio <- trio; female.young.trio.res <- trio.res
rownames(female.young.trio) <- apply(female.young.trio, 1, paste, collapse = "-")
rownames(female.young.trio.res) <- apply(female.young.trio, 1, paste, collapse = "-")

# female old
load(args[2])
female.old.trio <- trio; female.old.trio.res <- trio.res
rownames(female.old.trio) <- apply(female.old.trio, 1, paste, collapse = "-")
rownames(female.old.trio.res) <- apply(female.old.trio, 1, paste, collapse = "-")

# male young
load(args[3])
male.young.trio <- trio; male.young.trio.res <- trio.res
rownames(male.young.trio) <- apply(male.young.trio, 1, paste, collapse = "-")
rownames(male.young.trio.res) <- apply(male.young.trio, 1, paste, collapse = "-")

# male old
load(args[4])
male.old.trio <- trio; male.old.trio.res <- trio.res
rownames(male.old.trio) <- apply(male.old.trio, 1, paste, collapse = "-")
rownames(male.old.trio.res) <- apply(male.old.trio, 1, paste, collapse = "-")

# gene info
gene.info <- read.table(args[5], header = FALSE, as.is = TRUE)
rownames(gene.info) <- gene.info[, 1]

# female young
cat("Young females: ", sum(p.adjust(female.young.trio.res[, 4], method = "BH") <= 0.05), "\n")
cat("Aged females: ", sum(p.adjust(female.old.trio.res[, 4], method = "BH") <= 0.05), "\n")
cat("Young males: ", sum(p.adjust(male.young.trio.res[, 4], method = "BH") <= 0.05), "\n")
cat("Aged males: ", sum(p.adjust(male.old.trio.res[, 4], method = "BH") <= 0.05), "\n")



# female
# ============================================================

final.table <- rbind(as.matrix(cbind("Young females", female.young.trio[, c(1, 2, 5)], 
																					  female.young.trio.res[, c(1, 4)], p.adjust(female.young.trio.res[, 4], method = "BH"),
																						female.young.trio.res[, c(5, 8)], p.adjust(female.young.trio.res[, 8], method = "BH"),
																						female.young.trio.res[, 9],
																						female.young.trio.res[, c(13, 16)], p.adjust(female.young.trio.res[, 16], method = "BH"))),
										 as.matrix(cbind("Aged females", female.old.trio[, c(1, 2, 5)], 
																					female.old.trio.res[, c(1, 4)], p.adjust(female.old.trio.res[, 4], method = "BH"),
																					female.old.trio.res[, c(5, 8)], p.adjust(female.old.trio.res[, 8], method = "BH"),
																					female.old.trio.res[, 9],
																					female.old.trio.res[, c(13, 16)], p.adjust(female.old.trio.res[, 16], method = "BH"))),
										 as.matrix(cbind("Young males", male.young.trio[, c(1, 2, 5)], 
																						male.young.trio.res[, c(1, 4)], p.adjust(male.young.trio.res[, 4], method = "BH"),
																						male.young.trio.res[, c(5, 8)], p.adjust(male.young.trio.res[, 8], method = "BH"),
																						male.young.trio.res[, 9],
																						male.young.trio.res[, c(13, 16)], p.adjust(male.young.trio.res[, 16], method = "BH"))),
										 as.matrix(cbind("Aged males", male.old.trio[, c(1, 2, 5)], 
																					male.old.trio.res[, c(1, 4)], p.adjust(male.old.trio.res[, 4], method = "BH"),
																					male.old.trio.res[, c(5, 8)], p.adjust(male.old.trio.res[, 8], method = "BH"),
																					male.old.trio.res[, 9],
																					male.old.trio.res[, c(13, 16)], p.adjust(male.old.trio.res[, 16], method = "BH"))) )

final.table <- cbind(final.table[, 1:2], gene.info[final.table[, 3], ], gene.info[final.table[, 4], ], final.table[, 5:ncol(final.table)])

# save table
# ============================================================

write.xlsx(list(table = final.table), file = args[6], colNames = FALSE, rowNames = FALSE)


