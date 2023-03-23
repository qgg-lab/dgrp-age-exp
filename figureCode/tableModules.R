# ===============================
# = table for module identities =
# ===============================

args <- commandArgs(TRUE) # args <- c("../figureData/female.wgcna.RData", "../figureData/male.wgcna.RData", "../figureData/gene.info")
library("openxlsx")

# load data
# ============================================================

gene.info <- read.table(args[3], header = FALSE, as.is = TRUE, quote = "\"", row.names = 1)

load(args[1])
female.young.module <- sort(factor(blup.25c.tree, levels = c(1:max(unique(blup.25c.tree)),0)))
female.young.module <- cbind("Female", "Young", names(female.young.module), gene.info[names(female.young.module), ], as.character(female.young.module))

names(blup.18c.tree) <- names(blup.25c.tree)
female.old.module <- sort(factor(blup.18c.tree, levels = c(1:max(unique(blup.18c.tree)),0)))
female.old.module <- cbind("Female", "Old", names(female.old.module), gene.info[names(female.old.module), ], as.character(female.old.module))

load(args[2])
male.young.module <- sort(factor(blup.25c.tree, levels = c(1:max(unique(blup.25c.tree)),0)))
male.young.module <- cbind("Male", "Young", names(male.young.module), gene.info[names(male.young.module), ], as.character(male.young.module))

names(blup.18c.tree) <- names(blup.25c.tree)
male.old.module <- sort(factor(blup.18c.tree, levels = c(1:max(unique(blup.18c.tree)),0)))
male.old.module <- cbind("Male", "Old", names(male.old.module), gene.info[names(male.old.module), ], as.character(male.old.module))

final.table = mapply(c, female.young.module, female.old.module, male.young.module, male.old.module)
colnames(final.table) <- c("sex", "age", "FBgn", "symbol", "type", "pos", "module")

write.xlsx(list('module' = final.table), file = args[4], colNames = TRUE, rowNames = FALSE)

