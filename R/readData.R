# =============
# = read data =
# =============

args <- commandArgs(TRUE) # args <- c("/mnt/research/qgg/dgrp-age-exp/joinExp/female.exp.tpm", "/mnt/research/qgg/dgrp-age-exp/joinExp/male.exp.tpm", "/mnt/research/qgg/dgrp-age-exp/joinExp/all.sample.flow.cell.info")

# flow cell information
# ============================================================

all.flow.cell <- read.table(args[3], header = FALSE, as.is = TRUE)
rownames(all.flow.cell) <- all.flow.cell[, 1]

# female gene expression
# ============================================================

gene.exp <- read.table(args[1], header = T, as.is = TRUE)
rownames(gene.exp) <- gene.exp[, 1]
gene.exp <- gene.exp[, -1]

sample.name <- gsub("^X", "", colnames(gene.exp))
flow.cell <- all.flow.cell[sample.name, ]
age <- factor(gsub("_.*", "", sample.name), levels = c("baseline", "3wk"))

# filter based on gene expression, must have median >= 1 in at least one age
max.median <- apply(gene.exp, 1, function(x){ return(max(sapply(split(unlist(x), age), median, na.rm = T), na.rm = T)) })
gene.mad <- apply(gene.exp, 1, mad, na.rm = T)
gene.exp <- gene.exp[max.median >= 1 & gene.mad > 0.1, ]

save(age, flow.cell, gene.exp, sample.name, file = args[4])

# male gene expression
# ============================================================

gene.exp <- read.table(args[2], header = T, as.is = TRUE)
rownames(gene.exp) <- gene.exp[, 1]
gene.exp <- gene.exp[, -1]

sample.name <- gsub("^X", "", colnames(gene.exp))
flow.cell <- all.flow.cell[sample.name, ]
age <- factor(gsub("_.*", "", sample.name), levels = c("baseline", "3wk"))

# filter based on gene expression, must have median >= 1 in at least one age
max.median <- apply(gene.exp, 1, function(x){ return(max(sapply(split(unlist(x), age), median, na.rm = T), na.rm = T)) })
gene.mad <- apply(gene.exp, 1, mad, na.rm = T)
gene.exp <- gene.exp[max.median >= 1 & gene.mad > 0.1, ]

save(age, flow.cell, gene.exp, sample.name, file = args[5])
