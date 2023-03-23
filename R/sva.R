# ==================
# = sva correction =
# ==================

args <- commandArgs(TRUE) # args <- c("/mnt/research/qgg/dgrp-age-exp/normalize/male.gene.exp.filtered.RData", "/mnt/research/qgg/resource/dgrp/dgrp2.r6.adjustData.RData")

library("sva")

# nqt function
# ============================================================

nqt <- function(x) {
  
  return(qnorm(rank(x)/(length(x) + 1), sd = mad(x)) + median(x))
  
}


# load data
# ============================================================

load(args[1])
load(args[2])

# convert age to temp
# ============================================================

temp <- as.character(age)
temp[temp == "baseline"] <- "25C"
temp[temp == "3wk"] <- "18C"

# find line id and temp
# ============================================================

temp <- factor(temp, levels = c("25C", "18C"))
line.id <- factor(paste("line_", unlist(strsplit(sample.name, split = "_"))[seq(from = 2, length = length(sample.name), by = 3)], sep = ""))
flow.cell <- flow.cell[, 2]
wolba <- factor(wolba[as.character(line.id), 1])

# there are a few missing data, fill in those with mean expression across all samples
# ============================================================

gene.exp <- as.matrix(gene.exp)
missing.gene <- which(apply(gene.exp, 1, function(x) { sum(is.na(x)) }) > 0)
if (length(missing.gene) > 0) {
  for (i in missing.gene) {
    gene.exp[i, ] <- ifelse(is.na(gene.exp[i, ]), mean(gene.exp[i, ], na.rm = TRUE), gene.exp[i, ])
  }
}

# estimate sv number
# ============================================================

gene.exp <- as.matrix(log2(gene.exp + 1/32))
gene.exp.nqt <- gene.exp
for (i in 1:nrow(gene.exp)) { gene.exp.nqt[i, ] <- nqt(gene.exp[i, ]) }

mod = model.matrix( ~ temp + wolba)
mod0 = matrix(1, nrow = ncol(gene.exp.nqt), ncol = 1)
colnames(mod0) <- "Intercept"
nsv <- num.sv(gene.exp.nqt, mod = mod, method = "be", B = 20, seed = 1)

sva.adj <- sva(gene.exp.nqt, mod, mod0, n.sv = nsv)

gene.exp.adj <- gene.exp.nqt
for (i in 1:nrow(gene.exp.nqt)) {
  this.fit <- lm(gene.exp.nqt[i, ] ~ sva.adj$sv)
  gene.exp.adj[i, ] <- coefficients(this.fit)[1] + residuals(this.fit)
}

save(sva.adj, flow.cell, file = args[3])
save(gene.exp.adj, wolba, temp, line.id, file = args[4])

# session info
# ============================================================

sessionInfo()
