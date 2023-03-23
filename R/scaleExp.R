# =========================================
# = scale expression to identify outliers =
# =========================================

args <- commandArgs(TRUE) # args <- ("/mnt/research/qgg/dgrp-age-exp/normalize/female.gene.exp.RData")

load(args[1])

# scale expression
# ============================================================

# split by age
gene.exp.baseline <- gene.exp[, age == "baseline"]
gene.exp.old <- gene.exp[, age == "3wk"]

baseline.scaled <- t(scale(t(log2(gene.exp.baseline + 1/32))))
old.scaled <- t(scale(t(log2(gene.exp.old + 1/32))))

# scale within age
baseline.quantile <- apply(baseline.scaled, 2, function(x) { quantile(x, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), na.rm = T) })
old.quantile <- apply(old.scaled, 2, function(x) { quantile(x, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), na.rm = T) })


save(age, sample.name, flow.cell, baseline.quantile, old.quantile, file = args[2])

sessionInfo()
