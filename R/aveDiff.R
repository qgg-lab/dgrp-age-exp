# =================================================================
# = calculate average difference between randomly sampled samples =
# =================================================================

args <- commandArgs(TRUE) # args <- c("normalize/female.gene.exp.adjust.RData", "normalize/male.gene.exp.adjust.RData")

load(args[1])
female.young <- gene.exp.adj[, temp == "25C"]
female.old <- gene.exp.adj[, temp == "18C"]
load(args[2])
male.young <- gene.exp.adj[, temp == "25C"]
male.old <- gene.exp.adj[, temp == "18C"]

# calculate mean
# ============================================================

female.young.diff <- numeric(1000)
set.seed(1)

for (i in 1:1000) {
  
  female.young.diff[i] <- sum(abs(apply(female.young[, sample(1:ncol(female.young), 2)], 1, function(x) { x[1] - x[2]})) >= 1)
  
}

# session info
# ============================================================

sessionInfo()
