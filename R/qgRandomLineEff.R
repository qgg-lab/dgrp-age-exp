# ============================================
# = quantitative genetics of gene expression =
# ============================================

# randomly assign age groups
# ============================================================

args <- commandArgs(TRUE) # args <- c("normalize/female.gene.exp.adjust.RData", 8)
rdata.file <- args[1]
library("lmerTest")

# load data
# ============================================================

load(rdata.file)

# test for aging effect
# 1. just age in the model
# 2. age and wolbachia but get age effect
# 3. stratify by wolbachia status, test for age
# ============================================================

gene.name <- rownames(gene.exp.adj)

old.gene.exp.adj <- gene.exp.adj[, temp == "18C"]
young.gene.exp.adj <- gene.exp.adj[, temp == "25C"]
old.line.id <- line.id[temp == "18C"]
young.line.id <- line.id[temp == "25C"]
old.wolba <- wolba[temp == "18C"]
young.wolba <- wolba[temp == "25C"]

young.gene.diff.num <- numeric(10)
old.gene.diff.num <- numeric(10)

for (i in 1:10) {
  
  set.seed(i)
  
  # old
  old.temp <- factor(ifelse(old.line.id %in% sample(levels(old.line.id), 100), "young", "old"), levels = c("young", "old"))
  
  res <- matrix(nrow = nrow(old.gene.exp.adj), ncol = 3)
  
  for (j in 1:nrow(old.gene.exp.adj)) {
    y <- unlist(old.gene.exp.adj[j, ])
    lmer.pool <- lmer(y ~ old.temp + (1|old.line.id))
    res[j, ] <- summary(lmer.pool)$coefficients[2, c(1, 2, 5)]
    cat(j, "\n")
  }
  
  old.gene.diff.num[i] <- sum(abs(res[, 1]) >= 1 & p.adjust(res[, 3], method = "BH") <= 0.05)
  
  # young
  young.temp <- factor(ifelse(young.line.id %in% sample(levels(young.line.id), 100), "young", "old"), levels = c("young", "old"))
  
  res <- matrix(nrow = nrow(young.gene.exp.adj), ncol = 3)
  
  for (j in 1:nrow(young.gene.exp.adj)) {
    y <- unlist(young.gene.exp.adj[j, ])
    lmer.pool <- lmer(y ~ young.temp + (1|young.line.id))
    res[j, ] <- summary(lmer.pool)$coefficients[2, c(1, 2, 5)]
    cat(j, "\n")
  }
  
  young.gene.diff.num[i] <- sum(abs(res[, 1]) >= 1 & p.adjust(res[, 3], method = "BH") <= 0.05)
  
}

cat("among old, randomly splitting:", old.gene.diff.num, "\n")
cat("among young, randomly splitting:", young.gene.diff.num, "\n")

# session info
# ============================================================

sessionInfo()
