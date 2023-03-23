# ============================================
# = quantitative genetics of gene expression =
# ============================================

# use normal quantiles to transform gene expression
# within each temperature
# ============================================================

args <- commandArgs(TRUE) # args <- c("normalize/female.gene.exp.adjust.RData", 8)
rdata.file <- args[1]
n.cpu <- as.numeric(args[2])

library("lmerTest")
library("doMC")
registerDoMC(n.cpu)

# load data
# ============================================================

load(rdata.file)

# test for aging effect
# 1. just age in the model
# 2. age and wolbachia but get age effect
# 3. stratify by wolbachia status, test for age
# ============================================================

gene.name <- rownames(gene.exp.adj)
temp.wolba <- temp[wolba == "y"]
temp.nowolba <- temp[wolba == "n"]
line.id.wolba <- line.id[wolba == "y"]
line.id.nowolba <- line.id[wolba == "n"]

age.effect <- foreach (i = 1:nrow(gene.exp.adj), .combine = rbind) %dopar% {
  
  cat(gene.name[i], "\n")
  
  y <- unlist(gene.exp.adj[i, ])
  
  # result vector
  res <- rep(NA, 12)
  
  # no wolbachia
  # ============================================================
  
  lmer.pool <- lmer(y ~ temp + (1|line.id) + (1|line.id:temp))
  res[1:3] <- summary(lmer.pool)$coefficients[2, c(1, 2, 5)] # fixed effect est, ste, p

  # in the presence of wolbachia
  # ============================================================
  
  lmer.pool <- lmer(y ~ wolba + temp + (1|line.id) + (1|line.id:temp))
  res[4:6] <- summary(lmer.pool)$coefficients[3, c(1, 2, 5)]
  
  # wolbachia "y"
  # ============================================================
  
  y.wolba <- y[wolba == "y"]
  lmer.pool <- lmer(y.wolba ~ temp.wolba + (1|line.id.wolba) + (1|line.id.wolba:temp.wolba))
  res[7:9] <- summary(lmer.pool)$coefficients[2, c(1, 2, 5)] # fixed effect est, ste, p
  
  # wolbachia "n"
  # ============================================================
  
  y.nowolba <- y[wolba == "n"]
  lmer.pool <- lmer(y.nowolba ~ temp.nowolba + (1|line.id.nowolba) + (1|line.id.nowolba:temp.nowolba))
  res[10:12] <- summary(lmer.pool)$coefficients[2, c(1, 2, 5)] # fixed effect est, ste, p
  
  res
  
}

colnames(age.effect) <- c("fixed.est", "fixed.ste", "fixed.p", "add.wolba.fixed.est", "add.wolba.fixed.ste", "add.wolba.fixed.p", "wolba.fixed.est", "wolba.fixed.ste", "wolba.fixed.p", "nowolba.fixed.est", "nowolba.fixed.ste", "nowolba.fixed.p")
save(age.effect, gene.name, file = args[3])

# session info
# ============================================================

sessionInfo()
