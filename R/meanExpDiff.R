# get mean difference between ages for females and males
# ============================================================

args <- commandArgs(TRUE) # args <- c("normalize/female.gene.exp.adjust.RData", "normalize/male.gene.exp.adjust.RData")

load(rdata.file)

# transformation function
# ============================================================

nqt <- function(x) {
  
  return(qnorm(rank(x)/(length(x) + 1), sd = mad(x)) + median(x))
  
}

# find line id and temperature
# ============================================================

unique.id <- factor(paste(line.id, temp, sep = "-"))

# perform QG analysis
# need
# 1. fit a linear mixed model: temp + (1|line) + (1|line:temp)
# 2. test for significance of line by temp interaction
# 3. get blup from the full model
# 4. get variance components from single temperature, not blup
# ============================================================

gene.name <- rownames(gene.exp.adj)

gxe <- foreach (i = 1:nrow(gene.exp.adj), .combine = rbind) %dopar% {
  
  cat(gene.name[i], "\n")
  
  y <- unlist(gene.exp.adj[i, ])
  
  # result vector
  res <- rep(NA, 13)
  
  # no wolbachia
  # ============================================================
  
  lmer.pool <- lmer(y ~ temp + (1|line.id) + (1|line.id:temp))
  
  res[1:3] <- summary(lmer.pool)$coefficients[2, c(1, 2, 5)] # fixed effect est, ste, p
  res[4:6] <- as.data.frame(VarCorr(lmer.pool))[, 4] # var line x temp, line, res
  res[7:8] <- rand(lmer.pool)$Pr[2:3] # p value for line, line x temp
  
  # in the presence of wolbachia
  # ============================================================
  
  lmer.pool.wolba <- lmer(y ~ temp*wolba + (1|line.id) + (1|line.id:temp))
  res[9:11] <- as.data.frame(VarCorr(lmer.pool.wolba))[, 4]
  res[12:13] <- rand(lmer.pool.wolba)$Pr[2:3]
  
  # return result
  # ============================================================
  
  res
  
}

colnames(gxe) <- c("fixed.est", "fixed.ste", "fixed.p", "var.gxe", "var.g", "var.e", "p.g", "p.gxe", "wolba.var.gxe", "wolba.var.g", "wolba.var.e", "wolba.p.g", "wolba.p.gxe")
save(gxe, gene.name, file = args[3])

# session info
# ============================================================

sessionInfo()
