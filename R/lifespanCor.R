# =====================================================
# = correlation between life span and expression diff =
# =====================================================

args <- commandArgs(TRUE) # args <- c("qg/female.adj.qgBLUP.RData", 8)
rdata.file <- args[1]
n.cpu <- as.numeric(args[2])

library("doMC")
registerDoMC(n.cpu)

# load data
# ============================================================

load(rdata.file)
female.25c.mean.pheno <- read.table("female.25c.mean.pheno", as.is = TRUE)
female.25c.mean.pheno$V1 <- paste("line_", female.25c.mean.pheno$V1, sep = "")

# correlation
# ============================================================

female.25c.life <- female.25c.mean.pheno[which(!is.na(match(female.25c.mean.pheno$V1, line.order))), 2]

blup.diff <- blup.18c - blup.25c

blup.diff.match <- blup.diff[,na.omit(match(female.25c.mean.pheno$V1, line.order))]

correlation <- foreach (i = 1:nrow(blup.diff.match), .combine = rbind) %dopar% {
  res <- cor(blup.diff.match[i, ], female.25c.life)
  res
}

rownames(correlation) <- gene.name

#head(correlation[order(correlation[,1]),])
#tail(correlation[order(correlation[,1]),])

save(correlation, file = args[3])

# session info
# ============================================================

sessionInfo()