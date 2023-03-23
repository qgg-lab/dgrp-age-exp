# ================================
# = mediation analytis for trios =
# ================================

args = commandArgs(TRUE) # args <- c("../eqtl/eqtl.fdr05.snp.geno.tped", "../eqtl/eqtl.fdr05.snp.geno.tfam", "../eqtl/female.25c.pheno", "female.young.eqtl.trio")
# takes seven arguments
# 1. genotype file
# 2. tfam fle
# 3. phenotype file
# 4. trio file
library("mediation")

tped.file <- args[1]
tfam.file <- args[2]
pheno.file <- args[3]
trio.file <- args[4]

# read genotype
tfam <- read.table(tfam.file, header = FALSE, as.is = TRUE)
n.line <- nrow(tfam)
tped <- matrix(scan(tped.file, what = ""), ncol = 4 + 2*n.line, byrow = TRUE)

snp.name <- paste("VAR_", tped[, 2], sep = "")
snp.loc <- tped[, c(1, 4)]
#chr.name <- c("2L", "2R", "3L", "3R", "X", "4")
#snp.loc[, 1] <- chr.name[as.numeric(snp.loc[, 1])]
snp.loc <- data.frame(snp.loc)
snp.loc[,1] <- gsub("23","X",snp.loc[,1])
rownames(snp.loc) <- snp.name

geno.code <- matrix(4 - (as.numeric(tped[, seq(from = 5, length = n.line, by = 2)]) + as.numeric(tped[, seq(from = 6, length = n.line, by = 2)])), ncol = n.line)
geno.code[geno.code == 4] <- NA
rownames(geno.code) <- snp.name

# read phenotype
pheno <- read.table(pheno.file, header = TRUE, as.is = TRUE, row.names = 1)
pheno <- pheno[tfam[, 1], ]

# read trios
trio <- read.table(args[4], as.is = TRUE, header = FALSE)
trio.res <- matrix(NA, nrow = nrow(trio), ncol = 16)

for (i in 1:nrow(trio)) {
	
	cis.gene <- trio[i, 2]
	trans.gene <- trio[i, 5]
	snp <- trio[i, 1]
	
	cis.gene.exp <- pheno[, cis.gene]
	trans.gene.exp <- pheno[, trans.gene]
	snp.geno <- geno.code[paste("VAR_", snp, sep = ""), ]
	
	this.data <- na.omit(data.frame(cis = cis.gene.exp, trans = trans.gene.exp, snp = snp.geno))
	
  m <- lm(cis ~ snp, data = this.data)
  y <- lm(trans ~ cis + snp, data = this.data)
  res <- mediate(m, y, sims = 10000, treat = "snp", mediator = "cis")
	
	trio.res[i, ] <- c(res$d0, res$d0.ci, res$d0.p, res$z0, res$z0.ci, res$z0.p, res$n0, res$n0.ci, res$n0.p, res$tau.coef, res$tau.ci, res$tau.p)
	
}

save(trio, trio.res, file = args[5])
