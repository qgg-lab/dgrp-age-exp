# ==========================
# = look for special trios =
# ==========================

args <- commandArgs(TRUE) # args <- c("female.young.trio.mediation.RData", "../eqtl/eqtl.fdr05.snp.geno.tped", "../eqtl/eqtl.fdr05.snp.geno.tfam", "../eqtl/female.25c.pheno", "../eqtl/female.18c.pheno")

# read data
# ============================================================

load(args[1])

tped.file <- args[2]
tfam.file <- args[3]

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

# phenotypes
pheno1 <- read.table(args[4], header = TRUE, as.is = TRUE, row.names = 1)
pheno1 <- pheno1[tfam[, 1], ]

pheno2 <- read.table(args[5], header = TRUE, as.is = TRUE, row.names = 1)
pheno2 <- pheno2[tfam[, 1], ]

# add four more columns
# p value for cis pheno1, p value for cis pheno2
# p value for trans pheno1, p value for trans pheno2
trio.info <- cbind(trio, trio.res, NA, NA, NA, NA, NA, NA)

# loop through trios
# ============================================================

for (i in 1:nrow(trio)) {
	
	cis.gene <- trio[i, 2]
	trans.gene <- trio[i, 5]
	snp <- trio[i, 1]
	
	cis.gene.exp1 <- pheno1[, cis.gene]
	trans.gene.exp1 <- pheno1[, trans.gene]
	snp.geno <- geno.code[paste("VAR_", snp, sep = ""), ]
	
	if (length(intersect(colnames(pheno2), c(cis.gene, trans.gene))) == 2) {
		
		cis.gene.exp2 <- pheno2[, cis.gene]
		trans.gene.exp2 <- pheno2[, trans.gene]
	
		# pheno1 cis
		trio.info[i, 24] <- summary(lm(cis.gene.exp1 ~ snp.geno))$coefficients[2, 4]
		
		# pheno2 cis
		trio.info[i, 25] <- summary(lm(cis.gene.exp2 ~ snp.geno))$coefficients[2, 4]
		
		# pheno1 trans
		trio.info[i, 26] <- summary(lm(trans.gene.exp1 ~ snp.geno))$coefficients[2, 4]
		
		# pheno2 trans
		trio.info[i, 27] <- summary(lm(trans.gene.exp2 ~ snp.geno))$coefficients[2, 4]
		
		# cor pheno1
		trio.info[i, 28] <- cor(cis.gene.exp1, trans.gene.exp1, use = "pairwise")
		
		# cor pheno2
		trio.info[i, 29] <- cor(cis.gene.exp2, trans.gene.exp2, use = "pairwise")
		
	}

}
	
trio.info[trio.info[, 24] < 1e-5 & trio.info[, 25] < 1e-5 & trio.info[, 26] < 1e-5 & trio.info[, 27] > 1e-5 & trio.info[, 11] < 0.001 & trio.info[, 28] > 0.4 & trio.info[, 29] < 0.3, ]

i = 91

cis.gene <- trio[i, 2]
trans.gene <- trio[i, 5]
snp <- trio[i, 1]

cis.gene.exp1 <- pheno1[, cis.gene]
trans.gene.exp1 <- pheno1[, trans.gene]
snp.geno <- geno.code[paste("VAR_", snp, sep = ""), ]

cis.gene.exp2 <- pheno2[, cis.gene]
trans.gene.exp2 <- pheno2[, trans.gene]

save(cis.gene, trans.gene, snp, cis.gene.exp1, trans.gene.exp1, snp.geno, cis.gene.exp2, trans.gene.exp2, file = args[6])

# 91 3R_4680760 FBgn0263346 select cis FBgn0026620 map trans



