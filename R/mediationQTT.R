# =====================================================
# = mediation analysis for trios of eqtl, gene, trait =
# =====================================================

args = commandArgs(TRUE) # args <- c("../eqtl/eqtl.fdr05.snp.geno.tped", "../eqtl/eqtl.fdr05.snp.geno.tfam", "../qtt/male.phototaxis.decline.pheno,../qtt/male.phototaxis.decline.pheno,../qtt/male.speed.decline.pheno,../qtt/male.endurance.decline.pheno,../qtt/male.lifespan.pheno", "male.eqtl.gwas.trio", "../eqtl/male.25c.pheno", "../eqtl/male.18c.pheno")
# takes seven arguments
# 1. genotype file
# 2. tfam fle
# 3. phenotype files
# 4. trio file
library("mediation")

tped.file <- args[1]
tfam.file <- args[2]
pheno.file <- args[3]
trio.file <- args[4]
young.exp.file <- args[5]
old.exp.file <- args[6]

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
pheno.list <- unlist(strsplit(pheno.file, split = ","))
pheno.names <- gsub("\\.pheno", "", gsub(".*/", "", pheno.list))

pheno <- read.table(pheno.list[1], header = TRUE, as.is = TRUE, row.names = 1)

for (i in 2:length(pheno.list)) {
  pheno <- merge(pheno, read.table(pheno.list[i], header = TRUE, as.is = TRUE, row.names = 1), by.x = 1, by.y = 1, all = TRUE)
}
rownames(pheno) <- pheno[, 1]
common.lines <- intersect(tfam[, 1], rownames(pheno))

pheno <- pheno[common.lines, ]

colnames(geno.code) <- tfam[, 1]
geno.code <- geno.code[, common.lines]

# read trios
trio <- read.table(trio.file, as.is = TRUE, header = FALSE)

# process trio
snp.gene.pair <- paste(paste(trio[, 1], trio[, 2], sep = "-"), trio[, 3], sep = "-")
snp.gene.pair.count <- table(snp.gene.pair)
multi.pair <- which(snp.gene.pair.count > 1)
if (length(multi.pair) > 0) {
	trio[match(names(multi.pair), snp.gene.pair), 5] <- "both"
	trio <- trio[-which(snp.gene.pair %in% names(multi.pair) & trio[, 5] != "both"), ]
}

# read expression
young.exp <- read.table(young.exp.file, header = TRUE, as.is = TRUE, row.names = 1)
young.exp <- young.exp[common.lines, ]

old.exp <- read.table(old.exp.file, header = TRUE, as.is = TRUE, row.names = 1)
old.exp <- old.exp[common.lines, ]

# for each trait and gene

count = 0
trio.res <- matrix(NA, nrow = nrow(trio), ncol = 68)

for (trait in unique(trio[, 3])) {
	
	this.trait.trio <- trio[trio[, 3] == trait, ]
	this.pheno <- pheno[, match(trait, pheno.names) + 1]
	cat(trait, "\n")
	
	for (gene in unique(this.trait.trio[, 2])) {
		cat(gene, "\n")
		
		this.trait.gene.trio <- this.trait.trio[this.trait.trio[, 2] == gene, ]
		young.exp.present <- 0
		old.exp.present <- 0
		this.young.exp <- rep(NA, length(this.pheno))
		this.old.exp <- rep(NA, length(this.pheno))
		
		if (sum(colnames(young.exp) == gene) > 0) {
			this.young.exp <- young.exp[, gene]
			young.exp.present <- 1
		}
		
		if (sum(colnames(old.exp) == gene) > 0) {
			this.old.exp <- old.exp[, gene]
			old.exp.present <- 1
		}
		
		if (nrow(this.trait.gene.trio) > 1) {
			
			this.gene.snp.geno <- t(geno.code[paste("VAR", this.trait.gene.trio[, 1], sep = "_"), ])
			
		}
		
		for (i in 1:nrow(this.trait.gene.trio)) {

			count = count + 1
			this.snp <- this.trait.gene.trio[i, 1]
			trio.res[count, 1:3] <- c(this.snp, gene, trait)
			this.snp.geno <- geno.code[paste("VAR_", this.snp, sep = ""), ]
			
			# r2 with other snps
			if (nrow(this.trait.gene.trio) > 1) {
			
				this.snp.r2 <- paste(cor(this.gene.snp.geno, use = "pairwise")[i, -i], collapse = ",")
				this.snp.r2.snps <- paste(this.trait.gene.trio[-i, 1], collapse = ",")
				trio.res[count, 4:5] <- c(this.snp.r2.snps, this.snp.r2)
			
			}
			
			trio.res[count, 6] <- this.trait.gene.trio[i, 5]
			
			this.data <- na.omit(data.frame(young = this.young.exp, old = this.old.exp, diff = this.old.exp - this.young.exp, snp = this.snp.geno, pheno = this.pheno))
			
			young.cor.test <- cor.test(this.young.exp, this.pheno, use = "pairwise")
			old.cor.test <- cor.test(this.old.exp, this.pheno, use = "pairwise")
			diff.cor.test <- cor.test(this.old.exp - this.young.exp, this.pheno, use = "pairwise")
			
			pheno.qtl <- lm(this.pheno ~ this.snp.geno)
			young.eqtl <- lm(this.young.exp ~ this.snp.geno)
			old.eqtl <- lm(this.old.exp ~ this.snp.geno)
			diff.eqtl <- lm(this.old.exp - this.young.exp ~ this.snp.geno)
			
			trio.res[count, 7:14] <- c(summary(pheno.qtl)$coefficients[2, c(1, 4)], summary(young.eqtl)$coefficients[2, c(1, 4)], summary(old.eqtl)$coefficients[2, c(1, 4)], summary(diff.eqtl)$coefficients[2, c(1, 4)])
			
			trio.res[count, 15:20]<- c(young.cor.test$estimate, young.cor.test$p.value, old.cor.test$estimate, old.cor.test$p.value, diff.cor.test$estimate, diff.cor.test$p.value)
			
			
			if (young.exp.present) {
			  m <- lm(young ~ snp, data = this.data)
			  y <- lm(pheno ~ young + snp, data = this.data)
			  young.res <- mediate(m, y, sims = 10000, treat = "snp", mediator = "young")
				trio.res[count, 21:36] <- c(young.res$d0, young.res$d0.ci, young.res$d0.p, young.res$z0, young.res$z0.ci, young.res$z0.p, young.res$n0, young.res$n0.ci, young.res$n0.p, young.res$tau.coef, young.res$tau.ci, young.res$tau.p)
				
			}
			
			if (old.exp.present) {
			  m <- lm(young ~ snp, data = this.data)
			  y <- lm(pheno ~ old + snp, data = this.data)
			  old.res <- mediate(m, y, sims = 10000, treat = "snp", mediator = "old")
				trio.res[count, 37:52] <- c(old.res$d0, old.res$d0.ci, old.res$d0.p, old.res$z0, old.res$z0.ci, old.res$z0.p, old.res$n0, old.res$n0.ci, old.res$n0.p, old.res$tau.coef, old.res$tau.ci, old.res$tau.p)
				
			}
			
			if (young.exp.present & old.exp.present) {
			  m <- lm(young ~ snp, data = this.data)
			  y <- lm(pheno ~ diff + snp, data = this.data)
			  diff.res <- mediate(m, y, sims = 10000, treat = "snp", mediator = "diff")
				trio.res[count, 53:68] <- c(diff.res$d0, diff.res$d0.ci, diff.res$d0.p, diff.res$z0, diff.res$z0.ci, diff.res$z0.p, diff.res$n0, diff.res$n0.ci, diff.res$n0.p, diff.res$tau.coef, diff.res$tau.ci, diff.res$tau.p)
				
			}
		
			
		}
		
	}
	
}



save(trio.res, file = args[7])
