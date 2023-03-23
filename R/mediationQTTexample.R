# ========================
# = extract example data =
# ========================

args = commandArgs(TRUE) # args <- c("../eqtl/eqtl.fdr05.snp.geno.tped", "../eqtl/eqtl.fdr05.snp.geno.tfam", "../qtt/male.phototaxis.decline.pheno,../qtt/male.phototaxis.decline.pheno,../qtt/male.speed.decline.pheno,../qtt/male.endurance.decline.pheno,../qtt/male.lifespan.pheno", "male.eqtl.gwas.trio", "../eqtl/male.25c.pheno", "../eqtl/male.18c.pheno", "male.speed.decline", "3R_12642725", "FBgn0026207")
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
select.traits <- unlist(strsplit(args[7], split = ","))
select.snps <- unlist(strsplit(args[8], split = ","))
select.genes <- unlist(strsplit(args[9], split = ","))

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

# expression
# read expression
young.exp <- read.table(young.exp.file, header = TRUE, as.is = TRUE, row.names = 1)
young.exp <- young.exp[common.lines, ]

old.exp <- read.table(old.exp.file, header = TRUE, as.is = TRUE, row.names = 1)
old.exp <- old.exp[common.lines, ]

# get phenotype
select.pheno <- as.matrix(pheno[, match(select.traits, pheno.names) + 1])
colnames(select.pheno) <- select.traits

select.young.exp <- as.matrix(young.exp[, select.genes])
colnames(select.young.exp) <- select.genes

select.old.exp <- as.matrix(old.exp[, select.genes])
colnames(select.old.exp) <- select.genes

select.geno <- matrix(geno.code[paste("VAR_", select.snps, sep = ""), ], nrow = length(select.snps))
rownames(select.geno) <- select.snps

save(select.pheno, select.young.exp, select.old.exp, select.geno, common.lines, file = args[10])
