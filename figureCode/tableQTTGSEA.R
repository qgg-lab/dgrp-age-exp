# =====================================================
# = make a table for the QTT GSEA analysis =
# =====================================================

args <- commandArgs(TRUE) # args <- c("phototaxis.decline,fecundity.decline,lifespan", "phototaxis.decline,speed.decline,endurance.decline,lifespan", "../figureData/kegg.id")
library("openxlsx")
library("GO.db")

# read and process female data
# ============================================================
kegg.id <- read.table(args[3], header = FALSE, as.is = TRUE, sep = "\t")      
rownames(kegg.id) <- kegg.id[, 1]

female.traits <- unlist(strsplit(args[1], split = ","))
male.traits <- unlist(strsplit(args[2], split = ","))

all.go.gsea <- NULL

for (trait in female.traits) {
	
	for (age in c("young", "old", "diff")) {
	
		load(paste("../figureData/female.", trait, ".qtt.", age, ".gsea.fdr.RData", sep = ""))

		bp.pos.fdr <- bp.pos.fdr[!is.na(as.vector(Term(as.character(bp.pos.fdr[, 1])))), ]
		bp.neg.fdr <- bp.neg.fdr[!is.na(as.vector(Term(as.character(bp.neg.fdr[, 1])))), ]
		mf.pos.fdr <- mf.pos.fdr[!is.na(as.vector(Term(as.character(mf.pos.fdr[, 1])))), ]
		mf.neg.fdr <- mf.neg.fdr[!is.na(as.vector(Term(as.character(mf.neg.fdr[, 1])))), ]
		this.bp <- unname(cbind(trait, "Female", age, "BP", as.character(bp.pos.fdr[, 1]), as.vector(Term(as.character(bp.pos.fdr[, 1]))), bp.pos.fdr[, 2:3], trait, "Female", age, "BP",
                          as.character(bp.neg.fdr[, 1]), as.vector(Term(as.character(bp.neg.fdr[, 1]))), bp.neg.fdr[, 2:3]))
                          
		this.mf <- unname(cbind(trait, "Female", age, "MF", as.character(mf.pos.fdr[, 1]), as.vector(Term(as.character(mf.pos.fdr[, 1]))), mf.pos.fdr[, 2:3], trait, "Female", age, "MF",
                          as.character(mf.neg.fdr[, 1]), as.vector(Term(as.character(mf.neg.fdr[, 1]))), mf.neg.fdr[, 2:3]))
		
		all.go.gsea <- rbind(all.go.gsea, as.matrix(this.bp))
		all.go.gsea <- rbind(all.go.gsea, as.matrix(this.mf))
	
	}
}
                        
for (trait in male.traits) {
	
	for (age in c("young", "old", "diff")) {
	
		load(paste("../figureData/male.", trait, ".qtt.", age, ".gsea.fdr.RData", sep = ""))

		bp.pos.fdr <- bp.pos.fdr[!is.na(as.vector(Term(as.character(bp.pos.fdr[, 1])))), ]
		bp.neg.fdr <- bp.neg.fdr[!is.na(as.vector(Term(as.character(bp.neg.fdr[, 1])))), ]
		mf.pos.fdr <- mf.pos.fdr[!is.na(as.vector(Term(as.character(mf.pos.fdr[, 1])))), ]
		mf.neg.fdr <- mf.neg.fdr[!is.na(as.vector(Term(as.character(mf.neg.fdr[, 1])))), ]
		this.bp <- unname(cbind(trait, "Male", age, "BP", as.character(bp.pos.fdr[, 1]), as.vector(Term(as.character(bp.pos.fdr[, 1]))), bp.pos.fdr[, 2:3], trait, "Male", age, "BP", 
                          as.character(bp.neg.fdr[, 1]), as.vector(Term(as.character(bp.neg.fdr[, 1]))), bp.neg.fdr[, 2:3]))
                          
		this.mf <- unname(cbind(trait, "Male", age, "MF", as.character(mf.pos.fdr[, 1]), as.vector(Term(as.character(mf.pos.fdr[, 1]))), mf.pos.fdr[, 2:3], trait, "Male", age, "MF", 
                          as.character(mf.neg.fdr[, 1]), as.vector(Term(as.character(mf.neg.fdr[, 1]))), mf.neg.fdr[, 2:3]))
		
		all.go.gsea <- rbind(all.go.gsea, as.matrix(this.bp))
		all.go.gsea <- rbind(all.go.gsea, as.matrix(this.mf))
	
	}
}

all.kegg.gsea <- NULL

for (trait in female.traits) {
	
	for (age in c("young", "old", "diff")) {
	
		load(paste("../figureData/female.", trait, ".qtt.", age, ".kegg.fdr.RData", sep = ""))

		this.kegg <- unname(cbind(trait, "Female", age, "KEGG", as.character(kegg.pos.fdr[, 1]), kegg.id[as.character(kegg.pos.fdr[, 1]), 2], kegg.pos.fdr[, 2:3], trait, "Female", age, "KEGG", as.character(kegg.neg.fdr[, 1]), kegg.id[as.character(kegg.neg.fdr[, 1]), 2], kegg.neg.fdr[, 2:3]))         
		
		all.kegg.gsea <- rbind(all.kegg.gsea, as.matrix(this.kegg))
	
	}
}

for (trait in male.traits) {
	
	for (age in c("young", "old", "diff")) {
	
		load(paste("../figureData/male.", trait, ".qtt.", age, ".kegg.fdr.RData", sep = ""))

		this.kegg <- unname(cbind(trait, "Male", age, "KEGG", as.character(kegg.pos.fdr[, 1]), kegg.id[as.character(kegg.pos.fdr[, 1]), 2], kegg.pos.fdr[, 2:3], trait, "Male", age, "KEGG", as.character(kegg.neg.fdr[, 1]), kegg.id[as.character(kegg.neg.fdr[, 1]), 2], kegg.neg.fdr[, 2:3]))         
		
		all.kegg.gsea <- rbind(all.kegg.gsea, as.matrix(this.kegg))
	
	}
}

# write excel file
# ============================================================

write.xlsx(list('go' = all.go.gsea, 'kegg' = all.kegg.gsea), file = args[4], colNames = TRUE, rowNames = FALSE)
