# =====================================================
# = make a table for the GxA GSEA analysis =
# =====================================================

args <- commandArgs(TRUE) # args <- c("../figureData/female.gei.gsea.fdr.RData", "../figureData/male.gei.gsea.fdr.RData", "../figureData/female.gei.kegg.fdr.RData", "../figureData/male.gei.kegg.fdr.RData", "../figureData/kegg.id")
library("openxlsx")
library("GO.db")

# read and process female data
# ============================================================
load(args[1])

bp.pos.fdr <- bp.pos.fdr[!is.na(as.vector(Term(as.character(bp.pos.fdr[, 1])))), ]
bp.neg.fdr <- bp.neg.fdr[!is.na(as.vector(Term(as.character(bp.neg.fdr[, 1])))), ]
mf.pos.fdr <- mf.pos.fdr[!is.na(as.vector(Term(as.character(mf.pos.fdr[, 1])))), ]
mf.neg.fdr <- mf.neg.fdr[!is.na(as.vector(Term(as.character(mf.neg.fdr[, 1])))), ]
female.bp <- unname(cbind("BP", as.character(bp.pos.fdr[, 1]), as.vector(Term(as.character(bp.pos.fdr[, 1]))), bp.pos.fdr[, 2:3],
                          as.character(bp.neg.fdr[, 1]), as.vector(Term(as.character(bp.neg.fdr[, 1]))), bp.neg.fdr[, 2:3]))
                          
female.mf <- unname(cbind("MF", as.character(mf.pos.fdr[, 1]), as.vector(Term(as.character(mf.pos.fdr[, 1]))), mf.pos.fdr[, 2:3],
                          as.character(mf.neg.fdr[, 1]), as.vector(Term(as.character(mf.neg.fdr[, 1]))), mf.neg.fdr[, 2:3]))
                        
load(args[2])
bp.pos.fdr <- bp.pos.fdr[!is.na(as.vector(Term(as.character(bp.pos.fdr[, 1])))), ]
bp.neg.fdr <- bp.neg.fdr[!is.na(as.vector(Term(as.character(bp.neg.fdr[, 1])))), ]
mf.pos.fdr <- mf.pos.fdr[!is.na(as.vector(Term(as.character(mf.pos.fdr[, 1])))), ]
mf.neg.fdr <- mf.neg.fdr[!is.na(as.vector(Term(as.character(mf.neg.fdr[, 1])))), ]
male.bp <- unname(cbind("BP", as.character(bp.pos.fdr[, 1]), as.vector(Term(as.character(bp.pos.fdr[, 1]))), bp.pos.fdr[, 2:3],
                        as.character(bp.neg.fdr[, 1]), as.vector(Term(as.character(bp.neg.fdr[, 1]))), bp.neg.fdr[, 2:3]))
                          
male.mf <- unname(cbind("MF", as.character(mf.pos.fdr[, 1]), as.vector(Term(as.character(mf.pos.fdr[, 1]))), mf.pos.fdr[, 2:3],
                        as.character(mf.neg.fdr[, 1]), as.vector(Term(as.character(mf.neg.fdr[, 1]))), mf.neg.fdr[, 2:3]))
                        
                        
												kegg.id <- read.table(args[5], header = FALSE, as.is = TRUE, sep = "\t")                        
												rownames(kegg.id) <- kegg.id[, 1]
												
load(args[3])
female.kegg <- unname(cbind("KEGG", as.character(kegg.pos.fdr[, 1]), kegg.id[as.character(kegg.pos.fdr[, 1]), 2], kegg.pos.fdr[, 2:3], as.character(kegg.neg.fdr[, 1]), kegg.id[as.character(kegg.neg.fdr[, 1]), 2], kegg.neg.fdr[, 2:3]))

load(args[4])
male.kegg <- unname(cbind("KEGG", as.character(kegg.pos.fdr[, 1]), kegg.id[as.character(kegg.pos.fdr[, 1]), 2], kegg.pos.fdr[, 2:3],	 as.character(kegg.neg.fdr[, 1]), kegg.id[as.character(kegg.neg.fdr[, 1]), 2], kegg.neg.fdr[, 2:3]))

# combine data
# ============================================================


final.table <- rbind(cbind("Female", rbind(as.matrix(female.bp), as.matrix(female.mf))),
                     cbind("Male", rbind(as.matrix(male.bp), as.matrix(male.mf))),
										 cbind("Female", as.matrix(female.kegg)),
										 cbind("Male", as.matrix(male.kegg)))


# write excel file
# ============================================================

write.xlsx(list('final' = final.table), file = args[6], colNames = TRUE, rowNames = FALSE)
