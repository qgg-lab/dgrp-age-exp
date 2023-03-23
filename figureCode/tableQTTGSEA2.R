# =====================================================
# = make a table for the QTT GSEA analysis =
# =====================================================

args <- commandArgs(TRUE) # args <- c("../figureData/female.qtt.young.gsea.fdr.RData", "../figureData/female.qtt.old.gsea.fdr.RData", "../figureData/female.qtt.diff.gsea.fdr.RData", "../figureData/male.qtt.young.gsea.fdr.RData", "../figureData/male.qtt.old.gsea.fdr.RData", "../figureData/male.qtt.diff.gsea.fdr.RData", "../figureData/female.qtt.young.kegg.fdr.RData", "../figureData/female.qtt.old.kegg.fdr.RData", "../figureData/female.qtt.diff.kegg.fdr.RData", "../figureData/male.qtt.young.kegg.fdr.RData", "../figureData/male.qtt.old.kegg.fdr.RData", "../figureData/male.qtt.diff.kegg.fdr.RData", "../figureData/kegg.id")
library("openxlsx")
library("GO.db")

# read and process female data
# ============================================================
kegg.id <- read.table(args[13], header = FALSE, as.is = TRUE, sep = "\t")                        
rownames(kegg.id) <- kegg.id[, 1]


load(args[1])

bp.pos.fdr <- bp.pos.fdr[!is.na(as.vector(Term(as.character(bp.pos.fdr[, 1])))), ]
bp.neg.fdr <- bp.neg.fdr[!is.na(as.vector(Term(as.character(bp.neg.fdr[, 1])))), ]
mf.pos.fdr <- mf.pos.fdr[!is.na(as.vector(Term(as.character(mf.pos.fdr[, 1])))), ]
mf.neg.fdr <- mf.neg.fdr[!is.na(as.vector(Term(as.character(mf.neg.fdr[, 1])))), ]
female.young.bp <- unname(cbind("Female", "Young", "BP", as.character(bp.pos.fdr[, 1]), as.vector(Term(as.character(bp.pos.fdr[, 1]))), bp.pos.fdr[, 2:3],
                          as.character(bp.neg.fdr[, 1]), as.vector(Term(as.character(bp.neg.fdr[, 1]))), bp.neg.fdr[, 2:3]))
                          
female.young.mf <- unname(cbind("Female", "Young", "MF", as.character(mf.pos.fdr[, 1]), as.vector(Term(as.character(mf.pos.fdr[, 1]))), mf.pos.fdr[, 2:3],
                          as.character(mf.neg.fdr[, 1]), as.vector(Term(as.character(mf.neg.fdr[, 1]))), mf.neg.fdr[, 2:3]))
                        
load(args[2])

bp.pos.fdr <- bp.pos.fdr[!is.na(as.vector(Term(as.character(bp.pos.fdr[, 1])))), ]
bp.neg.fdr <- bp.neg.fdr[!is.na(as.vector(Term(as.character(bp.neg.fdr[, 1])))), ]
mf.pos.fdr <- mf.pos.fdr[!is.na(as.vector(Term(as.character(mf.pos.fdr[, 1])))), ]
mf.neg.fdr <- mf.neg.fdr[!is.na(as.vector(Term(as.character(mf.neg.fdr[, 1])))), ]
female.old.bp <- unname(cbind("Female", "Old", "BP", as.character(bp.pos.fdr[, 1]), as.vector(Term(as.character(bp.pos.fdr[, 1]))), bp.pos.fdr[, 2:3],
                        as.character(bp.neg.fdr[, 1]), as.vector(Term(as.character(bp.neg.fdr[, 1]))), bp.neg.fdr[, 2:3]))
                          
female.old.mf <- unname(cbind("Female", "Old", "MF", as.character(mf.pos.fdr[, 1]), as.vector(Term(as.character(mf.pos.fdr[, 1]))), mf.pos.fdr[, 2:3],
                        as.character(mf.neg.fdr[, 1]), as.vector(Term(as.character(mf.neg.fdr[, 1]))), mf.neg.fdr[, 2:3]))
                        
load(args[3])

bp.pos.fdr <- bp.pos.fdr[!is.na(as.vector(Term(as.character(bp.pos.fdr[, 1])))), ]
bp.neg.fdr <- bp.neg.fdr[!is.na(as.vector(Term(as.character(bp.neg.fdr[, 1])))), ]
mf.pos.fdr <- mf.pos.fdr[!is.na(as.vector(Term(as.character(mf.pos.fdr[, 1])))), ]
mf.neg.fdr <- mf.neg.fdr[!is.na(as.vector(Term(as.character(mf.neg.fdr[, 1])))), ]
female.diff.bp <- unname(cbind("Female", "Old - Young", "BP", as.character(bp.pos.fdr[, 1]), as.vector(Term(as.character(bp.pos.fdr[, 1]))), bp.pos.fdr[, 2:3],
                        as.character(bp.neg.fdr[, 1]), as.vector(Term(as.character(bp.neg.fdr[, 1]))), bp.neg.fdr[, 2:3]))
                          
female.diff.mf <- unname(cbind("Female", "Old - Young", "MF", as.character(mf.pos.fdr[, 1]), as.vector(Term(as.character(mf.pos.fdr[, 1]))), mf.pos.fdr[, 2:3],
                        as.character(mf.neg.fdr[, 1]), as.vector(Term(as.character(mf.neg.fdr[, 1]))), mf.neg.fdr[, 2:3]))

load(args[4])

bp.pos.fdr <- bp.pos.fdr[!is.na(as.vector(Term(as.character(bp.pos.fdr[, 1])))), ]
bp.neg.fdr <- bp.neg.fdr[!is.na(as.vector(Term(as.character(bp.neg.fdr[, 1])))), ]
mf.pos.fdr <- mf.pos.fdr[!is.na(as.vector(Term(as.character(mf.pos.fdr[, 1])))), ]
mf.neg.fdr <- mf.neg.fdr[!is.na(as.vector(Term(as.character(mf.neg.fdr[, 1])))), ]
male.young.bp <- unname(cbind("Male", "Young", "BP", as.character(bp.pos.fdr[, 1]), as.vector(Term(as.character(bp.pos.fdr[, 1]))), bp.pos.fdr[, 2:3],
                          as.character(bp.neg.fdr[, 1]), as.vector(Term(as.character(bp.neg.fdr[, 1]))), bp.neg.fdr[, 2:3]))
                          
male.young.mf <- unname(cbind("Male", "Young", "MF", as.character(mf.pos.fdr[, 1]), as.vector(Term(as.character(mf.pos.fdr[, 1]))), mf.pos.fdr[, 2:3],
                          as.character(mf.neg.fdr[, 1]), as.vector(Term(as.character(mf.neg.fdr[, 1]))), mf.neg.fdr[, 2:3]))
                        
load(args[5])

bp.pos.fdr <- bp.pos.fdr[!is.na(as.vector(Term(as.character(bp.pos.fdr[, 1])))), ]
bp.neg.fdr <- bp.neg.fdr[!is.na(as.vector(Term(as.character(bp.neg.fdr[, 1])))), ]
mf.pos.fdr <- mf.pos.fdr[!is.na(as.vector(Term(as.character(mf.pos.fdr[, 1])))), ]
mf.neg.fdr <- mf.neg.fdr[!is.na(as.vector(Term(as.character(mf.neg.fdr[, 1])))), ]
male.old.bp <- unname(cbind("Male", "Old", "BP", as.character(bp.pos.fdr[, 1]), as.vector(Term(as.character(bp.pos.fdr[, 1]))), bp.pos.fdr[, 2:3],
                        as.character(bp.neg.fdr[, 1]), as.vector(Term(as.character(bp.neg.fdr[, 1]))), bp.neg.fdr[, 2:3]))
                          
male.old.mf <- unname(cbind("Male", "Old", "MF", as.character(mf.pos.fdr[, 1]), as.vector(Term(as.character(mf.pos.fdr[, 1]))), mf.pos.fdr[, 2:3],
                        as.character(mf.neg.fdr[, 1]), as.vector(Term(as.character(mf.neg.fdr[, 1]))), mf.neg.fdr[, 2:3]))
                        
load(args[6])

bp.pos.fdr <- bp.pos.fdr[!is.na(as.vector(Term(as.character(bp.pos.fdr[, 1])))), ]
bp.neg.fdr <- bp.neg.fdr[!is.na(as.vector(Term(as.character(bp.neg.fdr[, 1])))), ]
mf.pos.fdr <- mf.pos.fdr[!is.na(as.vector(Term(as.character(mf.pos.fdr[, 1])))), ]
mf.neg.fdr <- mf.neg.fdr[!is.na(as.vector(Term(as.character(mf.neg.fdr[, 1])))), ]
male.diff.bp <- unname(cbind("Male", "Old - Young", "BP", as.character(bp.pos.fdr[, 1]), as.vector(Term(as.character(bp.pos.fdr[, 1]))), bp.pos.fdr[, 2:3],
                        as.character(bp.neg.fdr[, 1]), as.vector(Term(as.character(bp.neg.fdr[, 1]))), bp.neg.fdr[, 2:3]))
                          
male.diff.mf <- unname(cbind("Male", "Old - Young", "MF", as.character(mf.pos.fdr[, 1]), as.vector(Term(as.character(mf.pos.fdr[, 1]))), mf.pos.fdr[, 2:3],
                        as.character(mf.neg.fdr[, 1]), as.vector(Term(as.character(mf.neg.fdr[, 1]))), mf.neg.fdr[, 2:3]))


load(args[7])
female.young.kegg <- unname(cbind("Female", "Young", "KEGG", as.character(kegg.pos.fdr[, 1]), kegg.id[as.character(kegg.pos.fdr[, 1]), 2], kegg.pos.fdr[, 2:3], as.character(kegg.neg.fdr[, 1]), kegg.id[as.character(kegg.neg.fdr[, 1]), 2], kegg.neg.fdr[, 2:3]))
  
load(args[8])
female.old.kegg <- unname(cbind("Female", "Old", "KEGG", as.character(kegg.pos.fdr[, 1]), kegg.id[as.character(kegg.pos.fdr[, 1]), 2], kegg.pos.fdr[, 2:3], as.character(kegg.neg.fdr[, 1]), kegg.id[as.character(kegg.neg.fdr[, 1]), 2], kegg.neg.fdr[, 2:3]))


load(args[9])
female.diff.kegg <- unname(cbind("Female", "Old - Young", "KEGG", as.character(kegg.pos.fdr[, 1]), kegg.id[as.character(kegg.pos.fdr[, 1]), 2], kegg.pos.fdr[, 2:3], as.character(kegg.neg.fdr[, 1]), kegg.id[as.character(kegg.neg.fdr[, 1]), 2], kegg.neg.fdr[, 2:3]))
  

load(args[10])
male.young.kegg <- unname(cbind("Male", "Young", "KEGG", as.character(kegg.pos.fdr[, 1]), kegg.id[as.character(kegg.pos.fdr[, 1]), 2], kegg.pos.fdr[, 2:3], as.character(kegg.neg.fdr[, 1]), kegg.id[as.character(kegg.neg.fdr[, 1]), 2], kegg.neg.fdr[, 2:3]))
  
load(args[11])
male.old.kegg <- unname(cbind("Male", "Old", "KEGG", as.character(kegg.pos.fdr[, 1]), kegg.id[as.character(kegg.pos.fdr[, 1]), 2], kegg.pos.fdr[, 2:3], as.character(kegg.neg.fdr[, 1]), kegg.id[as.character(kegg.neg.fdr[, 1]), 2], kegg.neg.fdr[, 2:3]))


load(args[12])
male.diff.kegg <- unname(cbind("Male", "Old - Young", "KEGG", as.character(kegg.pos.fdr[, 1]), kegg.id[as.character(kegg.pos.fdr[, 1]), 2], kegg.pos.fdr[, 2:3], as.character(kegg.neg.fdr[, 1]), kegg.id[as.character(kegg.neg.fdr[, 1]), 2], kegg.neg.fdr[, 2:3]))
  
  
# combine data
# ============================================================

final.table <- rbind(as.matrix(female.young.bp), as.matrix(female.young.mf), as.matrix(female.old.bp), as.matrix(female.old.mf), as.matrix(female.diff.bp), as.matrix(female.diff.mf), as.matrix(male.young.bp), as.matrix(male.young.mf), as.matrix(male.old.bp), as.matrix(male.old.mf), as.matrix(male.diff.bp), as.matrix(male.diff.mf), as.matrix(female.young.kegg), as.matrix(female.old.kegg), as.matrix(female.diff.kegg), as.matrix(male.young.kegg), as.matrix(male.old.kegg), as.matrix(male.diff.kegg))

# write excel file
# ============================================================

write.xlsx(list('final' = final.table), file = args[14], colNames = TRUE, rowNames = FALSE)
