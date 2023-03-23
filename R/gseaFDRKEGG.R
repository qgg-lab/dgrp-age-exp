# ==========================
# = calculate FDR for GSEA =
# ==========================

args <- commandArgs(TRUE) # args <- c("gsea/female.gei.gsea.RData")
load(args[1])

# 1. normalize es score using means of the permuted data sets
# ============================================================

kegg.mean.es.pos <- rowMeans(kegg.es.perm[, seq(1, ncol(kegg.es.perm), 2)])
kegg.mean.es.neg <- rowMeans(kegg.es.perm[, seq(2, ncol(kegg.es.perm), 2)])

kegg.es.pos <- kegg.es[, 1]/kegg.mean.es.pos
kegg.es.neg <- kegg.es[, 2]/kegg.mean.es.neg

kegg.mean.es.pos.perm <- apply(kegg.es.perm[, seq(1, ncol(kegg.es.perm), 2)], 2, function(x) { return(x/kegg.mean.es.pos) })
kegg.mean.es.neg.perm <- apply(kegg.es.perm[, seq(2, ncol(kegg.es.perm), 2)], 2, function(x) { return(x/kegg.mean.es.neg) })


# 2. calculate FDR
# ============================================================

kegg.es.pos.sorted <- sort(kegg.es.pos, decreasing = TRUE)
kegg.es.neg.sorted <- sort(kegg.es.neg, decreasing = TRUE)

kegg.pos.fdr <- data.frame(go = names(kegg.es.pos.sorted), nes = kegg.es.pos.sorted, fdr = NA)
kegg.neg.fdr <- data.frame(go = names(kegg.es.neg.sorted), nes = kegg.es.neg.sorted, fdr = NA)

for (i in 1:length(kegg.es.pos.sorted)) {
  
  kegg.pos.fdr[i, 3] <- (sum(kegg.mean.es.pos.perm >= kegg.es.pos.sorted[i]) + 1)/(ncol(kegg.es.perm)/2 + 1)/i
  
}

for (i in 1:length(kegg.es.neg.sorted)) {
  
  kegg.neg.fdr[i, 3] <- (sum(kegg.mean.es.neg.perm >= kegg.es.neg.sorted[i]) + 1)/(ncol(kegg.es.perm)/2 + 1)/i
  
}

kegg.pos.fdr[, 3] <- ifelse(cummax(kegg.pos.fdr[, 3]) > 1, 1, cummax(kegg.pos.fdr[, 3]))
kegg.neg.fdr[, 3] <- ifelse(cummax(kegg.neg.fdr[, 3]) > 1, 1, cummax(kegg.neg.fdr[, 3]))



save(kegg.pos.fdr, kegg.neg.fdr, file = args[2])
