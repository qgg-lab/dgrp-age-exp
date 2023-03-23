# =======================
# = GSEA main R program =
# =======================

args <- commandArgs(TRUE)
load(args[1]) # load GSEA annotation data
score.list <- read.table(args[2], header = FALSE, as.is = TRUE)

# rank the score list
# ============================================================

score.list <- score.list[order(score.list[, 2], decreasing = TRUE), ]
rownames(score.list) <- score.list[, 1]

# function to calculate ES
# takes 1) gene.set (set gene IDs), which has been processed
# to contain only gene IDs present in the score list
# this is not exactly necessary but will improve efficiency
# if pre-computed
# 2) score list
# ============================================================

gsea <- function(gene.set, score) {
  
  N <- nrow(score)
  
  p <- rep(-1, N)
  S <- sort(match(gene.set, score[, 1]))
  p[S] <- abs(score[S, 2])/sum(abs(score[S, 2]))
  p[-S] <- p[-S]/(N - length(S))
  
  p.cumsum <- cumsum(p)
  return(c(
    max(c(max(p.cumsum), 0)),
    min(c(0, min(p.cumsum)))
  ))
  
}

# loop the main program
# ============================================================

kegg.es <- matrix(NA, nrow = length(kegg.gene.set), ncol = 2)

kegg.score <- score.list[score.list[, 1] %in% kegg.gene.id, ]

rownames(kegg.es) <- names(kegg.gene.set)

for (i in 1:length(kegg.gene.set)) {
  
  kegg.es[i, ] <- gsea(kegg.gene.set[[i]], kegg.score)
  
}

# permutation results
# ============================================================

kegg.es.perm <- matrix(NA, nrow = length(kegg.gene.set), ncol = 2*as.numeric(args[3]))

set.seed(as.numeric(args[3]))

for (i in 1:as.numeric(args[3])) {
  
  score.permute <- data.frame(gene = sample(score.list[, 1]), score = score.list[, 2])
  
  kegg.score.permute <- score.permute[score.permute[, 1] %in% kegg.gene.id, ]
 
  for (j in 1:length(kegg.gene.set)) {
    kegg.es.perm[j, c(2*i - 1, 2*i)] <- gsea(kegg.gene.set[[j]], kegg.score.permute)
  }
  
  cat(i, "\n")
  
}

# save results
# ============================================================

save(kegg.es, kegg.es.perm, file = args[4])
