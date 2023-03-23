# ================================================
# = prepare annotation data sets for KEGG to use =
# ================================================

args <- commandArgs(TRUE)

gene.id <- read.table(args[1], header = FALSE, as.is = TRUE)[, 1]
annot <- read.table(args[2], header = FALSE, as.is = TRUE, sep = "\t") # annotation table
min.kegg <- as.numeric(args[3])

# subset annotations
# ============================================================

kegg.annot <- annot[annot[, 1] %in% gene.id, ]
kegg.gene.id <- unique(kegg.annot[, 1])
unique.kegg <- unique(kegg.annot[, 2])

# prepare list of GO genes
# ============================================================

kegg.gene.set <- list()
i = 1
for (kegg in unique.kegg) {
  kegg.genes <- kegg.annot[kegg.annot[, 2] == kegg, 1]
  if (length(kegg.genes) >= min.kegg) {
    kegg.gene.set[[i]] <- kegg.genes
    names(kegg.gene.set)[i] <- kegg
    i = i + 1
  }
}

# save data set
save(kegg.gene.set, kegg.gene.id, file = args[4])
