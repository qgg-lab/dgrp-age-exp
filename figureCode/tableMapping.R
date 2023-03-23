# ==============================================
# = make summary tables for mapping statistics =
# ==============================================

args <- commandArgs(TRUE) # args <- c("../figureData/BioSampleObjects.txt", "../figureData/BioSampleObjects2.txt", "../figureData/baseline.biosample.txt", "../figureData/bad.sample.id", "../figureData/sample.map.summary")
library(xlsx)

# read data
# ============================================================

old.biosamples <- read.table(args[1], header = TRUE, as.is = TRUE, sep = "\t")
new.young.biosamples <- read.table(args[2], header = TRUE, as.is = TRUE, sep = "\t")
young.biosamples <- read.table(args[3], header = FALSE, as.is = TRUE, sep = "\t")
bad.sample <- read.table(args[4], header = FALSE, as.is = TRUE)[, 1]
map.summary <- read.table(args[5], header = FALSE, as.is = TRUE, sep = "\t")

# combine data, get unique ID
# ============================================================

map.summary$id <- paste("DGRP", "Line", map.summary[, 7], map.summary[, 6], map.summary[, 8], "Rep", map.summary[, 9], sep = "_")
map.summary.reads <- sapply(split(map.summary[, 1], map.summary$id), sum)
map.summary.non.mapped <- sapply(split(map.summary[, 2], map.summary$id), sum)
map.summary.unique.mapped <- sapply(split(map.summary[, 3], map.summary$id), sum)
map.summary.multi.mapped <- sapply(split(map.summary[, 4], map.summary$id), sum)

rownames(old.biosamples) <- old.biosamples[, 2]
rownames(new.young.biosamples) <- new.young.biosamples[, 2]
rownames(young.biosamples) <- paste("DGRP_Line", gsub("Male", "baseline_M", gsub("Female", "baseline_F", young.biosamples[, 1])), sep = "_")

biosamples <- rbind(as.matrix(old.biosamples[, c(1, 2)]), as.matrix(new.young.biosamples[, c(1, 2)]), as.matrix(young.biosamples[, c(2, 1)]))

bad.sample <- gsub("M", "M_Rep_", gsub("F", "F_Rep_", apply(matrix(unlist(strsplit(bad.sample, split = "_")), ncol = 3, byrow = T), 1, function(x) { return(paste("DGRP_Line", x[2], x[1], x[3], sep = "_")) })))

biosamples.summary <- biosamples[rownames(biosamples) %in% map.summary$id & !(rownames(biosamples) %in% bad.sample), ]

final.table <- cbind(biosamples.summary[, 1], t(apply(matrix(unlist(strsplit(rownames(biosamples.summary), split = "_")), ncol = 7, byrow = T), 1, function(x) { return( c(paste("R", x[3], sep = ""), ifelse(x[4] == "3wk", "3 weeks", "3-5 days"), x[5], x[7]) )})), map.summary.reads[rownames(biosamples.summary)], map.summary.unique.mapped[rownames(biosamples.summary)], map.summary.multi.mapped[rownames(biosamples.summary)], map.summary.non.mapped[rownames(biosamples.summary)])

colnames(final.table) <- c("Biosamples", "DGRP", "Age", "Sex", "Replicate", "Sequenced_reads", "Uniquely_mapped_reads", "Multi_mapped_reads", "Non_mapped_reads")

# summary statistics
cat("on average", mean(as.numeric(final.table[, 6])), "reads were sequenced.\n")
cat("that range", range(as.numeric(final.table[, 6])), "\n")

# write excel file
# ============================================================
write.xlsx(final.table, file = args[6], col.names = TRUE, row.names = FALSE)
