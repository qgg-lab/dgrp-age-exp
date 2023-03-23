# 1. run GSEA for canalization/decanalization score
# ============================================================
cd /mnt/research/qgg/dgrp-age-exp/joinExp

Rscript -e 'load("qg/female.adj.qgVarHet.RData"); load("qg/female.adj.qgSingleAge.RData"); sig.h2 <- which(p.adjust(single.age[, 9], method = "BH") < 0.05 | p.adjust(single.age[, 12], method = "BH") < 0.05); h2.3wk <- as.numeric(var.het[, 8])[sig.h2]; h2.baseline <- as.numeric(var.het[, 7])[sig.h2]; score <- (h2.3wk - h2.baseline)/(h2.3wk + h2.baseline); write.table(na.omit(cbind(gene.name[sig.h2], score)), file = "gsea/female.decan.score", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE, na = "NA");'

Rscript -e 'load("qg/male.adj.qgVarHet.RData"); load("qg/male.adj.qgSingleAge.RData"); sig.h2 <- which(p.adjust(single.age[, 9], method = "BH") < 0.05 | p.adjust(single.age[, 12], method = "BH") < 0.05); h2.3wk <- as.numeric(var.het[, 8])[sig.h2]; h2.baseline <- as.numeric(var.het[, 7])[sig.h2]; score <- (h2.3wk - h2.baseline)/(h2.3wk + h2.baseline); write.table(na.omit(cbind(gene.name[sig.h2], score)), file = "gsea/male.decan.score", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE, na = "NA");'

Rscript prepGSEA.R gsea/female.decan.score /mnt/research/qgg/flybase/fb-r5.57/gene_go.table 20 gsea/female.decan.gsea.data.RData > log/female.decan.prepGSEA.Rout 2>&1 &
Rscript prepGSEA.R gsea/male.decan.score /mnt/research/qgg/flybase/fb-r5.57/gene_go.table 20 gsea/male.decan.gsea.data.RData > log/male.decan.prepGSEA.Rout 2>&1 &

Rscript gsea.R gsea/female.decan.gsea.data.RData gsea/female.decan.score 1000 gsea/female.decan.gsea.RData > log/female.decan.gsea.Rout 2>&1 &
Rscript gsea.R gsea/male.decan.gsea.data.RData gsea/male.decan.score 1000 gsea/male.decan.gsea.RData > log/male.decan.gsea.Rout 2>&1 &

Rscript gseaFDR.R gsea/female.decan.gsea.RData gsea/female.decan.gsea.fdr.RData > log/female.decan.gsea.fdr.Rout 2>&1 &
Rscript gseaFDR.R gsea/male.decan.gsea.RData gsea/male.decan.gsea.fdr.RData > log/male.decan.gsea.fdr.Rout 2>&1 &

# extract examples
# ============================================================

#Rscript gseaExample.R gsea/female.decan.gsea.data.RData gsea/female.decan.score "MF>GO:0005214" gsea/female.decan.gsea.example.RData > log/female.decan.gsea.example.Rout 2>&1 &
#Rscript gseaExample.R gsea/gsea.RData gsea/male.decan.score "MF>GO:0005549,MF>GO:0003700" gsea/male.decan.gsea.example.RData > log/male.decan.gsea.example.Rout 2>&1 &

# 3. prepare data for GEI GSEA
# ============================================================

Rscript -e 'load("qg/female.adj.qgGxE.RData"); gei <- as.numeric(gxe[, 6])/(as.numeric(gxe[, 6]) + as.numeric(gxe[, 7])); int.fdr <- p.adjust(as.numeric(gxe[, 10]), method = "BH"); line.fdr <- p.adjust(as.numeric(gxe[, 9]), method = "BH"); write.table(na.omit(cbind(gene.name[int.fdr < 0.05 | line.fdr < 0.05], 2*gei[int.fdr < 0.05 | line.fdr < 0.05] - 1)), file = "gsea/female.gei.score", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE, na = "NA");'

Rscript -e 'load("qg/male.adj.qgGxE.RData"); gei <- as.numeric(gxe[, 6])/(as.numeric(gxe[, 6]) + as.numeric(gxe[, 7])); int.fdr <- p.adjust(as.numeric(gxe[, 10]), method = "BH"); line.fdr <- p.adjust(as.numeric(gxe[, 9]), method = "BH"); write.table(na.omit(cbind(gene.name[int.fdr < 0.05 | line.fdr < 0.05], 2*gei[int.fdr < 0.05 | line.fdr < 0.05] - 1)), file = "gsea/male.gei.score", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE, na = "NA");'

Rscript prepGSEA.R gsea/female.gei.score /mnt/research/qgg/flybase/fb-r5.57/gene_go.table 20 gsea/female.gei.gsea.data.RData > log/female.gei.prepGSEA.Rout 2>&1 &
Rscript prepGSEA.R gsea/male.gei.score /mnt/research/qgg/flybase/fb-r5.57/gene_go.table 20 gsea/male.gei.gsea.data.RData > log/male.gei.prepGSEA.Rout 2>&1 &

Rscript gsea.R gsea/female.gei.gsea.data.RData gsea/male.gei.score 1000 gsea/male.gei.gsea.RData > log/male.gei.gsea.Rout 2>&1 &
Rscript gsea.R gsea/male.gei.gsea.data.RData gsea/female.gei.score 1000 gsea/female.gei.gsea.RData > log/female.gei.gsea.Rout 2>&1 &

Rscript gseaFDR.R gsea/female.gei.gsea.RData gsea/female.gei.gsea.fdr.RData > log/female.gei.gsea.fdr.Rout 2>&1 &
Rscript gseaFDR.R gsea/male.gei.gsea.RData gsea/male.gei.gsea.fdr.RData > log/male.gei.gsea.fdr.Rout 2>&1 &

