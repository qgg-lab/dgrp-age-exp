# 1. run GSEA for aging effects, convert t statistic to correlation
# ============================================================

~/qgg/software/R-3.6.0/bin/Rscript -e 'load("figureData/female.adj.qgAgeEff.RData"); t <- age.effect[, 4]/age.effect[, 5]; r <- sign(t)*sqrt(t^2/(198 + t^2)); write.table(cbind(gene.name, r), file = "gsea/female.age.effect.score", sep = " ", col.names = F, row.names = F, quote = F, na = "NA"); load("figureData/male.adj.qgAgeEff.RData"); t <- age.effect[, 4]/age.effect[, 5]; r <- sign(t)*sqrt(t^2/(198 + t^2)); write.table(cbind(gene.name, r), file = "gsea/male.age.effect.score", sep = " ", col.names = F, row.names = F, quote = F, na = "NA");'

~/qgg/software/R-3.6.0/bin/Rscript R/prepGSEA.R gsea/female.age.effect.score /mnt/research/qgg/resource/flybase/r6.36/annot/gene_go.table 20 gsea/female.age.effect.gsea.data.RData > log/female.age.effect.prepGSEA.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/prepGSEA.R gsea/male.age.effect.score /mnt/research/qgg/resource/flybase/r6.36/annot/gene_go.table 20 gsea/male.age.effect.gsea.data.RData > log/male.age.effect.prepGSEA.Rout 2>&1 &


~/qgg/software/R-3.6.0/bin/Rscript R/prepKEGG.R gsea/female.age.effect.score /mnt/research/qgg/resource/flybase/r6.36/annot/fbgn.kegg 20 gsea/female.age.effect.kegg.data.RData > log/female.age.effect.prepKEGG.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/prepKEGG.R gsea/male.age.effect.score /mnt/research/qgg/resource/flybase/r6.36/annot/fbgn.kegg 20 gsea/male.age.effect.kegg.data.RData > log/male.age.effect.prepKEGG.Rout 2>&1 &


~/qgg/software/R-3.6.0/bin/Rscript R/gsea.R gsea/female.age.effect.gsea.data.RData gsea/female.age.effect.score 1000 gsea/female.age.effect.gsea.RData > log/female.age.effect.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/gsea.R gsea/male.age.effect.gsea.data.RData gsea/male.age.effect.score 1000 gsea/male.age.effect.gsea.RData > log/male.age.effect.gsea.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript R/gseaKEGG.R gsea/female.age.effect.kegg.data.RData gsea/female.age.effect.score 1000 gsea/female.age.effect.kegg.RData > log/female.age.effect.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/gseaKEGG.R gsea/male.age.effect.kegg.data.RData gsea/male.age.effect.score 1000 gsea/male.age.effect.kegg.RData > log/male.age.effect.kegg.Rout 2>&1 &


~/qgg/software/R-3.6.0/bin/Rscript R/gseaFDR.R gsea/female.age.effect.gsea.RData figureData/female.age.effect.gsea.fdr.RData > log/female.age.effect.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/gseaFDR.R gsea/male.age.effect.gsea.RData figureData/male.age.effect.gsea.fdr.RData > log/male.age.effect.gsea.fdr.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript R/gseaFDRKEGG.R gsea/female.age.effect.kegg.RData figureData/female.age.effect.kegg.fdr.RData > log/female.age.effect.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/gseaFDRKEGG.R gsea/male.age.effect.kegg.RData figureData/male.age.effect.kegg.fdr.RData > log/male.age.effect.kegg.fdr.Rout 2>&1 &


# gsea example
~/qgg/software/R-3.6.0/bin/Rscript R/gseaExample.R gsea/female.age.effect.gsea.data.RData gsea/female.age.effect.score "BP>GO:0006099" gsea/female.age.effect.kegg.data.RData "dme04213,dme00020" figureData/female.age.effect.gsea.example.RData > log/female.age.effect.gsea.example.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/gseaExample.R gsea/male.age.effect.gsea.data.RData gsea/male.age.effect.score "BP>GO:0006099" gsea/male.age.effect.kegg.data.RData "dme04213,dme00020" figureData/male.age.effect.gsea.example.RData > log/male.age.effect.gsea.example.Rout 2>&1 &

# 2. run GSEA for GxE
# ============================================================

~/qgg/software/R-3.6.0/bin/Rscript -e 'load("figureData/female.adj.qgGxE.RData"); gei <- as.numeric(gxe[, 9])/(as.numeric(gxe[, 9]) + as.numeric(gxe[, 10])); int.fdr <- p.adjust(as.numeric(gxe[, 13]), method = "BH"); line.fdr <- p.adjust(as.numeric(gxe[, 12]), method = "BH"); write.table(na.omit(cbind(gene.name[int.fdr < 0.05 | line.fdr < 0.05], 2*gei[int.fdr < 0.05 | line.fdr < 0.05] - 1)), file = "gsea/female.gei.score", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE, na = "NA");'

~/qgg/software/R-3.6.0/bin/Rscript -e 'load("figureData/male.adj.qgGxE.RData"); gei <- as.numeric(gxe[, 9])/(as.numeric(gxe[, 9]) + as.numeric(gxe[, 10])); int.fdr <- p.adjust(as.numeric(gxe[, 13]), method = "BH"); line.fdr <- p.adjust(as.numeric(gxe[, 12]), method = "BH"); write.table(na.omit(cbind(gene.name[int.fdr < 0.05 | line.fdr < 0.05], 2*gei[int.fdr < 0.05 | line.fdr < 0.05] - 1)), file = "gsea/male.gei.score", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE, na = "NA");'

~/qgg/software/R-3.6.0/bin/Rscript R/prepGSEA.R gsea/female.gei.score /mnt/research/qgg/resource/flybase/r6.36/annot/gene_go.table 20 gsea/female.gei.gsea.data.RData > log/female.gei.prepGSEA.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/prepGSEA.R gsea/male.gei.score /mnt/research/qgg/resource/flybase/r6.36/annot/gene_go.table 20 gsea/male.gei.gsea.data.RData > log/male.gei.prepGSEA.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript R/prepKEGG.R gsea/female.gei.score /mnt/research/qgg/resource/flybase/r6.36/annot/fbgn.kegg 20 gsea/female.gei.kegg.data.RData > log/female.gei.prepKEGG.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/prepKEGG.R gsea/male.gei.score /mnt/research/qgg/resource/flybase/r6.36/annot/fbgn.kegg 20 gsea/male.gei.kegg.data.RData > log/male.gei.prepKEGG.Rout 2>&1 &


~/qgg/software/R-3.6.0/bin/Rscript R/gsea.R gsea/female.gei.gsea.data.RData gsea/female.gei.score 1000 gsea/female.gei.gsea.RData > log/female.gei.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/gsea.R gsea/male.gei.gsea.data.RData gsea/male.gei.score 1000 gsea/male.gei.gsea.RData > log/male.gei.gsea.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript R/gseaKEGG.R gsea/female.gei.kegg.data.RData gsea/female.gei.score 1000 gsea/female.gei.kegg.RData > log/female.gei.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/gseaKEGG.R gsea/male.gei.kegg.data.RData gsea/male.gei.score 1000 gsea/male.gei.kegg.RData > log/male.gei.kegg.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript R/gseaFDR.R gsea/female.gei.gsea.RData figureData/female.gei.gsea.fdr.RData > log/female.gei.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/gseaFDR.R gsea/male.gei.gsea.RData figureData/male.gei.gsea.fdr.RData > log/male.gei.gsea.fdr.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript R/gseaFDRKEGG.R gsea/female.gei.kegg.RData figureData/female.gei.kegg.fdr.RData > log/female.gei.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/gseaFDRKEGG.R gsea/male.gei.kegg.RData figureData/male.gei.kegg.fdr.RData > log/male.gei.kegg.fdr.Rout 2>&1 &

# 3. run GSEA for canalization/decanalization score
# ============================================================

~/qgg/software/R-3.6.0/bin/Rscript -e 'load("qg/female.adj.qgVarHet.RData"); load("qg/female.adj.qgSingleTemp.RData"); sig.h2 <- which(p.adjust(single.temp[, 9], method = "BH") < 0.05 | p.adjust(single.temp[, 12], method = "BH") < 0.05); h2.18c <- as.numeric(var.het[, 8])[sig.h2]; h2.25c <- as.numeric(var.het[, 7])[sig.h2]; score <- (h2.18c - h2.25c)/(h2.18c + h2.25c); write.table(na.omit(cbind(gene.name[sig.h2], score)), file = "gsea/female.decan.score", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE, na = "NA");'

~/qgg/software/R-3.6.0/bin/Rscript -e 'load("qg/male.adj.qgVarHet.RData"); load("qg/male.adj.qgSingleTemp.RData"); sig.h2 <- which(p.adjust(single.temp[, 9], method = "BH") < 0.05 | p.adjust(single.temp[, 12], method = "BH") < 0.05); h2.18c <- as.numeric(var.het[, 8])[sig.h2]; h2.25c <- as.numeric(var.het[, 7])[sig.h2]; score <- (h2.18c - h2.25c)/(h2.18c + h2.25c); write.table(na.omit(cbind(gene.name[sig.h2], score)), file = "gsea/male.decan.score", col.names = FALSE, row.names = FALSE, sep = " ", quote = FALSE, na = "NA");'

~/qgg/software/R-3.6.0/bin/Rscript R/prepGSEA.R gsea/female.decan.score /mnt/research/qgg/resource/flybase/r6.36/annot/gene_go.table 20 gsea/female.decan.gsea.data.RData > log/female.decan.prepGSEA.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/prepGSEA.R gsea/male.decan.score /mnt/research/qgg/resource/flybase/r6.36/annot/gene_go.table 20 gsea/male.decan.gsea.data.RData > log/male.decan.prepGSEA.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript R/prepKEGG.R gsea/female.decan.score /mnt/research/qgg/resource/flybase/r6.36/annot/fbgn.kegg 20 gsea/female.decan.kegg.data.RData > log/female.decan.prepKEGG.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/prepKEGG.R gsea/male.decan.score /mnt/research/qgg/resource/flybase/r6.36/annot/fbgn.kegg 20 gsea/male.decan.kegg.data.RData > log/male.decan.prepKEGG.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript R/gsea.R gsea/female.decan.gsea.data.RData gsea/female.decan.score 1000 gsea/female.decan.gsea.RData > log/female.decan.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/gsea.R gsea/male.decan.gsea.data.RData gsea/male.decan.score 1000 gsea/male.decan.gsea.RData > log/male.decan.gsea.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript R/gseaKEGG.R gsea/female.decan.kegg.data.RData gsea/female.decan.score 1000 gsea/female.decan.kegg.RData > log/female.decan.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/gseaKEGG.R gsea/male.decan.kegg.data.RData gsea/male.decan.score 1000 gsea/male.decan.kegg.RData > log/male.decan.kegg.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript R/gseaFDR.R gsea/female.decan.gsea.RData figureData/female.decan.gsea.fdr.RData > log/female.decan.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/gseaFDR.R gsea/male.decan.gsea.RData figureData/male.decan.gsea.fdr.RData > log/male.decan.gsea.fdr.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript R/gseaFDRKEGG.R gsea/female.decan.kegg.RData figureData/female.decan.kegg.fdr.RData > log/female.decan.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript R/gseaFDRKEGG.R gsea/male.decan.kegg.RData figureData/male.decan.kegg.fdr.RData > log/male.decan.kegg.fdr.Rout 2>&1 &
