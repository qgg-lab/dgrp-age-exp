# 0. prepare phenotype, including lifespan, fecundity, 
#    visual decline, activity decline
# ============================================================

# run qttPrepPheno.R

# 1. QTT detection
# ============================================================

~/qgg/software/R-3.6.0/bin/Rscript ../R/qtt.R ../qg/female.sig.gene.exp.RData female.fecundity.decline.pheno ../eqtl/dgrp.common.fam female.fecundity.decline.qtt.RData > female.fecundity.decline.qtt.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/qtt.R ../qg/female.sig.gene.exp.RData female.phototaxis.decline.pheno ../eqtl/dgrp.common.fam female.phototaxis.decline.qtt.RData > female.phototaxis.decline.qtt.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/qtt.R ../qg/female.sig.gene.exp.RData female.lifespan.pheno ../eqtl/dgrp.common.fam female.lifespan.qtt.RData > female.lifespan.qtt.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript ../R/qtt.R ../qg/male.sig.gene.exp.RData male.phototaxis.decline.pheno ../eqtl/dgrp.common.fam male.phototaxis.decline.qtt.RData > male.phototaxis.decline.qtt.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/qtt.R ../qg/male.sig.gene.exp.RData male.speed.decline.pheno ../eqtl/dgrp.common.fam male.speed.decline.qtt.RData > male.speed.decline.qtt.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/qtt.R ../qg/male.sig.gene.exp.RData male.endurance.decline.pheno ../eqtl/dgrp.common.fam male.endurance.decline.qtt.RData > male.endurance.decline.qtt.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/qtt.R ../qg/male.sig.gene.exp.RData male.lifespan.pheno ../eqtl/dgrp.common.fam male.lifespan.qtt.RData > male.lifespan.qtt.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript -e 'load("female.fecundity.decline.qtt.RData"); write.table(cbind(gene.name, qtt.cor.young), file = "female.fecundity.decline.qtt.young.score", sep = " ", row.names = F, col.names = F, quote = F); write.table(cbind(gene.name, qtt.cor.old), file = "female.fecundity.decline.qtt.old.score", sep = " ", row.names = F, col.names = F, quote = F); write.table(cbind(gene.name, qtt.cor.diff), file = "female.fecundity.decline.qtt.diff.score", sep = " ", row.names = F, col.names = F, quote = F);'
~/qgg/software/R-3.6.0/bin/Rscript -e 'load("female.phototaxis.decline.qtt.RData"); write.table(cbind(gene.name, qtt.cor.young), file = "female.phototaxis.decline.qtt.young.score", sep = " ", row.names = F, col.names = F, quote = F); write.table(cbind(gene.name, qtt.cor.old), file = "female.phototaxis.decline.qtt.old.score", sep = " ", row.names = F, col.names = F, quote = F); write.table(cbind(gene.name, qtt.cor.diff), file = "female.phototaxis.decline.qtt.diff.score", sep = " ", row.names = F, col.names = F, quote = F);'
~/qgg/software/R-3.6.0/bin/Rscript -e 'load("female.lifespan.qtt.RData"); write.table(cbind(gene.name, qtt.cor.young), file = "female.lifespan.qtt.young.score", sep = " ", row.names = F, col.names = F, quote = F); write.table(cbind(gene.name, qtt.cor.old), file = "female.lifespan.qtt.old.score", sep = " ", row.names = F, col.names = F, quote = F); write.table(cbind(gene.name, qtt.cor.diff), file = "female.lifespan.qtt.diff.score", sep = " ", row.names = F, col.names = F, quote = F);'


~/qgg/software/R-3.6.0/bin/Rscript -e 'load("male.phototaxis.decline.qtt.RData"); write.table(cbind(gene.name, qtt.cor.young), file = "male.phototaxis.decline.qtt.young.score", sep = " ", row.names = F, col.names = F, quote = F); write.table(cbind(gene.name, qtt.cor.old), file = "male.phototaxis.decline.qtt.old.score", sep = " ", row.names = F, col.names = F, quote = F); write.table(cbind(gene.name, qtt.cor.diff), file = "male.phototaxis.decline.qtt.diff.score", sep = " ", row.names = F, col.names = F, quote = F);'
~/qgg/software/R-3.6.0/bin/Rscript -e 'load("male.speed.decline.qtt.RData"); write.table(cbind(gene.name, qtt.cor.young), file = "male.speed.decline.qtt.young.score", sep = " ", row.names = F, col.names = F, quote = F); write.table(cbind(gene.name, qtt.cor.old), file = "male.speed.decline.qtt.old.score", sep = " ", row.names = F, col.names = F, quote = F); write.table(cbind(gene.name, qtt.cor.diff), file = "male.speed.decline.qtt.diff.score", sep = " ", row.names = F, col.names = F, quote = F);'
~/qgg/software/R-3.6.0/bin/Rscript -e 'load("male.endurance.decline.qtt.RData"); write.table(cbind(gene.name, qtt.cor.young), file = "male.endurance.decline.qtt.young.score", sep = " ", row.names = F, col.names = F, quote = F); write.table(cbind(gene.name, qtt.cor.old), file = "male.endurance.decline.qtt.old.score", sep = " ", row.names = F, col.names = F, quote = F); write.table(cbind(gene.name, qtt.cor.diff), file = "male.endurance.decline.qtt.diff.score", sep = " ", row.names = F, col.names = F, quote = F);'
~/qgg/software/R-3.6.0/bin/Rscript -e 'load("male.lifespan.qtt.RData"); write.table(cbind(gene.name, qtt.cor.young), file = "male.lifespan.qtt.young.score", sep = " ", row.names = F, col.names = F, quote = F); write.table(cbind(gene.name, qtt.cor.old), file = "male.lifespan.qtt.old.score", sep = " ", row.names = F, col.names = F, quote = F); write.table(cbind(gene.name, qtt.cor.diff), file = "male.lifespan.qtt.diff.score", sep = " ", row.names = F, col.names = F, quote = F);'

# 2. gsea
# ============================================================

# all of old, young, diff, and all traits have the same GSEA sets
~/qgg/software/R-3.6.0/bin/Rscript ../R/prepGSEA.R female.fecundity.decline.qtt.young.score /mnt/research/qgg/resource/flybase/r6.36/annot/gene_go.table 20 female.qtt.gsea.data.RData > female.qtt.prepGSEA.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/prepKEGG.R female.fecundity.decline.qtt.young.score /mnt/research/qgg/resource/flybase/r6.36/annot/fbgn.kegg 20 female.qtt.kegg.data.RData > female.qtt.prepKEGG.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript ../R/prepGSEA.R male.phototaxis.decline.qtt.young.score /mnt/research/qgg/resource/flybase/r6.36/annot/gene_go.table 20 male.qtt.gsea.data.RData > male.qtt.prepGSEA.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/prepKEGG.R male.phototaxis.decline.qtt.young.score /mnt/research/qgg/resource/flybase/r6.36/annot/fbgn.kegg 20 male.qtt.kegg.data.RData > male.qtt.prepKEGG.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R female.qtt.gsea.data.RData female.fecundity.decline.qtt.young.score 1000 female.fecundity.decline.qtt.young.gsea.RData > female.fecundity.decline.qtt.young.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R female.qtt.gsea.data.RData female.fecundity.decline.qtt.old.score 1000 female.fecundity.decline.qtt.old.gsea.RData > female.fecundity.decline.qtt.old.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R female.qtt.gsea.data.RData female.fecundity.decline.qtt.diff.score 1000 female.fecundity.decline.qtt.diff.gsea.RData > female.fecundity.decline.qtt.diff.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R female.qtt.kegg.data.RData female.fecundity.decline.qtt.young.score 1000 female.fecundity.decline.qtt.young.kegg.RData > female.fecundity.decline.qtt.young.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R female.qtt.kegg.data.RData female.fecundity.decline.qtt.old.score 1000 female.fecundity.decline.qtt.old.kegg.RData > female.fecundity.decline.qtt.old.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R female.qtt.kegg.data.RData female.fecundity.decline.qtt.diff.score 1000 female.fecundity.decline.qtt.diff.kegg.RData > female.fecundity.decline.qtt.diff.kegg.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R female.qtt.gsea.data.RData female.phototaxis.decline.qtt.young.score 1000 female.phototaxis.decline.qtt.young.gsea.RData > female.phototaxis.decline.qtt.young.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R female.qtt.gsea.data.RData female.phototaxis.decline.qtt.old.score 1000 female.phototaxis.decline.qtt.old.gsea.RData > female.phototaxis.decline.qtt.old.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R female.qtt.gsea.data.RData female.phototaxis.decline.qtt.diff.score 1000 female.phototaxis.decline.qtt.diff.gsea.RData > female.phototaxis.decline.qtt.diff.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R female.qtt.kegg.data.RData female.phototaxis.decline.qtt.young.score 1000 female.phototaxis.decline.qtt.young.kegg.RData > female.phototaxis.decline.qtt.young.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R female.qtt.kegg.data.RData female.phototaxis.decline.qtt.old.score 1000 female.phototaxis.decline.qtt.old.kegg.RData > female.phototaxis.decline.qtt.old.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R female.qtt.kegg.data.RData female.phototaxis.decline.qtt.diff.score 1000 female.phototaxis.decline.qtt.diff.kegg.RData > female.phototaxis.decline.qtt.diff.kegg.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R female.qtt.gsea.data.RData female.lifespan.qtt.young.score 1000 female.lifespan.qtt.young.gsea.RData > female.lifespan.qtt.young.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R female.qtt.gsea.data.RData female.lifespan.qtt.old.score 1000 female.lifespan.qtt.old.gsea.RData > female.lifespan.qtt.old.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R female.qtt.gsea.data.RData female.lifespan.qtt.diff.score 1000 female.lifespan.qtt.diff.gsea.RData > female.lifespan.qtt.diff.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R female.qtt.kegg.data.RData female.lifespan.qtt.young.score 1000 female.lifespan.qtt.young.kegg.RData > female.lifespan.qtt.young.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R female.qtt.kegg.data.RData female.lifespan.qtt.old.score 1000 female.lifespan.qtt.old.kegg.RData > female.lifespan.qtt.old.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R female.qtt.kegg.data.RData female.lifespan.qtt.diff.score 1000 female.lifespan.qtt.diff.kegg.RData > female.lifespan.qtt.diff.kegg.Rout 2>&1 &

# males
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R male.qtt.gsea.data.RData male.phototaxis.decline.qtt.young.score 1000 male.phototaxis.decline.qtt.young.gsea.RData > male.phototaxis.decline.qtt.young.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R male.qtt.gsea.data.RData male.phototaxis.decline.qtt.old.score 1000 male.phototaxis.decline.qtt.old.gsea.RData > male.phototaxis.decline.qtt.old.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R male.qtt.gsea.data.RData male.phototaxis.decline.qtt.diff.score 1000 male.phototaxis.decline.qtt.diff.gsea.RData > male.phototaxis.decline.qtt.diff.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R male.qtt.kegg.data.RData male.phototaxis.decline.qtt.young.score 1000 male.phototaxis.decline.qtt.young.kegg.RData > male.phototaxis.decline.qtt.young.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R male.qtt.kegg.data.RData male.phototaxis.decline.qtt.old.score 1000 male.phototaxis.decline.qtt.old.kegg.RData > male.phototaxis.decline.qtt.old.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R male.qtt.kegg.data.RData male.phototaxis.decline.qtt.diff.score 1000 male.phototaxis.decline.qtt.diff.kegg.RData > male.phototaxis.decline.qtt.diff.kegg.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R male.qtt.gsea.data.RData male.speed.decline.qtt.young.score 1000 male.speed.decline.qtt.young.gsea.RData > male.speed.decline.qtt.young.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R male.qtt.gsea.data.RData male.speed.decline.qtt.old.score 1000 male.speed.decline.qtt.old.gsea.RData > male.speed.decline.qtt.old.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R male.qtt.gsea.data.RData male.speed.decline.qtt.diff.score 1000 male.speed.decline.qtt.diff.gsea.RData > male.speed.decline.qtt.diff.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R male.qtt.kegg.data.RData male.speed.decline.qtt.young.score 1000 male.speed.decline.qtt.young.kegg.RData > male.speed.decline.qtt.young.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R male.qtt.kegg.data.RData male.speed.decline.qtt.old.score 1000 male.speed.decline.qtt.old.kegg.RData > male.speed.decline.qtt.old.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R male.qtt.kegg.data.RData male.speed.decline.qtt.diff.score 1000 male.speed.decline.qtt.diff.kegg.RData > male.speed.decline.qtt.diff.kegg.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R male.qtt.gsea.data.RData male.endurance.decline.qtt.young.score 1000 male.endurance.decline.qtt.young.gsea.RData > male.endurance.decline.qtt.young.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R male.qtt.gsea.data.RData male.endurance.decline.qtt.old.score 1000 male.endurance.decline.qtt.old.gsea.RData > male.endurance.decline.qtt.old.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R male.qtt.gsea.data.RData male.endurance.decline.qtt.diff.score 1000 male.endurance.decline.qtt.diff.gsea.RData > male.endurance.decline.qtt.diff.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R male.qtt.kegg.data.RData male.endurance.decline.qtt.young.score 1000 male.endurance.decline.qtt.young.kegg.RData > male.endurance.decline.qtt.young.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R male.qtt.kegg.data.RData male.endurance.decline.qtt.old.score 1000 male.endurance.decline.qtt.old.kegg.RData > male.endurance.decline.qtt.old.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R male.qtt.kegg.data.RData male.endurance.decline.qtt.diff.score 1000 male.endurance.decline.qtt.diff.kegg.RData > male.endurance.decline.qtt.diff.kegg.Rout 2>&1 &


~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R male.qtt.gsea.data.RData male.lifespan.qtt.young.score 1000 male.lifespan.qtt.young.gsea.RData > male.lifespan.qtt.young.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R male.qtt.gsea.data.RData male.lifespan.qtt.old.score 1000 male.lifespan.qtt.old.gsea.RData > male.lifespan.qtt.old.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R male.qtt.gsea.data.RData male.lifespan.qtt.diff.score 1000 male.lifespan.qtt.diff.gsea.RData > male.lifespan.qtt.diff.gsea.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R male.qtt.kegg.data.RData male.lifespan.qtt.young.score 1000 male.lifespan.qtt.young.kegg.RData > male.lifespan.qtt.young.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R male.qtt.kegg.data.RData male.lifespan.qtt.old.score 1000 male.lifespan.qtt.old.kegg.RData > male.lifespan.qtt.old.kegg.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaKEGG.R male.qtt.kegg.data.RData male.lifespan.qtt.diff.score 1000 male.lifespan.qtt.diff.kegg.RData > male.lifespan.qtt.diff.kegg.Rout 2>&1 &

# calculate FDR
# ============================================================

~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R female.fecundity.decline.qtt.young.gsea.RData female.fecundity.decline.qtt.young.gsea.fdr.RData > female.fecundity.decline.qtt.young.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R female.fecundity.decline.qtt.old.gsea.RData female.fecundity.decline.qtt.old.gsea.fdr.RData > female.fecundity.decline.qtt.old.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R female.fecundity.decline.qtt.diff.gsea.RData female.fecundity.decline.qtt.diff.gsea.fdr.RData > female.fecundity.decline.qtt.diff.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R female.fecundity.decline.qtt.young.kegg.RData female.fecundity.decline.qtt.young.kegg.fdr.RData > female.fecundity.decline.qtt.young.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R female.fecundity.decline.qtt.old.kegg.RData female.fecundity.decline.qtt.old.kegg.fdr.RData > female.fecundity.decline.qtt.old.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R female.fecundity.decline.qtt.diff.kegg.RData female.fecundity.decline.qtt.diff.kegg.fdr.RData > female.fecundity.decline.qtt.diff.kegg.fdr.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R female.phototaxis.decline.qtt.young.gsea.RData female.phototaxis.decline.qtt.young.gsea.fdr.RData > female.phototaxis.decline.qtt.young.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R female.phototaxis.decline.qtt.old.gsea.RData female.phototaxis.decline.qtt.old.gsea.fdr.RData > female.phototaxis.decline.qtt.old.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R female.phototaxis.decline.qtt.diff.gsea.RData female.phototaxis.decline.qtt.diff.gsea.fdr.RData > female.phototaxis.decline.qtt.diff.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R female.phototaxis.decline.qtt.young.kegg.RData female.phototaxis.decline.qtt.young.kegg.fdr.RData > female.phototaxis.decline.qtt.young.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R female.phototaxis.decline.qtt.old.kegg.RData female.phototaxis.decline.qtt.old.kegg.fdr.RData > female.phototaxis.decline.qtt.old.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R female.phototaxis.decline.qtt.diff.kegg.RData female.phototaxis.decline.qtt.diff.kegg.fdr.RData > female.phototaxis.decline.qtt.diff.kegg.fdr.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R female.lifespan.qtt.young.gsea.RData female.lifespan.qtt.young.gsea.fdr.RData > female.lifespan.qtt.young.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R female.lifespan.qtt.old.gsea.RData female.lifespan.qtt.old.gsea.fdr.RData > female.lifespan.qtt.old.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R female.lifespan.qtt.diff.gsea.RData female.lifespan.qtt.diff.gsea.fdr.RData > female.lifespan.qtt.diff.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R female.lifespan.qtt.young.kegg.RData female.lifespan.qtt.young.kegg.fdr.RData > female.lifespan.qtt.young.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R female.lifespan.qtt.old.kegg.RData female.lifespan.qtt.old.kegg.fdr.RData > female.lifespan.qtt.old.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R female.lifespan.qtt.diff.kegg.RData female.lifespan.qtt.diff.kegg.fdr.RData > female.lifespan.qtt.diff.kegg.fdr.Rout 2>&1 &

# males
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R male.phototaxis.decline.qtt.young.gsea.RData male.phototaxis.decline.qtt.young.gsea.fdr.RData > male.phototaxis.decline.qtt.young.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R male.phototaxis.decline.qtt.old.gsea.RData male.phototaxis.decline.qtt.old.gsea.fdr.RData > male.phototaxis.decline.qtt.old.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R male.phototaxis.decline.qtt.diff.gsea.RData male.phototaxis.decline.qtt.diff.gsea.fdr.RData > male.phototaxis.decline.qtt.diff.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R male.phototaxis.decline.qtt.young.kegg.RData male.phototaxis.decline.qtt.young.kegg.fdr.RData > male.phototaxis.decline.qtt.young.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R male.phototaxis.decline.qtt.old.kegg.RData male.phototaxis.decline.qtt.old.kegg.fdr.RData > male.phototaxis.decline.qtt.old.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R male.phototaxis.decline.qtt.diff.kegg.RData male.phototaxis.decline.qtt.diff.kegg.fdr.RData > male.phototaxis.decline.qtt.diff.kegg.fdr.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R male.speed.decline.qtt.young.gsea.RData male.speed.decline.qtt.young.gsea.fdr.RData > male.speed.decline.qtt.young.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R male.speed.decline.qtt.old.gsea.RData male.speed.decline.qtt.old.gsea.fdr.RData > male.speed.decline.qtt.old.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R male.speed.decline.qtt.diff.gsea.RData male.speed.decline.qtt.diff.gsea.fdr.RData > male.speed.decline.qtt.diff.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R male.speed.decline.qtt.young.kegg.RData male.speed.decline.qtt.young.kegg.fdr.RData > male.speed.decline.qtt.young.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R male.speed.decline.qtt.old.kegg.RData male.speed.decline.qtt.old.kegg.fdr.RData > male.speed.decline.qtt.old.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R male.speed.decline.qtt.diff.kegg.RData male.speed.decline.qtt.diff.kegg.fdr.RData > male.speed.decline.qtt.diff.kegg.fdr.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R male.endurance.decline.qtt.young.gsea.RData male.endurance.decline.qtt.young.gsea.fdr.RData > male.endurance.decline.qtt.young.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R male.endurance.decline.qtt.old.gsea.RData male.endurance.decline.qtt.old.gsea.fdr.RData > male.endurance.decline.qtt.old.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R male.endurance.decline.qtt.diff.gsea.RData male.endurance.decline.qtt.diff.gsea.fdr.RData > male.endurance.decline.qtt.diff.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R male.endurance.decline.qtt.young.kegg.RData male.endurance.decline.qtt.young.kegg.fdr.RData > male.endurance.decline.qtt.young.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R male.endurance.decline.qtt.old.kegg.RData male.endurance.decline.qtt.old.kegg.fdr.RData > male.endurance.decline.qtt.old.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R male.endurance.decline.qtt.diff.kegg.RData male.endurance.decline.qtt.diff.kegg.fdr.RData > male.endurance.decline.qtt.diff.kegg.fdr.Rout 2>&1 &

~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R male.lifespan.qtt.young.gsea.RData male.lifespan.qtt.young.gsea.fdr.RData > male.lifespan.qtt.young.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R male.lifespan.qtt.old.gsea.RData male.lifespan.qtt.old.gsea.fdr.RData > male.lifespan.qtt.old.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDR.R male.lifespan.qtt.diff.gsea.RData male.lifespan.qtt.diff.gsea.fdr.RData > male.lifespan.qtt.diff.gsea.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R male.lifespan.qtt.young.kegg.RData male.lifespan.qtt.young.kegg.fdr.RData > male.lifespan.qtt.young.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R male.lifespan.qtt.old.kegg.RData male.lifespan.qtt.old.kegg.fdr.RData > male.lifespan.qtt.old.kegg.fdr.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaFDRKEGG.R male.lifespan.qtt.diff.kegg.RData male.lifespan.qtt.diff.kegg.fdr.RData > male.lifespan.qtt.diff.kegg.fdr.Rout 2>&1 &

# # 3. tx mixed model
# # ============================================================
#
# module load GCC/11.2.0 OpenMPI/4.1.1 R/4.1.2 && Rscript ../R/txMixed.R ../qg/female.sig.gene.exp.RData female.fecundity.decline.pheno ../eqtl/dgrp.common.fam female.fecundity.decline.qtt.mixed.RData > female.fecundity.decline.qtt.mixed.Rout 2>&1 &
# module load GCC/11.2.0 OpenMPI/4.1.1 R/4.1.2 && Rscript ../R/txMixed.R ../qg/female.sig.gene.exp.RData female.phototaxis.decline.pheno ../eqtl/dgrp.common.fam female.phototaxis.decline.qtt.mixed.RData > female.phototaxis.decline.qtt.mixed.Rout 2>&1 &
# module load GCC/11.2.0 OpenMPI/4.1.1 R/4.1.2 && Rscript ../R/txMixed.R ../qg/female.sig.gene.exp.RData female.lifespan.pheno ../eqtl/dgrp.common.fam female.lifespan.qtt.mixed.RData > female.lifespan.qtt.mixed.Rout 2>&1 &
#
# module load GCC/11.2.0 OpenMPI/4.1.1 R/4.1.2 && Rscript ../R/txMixed.R ../qg/female.sig.gene.exp.RData female.lifespan.raw.pheno ../eqtl/dgrp.common.fam female.lifespan.raw.qtt.mixed.RData > female.lifespan.raw.qtt.mixed.Rout 2>&1 &
#
#
#
# module load GCC/11.2.0 OpenMPI/4.1.1 R/4.1.2 && Rscript ../R/txMixed.R ../qg/male.sig.gene.exp.RData male.phototaxis.decline.pheno ../eqtl/dgrp.common.fam male.phototaxis.decline.qtt.mixed.RData > male.phototaxis.decline.qtt.mixed.Rout 2>&1 &
# module load GCC/11.2.0 OpenMPI/4.1.1 R/4.1.2 && Rscript ../R/txMixed.R ../qg/male.sig.gene.exp.RData male.speed.decline.pheno ../eqtl/dgrp.common.fam male.speed.decline.qtt.mixed.RData > male.speed.decline.qtt.mixed.Rout 2>&1 &
# module load GCC/11.2.0 OpenMPI/4.1.1 R/4.1.2 && Rscript ../R/txMixed.R ../qg/male.sig.gene.exp.RData male.endurance.decline.pheno ../eqtl/dgrp.common.fam male.endurance.decline.qtt.mixed.RData > male.endurance.decline.qtt.mixed.Rout 2>&1 &
# module load GCC/11.2.0 OpenMPI/4.1.1 R/4.1.2 && Rscript ../R/txMixed.R ../qg/male.sig.gene.exp.RData male.lifespan.pheno ../eqtl/dgrp.common.fam male.lifespan.qtt.mixed.RData > male.lifespan.qtt.mixed.Rout 2>&1 &
#

# 4. eigen gene association
# ============================================================

for trait in phototaxis.decline fecundity.decline lifespan
do
  module load R && Rscript ../R/eigenGeneAssoc.R ../qg/female.sig.gene.exp.RData ../figureData/female.wgcna.RData ../eqtl/dgrp.common.fam female.$trait.pheno female.$trait.qtt.eigen.assoc.RData > female.$trait.qtt.eigen.assoc.Rout 2>&1 &
done

for trait in phototaxis.decline speed.decline endurance.decline lifespan
do
  module load R && Rscript ../R/eigenGeneAssoc.R ../qg/male.sig.gene.exp.RData ../figureData/male.wgcna.RData ../eqtl/dgrp.common.fam male.$trait.pheno male.$trait.qtt.eigen.assoc.RData > male.$trait.qtt.eigen.assoc.Rout 2>&1 &
done

# 5. gsea example
# ============================================================
~/qgg/software/R-3.6.0/bin/Rscript ../R/gsea.R female.qtt.gsea.data.RData female.lifespan.qtt.young.score 1000 female.lifespan.qtt.young.gsea.RData > female.lifespan.qtt.young.gsea.Rout 2>&1 &
"gsea/female.age.effect.kegg.data.RData", "dme04213,dme00020", "figureData/female.age.effect.gsea.example.RData")
  
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaExample.R female.qtt.gsea.data.RData female.lifespan.qtt.young.score "BP>GO:0006099" female.qtt.kegg.data.RData "dme00020" female.lifespan.qtt.young.gsea.example.RData > female.lifespan.qtt.young.gsea.example.Rout 2>&1 &
~/qgg/software/R-3.6.0/bin/Rscript ../R/gseaExample.R female.qtt.gsea.data.RData female.lifespan.qtt.old.score "BP>GO:0006099" female.qtt.kegg.data.RData "dme00020" female.lifespan.qtt.old.gsea.example.RData > female.lifespan.qtt.old.gsea.example.Rout 2>&1 &

