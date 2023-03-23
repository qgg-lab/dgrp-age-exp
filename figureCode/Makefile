# ===============================================
# = make file for generating tables and figures =
# ===============================================

all: ../figure/figureScaledExpOutlier.pdf ../figure/figureAgingSummary.pdf ../figure/tableAgeEffGSEA.xlsx ../figure/figureGxE.pdf ../figure/tableGxAGSEA.xlsx ../figure/tableMapping.xlsx ../figure/figureVarHetSigmaE.pdf ../figure/figureGxESexDiff.pdf ../figure/tableDecanGSEA.xlsx ../figure/tableVarComp.xlsx ../figure/figureH2.pdf ../figure/figureModules.pdf ../figure/figureTCAanalysisFemale.pdf ../figure/figureTCAanalysisMale.pdf ../figure/tableeQTL.xlsx ../figure/figureSummaryeQTL.pdf ../figure/figureEffeQTL.pdf ../figure/tableMediation.xlsx ../figure/figureMediation.pdf ../figure/tableNetRep.xlsx ../figure/tableTxMixed.xlsx ../figure/figureQTTsingle.pdf ../figure/tableEigenGene.xlsx ../figure/tableQTTGSEA.xlsx ../figure/figureQTTGSEA.pdf ../figure/tableModules.xlsx ../figure/tableTCAgenes.xlsx ../figure/tableNetworkChangeFemale.xlsx ../figure/tableNetworkChangeMale.xlsx ../figure/tableTCAQTT.xlsx ../figure/figurePhenoCorr.pdf ../figure/tableQTTsingle.xlsx ../figure/tableQTTmediation.xlsx ../figure/figureQTTTCA.pdf ../figure/figureQTTMediation.pdf
	
FONT="Myriad Pro"

# figure for outlier detection
# ============================================================

../figure/figureScaledExpOutlier.pdf: figureScaledExpOutlier.R ../figureData/female.gene.exp.scaled.quantile.RData ../figureData/male.gene.exp.scaled.quantile.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureScaledExpOutlier.R ../figureData/female.gene.exp.scaled.quantile.RData ../figureData/male.gene.exp.scaled.quantile.RData ../figure/figureScaledExpOutlier.pdf $(FONT)

# table for mapping summary
#	cat ../figureData/biosample_result.txt | paste - - - - - - | perl -wne 'chomp $_; if ($_ =~ m/\d+\. DGRP Line (.*Rep \d).*BioSample: (SAMN.*?);/) { print $1, "\t", $2, "\n"; }' | sed 's/ /_/g' > ../figureData/baseline.biosample.txt
# ============================================================

../figure/tableMapping.xlsx: tableMapping.R ../figureData/baseline.biosample.txt ../figureData/BioSampleObjects.txt ../figureData/BioSampleObjects2.txt ../figureData/baseline.biosample.txt ../figureData/bad.sample.id ../figureData/sample.map.summary
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableMapping.R ../figureData/BioSampleObjects.txt ../figureData/BioSampleObjects2.txt ../figureData/baseline.biosample.txt ../figureData/bad.sample.id ../figureData/sample.map.summary ../figure/tableMapping.xlsx
# figure for 1st summary
# ============================================================

../figure/figureAgingSummary.pdf: figureAgingSummary.R ../figureData/female.adj.qgAgeEff.RData ../figureData/male.adj.qgAgeEff.RData ../figureData/gene.exp.mean.RData ../figureData/female.age.effect.gsea.example.RData ../figureData/male.age.effect.gsea.example.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureAgingSummary.R ../figureData/female.adj.qgAgeEff.RData ../figureData/male.adj.qgAgeEff.RData ../figureData/gene.exp.mean.RData ../figureData/female.age.effect.gsea.example.RData ../figureData/male.age.effect.gsea.example.RData ../figure/figureAgingSummary.pdf $(FONT)

# table for age effect GSEA
# ============================================================

../figure/tableAgeEffGSEA.xlsx: tableAgeEffGSEA.R ../figureData/female.age.effect.gsea.fdr.RData ../figureData/male.age.effect.gsea.fdr.RData ../figureData/female.age.effect.kegg.fdr.RData ../figureData/male.age.effect.kegg.fdr.RData ../figureData/kegg.id
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableAgeEffGSEA.R ../figureData/female.age.effect.gsea.fdr.RData ../figureData/male.age.effect.gsea.fdr.RData ../figureData/female.age.effect.kegg.fdr.RData ../figureData/male.age.effect.kegg.fdr.RData ../figureData/kegg.id ../figure/tableAgeEffGSEA.xlsx

# figure for gxe
# ============================================================

../figure/figureGxE.pdf: figureGxE.R ../figureData/female.adj.qgGxE.RData ../figureData/male.adj.qgGxE.RData ../figureData/female.adj.qgVarHet.RData ../figureData/female.adj.qgSingleTemp.RData ../figureData/male.adj.qgVarHet.RData ../figureData/male.adj.qgSingleTemp.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureGxE.R ../figureData/female.adj.qgGxE.RData ../figureData/male.adj.qgGxE.RData ../figureData/female.adj.qgVarHet.RData ../figureData/female.adj.qgSingleTemp.RData ../figureData/male.adj.qgVarHet.RData ../figureData/male.adj.qgSingleTemp.RData ../figure/figureGxE.pdf $(FONT)


# figure for H2
# ============================================================

../figure/figureH2.pdf: figureH2.R ../figureData/female.adj.qgSingleTemp.RData ../figureData/male.adj.qgSingleTemp.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureH2.R ../figureData/female.adj.qgSingleTemp.RData ../figureData/male.adj.qgSingleTemp.RData ../figure/figureH2.pdf $(FONT)


# figure for environmental variance 
# ============================================================

../figure/figureVarHetSigmaE.pdf: figureVarHetSigmaE.R ../figureData/female.adj.qgVarHet.RData ../figureData/male.adj.qgVarHet.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureVarHetSigmaE.R ../figureData/female.adj.qgVarHet.RData ../figureData/male.adj.qgVarHet.RData ../figure/figureVarHetSigmaE.pdf $(FONT)

# figure for GxE in two sexes
# ============================================================

../figure/figureGxESexDiff.pdf: figureGxESexDiff.R ../figureData/female.adj.qgGxE.RData ../figureData/male.adj.qgGxE.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureGxESexDiff.R ../figureData/female.adj.qgGxE.RData ../figureData/male.adj.qgGxE.RData ../figure/figureGxESexDiff.pdf $(FONT)

# table for variance components
# ============================================================

../figure/tableVarComp.xlsx: tableVarComponents.R ../figureData/female.adj.qgSingleTemp.RData ../figureData/male.adj.qgSingleTemp.RData ../figureData/female.adj.qgVarHet.RData ../figureData/male.adj.qgVarHet.RData ../figureData/female.adj.qgGxE.RData ../figureData/male.adj.qgGxE.RData ../figureData/gene.info ../figureData/female.adj.qgAgeEff.RData ../figureData/male.adj.qgAgeEff.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableVarComponents.R ../figureData/female.adj.qgSingleTemp.RData ../figureData/male.adj.qgSingleTemp.RData ../figureData/female.adj.qgVarHet.RData ../figureData/male.adj.qgVarHet.RData ../figureData/female.adj.qgGxE.RData ../figureData/male.adj.qgGxE.RData ../figureData/gene.info ../figureData/female.adj.qgAgeEff.RData ../figureData/male.adj.qgAgeEff.RData ../figure/tableVarComp.xlsx

# table for GxA GSEA
# ============================================================

../figure/tableGxAGSEA.xlsx: tableGxAGSEA.R ../figureData/female.gei.gsea.fdr.RData ../figureData/male.gei.gsea.fdr.RData ../figureData/female.gei.kegg.fdr.RData ../figureData/male.gei.kegg.fdr.RData ../figureData/kegg.id
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableGxAGSEA.R ../figureData/female.gei.gsea.fdr.RData  ../figureData/male.gei.gsea.fdr.RData ../figureData/female.gei.kegg.fdr.RData ../figureData/male.gei.kegg.fdr.RData ../figureData/kegg.id ../figure/tableGxAGSEA.xlsx

# table for decanalization GSEA
# ============================================================

../figure/tableDecanGSEA.xlsx: tableDecanGSEA.R ../figureData/female.decan.gsea.fdr.RData ../figureData/male.decan.gsea.fdr.RData ../figureData/female.decan.kegg.fdr.RData ../figureData/male.decan.kegg.fdr.RData ../figureData/kegg.id
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableDecanGSEA.R ../figureData/female.decan.gsea.fdr.RData ../figureData/male.decan.gsea.fdr.RData ../figureData/female.decan.kegg.fdr.RData ../figureData/male.decan.kegg.fdr.RData ../figureData/kegg.id ../figure/tableDecanGSEA.xlsx

# figure for modules
# ============================================================

../figure/figureModules.pdf: figureModules.R ../figureData/female.wgcna.RData  ../figureData/male.wgcna.RData 
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureModules.R ../figureData/female.wgcna.RData ../figureData/male.wgcna.RData  ../figure/figureModules.pdf $(FONT)


# table for modules
# ============================================================

../figure/tableModules.xlsx: tableModules.R ../figureData/female.wgcna.RData ../figureData/male.wgcna.RData ../figureData/gene.info
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableModules.R ../figureData/female.wgcna.RData ../figureData/male.wgcna.RData ../figureData/gene.info ../figure/tableModules.xlsx

# table for TCA genes
# ============================================================

../figure/tableTCAgenes.xlsx: tableTCAgenes.R ../figureData/female.wgcna.RData ../figureData/male.wgcna.RData ../figureData/gene_go.table ../figureData/fbgn.kegg ../figureData/gene.info
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableTCAgenes.R ../figureData/female.wgcna.RData ../figureData/male.wgcna.RData ../figureData/gene_go.table ../figureData/fbgn.kegg ../figureData/gene.info ../figure/tableTCAgenes.xlsx


# figure for TCA
# ============================================================

../figure/figureTCAanalysisFemale.pdf: figureTCAanalysisFemale.R ../figureData/female.wgcna.RData  ../figureData/gene_go.table ../figureData/fbgn.kegg
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureTCAanalysisFemale.R ../figureData/female.wgcna.RData ../figureData/gene_go.table ../figureData/fbgn.kegg ../figureData/gene.info ../figure/figureTCAanalysisFemale.pdf $(FONT)

# table for TCA difference
# ============================================================

../figure/tableNetworkChangeFemale.xlsx: tableNetworkChangeAnalysis.R ../figureData/female.wgcna.RData ../figureData/gene_go.table ../figureData/fbgn.kegg ../figureData/gene.info ../figureData/kegg.id
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableNetworkChangeAnalysis.R ../figureData/female.wgcna.RData ../figureData/gene_go.table ../figureData/fbgn.kegg ../figureData/gene.info ../figureData/kegg.id ../figure/tableNetworkChangeFemale.xlsx
	
../figure/tableNetworkChangeMale.xlsx: tableNetworkChangeAnalysis.R ../figureData/male.wgcna.RData ../figureData/gene_go.table ../figureData/fbgn.kegg ../figureData/gene.info ../figureData/kegg.id
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableNetworkChangeAnalysis.R ../figureData/male.wgcna.RData ../figureData/gene_go.table ../figureData/fbgn.kegg ../figureData/gene.info ../figureData/kegg.id ../figure/tableNetworkChangeMale.xlsx
	
# eqtl table
# ============================================================

../figure/tableeQTL.xlsx:	tableModelSelectTable.R ../figureData/female.18c.fdr05.eqtl.table.txt ../figureData/female.25c.fdr05.eqtl.table.txt ../figureData/male.18c.fdr05.eqtl.table.txt ../figureData/male.25c.fdr05.eqtl.table.txt ../figureData/gene.info ../figureData/female.18c.genvar.id ../figureData/female.25c.genvar.id ../figureData/male.18c.genvar.id ../figureData/male.25c.genvar.id
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableModelSelectTable.R ../figureData/female.18c.fdr05.eqtl.table.txt ../figureData/female.25c.fdr05.eqtl.table.txt ../figureData/male.18c.fdr05.eqtl.table.txt ../figureData/male.25c.fdr05.eqtl.table.txt ../figureData/gene.info ../figureData/female.18c.genvar.id ../figureData/female.25c.genvar.id ../figureData/male.18c.genvar.id ../figureData/male.25c.genvar.id ../figureData/female.adj.qgGxE.RData ../figureData/male.adj.qgGxE.RData ../figure/tableeQTL.xlsx

# eqtl figure
# ============================================================

../figure/figureSummaryeQTL.pdf: figureSummaryeQTL.R ../figureData/female.18c.fdr05.eqtl.table.txt ../figureData/female.25c.fdr05.eqtl.table.txt ../figureData/male.18c.fdr05.eqtl.table.txt ../figureData/male.25c.fdr05.eqtl.table.txt ../figureData/gene.info ../figureData/female.18c.genvar.id ../figureData/female.25c.genvar.id ../figureData/male.18c.genvar.id ../figureData/male.25c.genvar.id ../figureData/female.tf.summary.RData ../figureData/male.tf.summary.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureSummaryeQTL.R ../figureData/female.18c.fdr05.eqtl.table.txt ../figureData/female.25c.fdr05.eqtl.table.txt ../figureData/male.18c.fdr05.eqtl.table.txt ../figureData/male.25c.fdr05.eqtl.table.txt ../figureData/gene.info ../figureData/female.18c.genvar.id ../figureData/female.25c.genvar.id ../figureData/male.18c.genvar.id ../figureData/male.25c.genvar.id ../figureData/female.adj.qgGxE.RData ../figureData/male.adj.qgGxE.RData ../figureData/female.tf.summary.RData ../figureData/male.tf.summary.RData ../figure/figureSummaryeQTL.pdf $(FONT)
	

# eqtl eff
# ============================================================

../figure/figureEffeQTL.pdf: figureEffeQTL.R ../figureData/female.eqtl.eff.RData ../figureData/male.eqtl.eff.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureEffeQTL.R ../figureData/female.eqtl.eff.RData ../figureData/male.eqtl.eff.RData ../figure/figureEffeQTL.pdf  $(FONT)

# mediation table
# ============================================================

../figure/tableMediation.xlsx: tableMediation.R ../figureData/female.young.trio.mediation.RData ../figureData/female.old.trio.mediation.RData ../figureData/male.young.trio.mediation.RData ../figureData/male.old.trio.mediation.RData ../figureData/gene.info
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableMediation.R ../figureData/female.young.trio.mediation.RData ../figureData/female.old.trio.mediation.RData ../figureData/male.young.trio.mediation.RData ../figureData/male.old.trio.mediation.RData ../figureData/gene.info ../figure/tableMediation.xlsx


# mediation figure
# ============================================================

../figure/figureMediation.pdf: figureMediation.R ../figureData/female.example.mediation.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureMediation.R ../figureData/female.example.mediation.RData ../figure/figureMediation.pdf $(FONT)

# netrep
# ============================================================

../figure/tableNetRep.xlsx: tableNetRep.R ../figureData/female.wgcna.RData ../figureData/male.wgcna.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableNetRep.R ../figureData/female.wgcna.RData ../figureData/male.wgcna.RData ../figure/tableNetRep.xlsx
	
# tx mixed model
# ============================================================

../figure/tableTxMixed.xlsx: tableTxMixed.R ../figureData/female.qtt.mixed.RData ../figureData/male.qtt.mixed.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableTxMixed.R ../figureData/female.qtt.mixed.RData ../figureData/male.qtt.mixed.RData ../figure/tableTxMixed.xlsx


# QTT tca genes
# ============================================================

../figure/tableTCAQTT.xlsx: tableTCAQTT.R ../figureData/female.qtt.RData ../figureData/male.qtt.RData ../figureData/gene_go.table ../figureData/fbgn.kegg ../figureData/gene.info
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableTCAQTT.R ../figureData/female.qtt.RData ../figureData/male.qtt.RData ../figureData/gene_go.table ../figureData/fbgn.kegg ../figureData/gene.info ../figure/tableTCAQTT.xlsx
		
# QTT figure
# ============================================================

../figure/figureQTTsingle.pdf: figureQTTsingle.R ../figureData/female.phototaxis.decline.qtt.RData ../figureData/female.fecundity.decline.qtt.RData ../figureData/female.lifespan.qtt.RData ../figureData/male.phototaxis.decline.qtt.RData ../figureData/male.speed.decline.qtt.RData ../figureData/male.endurance.decline.qtt.RData ../figureData/male.lifespan.qtt.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureQTTsingle.R ../figureData/female.phototaxis.decline.qtt.RData ../figureData/female.fecundity.decline.qtt.RData ../figureData/female.lifespan.qtt.RData ../figureData/male.phototaxis.decline.qtt.RData ../figureData/male.speed.decline.qtt.RData ../figureData/male.endurance.decline.qtt.RData ../figureData/male.lifespan.qtt.RData ../figure/figureQTTsingle.pdf $(FONT)

# QTT table
# ============================================================

../figure/tableQTTsingle.xlsx: tableQTTsingle.R ../figureData/female.phototaxis.decline.qtt.RData ../figureData/female.fecundity.decline.qtt.RData ../figureData/female.lifespan.qtt.RData ../figureData/male.phototaxis.decline.qtt.RData ../figureData/male.speed.decline.qtt.RData ../figureData/male.endurance.decline.qtt.RData ../figureData/male.lifespan.qtt.RData ../figureData/gene.info
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableQTTsingle.R ../figureData/female.phototaxis.decline.qtt.RData ../figureData/female.fecundity.decline.qtt.RData ../figureData/female.lifespan.qtt.RData ../figureData/male.phototaxis.decline.qtt.RData ../figureData/male.speed.decline.qtt.RData ../figureData/male.endurance.decline.qtt.RData ../figureData/male.lifespan.qtt.RData ../figureData/gene.info ../figure/tableQTTsingle.xlsx
	
# eigengene
# ============================================================

../figure/tableEigenGene.xlsx: tableEigenGene.R ../figureData/female.phototaxis.decline.qtt.eigen.assoc.RData ../figureData/female.fecundity.decline.qtt.eigen.assoc.RData ../figureData/female.lifespan.qtt.eigen.assoc.RData ../figureData/male.phototaxis.decline.qtt.eigen.assoc.RData ../figureData/male.speed.decline.qtt.eigen.assoc.RData ../figureData/male.endurance.decline.qtt.eigen.assoc.RData ../figureData/male.lifespan.qtt.eigen.assoc.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableEigenGene.R ../figureData/female.phototaxis.decline.qtt.eigen.assoc.RData ../figureData/female.fecundity.decline.qtt.eigen.assoc.RData ../figureData/female.lifespan.qtt.eigen.assoc.RData ../figureData/male.phototaxis.decline.qtt.eigen.assoc.RData ../figureData/male.speed.decline.qtt.eigen.assoc.RData ../figureData/male.endurance.decline.qtt.eigen.assoc.RData ../figureData/male.lifespan.qtt.eigen.assoc.RData ../figure/tableEigenGene.xlsx
	
# QTT GSEA
# ============================================================

../figure/tableQTTGSEA.xlsx: tableQTTGSEA.R ../figureData/kegg.id
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableQTTGSEA.R phototaxis.decline,fecundity.decline,lifespan phototaxis.decline,speed.decline,endurance.decline,lifespan ../figureData/kegg.id ../figure/tableQTTGSEA.xlsx
	
	
# QTT mediation
# ============================================================

../figure/tableQTTmediation.xlsx: tableQTTmediation.R ../figureData/female.eqtl.gwas.trio.mediation.RData ../figureData/male.eqtl.gwas.trio.mediation.RData ../figureData/gene.info
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript tableQTTmediation.R ../figureData/female.eqtl.gwas.trio.mediation.RData ../figureData/male.eqtl.gwas.trio.mediation.RData ../figureData/gene.info ../figure/tableQTTmediation.xlsx

# QTT TCA figure
# ============================================================

../figure/figureQTTGSEA.pdf: figureQTTGSEA.R ../figureData/female.young.qtt.gsea.example.RData ../figureData/female.old.qtt.gsea.example.RData ../figureData/female.qtt.RData ../figureData/gene_go.table ../figureData/fbgn.kegg
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureQTTGSEA.R ../figureData/female.young.qtt.gsea.example.RData ../figureData/female.old.qtt.gsea.example.RData ../figureData/female.qtt.RData ../figureData/gene_go.table ../figureData/fbgn.kegg ../figure/figureQTTGSEA.pdf $(FONT)
	
# QTT TCA figure
# ============================================================

../figure/figureQTTTCA.pdf: figureQTTTCA.R ../figureData/female.phototaxis.decline.qtt.RData ../figureData/female.fecundity.decline.qtt.RData ../figureData/female.lifespan.qtt.RData ../figureData/male.phototaxis.decline.qtt.RData ../figureData/male.speed.decline.qtt.RData ../figureData/male.endurance.decline.qtt.RData ../figureData/male.lifespan.qtt.RData ../figureData/tca.genes.txt ../figureData/female.lifespan.qtt.young.gsea.example.RData ../figureData/female.lifespan.qtt.old.gsea.example.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureQTTTCA.R ../figureData/female.phototaxis.decline.qtt.RData ../figureData/female.fecundity.decline.qtt.RData ../figureData/female.lifespan.qtt.RData ../figureData/male.phototaxis.decline.qtt.RData ../figureData/male.speed.decline.qtt.RData ../figureData/male.endurance.decline.qtt.RData ../figureData/male.lifespan.qtt.RData ../figureData/tca.genes.txt ../figureData/female.lifespan.qtt.young.gsea.example.RData ../figureData/female.lifespan.qtt.old.gsea.example.RData ../figure/figureQTTTCA.pdf $(FONT)


# phenotypic correlation
# ============================================================

../figure/figurePhenoCorr.pdf: figurePhenoCorr.R ../figureData/pheno.merge.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figurePhenoCorr.R ../figureData/pheno.merge.RData ../figure/figurePhenoCorr.pdf $(FONT)


# QTT mediation figure
# ============================================================

../figure/figureQTTMediation.pdf: figureQTTMediation.R ../figureData/male.eqtl.gwas.trio.mediation.example.RData
	/Library/Frameworks/R.framework/Versions/Current/Resources/bin/Rscript figureQTTMediation.R ../figureData/male.eqtl.gwas.trio.mediation.example.RData ../figure/figureQTTMediation.pdf $(FONT)

	