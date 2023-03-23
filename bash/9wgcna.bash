# =====================================
# = co-expression network using wgcna =
# =====================================

# 1. get expression data
# ============================================================

Rscript -e 'load("figureData/female.adj.qgSingleTemp.RData"); sig.18c.gene <- which(p.adjust(single.temp[, 9], method = "BH") < 0.05); sig.25c.gene <- which(p.adjust(single.temp[, 12], method = "BH") < 0.05); sig.gene <- sort(unique(c(sig.18c.gene, sig.25c.gene))); geno.line.order <- read.table("eqtl/dgrp.common.fam", as.is = TRUE)[, 1]; load("qg/female.adj.qgBLUP.RData"); blup.18c.data <- blup.18c[sig.gene, match(geno.line.order, line.order)]; blup.25c.data <- blup.25c[sig.gene, match(geno.line.order, line.order)]; sig.18c.gene <- gene.name[sig.18c.gene]; sig.25c.gene <- gene.name[sig.25c.gene]; sig.gene.name <- gene.name[sig.gene]; save(blup.18c.data, blup.25c.data, sig.18c.gene, sig.25c.gene, sig.gene.name, file = "qg/female.sig.gene.exp.RData")' > log/female.sig.gene.exp.Rout 2>&1 &

Rscript -e 'load("figureData/male.adj.qgSingleTemp.RData"); sig.18c.gene <- which(p.adjust(single.temp[, 9], method = "BH") < 0.05); sig.25c.gene <- which(p.adjust(single.temp[, 12], method = "BH") < 0.05); sig.gene <- sort(unique(c(sig.18c.gene, sig.25c.gene))); geno.line.order <- read.table("eqtl/dgrp.common.fam", as.is = TRUE)[, 1]; load("qg/male.adj.qgBLUP.RData"); blup.18c.data <- blup.18c[sig.gene, match(geno.line.order, line.order)]; blup.25c.data <- blup.25c[sig.gene, match(geno.line.order, line.order)]; sig.18c.gene <- gene.name[sig.18c.gene]; sig.25c.gene <- gene.name[sig.25c.gene]; sig.gene.name <- gene.name[sig.gene]; save(blup.18c.data, blup.25c.data, sig.18c.gene, sig.25c.gene, sig.gene.name, file = "qg/male.sig.gene.exp.RData")' > log/male.sig.gene.exp.Rout 2>&1 &

# 2. simulation
# ============================================================

sbatch sbatch/femaleSimCorr18CJobArray.sbatch
sbatch sbatch/maleSimCorr18CJobArray.sbatch

# 3. summarize simulation results
# ============================================================

Rscript R/summarizeSimCorr18C.R qg/female.sig.gene.exp.RData figureData/female.adj.qgGxE.RData 10 simCorr/female.sim.corr.18c.perm 1000 simCorr/female.simCorr18C.RData > log/female.simCorr18C.Rout 2>&1 &
Rscript R/summarizeSimCorr18C.R qg/male.sig.gene.exp.RData figureData/male.adj.qgGxE.RData 10 simCorr/male.sim.corr.18c.perm 1000 simCorr/male.simCorr18C.RData > log/male.simCorr18C.Rout 2>&1 &

# 4. network analysis
# ============================================================

module load GCC/11.2.0 OpenMPI/4.1.1 R/4.1.2
Rscript R/runWGCNA.R qg/female.sig.gene.exp.RData figureData/female.wgcna.RData > log/female.wgcnaOnly.Rout 2>&1 &
Rscript R/runWGCNA.R qg/male.sig.gene.exp.RData figureData/male.wgcna.RData > log/male.wgcnaOnly.Rout 2>&1 &


sbatch --output log/female.runNetwork.out --error log/female.runNetwork.err sbatch/female.runNetwork.sbatch
sbatch --output log/male.runNetwork.out --error log/male.runNetwork.err sbatch/male.runNetwork.sbatch
