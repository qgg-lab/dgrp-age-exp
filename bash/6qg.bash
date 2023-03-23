# 1. QG and blup
# ============================================================
cd /mnt/research/qgg/dgrp-age-exp
mkdir qg

module load R && srun --cpus-per-task=16 --ntasks-per-node=1 --mem=8G --time=4:00:00 Rscript R/qgAgeEff.R normalize/female.gene.exp.adjust.RData 16 figureData/female.adj.qgAgeEff.RData > log/female.adj.qgAgeEff.Rout 2>&1 &
module load R && srun --cpus-per-task=16 --ntasks-per-node=1 --mem=8G --time=4:00:00 Rscript R/qgAgeEff.R normalize/male.gene.exp.adjust.RData 16 figureData/male.adj.qgAgeEff.RData > log/male.adj.qgAgeEff.Rout 2>&1 &

module load R && srun --cpus-per-task=16 --ntasks-per-node=1 --mem=8G --time=8:00:00 Rscript R/qgRandomLineEff.R normalize/female.gene.exp.adjust.RData > log2/female.qgRandomLineEff.Rout 2>&1 &
module load R && srun --cpus-per-task=16 --ntasks-per-node=1 --mem=8G --time=8:00:00 Rscript R/qgRandomLineEff.R normalize/male.gene.exp.adjust.RData > log2/male.qgRandomLineEff.Rout 2>&1 &

Rscript R/qgGxE.R normalize/female.gene.exp.adjust.RData 16 figureData/female.adj.qgGxE.RData > log2/female.adj.qgGxE.Rout 2>&1 &
Rscript R/qgGxE.R normalize/male.gene.exp.adjust.RData 16 figureData/male.adj.qgGxE.RData > log2/male.adj.qgGxE.Rout 2>&1 &

module load R && srun --cpus-per-task=16 --ntasks-per-node=1 --mem=8G --time=4:00:00 Rscript R/qgSingleTemp.R normalize/female.gene.exp.adjust.RData 16 figureData/female.adj.qgSingleTemp.RData > log2/female.adj.qgSingleTemp.Rout 2>&1 &
module load R && srun --cpus-per-task=16 --ntasks-per-node=1 --mem=8G --time=4:00:00 Rscript R/qgSingleTemp.R normalize/male.gene.exp.adjust.RData 16 figureData/male.adj.qgSingleTemp.RData > log2/male.adj.qgSingleTemp.Rout 2>&1 &

module load R && srun --cpus-per-task=16 --ntasks-per-node=1 --mem=8G --time=4:00:00 Rscript R/qgVarHet.R normalize/female.gene.exp.adjust.RData 16 figureData/female.adj.qgVarHet.RData > log/female.adj.qgVarHet.Rout 2>&1 &
module load R && srun --cpus-per-task=16 --ntasks-per-node=1 --mem=8G --time=4:00:00 Rscript R/qgVarHet.R normalize/male.gene.exp.adjust.RData 16 figureData/male.adj.qgVarHet.RData > log/male.adj.qgVarHet.Rout 2>&1 &

Rscript R/qgBLUP.R normalize/female.gene.exp.adjust.RData /mnt/research/qgg/resource/dgrp/dgrp2.r6.adjustData.RData 16 qg/female.adj.qgBLUP.RData > log2/female.adj.qgBLUP.Rout 2>&1 &
Rscript R/qgBLUP.R normalize/male.gene.exp.adjust.RData /mnt/research/qgg/resource/dgrp/dgrp2.r6.adjustData.RData 16 qg/male.adj.qgBLUP.RData > log2/male.adj.qgBLUP.Rout 2>&1 &

# also get mean difference between the two ages
# ============================================================

~/qgg/software/R-3.6.0/bin/Rscript -e 'load("/mnt/research/qgg/dgrp-age-exp/normalize/female.gene.exp.adjust.RData"); female.old.mean <- rowMeans(gene.exp.adj[, temp == "18C"]); female.young.mean <- rowMeans(gene.exp.adj[, temp == "25C"]); load("/mnt/research/qgg/dgrp-age-exp/normalize/male.gene.exp.adjust.RData"); male.old.mean <- rowMeans(gene.exp.adj[, temp == "18C"]); male.young.mean <- rowMeans(gene.exp.adj[, temp == "25C"]); save(female.old.mean, female.young.mean, male.old.mean, male.young.mean, file = "/mnt/research/qgg/dgrp-age-exp/figureData/gene.exp.mean.RData")'

# examples of GxE small and large (find some Hsps)
# ============================================================

# ~/software/R-3.2.2/bin/Rscript R/qgGxEexample.R qg/female.adj.qgBLUP.RData qg/female.adj.qgGxE.RData cuff/gene.info qg/female.adj.qgGxEexample.RData > log/female.adj.qgGxEexample.Rout 2>&1 &
