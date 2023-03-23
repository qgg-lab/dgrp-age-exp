# =======================================================
# = analysis for quantitative genetics and eqtl mapping =
# =======================================================

cd /mnt/research/qgg/dgrp-age-exp/
mkdir eqtl

# get plink genotype data
# ============================================================

~/qgg/software/R-3.6.0/bin/Rscript -e 'load("normalize/female.gene.exp.adjust.RData"); female.line.id <- line.id; load("normalize/male.gene.exp.adjust.RData"); male.line.id <- line.id; union.id <- unique(c(as.character(female.line.id), as.character(male.line.id))); write.table(cbind(union.id, union.id), file = "eqtl/exp.line.id", sep = " ", quote = F, row.names = F, col.names = F)'

/mnt/research/qgg/software/plink-v1.90b6.18/plink --silent --bfile /mnt/research/qgg/resource/dgrp/dgrp2.biallelic.common --allow-extra-chr --keep eqtl/exp.line.id --geno 0.2 --maf 0.05 --make-bed --out eqtl/dgrp.common

# 1. find significant genes and map eQTL
# ============================================================

Rscript -e 'load("figureData/female.adj.qgSingleTemp.RData"); load("figureData/female.adj.qgGxE.RData"); sig.gene <- which(p.adjust(single.temp[, 9], method = "BH") < 0.05 | p.adjust(single.temp[, 12], method = "BH") < 0.05 | p.adjust(gxe[, 12], method = "BH") < 0.05 | p.adjust(gxe[, 13], method = "BH") < 0.05); sig.gene.name <- gene.name[sig.gene]; sig.gene.18c <- which(p.adjust(single.temp[, 9], method = "BH") < 0.05); sig.gene.25c <- which(p.adjust(single.temp[, 12], method = "BH") < 0.05); sig.gene.name.25c <- gene.name[sig.gene.25c]; sig.gene.name.18c <- gene.name[sig.gene.18c]; geno.line.order <- read.table("/mnt/research/qgg/dgrp-age-exp/eqtl/dgrp.common.fam", as.is = TRUE)[, 1]; load("qg/female.adj.qgBLUP.RData"); blup.18c.data <- blup.18c[sig.gene, match(geno.line.order, line.order)]; blup.25c.data <- blup.25c[sig.gene, match(geno.line.order, line.order)]; blup.diff.data <- blup.18c.data - blup.25c.data; save(blup.18c.data, blup.25c.data, blup.diff.data, sig.gene.name, sig.gene.name.25c, sig.gene.name.18c, file = "qg/female.sig.gene.exp.RData"); blup.25c.var <- apply(blup.25c.data, 1, var, na.rm=TRUE); blup.18c.var <- apply(blup.18c.data, 1, var, na.rm=TRUE); blup.diff.var <- apply(blup.diff.data, 1, var, na.rm=TRUE); write.table(rbind(c("FID", "IID", gene.name[sig.gene[blup.18c.var > 1e-8]]), cbind(geno.line.order, geno.line.order, t(blup.18c.data[blup.18c.var > 1e-8, ]))), file = "eqtl/female.18c.pheno", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE); write.table(rbind(c("FID", "IID", gene.name[sig.gene[blup.25c.var > 1e-8]]), cbind(geno.line.order, geno.line.order, t(blup.25c.data[blup.25c.var > 1e-8, ]))), file = "eqtl/female.25c.pheno", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE); write.table(rbind(c("FID", "IID", gene.name[sig.gene[blup.diff.var > 1e-8]]), cbind(geno.line.order, geno.line.order, t(blup.diff.data[blup.diff.var > 1e-8, ]))), file = "eqtl/female.diff.pheno", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE);' > log/female.prepPhenoeQTL.Rout 2>&1 &

Rscript -e 'load("figureData/male.adj.qgSingleTemp.RData"); load("figureData/male.adj.qgGxE.RData"); sig.gene <- which(p.adjust(single.temp[, 9], method = "BH") < 0.05 | p.adjust(single.temp[, 12], method = "BH") < 0.05 | p.adjust(gxe[, 12], method = "BH") < 0.05 | p.adjust(gxe[, 13], method = "BH") < 0.05); sig.gene.name <- gene.name[sig.gene]; sig.gene.18c <- which(p.adjust(single.temp[, 9], method = "BH") < 0.05); sig.gene.25c <- which(p.adjust(single.temp[, 12], method = "BH") < 0.05); sig.gene.name.25c <- gene.name[sig.gene.25c]; sig.gene.name.18c <- gene.name[sig.gene.18c]; geno.line.order <- read.table("/mnt/research/qgg/dgrp-age-exp/eqtl/dgrp.common.fam", as.is = TRUE)[, 1]; load("qg/male.adj.qgBLUP.RData"); blup.18c.data <- blup.18c[sig.gene, match(geno.line.order, line.order)]; blup.25c.data <- blup.25c[sig.gene, match(geno.line.order, line.order)]; blup.diff.data <- blup.18c.data - blup.25c.data; save(blup.18c.data, blup.25c.data, blup.diff.data, sig.gene.name, sig.gene.name.25c, sig.gene.name.18c, file = "qg/male.sig.gene.exp.RData"); blup.25c.var <- apply(blup.25c.data, 1, var, na.rm=TRUE); blup.18c.var <- apply(blup.18c.data, 1, var, na.rm=TRUE); blup.diff.var <- apply(blup.diff.data, 1, var, na.rm=TRUE); write.table(rbind(c("FID", "IID", gene.name[sig.gene[blup.18c.var > 1e-8]]), cbind(geno.line.order, geno.line.order, t(blup.18c.data[blup.18c.var > 1e-8, ]))), file = "eqtl/male.18c.pheno", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE); write.table(rbind(c("FID", "IID", gene.name[sig.gene[blup.25c.var > 1e-8]]), cbind(geno.line.order, geno.line.order, t(blup.25c.data[blup.25c.var > 1e-8, ]))), file = "eqtl/male.25c.pheno", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE); write.table(rbind(c("FID", "IID", gene.name[sig.gene[blup.diff.var > 1e-8]]), cbind(geno.line.order, geno.line.order, t(blup.diff.data[blup.diff.var > 1e-8, ]))), file = "eqtl/male.diff.pheno", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE);' > log/male.prepPhenoeQTL.Rout 2>&1 &

mkdir eqtl/female.18c/
sbatch --output=log/eqtl.female.18C.log --error=log/eqtl.female.18C.err sbatch/eQTLfemale18C.sbatch

mkdir eqtl/female.25c/
sbatch --output=log/eqtl.female.25C.log --error=log/eqtl.female.25C.err sbatch/eQTLfemale25C.sbatch

mkdir eqtl/male.18c/
sbatch --output=log/eqtl.male.18C.log --error=log/eqtl.male.18C.err sbatch/eQTLmale18C.sbatch

mkdir eqtl/male.25c/
sbatch --output=log/eqtl.male.25C.log --error=log/eqtl.male.25C.err sbatch/eQTLmale25C.sbatch

# 2. permutation 
# ============================================================

mkdir eqtl/permPheno
Rscript -e 'set.seed(1); perm.seq <- matrix(ncol = 200, nrow = 100); for (i in 1:100) { perm.seq[i, ] <- sample(1:200); }; write.table(perm.seq, file = "eqtl/permPheno/perm.seq", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ");'

Rscript -e 'perm.seq <- read.table("eqtl/permPheno/perm.seq", header = FALSE, as.is = TRUE); pheno <- read.table("eqtl/female.18c.pheno", header = TRUE, as.is = TRUE); na.line <- which(is.na(pheno[, 3])); for (i in 1:100) { for (j in na.line) { if (perm.seq[i, j] != j) { swap <- which(perm.seq[i, ] == j); perm.seq[i, swap] <- perm.seq[i, j]; perm.seq[i, j] <- j; } }; write.table(cbind(pheno[, 1:2], pheno[unlist(perm.seq[i, ]), -(1:2)]), file = paste("eqtl/permPheno/female.18c.perm", i, ".pheno", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = " "); cat(i, "\n"); }' > log/female.18c.perm.pheno.Rout 2>&1 &
Rscript -e 'perm.seq <- read.table("eqtl/permPheno/perm.seq", header = FALSE, as.is = TRUE); pheno <- read.table("eqtl/female.25c.pheno", header = TRUE, as.is = TRUE); na.line <- which(is.na(pheno[, 3])); for (i in 1:100) { for (j in na.line) { if (perm.seq[i, j] != j) { swap <- which(perm.seq[i, ] == j); perm.seq[i, swap] <- perm.seq[i, j]; perm.seq[i, j] <- j; } }; write.table(cbind(pheno[, 1:2], pheno[unlist(perm.seq[i, ]), -(1:2)]), file = paste("eqtl/permPheno/female.25c.perm", i, ".pheno", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = " "); cat(i, "\n"); }' > log/female.25c.perm.pheno.Rout 2>&1 &

Rscript -e 'perm.seq <- read.table("eqtl/permPheno/perm.seq", header = FALSE, as.is = TRUE); pheno <- read.table("eqtl/male.18c.pheno", header = TRUE, as.is = TRUE); na.line <- which(is.na(pheno[, 3])); for (i in 1:100) { for (j in na.line) { if (perm.seq[i, j] != j) { swap <- which(perm.seq[i, ] == j); perm.seq[i, swap] <- perm.seq[i, j]; perm.seq[i, j] <- j; } }; write.table(cbind(pheno[, 1:2], pheno[unlist(perm.seq[i, ]), -(1:2)]), file = paste("eqtl/permPheno/male.18c.perm", i, ".pheno", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = " "); cat(i, "\n"); }' > log/male.18c.perm.pheno.Rout 2>&1 &
Rscript -e 'perm.seq <- read.table("eqtl/permPheno/perm.seq", header = FALSE, as.is = TRUE); pheno <- read.table("eqtl/male.25c.pheno", header = TRUE, as.is = TRUE); na.line <- which(is.na(pheno[, 3])); for (i in 1:100) { for (j in na.line) { if (perm.seq[i, j] != j) { swap <- which(perm.seq[i, ] == j); perm.seq[i, swap] <- perm.seq[i, j]; perm.seq[i, j] <- j; } }; write.table(cbind(pheno[, 1:2], pheno[unlist(perm.seq[i, ]), -(1:2)]), file = paste("eqtl/permPheno/male.25c.perm", i, ".pheno", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = " "); cat(i, "\n"); }' > log/male.25c.perm.pheno.Rout 2>&1 &

# 4. perform eqtl mapping for the permuted data
# ============================================================

mkdir eqtl/perm

for i in $(seq 100)
do
  mkdir eqtl/perm/female.18c.perm"$i"
  mkdir eqtl/perm/female.25c.perm"$i"
  mkdir eqtl/perm/male.18c.perm"$i"
  mkdir eqtl/perm/male.25c.perm"$i"
done

# 5. run job arrays
# ============================================================

sbatch sbatch/eQTLfemale18CJobArray.sbatch
sbatch sbatch/eQTLfemale25CJobArray.sbatch

sbatch sbatch/eQTLmale18CJobArray.sbatch
sbatch sbatch/eQTLmale25CJobArray.sbatch

# summarize permutation results
# ============================================================

Rscript R/eQTLfdr.R eqtl/female.18c.pheno eqtl/female.18c/ eqtl/perm/female.18c.perm 0.05 100 10 eqtl/female.18c.eqtl.fdr05.out > log/eqtl.female.18c.fdr05.Rout 2>&1 &
Rscript R/eQTLfdr.R eqtl/female.25c.pheno eqtl/female.25c/ eqtl/perm/female.25c.perm 0.05 100 10 eqtl/female.25c.eqtl.fdr05.out > log/eqtl.female.25c.fdr05.Rout 2>&1 &

Rscript R/eQTLfdr.R eqtl/male.18c.pheno eqtl/male.18c/ eqtl/perm/male.18c.perm 0.05 100 10 eqtl/male.18c.eqtl.fdr05.out > log/eqtl.male.18c.fdr05.Rout 2>&1 &
Rscript R/eQTLfdr.R eqtl/male.25c.pheno eqtl/male.25c/ eqtl/perm/male.25c.perm 0.05 100 10 eqtl/male.25c.eqtl.fdr05.out > log/eqtl.male.25c.fdr05.Rout 2>&1 &

# model selection
# ============================================================

perl -ne 'chomp $_; @line = split /\t/, $_; if ($line[2] ne "exon") { next; } if ($line[8] =~ m/gene_id \"(.*?)\";/) { print $1, "\t", $line[0], "\t", $line[3], "\t", $line[4], "\t", $line[6], "\n"; }' /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/stringtie.merge/r6.36.plus.candidate.gtf  | sort -k1,1 | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2,3,4,5 -o distinct,min,max,distinct | awk '{printf $1"\t"$2"\t"$5"\t"; if ($5 == "+") print $3; else if ($5 == "-") print $4; else print int(($3 + $4)/2); }' > eqtl/gene.tss

# snp genotypes
cat eqtl/*.fdr05.out | cut -d " " -f 3 | grep -v NA | sed 's/,/\n/g' | sort | uniq > eqtl/eqtl.fdr05.snp
/mnt/research/qgg/software/plink-v1.90b6.18/plink --silent --bfile eqtl/dgrp.common --allow-extra-chr --extract eqtl/eqtl.fdr05.snp --recode12 --transpose --out eqtl/eqtl.fdr05.snp.geno

Rscript R/modelSelect.R eqtl/eqtl.fdr05.snp.geno.tped eqtl/eqtl.fdr05.snp.geno.tfam eqtl/female.18c.pheno eqtl/female.18c.eqtl.fdr05.out eqtl/gene.tss 10 0.00001 10 eqtl/female.18c.fdr05.model.select.RData > log/female.18c.fdr05.model.select.Rout 2>&1 &
Rscript R/modelSelect.R eqtl/eqtl.fdr05.snp.geno.tped eqtl/eqtl.fdr05.snp.geno.tfam eqtl/female.25c.pheno eqtl/female.25c.eqtl.fdr05.out eqtl/gene.tss 10 0.00001 10 eqtl/female.25c.fdr05.model.select.RData > log/female.25c.fdr05.model.select.Rout 2>&1 &

Rscript R/modelSelect.R eqtl/eqtl.fdr05.snp.geno.tped eqtl/eqtl.fdr05.snp.geno.tfam eqtl/male.18c.pheno eqtl/male.18c.eqtl.fdr05.out eqtl/gene.tss 10 0.00001 10 eqtl/male.18c.fdr05.model.select.RData > log/male.18c.fdr05.model.select.Rout 2>&1 &
Rscript R/modelSelect.R eqtl/eqtl.fdr05.snp.geno.tped eqtl/eqtl.fdr05.snp.geno.tfam eqtl/male.25c.pheno eqtl/male.25c.eqtl.fdr05.out eqtl/gene.tss 10 0.00001 10 eqtl/male.25c.fdr05.model.select.RData > log/male.25c.fdr05.model.select.Rout 2>&1 &

# summarize results
# ============================================================

Rscript -e 'load("figureData/female.adj.qgSingleTemp.RData"); sig.gene.18c <- which(p.adjust(single.temp[, 9], method = "BH") < 0.05); sig.gene.25c <- which(p.adjust(single.temp[, 12], method = "BH") < 0.05); write.table(gene.name[sig.gene.18c], file = "eqtl/female.18c.genvar.id", col.names = FALSE, row.names = FALSE, sep = "\t", quote = F); write.table(gene.name[sig.gene.25c], file = "eqtl/female.25c.genvar.id", col.names = FALSE, row.names = FALSE, sep = "\t", quote = F);'

Rscript -e 'load("figureData/male.adj.qgSingleTemp.RData"); sig.gene.18c <- which(p.adjust(single.temp[, 9], method = "BH") < 0.05); sig.gene.25c <- which(p.adjust(single.temp[, 12], method = "BH") < 0.05); write.table(gene.name[sig.gene.18c], file = "eqtl/male.18c.genvar.id", col.names = FALSE, row.names = FALSE, sep = "\t", quote = F); write.table(gene.name[sig.gene.25c], file = "eqtl/male.25c.genvar.id", col.names = FALSE, row.names = FALSE, sep = "\t", quote = F);'

echo "number of eGenes at 18C in females"
grep -v NA eqtl/female.18c.eqtl.fdr05.out | sed 's/ /\t/g' | sort -k1,1 | join -t $'\t' <(sort eqtl/female.18c.genvar.id) - | wc -l

echo "number of eGenes at 25C in females"
grep -v NA eqtl/female.25c.eqtl.fdr05.out | sed 's/ /\t/g' | sort -k1,1 | join -t $'\t' <(sort eqtl/female.25c.genvar.id) - | wc -l

echo "number of eGenes at 18C in males"
grep -v NA eqtl/male.18c.eqtl.fdr05.out | sed 's/ /\t/g' | sort -k1,1 | join -t $'\t' <(sort eqtl/male.18c.genvar.id) - | wc -l

echo "number of eGenes at 25C in males"
grep -v NA eqtl/male.25c.eqtl.fdr05.out | sed 's/ /\t/g' | sort -k1,1 | join -t $'\t' <(sort eqtl/male.25c.genvar.id) - | wc -l

Rscript -e 'load("eqtl/female.18c.fdr05.model.select.RData"); write.table(mc.select.eqtl, quote = F, row.names = F, col.names = F, sep = "\t")' | sort -k1,1 | join -t $'\t' - <(grep -v NA eqtl/female.18c.eqtl.fdr05.out | sed 's/ /\t/g' | sort -k1,1) | join -t $'\t' <(sort eqtl/female.18c.genvar.id) - | sed 's/VAR_//g' > eqtl/female.18c.fdr05.eqtl.table.txt
Rscript -e 'load("eqtl/female.25c.fdr05.model.select.RData"); write.table(mc.select.eqtl, quote = F, row.names = F, col.names = F, sep = "\t")' | sort -k1,1 | join -t $'\t' - <(grep -v NA eqtl/female.25c.eqtl.fdr05.out | sed 's/ /\t/g' | sort -k1,1) | join -t $'\t' <(sort eqtl/female.25c.genvar.id) - | sed 's/VAR_//g' > eqtl/female.25c.fdr05.eqtl.table.txt

Rscript -e 'load("eqtl/male.18c.fdr05.model.select.RData"); write.table(mc.select.eqtl, quote = F, row.names = F, col.names = F, sep = "\t")' | sort -k1,1 | join -t $'\t' - <(grep -v NA eqtl/male.18c.eqtl.fdr05.out | sed 's/ /\t/g' | sort -k1,1) | join -t $'\t' <(sort eqtl/male.18c.genvar.id) - | sed 's/VAR_//g' > eqtl/male.18c.fdr05.eqtl.table.txt
Rscript -e 'load("eqtl/male.25c.fdr05.model.select.RData"); write.table(mc.select.eqtl, quote = F, row.names = F, col.names = F, sep = "\t")' | sort -k1,1 | join -t $'\t' - <(grep -v NA eqtl/male.25c.eqtl.fdr05.out | sed 's/ /\t/g' | sort -k1,1) | join -t $'\t' <(sort eqtl/male.25c.genvar.id) - | sed 's/VAR_//g' > eqtl/male.25c.fdr05.eqtl.table.txt

# estimate qtl effects
# ============================================================

Rscript R/eQTLeff.R eqtl/female.18c.fdr05.model.select.RData eqtl/female.18c.eqtl.fdr05.out eqtl/female.18c.genvar.id eqtl/female.25c.fdr05.model.select.RData eqtl/female.25c.eqtl.fdr05.out eqtl/female.25c.genvar.id eqtl/eqtl.fdr05.snp.geno.tped eqtl/eqtl.fdr05.snp.geno.tfam eqtl/female.18c.pheno eqtl/female.25c.pheno eqtl/female.eqtl.eff.RData > log/female.eqtl.eff.Rout 2>&1 &
Rscript R/eQTLeff.R eqtl/male.18c.fdr05.model.select.RData eqtl/male.18c.eqtl.fdr05.out eqtl/male.18c.genvar.id eqtl/male.25c.fdr05.model.select.RData eqtl/male.25c.eqtl.fdr05.out eqtl/male.25c.genvar.id eqtl/eqtl.fdr05.snp.geno.tped eqtl/eqtl.fdr05.snp.geno.tfam eqtl/male.18c.pheno eqtl/male.25c.pheno eqtl/male.eqtl.eff.RData > log/male.eqtl.eff.Rout 2>&1 &

# prediction based on eQTLs
# ============================================================

Rscript R/eQTLpredict.R eqtl/female.18c.fdr05.model.select.RData eqtl/female.18c.genvar.id eqtl/female.25c.fdr05.model.select.RData eqtl/female.25c.genvar.id eqtl/eqtl.fdr05.snp.geno.tped eqtl/eqtl.fdr05.snp.geno.tfam eqtl/female.18c.pheno eqtl/female.25c.pheno eqtl/female.eqtl.pred.RData > log/female.eqtl.pred.Rout 2>&1 &
Rscript R/eQTLpredict.R eqtl/male.18c.fdr05.model.select.RData eqtl/male.18c.genvar.id eqtl/male.25c.fdr05.model.select.RData eqtl/male.25c.genvar.id eqtl/eqtl.fdr05.snp.geno.tped eqtl/eqtl.fdr05.snp.geno.tfam eqtl/male.18c.pheno eqtl/male.25c.pheno eqtl/male.eqtl.pred.RData > log/male.eqtl.pred.Rout 2>&1 &

# eQTL reg region enrichment
# ============================================================

awk -F "\t" '$3 == "enhancer" || $3 == "insulator" || $3 == "regulatory_region" || $3 == "silencer" || $3 == "TF_binding_site" || $3 == "transposable_element"' /mnt/research/qgg/resource/flybase/r6.36/annot/dmel-all-no-analysis.gff | perl -ne 'chomp $_; @line = split /\t/, $_; print $line[0], "\t", $line[3] - 1, "\t", $line[4], "\t"; if ($line[8] =~ m/ID=(.+?);/) { $id = $1; } $line[1] =~ s/ /_/g; print $line[2], "\t", $line[1], "\t", $id, "\n";' > flybase.track.bed

cut -f 2 /mnt/research/qgg/dgrp-age-exp/eqtl/dgrp.common.bim | awk -F "_" '{print $1"\t"$2-1"\t"$2"\t"$0}' | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools intersect -a - -b flybase.track.bed -wb | awk '{print $4"\t("$8"|"$9"|"$10")"}' | sort -k1,1 | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse -delim ";" > dgrp.common.all.reg.txt

Rscript -e 'load("eqtl/female.eqtl.eff.RData"); write.table(eqtl.pair[, 2:3], col.names = F, row.names = F, quote = F, sep = "\t");' | cut -f 1,2 | sort -k1,1 | join -t $'\t' -a 2 -e - -o 0,1.2,2.2 - dgrp.common.all.reg.txt > eqtl/female.eqtl.reg.overlap.txt &

Rscript -e 'load("eqtl/male.eqtl.eff.RData"); write.table(eqtl.pair[, 2:3], col.names = F, row.names = F, quote = F, sep = "\t");' | cut -f 1,2 | sort -k1,1 | join -t $'\t' -a 2 -e - -o 0,1.2,2.2 - dgrp.common.all.reg.txt > eqtl/male.eqtl.reg.overlap.txt &

# use a resampling based approach to get table
# ============================================================

awk '$4 == "TF_binding_site" {print $5}' flybase.track.bed | sort | uniq -c | awk '$1 >= 1000 {print $2}' | sort | uniq > eqtl/tf.site.list

# in females 1324 both, 581 18c, 1137 25c
female18cCount=`cut -f 2 eqtl/female.eqtl.reg.overlap.txt | sort | uniq -c | awk '$2 == "18c" {print $1}'`
female25cCount=`cut -f 2 eqtl/female.eqtl.reg.overlap.txt | sort | uniq -c | awk '$2 == "25c" {print $1}'`
femalebothCount=`cut -f 2 eqtl/female.eqtl.reg.overlap.txt | sort | uniq -c | awk '$2 == "both" {print $1}'`
femaletotalCount=`cut -f 2 eqtl/female.eqtl.reg.overlap.txt | awk '$1 != "-"' | wc -l`

cat <(seq 1 $femalebothCount | awk '{print "both"}') <(seq 1 $female18cCount | awk '{print "18c"}') <(seq 1 $female25cCount | awk '{print "25c"}') > eqtl/female.eqtl.count

awk '$2 != "-" {print $2"\tobs\t"$3}' eqtl/female.eqtl.reg.overlap.txt | awk '$3 ~ /TF_binding_site/' | perl -wne 'chomp $_; @line = split /\t/, $_; @sites = split /;/, $line[2]; for (my $i = 0; $i <= $#sites; $i++) { print $sites[$i], "\t", $line[0], "\t", $line[1], "\n"; }' | grep TF_binding_site | tr "\|" "\t" | awk '{print $2"\t"$4"\t"$5}' | sort -k1,1 | join -t $'\t' eqtl/tf.site.list - > eqtl/female.tf.perm.count

for i in $(seq 1 10000)
do
	# get random sites
	awk '{print "perm'$i'\t"$0}' eqtl/female.eqtl.reg.overlap.txt | shuf | head -n $femaletotalCount | awk '{print $1"\t"$4}' | paste eqtl/female.eqtl.count - | awk '$3 ~ /TF_binding_site/' | perl -wne 'chomp $_; @line = split /\t/, $_; @sites = split /;/, $line[2]; for (my $i = 0; $i <= $#sites; $i++) { print $sites[$i], "\t", $line[0], "\t", $line[1], "\n"; }' | grep TF_binding_site | tr "\|" "\t" | awk '{print $2"\t"$4"\t"$5}' | sort -k1,1 | join -t $'\t' eqtl/tf.site.list - >> eqtl/female.tf.perm.count
done &


# in males 820 both, 524 18c, 1248 25c
male18cCount=`cut -f 2 eqtl/male.eqtl.reg.overlap.txt | sort | uniq -c | awk '$2 == "18c" {print $1}'`
male25cCount=`cut -f 2 eqtl/male.eqtl.reg.overlap.txt | sort | uniq -c | awk '$2 == "25c" {print $1}'`
malebothCount=`cut -f 2 eqtl/male.eqtl.reg.overlap.txt | sort | uniq -c | awk '$2 == "both" {print $1}'`
maletotalCount=`cut -f 2 eqtl/male.eqtl.reg.overlap.txt | awk '$1 != "-"' | wc -l`

cat <(seq 1 $malebothCount | awk '{print "both"}') <(seq 1 $male18cCount | awk '{print "18c"}') <(seq 1 $male25cCount | awk '{print "25c"}') > eqtl/male.eqtl.count

awk '$2 != "-" {print $2"\tobs\t"$3}' eqtl/male.eqtl.reg.overlap.txt | awk '$3 ~ /TF_binding_site/' | perl -wne 'chomp $_; @line = split /\t/, $_; @sites = split /;/, $line[2]; for (my $i = 0; $i <= $#sites; $i++) { print $sites[$i], "\t", $line[0], "\t", $line[1], "\n"; }' | grep TF_binding_site | tr "\|" "\t" | awk '{print $2"\t"$4"\t"$5}' | sort -k1,1 | join -t $'\t' eqtl/tf.site.list - > eqtl/male.tf.perm.count

for i in $(seq 1 10000)
do
	# get random sites
	awk '{print "perm'$i'\t"$0}' eqtl/male.eqtl.reg.overlap.txt | shuf | head -n $maletotalCount | awk '{print $1"\t"$4}' | paste eqtl/male.eqtl.count - | awk '$3 ~ /TF_binding_site/' | perl -wne 'chomp $_; @line = split /\t/, $_; @sites = split /;/, $line[2]; for (my $i = 0; $i <= $#sites; $i++) { print $sites[$i], "\t", $line[0], "\t", $line[1], "\n"; }' | grep TF_binding_site | tr "\|" "\t" | awk '{print $2"\t"$4"\t"$5}' | sort -k1,1 | join -t $'\t' eqtl/tf.site.list - >> eqtl/male.tf.perm.count
done &

# summarize results
# ============================================================

Rscript R/summarizeTF.R eqtl/tf.site.list eqtl/female.tf.perm.count 10000 eqtl/female.tf.summary.RData > log/female.tf.summary.Rout 2>&1 &
Rscript R/summarizeTF.R eqtl/tf.site.list eqtl/male.tf.perm.count 10000 eqtl/male.tf.summary.RData > log/male.tf.summary.Rout 2>&1 &

