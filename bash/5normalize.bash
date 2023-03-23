# ========================================
# = normalize, filter, adjust expression =
# ========================================

# 1. prepare file with gene information
# ============================================================

mkdir /mnt/research/qgg/dgrp-age-exp/normalize
cd /mnt/research/qgg/dgrp-age-exp/normalize
perl -wne 'chomp $_; @line = split /\t/, $_; if ($line[8] =~ m/class_code \"j\"/) { next; }; if ($line[8] =~ m/gene_id \"(.*?)\"/) { $id = $1; }; $symbol = "-"; $type = "-"; if (!($id =~ /XLOC/)) { if ($line[8] =~ m/gene_symbol \"(.*?)\"/) { $symbol = $1; }; } else { $type = "NTR"; }; print $id, "\t", $line[0], "\t", $line[6], "\t", $line[3], "\t", $line[4], "\t", $symbol, "\t", $type, "\n"; ' /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/stringtie.merge/r6.36.plus.candidate.gtf | sort -k1,1 | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools groupby -i - -g 1 -c 2,3,4,5,6,7 -o distinct,distinct,min,max,distinct,distinct | awk '{print $1"\t"$6"\t"$7"\t"$2":"$3":"$4"-"$5}' | join -t $'\t' -a 1 -e - -o 0,1.2,1.3,1.4,2.2 - /mnt/research/qgg/resource/flybase/r6.36/annot/gene.biotype | awk '{ if ($3 == "NTR") { print $1"\t"$2"\t"$3"\t"$4 } else { print $1"\t"$2"\t"$5"\t"$4 } }'> /mnt/research/qgg/dgrp-age-exp/figureData/gene.info

# 2. flow cell informaiton
# ============================================================

awk '{print $2"_"$3"_"$4 $5"\t"$6}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis | sort -k1,1 | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse -i - > /mnt/research/qgg/dgrp-age-exp/joinExp/all.sample.flow.cell.info

# 3. read data into RData for subsequent processing
# separate by sex, in each sex
# need the data frame gene.exp.adj
# wolba
# and age for each sample's age
# ============================================================

/mnt/research/qgg/software/R-3.6.0/bin/Rscript /mnt/research/qgg/dgrp-age-exp/R/readData.R /mnt/research/qgg/dgrp-age-exp/joinExp/female.exp.tpm /mnt/research/qgg/dgrp-age-exp/joinExp/male.exp.tpm /mnt/research/qgg/dgrp-age-exp/joinExp/all.sample.flow.cell.info /mnt/research/qgg/dgrp-age-exp/normalize/female.gene.exp.RData /mnt/research/qgg/dgrp-age-exp/normalize/male.gene.exp.RData > /mnt/research/qgg/dgrp-age-exp/normalize/readData.Rout 2>&1 &

# 4. scale expression, check for outliers
# ============================================================

/mnt/research/qgg/software/R-3.6.0/bin/Rscript /mnt/research/qgg/dgrp-age-exp/R/scaleExp.R /mnt/research/qgg/dgrp-age-exp/normalize/female.gene.exp.RData /mnt/research/qgg/dgrp-age-exp/figureData/female.gene.exp.scaled.quantile.RData > /mnt/research/qgg/dgrp-age-exp/normalize/female.scaleExp.Rout 2>&1 &
/mnt/research/qgg/software/R-3.6.0/bin/Rscript /mnt/research/qgg/dgrp-age-exp/R/scaleExp.R /mnt/research/qgg/dgrp-age-exp/normalize/male.gene.exp.RData /mnt/research/qgg/dgrp-age-exp/figureData/male.gene.exp.scaled.quantile.RData > /mnt/research/qgg/dgrp-age-exp/normalize/male.scaleExp.Rout 2>&1 &

# identify outliers
/mnt/research/qgg/software/R-3.6.0/bin/Rscript -e 'load("/mnt/research/qgg/dgrp-age-exp/figureData/female.gene.exp.scaled.quantile.RData"); female.bad.sample <- c(gsub("^X", "", colnames(baseline.quantile))[which(baseline.quantile[1, ] < -5 | baseline.quantile[7, ] > 5)], gsub("^X", "", colnames(old.quantile))[which(old.quantile[1, ] < -5 | old.quantile[7, ] > 5)]); load("/mnt/research/qgg/dgrp-age-exp/figureData/male.gene.exp.scaled.quantile.RData"); male.bad.sample <- c(gsub("^X", "", colnames(baseline.quantile))[which(baseline.quantile[1, ] < -5 | baseline.quantile[7, ] > 5)], gsub("^X", "", colnames(old.quantile))[which(old.quantile[1, ] < -5 | old.quantile[7, ] > 5)]); write.table(c(female.bad.sample, male.bad.sample), file = "/mnt/research/qgg/dgrp-age-exp/normalize/bad.sample.id", sep = " ", col.names = F, row.names = F, quote = F);'

# remove outliers
/mnt/research/qgg/software/R-3.6.0/bin/Rscript -e 'bad.sample <- scan(file = "/mnt/research/qgg/dgrp-age-exp/normalize/bad.sample.id", what = ""); load("/mnt/research/qgg/dgrp-age-exp/normalize/female.gene.exp.RData"); sample.to.remove <- match(intersect(bad.sample, sample.name), sample.name); age <- age[-sample.to.remove]; flow.cell <- flow.cell[-sample.to.remove]; gene.exp <- gene.exp[, -sample.to.remove]; sample.name <- sample.name[-sample.to.remove]; flow.cell <- flow.cell[-sample.to.remove, ]; save(age, flow.cell, gene.exp, sample.name, file = "/mnt/research/qgg/dgrp-age-exp/normalize/female.gene.exp.filtered.RData");' > /mnt/research/qgg/dgrp-age-exp/normalize/female.gene.exp.filtered.Rout 2>&1 &

/mnt/research/qgg/software/R-3.6.0/bin/Rscript -e 'bad.sample <- scan(file = "/mnt/research/qgg/dgrp-age-exp/normalize/bad.sample.id", what = ""); load("/mnt/research/qgg/dgrp-age-exp/normalize/male.gene.exp.RData"); sample.to.remove <- match(intersect(bad.sample, sample.name), sample.name); age <- age[-sample.to.remove]; flow.cell <- flow.cell[-sample.to.remove]; gene.exp <- gene.exp[, -sample.to.remove]; sample.name <- sample.name[-sample.to.remove]; flow.cell <- flow.cell[-sample.to.remove, ]; save(age, flow.cell, gene.exp, sample.name, file = "/mnt/research/qgg/dgrp-age-exp/normalize/male.gene.exp.filtered.RData");' > /mnt/research/qgg/dgrp-age-exp/normalize/male.gene.exp.filtered.Rout 2>&1 &

# 5. perform PCA on expression
# ============================================================
cd /mnt/research/qgg/dgrp-age-exp/joinExp

awk '$4 == "F" {print $2"_"$3"_"$4 $5"\t"$6}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis|sort -k1,1 | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse -i -|awk -F , '{print $1}' > /mnt/research/qgg/dgrp-age-exp/joinExp/female.fc

awk '$4 == "M" {print $2"_"$3"_"$4 $5"\t"$6}' /mnt/research/qgg/dgrp-age-exp/checkID/final.samples.to.begin.analysis|sort -k1,1 | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2 -o collapse -i -|awk -F , '{print $1}' > /mnt/research/qgg/dgrp-age-exp/joinExp/male.fc

Rscript PCA.R joinedExp.female.txt female.fc female.pca.RData > log/female.pca.Rout 2>&1 &
Rscript PCA.R joinedExp.male.txt male.fc male.pca.RData > log/male.pca.Rout 2>&1 &

# known.gene 17612
# ntr.gene 696

# known.gene.filter 9878
# ntr.gene.filter 200

# 6. sva adjustment
# ============================================================

/mnt/research/qgg/software/R-3.6.0/bin/Rscript /mnt/research/qgg/dgrp-age-exp/R/sva.R /mnt/research/qgg/dgrp-age-exp/normalize/female.gene.exp.filtered.RData /mnt/research/qgg/resource/dgrp/dgrp2.r6.adjustData.RData /mnt/research/qgg/dgrp-age-exp/figureData/female.sva.adjust.RData /mnt/research/qgg/dgrp-age-exp/normalize/female.gene.exp.adjust.RData > /mnt/research/qgg/dgrp-age-exp/normalize/female.sva.Rout 2>&1 &

/mnt/research/qgg/software/R-3.6.0/bin/Rscript /mnt/research/qgg/dgrp-age-exp/R/sva.R /mnt/research/qgg/dgrp-age-exp/normalize/male.gene.exp.filtered.RData /mnt/research/qgg/resource/dgrp/dgrp2.r6.adjustData.RData /mnt/research/qgg/dgrp-age-exp/figureData/male.sva.adjust.RData /mnt/research/qgg/dgrp-age-exp/normalize/male.gene.exp.adjust.RData > /mnt/research/qgg/dgrp-age-exp/normalize/male.sva.Rout 2>&1 &

