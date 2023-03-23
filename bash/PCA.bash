# 1. prepare file with gene information
# ============================================================

perl -wne 'chomp $_; @line = split /\t/, $_; if ($line[8] =~ m/class_code \"j\"/) { next; }; if ($line[8] =~ m/gene_id \"(.*?)\"/) { $id = $1; }; $symbol = "-"; $type = "-"; if (!($id =~ /XLOC/)) { if ($line[8] =~ m/gene_symbol \"(.*?)\"/) { $symbol = $1; }; } else { $type = "NTR"; }; print $id, "\t", $line[0], "\t", $line[3], "\t", $line[4], "\t", $symbol, "\t", $type, "\n"; ' /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/stringtie.merge/r6.36.plus.candidate.gtf | sort -k1,1 | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools groupby -i - -g 1 -c 2,3,4,5,6 -o distinct,min,max,distinct,distinct | awk '{print $1"\t"$5"\t"$6"\t"$2":"$3"-"$4}' > /mnt/research/qgg/dgrp-age-exp/joinExp/gene.info

# 2. perform a pca on normal quantile transformed data
# first on two time combined, then within time
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

# 3. sva adjustment
# ============================================================

Rscript sva.R joinedExp.female.txt female.fc adjustData.RData female.sva.adjust.RData female.gene.exp.adjust.RData > log/female.sva.Rout 2>&1 &
Rscript sva.R joinedExp.male.txt male.fc adjustData.RData male.sva.adjust.RData male.gene.exp.adjust.RData > log/male.sva.Rout 2>&1 &

# figurePrelimNormOutlier
# ============================================================

cd /mnt/research/qgg/dgrp-age-exp/joinExp
FONT="Myriad Pro"
Rscript figurePrelimNormOutlier.R female.gene.exp.adjust.RData male.gene.exp.adjust.RData Figure_PreliminaryNormalizationOutlier.pdf $(FONT)

# figurePCA
# ============================================================
Rscript figurePCA.R female.pca.RData female.sva.adjust.RData male.pca.RData male.sva.adjust.RData figurePCA.pdf $(FONT)

