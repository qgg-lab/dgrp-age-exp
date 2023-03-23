# ======================
# = mediation analysis =
# ======================

cd /mnt/research/qgg/dgrp-age-exp
mkdir mediation

# use female 25c as an example

# ============================================
# prepare genotype:
awk '{for(x=1;x<=NF;x++)if(x % 2==0)printf "%s", $x (x == NF || x == (NF-1)?"\n":" ")}' eqtl/eqtl.fdr05.snp.geno.tped | cut -d " " -f 1,3- - > mediation/eqtl.fdr05.snp.geno

# ============================================
# phenotype:
# /mnt/research/qgg/dgrp-age-exp/eqtl/female.25c.pheno

# ============================================
# get cis-eqtl and trans-eqtl, cis as within the region of 500k upstream and dowmstream from the gene region, can be changed

# get gene range in bed
perl -wne 'chomp $_; @line = split /\t/, $_; if ($line[2] ne "exon") { next; } if ($line[8] =~ m/gene_id \"(.*?)\";/) { print $1, "\t", $line[0], "\t", $line[3], "\t", $line[4], "\t", $line[6], "\n"; }' /mnt/gs18/scratch/users/tansuxu/dgrp/prepare.GTF/stringtie.merge/r6.36.plus.candidate.gtf  | sort -k1,1 | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools groupby -g 1 -c 2,3,4,5 -o distinct,min,max,distinct | awk '{print $2"\t"$3-1"\t"$4"\t"$1}' > mediation/r6.36.plus.candidate.gene

# obtain gene-eqtl pair and add gene range to it
grep -v NA eqtl/female.25c.eqtl.fdr05.out | sort -k1,1 |awk '{print $1"\t"$3}'| join -t $'\t' <(awk '{print $4"\t"$1":"$2"-"$3}' mediation/r6.36.plus.candidate.gene | sort -k1,1) - > mediation/female.25c.eqtl

# list gene-eqtl pair
cd mediation

perl -wne 'chomp $_; @line = split /[\t,]/, $_; if ($#line == 2) { print join("\t", @line), "\n";} if ($#line > 2) {for (my $i=2; $i<=$#line; $i++) {print join("\t", @line[0..1]), "\t", $line[$i], "\n";} }' female.25c.eqtl > female.25c.eqtl.all

# obtain cis-eqtl for genes
awk '{print $1}' female.25c.eqtl.all > female.25c.eqtl.1
awk '{print $2}' female.25c.eqtl.all > female.25c.eqtl.2
awk '{print $3}' female.25c.eqtl.all > female.25c.eqtl.3

paste female.25c.eqtl.1 <(sed -E -e 's/(:|-)/\t/g' female.25c.eqtl.2) <(sed 's/_/\t/g' female.25c.eqtl.3) | awk '$5==$2 && $6>=$3-500000 && $6<=$4+500000 {print $0}' | awk '{print $5"_"$6"\t"$1}' > female.25c.eqtlCis.gene

# obtain trans-eqtl for genes
paste female.25c.eqtl.1 <(sed -E -e 's/(:|-)/\t/g' female.25c.eqtl.2) <(sed 's/_/\t/g' female.25c.eqtl.3) | awk '!($5==$2 && $6>=$3-500000 && $6<=$4+500000) {print $0}' | awk '{print $5"_"$6"\t"$1}' > female.25c.eqtlTrans.gene


# ============================================
# get the trios (L, C, T) in the genome showing both cis- and trans-eQTL associations
# i.e., Locus → Cis-gene and Locus → Trans-gene.

join -t $'\t' <(sort -k1,1 female.25c.eqtlCis.gene) <(sort -k1,1 female.25c.eqtlTrans.gene) > female.25c.eqtl.trio
# wc -l female.25c.eqtl.trio #19874

# ============================================
# For each trio, try to use gene from each module of WGCNA result
# For example, get the genes from module 1 (2607 genes) of female blup 25c 
Rscript -e 'load("/mnt/research/qgg/dgrp-age-exp/figureData/female.wgcna.RData"); module <- names(which(blup.25c.tree==1)); trios <- read.table("female.25c.eqtl.trio"); idx1 <- na.omit(match(module,trios$V2)); trios <- trios[idx1,]; idx2 <- na.omit(match(module,trios$V3)); trios <- trios[idx2,]' 
# Zero trios is kept after matching module 1 (2607 genes) against both cis-genes and trans-genes of 19874 trios

# ============================================
# run mediation analysis using all, not gene from WGCNA module

Rscript -e 'exp <- read.table("/mnt/research/qgg/dgrp-age-exp/eqtl/female.25c.pheno", header = T, row.names = 1, as.is = T); exp <- exp[,-1]; exp <- t(exp); fea.data <- exp[,-(match(colnames(exp)[colSums(is.na(exp)) > 0],colnames(exp)))]; fea.data <- as.matrix(fea.data); snp.data <- read.table("eqtl.fdr05.snp.geno",row.names = 1, as.is = T); snp.data <- snp.data[,-(match(colnames(exp)[colSums(is.na(exp)) > 0],colnames(exp)))]; snp.data <- as.matrix(snp.data); load("/mnt/research/qgg/dgrp-age-exp/figureData/female.wgcna.RData"); module <- names(which(blup.25c.tree==1)); trios <- read.table("female.25c.eqtl.trio"); trios$V1 <- match(trios$V1,rownames(snp.data)); trios$V2 <- match(trios$V2,rownames(fea.data)); trios$V3 <- match(trios$V3,rownames(fea.data)); trios <- as.matrix(trios); trios <- na.omit(trios); pc <- prcomp(t(fea.data), scale = T); cov.pool <- t(pc$x); known.confounder <- cov.pool[1:10,]; library(GMAC); cl = makeCluster(20); output <- gmac(cl = cl, known.conf = known.confounder, cov.pool = NULL, exp.dat = fea.data, snp.dat.cis = snp.data, trios.idx = trios, nperm = 1000, nominal.p = FALSE); save(fea.data, snp.data, trios, output, file = "female.25c.mediation.RData")' > log/female.25c.mediation.Rout 2>&1 &

# nrow(trios)
# 19874

# nrow(trios)after module match and na.omit
# 0

# sum(output$pvals[,1] < 0.05)
# 13746

# sum(p.adjust(output$pvals[,1], "BH") < 0.05)
# 12831