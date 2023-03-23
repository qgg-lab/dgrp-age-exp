# ======================
# = mediation analysis =
# ======================

# classify eQTLs into cis and trans
# for cis eQTL, only consider selected, for trans, also consider mapped
# ============================================================

sort -k1,1 gene.info | join -t $'\t' - <(sort -k1,1 ../eqtl/female.25c.fdr05.eqtl.table.txt) | perl -wne '
 chomp $_; @line = split /\t/, $_; @gene = split /:/, $line[3]; @bound = split /-/, $gene[2]; @select = split /,/, $line[4]; for (my $i = 0; $i <= $#select; $i++) { @thisinfo = split /_/, $select[$i]; if ($thisinfo[0] ne $gene[0]) { print $line[0], "\tselect\t", $select[$i], "\ttrans\n"; } else { if ($thisinfo[1] >= $bound[0] - 5000 && $thisinfo[1] <= $bound[1] + 5000) { print $line[0], "\tselect\t", $select[$i], "\tcis\n"; } else { print $line[0], "\tselect\t", $select[$i], "\ttrans\n"; } } } @map = split /,/, $line[6]; for (my $i = 0; $i <= $#map; $i++) { @thisinfo = split /_/, $map[$i]; if ($thisinfo[0] ne $gene[0]) { print $line[0], "\tmap\t", $map[$i], "\ttrans\n"; } else { if ($thisinfo[1] >= $bound[0] - 5000 && $thisinfo[1] <= $bound[1] + 5000) { print $line[0], "\tmap\t", $map[$i], "\tcis\n"; } else { print $line[0], "\tmap\t", $map[$i], "\ttrans\n"; } } }' > female.young.eqtl.class

awk '$2 == "select"&& $4 == "cis" { print $3"\t"$1"\t"$2"\t"$4 }' female.young.eqtl.class | sort -k1,1 | join -t $'\t' - <(awk '$2 == "map" && $4 == "trans" { print $3"\t"$1"\t"$2"\t"$4 }' female.young.eqtl.class | sort -k1,1) | awk '$2 != $5' > female.young.eqtl.trio


sort -k1,1 gene.info | join -t $'\t' - <(sort -k1,1 ../eqtl/female.18c.fdr05.eqtl.table.txt) | perl -wne '
 chomp $_; @line = split /\t/, $_; @gene = split /:/, $line[3]; @bound = split /-/, $gene[2]; @select = split /,/, $line[4]; for (my $i = 0; $i <= $#select; $i++) { @thisinfo = split /_/, $select[$i]; if ($thisinfo[0] ne $gene[0]) { print $line[0], "\tselect\t", $select[$i], "\ttrans\n"; } else { if ($thisinfo[1] >= $bound[0] - 5000 && $thisinfo[1] <= $bound[1] + 5000) { print $line[0], "\tselect\t", $select[$i], "\tcis\n"; } else { print $line[0], "\tselect\t", $select[$i], "\ttrans\n"; } } } @map = split /,/, $line[6]; for (my $i = 0; $i <= $#map; $i++) { @thisinfo = split /_/, $map[$i]; if ($thisinfo[0] ne $gene[0]) { print $line[0], "\tmap\t", $map[$i], "\ttrans\n"; } else { if ($thisinfo[1] >= $bound[0] - 5000 && $thisinfo[1] <= $bound[1] + 5000) { print $line[0], "\tmap\t", $map[$i], "\tcis\n"; } else { print $line[0], "\tmap\t", $map[$i], "\ttrans\n"; } } }' > female.old.eqtl.class

awk '$2 == "select"&& $4 == "cis" { print $3"\t"$1"\t"$2"\t"$4 }' female.old.eqtl.class | sort -k1,1 | join -t $'\t' - <(awk '$2 == "map" && $4 == "trans" { print $3"\t"$1"\t"$2"\t"$4 }' female.old.eqtl.class | sort -k1,1) | awk '$2 != $5' > female.old.eqtl.trio


sort -k1,1 gene.info | join -t $'\t' - <(sort -k1,1 ../eqtl/male.25c.fdr05.eqtl.table.txt) | perl -wne '
 chomp $_; @line = split /\t/, $_; @gene = split /:/, $line[3]; @bound = split /-/, $gene[2]; @select = split /,/, $line[4]; for (my $i = 0; $i <= $#select; $i++) { @thisinfo = split /_/, $select[$i]; if ($thisinfo[0] ne $gene[0]) { print $line[0], "\tselect\t", $select[$i], "\ttrans\n"; } else { if ($thisinfo[1] >= $bound[0] - 5000 && $thisinfo[1] <= $bound[1] + 5000) { print $line[0], "\tselect\t", $select[$i], "\tcis\n"; } else { print $line[0], "\tselect\t", $select[$i], "\ttrans\n"; } } } @map = split /,/, $line[6]; for (my $i = 0; $i <= $#map; $i++) { @thisinfo = split /_/, $map[$i]; if ($thisinfo[0] ne $gene[0]) { print $line[0], "\tmap\t", $map[$i], "\ttrans\n"; } else { if ($thisinfo[1] >= $bound[0] - 5000 && $thisinfo[1] <= $bound[1] + 5000) { print $line[0], "\tmap\t", $map[$i], "\tcis\n"; } else { print $line[0], "\tmap\t", $map[$i], "\ttrans\n"; } } }' > male.young.eqtl.class

awk '$2 == "select"&& $4 == "cis" { print $3"\t"$1"\t"$2"\t"$4 }' male.young.eqtl.class | sort -k1,1 | join -t $'\t' - <(awk '$2 == "map" && $4 == "trans" { print $3"\t"$1"\t"$2"\t"$4 }' male.young.eqtl.class | sort -k1,1) | awk '$2 != $5' > male.young.eqtl.trio



sort -k1,1 gene.info | join -t $'\t' - <(sort -k1,1 ../eqtl/male.18c.fdr05.eqtl.table.txt) | perl -wne '
 chomp $_; @line = split /\t/, $_; @gene = split /:/, $line[3]; @bound = split /-/, $gene[2]; @select = split /,/, $line[4]; for (my $i = 0; $i <= $#select; $i++) { @thisinfo = split /_/, $select[$i]; if ($thisinfo[0] ne $gene[0]) { print $line[0], "\tselect\t", $select[$i], "\ttrans\n"; } else { if ($thisinfo[1] >= $bound[0] - 5000 && $thisinfo[1] <= $bound[1] + 5000) { print $line[0], "\tselect\t", $select[$i], "\tcis\n"; } else { print $line[0], "\tselect\t", $select[$i], "\ttrans\n"; } } } @map = split /,/, $line[6]; for (my $i = 0; $i <= $#map; $i++) { @thisinfo = split /_/, $map[$i]; if ($thisinfo[0] ne $gene[0]) { print $line[0], "\tmap\t", $map[$i], "\ttrans\n"; } else { if ($thisinfo[1] >= $bound[0] - 5000 && $thisinfo[1] <= $bound[1] + 5000) { print $line[0], "\tmap\t", $map[$i], "\tcis\n"; } else { print $line[0], "\tmap\t", $map[$i], "\ttrans\n"; } } }' > male.old.eqtl.class

awk '$2 == "select"&& $4 == "cis" { print $3"\t"$1"\t"$2"\t"$4 }' male.old.eqtl.class | sort -k1,1 | join -t $'\t' - <(awk '$2 == "map" && $4 == "trans" { print $3"\t"$1"\t"$2"\t"$4 }' male.old.eqtl.class | sort -k1,1) | awk '$2 != $5' > male.old.eqtl.trio

file="male.old.eqtl.class"
# cis eQTLs
awk '$2 == "select" && $4 == "cis"' $file | wc -l
# cis eGenes
awk '$2 == "select" && $4 == "cis" {print $1}' $file | sort | uniq | wc -l

# trans eQTLs
awk '$2 == "select" && $4 == "trans"' $file | wc -l
# trans eGenes
awk '$2 == "select" && $4 == "trans" {print $1}' $file | sort | uniq | wc -l


# run mediation tests
# ============================================================
module load R && srun --cpus-per-task=1 --ntasks-per-node=1 --mem=8G --time=24:00:00 Rscript ../R/mediation.R ../eqtl/eqtl.fdr05.snp.geno.tped ../eqtl/eqtl.fdr05.snp.geno.tfam ../eqtl/female.25c.pheno female.young.eqtl.trio female.young.trio.mediation.RData > female.young.trio.mediation.Rout 2>&1 &

module load R && srun --cpus-per-task=1 --ntasks-per-node=1 --mem=8G --time=24:00:00 Rscript ../R/mediation.R ../eqtl/eqtl.fdr05.snp.geno.tped ../eqtl/eqtl.fdr05.snp.geno.tfam ../eqtl/female.18c.pheno female.old.eqtl.trio female.old.trio.mediation.RData > female.old.trio.mediation.Rout 2>&1 &

module load R && srun --cpus-per-task=1 --ntasks-per-node=1 --mem=8G --time=24:00:00 Rscript ../R/mediation.R ../eqtl/eqtl.fdr05.snp.geno.tped ../eqtl/eqtl.fdr05.snp.geno.tfam ../eqtl/male.25c.pheno male.young.eqtl.trio male.young.trio.mediation.RData > male.young.trio.mediation.Rout 2>&1 &

module load R && srun --cpus-per-task=1 --ntasks-per-node=1 --mem=8G --time=24:00:00 Rscript ../R/mediation.R ../eqtl/eqtl.fdr05.snp.geno.tped ../eqtl/eqtl.fdr05.snp.geno.tfam ../eqtl/male.18c.pheno male.old.eqtl.trio male.old.trio.mediation.RData > male.old.trio.mediation.Rout 2>&1 &

# get example
# ============================================================
module load R && Rscript ../R/screenTrio.R female.young.trio.mediation.RData ../eqtl/eqtl.fdr05.snp.geno.tped ../eqtl/eqtl.fdr05.snp.geno.tfam ../eqtl/female.25c.pheno ../eqtl/female.18c.pheno female.example.mediation.RData



