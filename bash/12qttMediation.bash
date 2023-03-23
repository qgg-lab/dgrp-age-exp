# GWAS for the different phenotypes
# ============================================================

for trait in female.phototaxis.decline female.fecundity.decline female.lifespan male.phototaxis.decline male.speed.decline male.endurance.decline male.lifespan
do
	/mnt/research/qgg/software/plink-v1.90b6.18/plink --silent --bfile /mnt/research/qgg/dgrp-age-exp/eqtl/dgrp.common --pheno ../qtt/"$trait".pheno --all-pheno --allow-extra-chr --allow-no-sex --assoc --out "$trait".gwas &
done

for trait in female.phototaxis.decline female.fecundity.decline female.lifespan male.phototaxis.decline male.speed.decline male.endurance.decline male.lifespan
do
	awk '$9 < 1e-4 { print "'$trait'\t" $0}' "$trait".gwas.pheno.qassoc
done > gwas.snp.qassoc



# get all eqtls
# ============================================================

grep -v NA ../eqtl/female.25c.eqtl.fdr05.out | sed 's/ /\t/g' | sort -k1,1 | join -t $'\t' <(sort ../eqtl/female.25c.genvar.id) - | perl -wne 'chomp $_; @line = split /\t/, $_; @snps = split /,/, $line[2]; for (my $i = 0; $i <= $#snps; $i++) { print $snps[$i], "\t", $line[0], "\n"; }' | sort -k1,1 | join -t $'\t' - <(awk '{print $3"\t"$1"\t"$10}' gwas.snp.qassoc | sort -k1,1) | grep female | awk '{print $0"\tyoung-egene"}' > female.eqtl.gwas.trio

grep -v NA ../eqtl/female.18c.eqtl.fdr05.out | sed 's/ /\t/g' | sort -k1,1 | join -t $'\t' <(sort ../eqtl/female.18c.genvar.id) - | perl -wne 'chomp $_; @line = split /\t/, $_; @snps = split /,/, $line[2]; for (my $i = 0; $i <= $#snps; $i++) { print $snps[$i], "\t", $line[0], "\n"; }' | sort -k1,1 | join -t $'\t' - <(awk '{print $3"\t"$1"\t"$10}' gwas.snp.qassoc | sort -k1,1) | grep female | awk '{print $0"\told-egene"}' >> female.eqtl.gwas.trio

grep -v NA ../eqtl/male.25c.eqtl.fdr05.out | sed 's/ /\t/g' | sort -k1,1 | join -t $'\t' <(sort ../eqtl/male.25c.genvar.id) - | perl -wne 'chomp $_; @line = split /\t/, $_; @snps = split /,/, $line[2]; for (my $i = 0; $i <= $#snps; $i++) { print $snps[$i], "\t", $line[0], "\n"; }' | sort -k1,1 | join -t $'\t' - <(awk '{print $3"\t"$1"\t"$10}' gwas.snp.qassoc | sort -k1,1) | grep male | grep -v female | awk '{print $0"\tyoung-egene"}' > male.eqtl.gwas.trio

grep -v NA ../eqtl/male.18c.eqtl.fdr05.out | sed 's/ /\t/g' | sort -k1,1 | join -t $'\t' <(sort ../eqtl/male.18c.genvar.id) - | perl -wne 'chomp $_; @line = split /\t/, $_; @snps = split /,/, $line[2]; for (my $i = 0; $i <= $#snps; $i++) { print $snps[$i], "\t", $line[0], "\n"; }' | sort -k1,1 | join -t $'\t' - <(awk '{print $3"\t"$1"\t"$10}' gwas.snp.qassoc | sort -k1,1) | grep male | grep -v female | awk '{print $0"\told-egene"}' >> male.eqtl.gwas.trio

# mediation analysis
# ============================================================

module load R && srun --cpus-per-task=1 --ntasks-per-node=1 --mem=16G --time=24:00:00 Rscript ../R/mediationQTT.R ../eqtl/eqtl.fdr05.snp.geno.tped ../eqtl/eqtl.fdr05.snp.geno.tfam ../qtt/female.phototaxis.decline.pheno,../qtt/female.fecundity.decline.pheno,../qtt/female.lifespan.pheno female.eqtl.gwas.trio ../eqtl/female.25c.pheno ../eqtl/female.18c.pheno female.eqtl.gwas.trio.mediation.RData > female.eqtl.gwas.trio.mediation.Rout 2>&1 &

module load R && srun --cpus-per-task=1 --ntasks-per-node=1 --mem=16G --time=24:00:00 Rscript ../R/mediationQTT.R ../eqtl/eqtl.fdr05.snp.geno.tped ../eqtl/eqtl.fdr05.snp.geno.tfam ../qtt/male.phototaxis.decline.pheno,../qtt/male.speed.decline.pheno,../qtt/male.endurance.decline.pheno,../qtt/male.lifespan.pheno male.eqtl.gwas.trio ../eqtl/male.25c.pheno ../eqtl/male.18c.pheno male.eqtl.gwas.trio.mediation.RData > male.eqtl.gwas.trio.mediation.Rout 2>&1 &

# get example data set
# ============================================================

module load R && Rscript ../R/mediationQTTexample.R ../eqtl/eqtl.fdr05.snp.geno.tped ../eqtl/eqtl.fdr05.snp.geno.tfam ../qtt/male.phototaxis.decline.pheno,../qtt/male.speed.decline.pheno,../qtt/male.endurance.decline.pheno,../qtt/male.lifespan.pheno male.eqtl.gwas.trio ../eqtl/male.25c.pheno ../eqtl/male.18c.pheno male.speed.decline,male.lifespan,male.phototaxis.decline 3R_12642719,3R_12642725,2R_4440030,2R_22212424 FBgn0026207,FBgn0025186,FBgn0260798 male.eqtl.gwas.trio.mediation.example.RData > male.eqtl.gwas.trio.mediation.example.Rout 2>&1 &

