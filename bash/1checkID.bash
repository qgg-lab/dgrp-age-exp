# =========================================================
# = check if RNA-Seq fastq's agree with DNA-Seq genotypes =
# =========================================================

# 1. filter DGRP genotypes to contain only
#    1) MAF > 0.2
#    2) only variants overlap with exons
#    3) only biallelic SNPs
#    retain imputed sequences
# ============================================================

# get exon list
awk '$3 == "exon" {print $1"\t"$4-1"\t"$5"\t.\t.\t"$7}' /mnt/research/qgg/resource/flybase/r6.36/annot/dmel-all.gtf | sort -k1,1 -k2,2n | /mnt/research/qgg/software/bedtools-2.29.2/bin/bedtools merge -i - -s > /mnt/research/qgg/dgrp/dmel.exon.merge.bed

# filter VCF based on the criteria above
/mnt/research/qgg/software/vcftools_0.1.13/bin/vcftools --vcf /mnt/gs18/scratch/users/tansuxu/dgrp/jgil/all.jgil.hpp.vcf --remove-indels --maf 0.2 --max-alleles 2 --min-alleles 2 --bed dmel.exon.merge.bed --recode --out /mnt/research/qgg/resource/dgrp/all.jgil.exon.common > /mnt/research/qgg/resource/dgrp/all.jgil.exon.common.vcf.log 2>&1

# index vcf
java -jar /mnt/research/qgg/software/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar IndexFeatureFile --input /mnt/research/qgg/resource/dgrp/all.jgil.exon.common.recode.vcf > /mnt/research/qgg/resource/dgrp/vcf.index.log 2>&1

# also prepare a file sorted by SNP ID with genotypes coded based on reference allele
# ref homo = 0, het = 0.5, alt homo = 1
tail -n+9 /mnt/research/qgg/resource/dgrp/all.jgil.exon.common.recode.vcf | perl -we 'while (<>) { chomp $_; @line = split /\t/, $_; for (my $i = 9; $i <= $#line; $i++) { @geno = split /:/, $line[$i]; if ($geno[0] eq "0/0") { $line[$i] = 0; } elsif ($geno[0] eq "0/1") { $line[$i] = 0.5 } elsif ($geno[0] eq "1/1") { $line[$i] = 1; } else { $line[$i] = "." } } print $line[2], "\t", join("\t", @line[9..$#line]), "\n"; }' | sort -k1,1 > /mnt/research/qgg/resource/dgrp/all.jgil.exon.common.sorted.tsv

# get ID list
sed -n '8p' /mnt/research/qgg/resource/dgrp/all.jgil.exon.common.recode.vcf | cut -f 10- | tr '\t' '\n' > /mnt/research/qgg/resource/dgrp/all.jgil.line.id.txt

# 2. loop through all fastq files
# calculate distance for each DGRP_BaseLine fastq
# ==========================================================================

mkdir /mnt/gs18/scratch/users/tansuxu/dgrp/chechID
cd /mnt/ufs18/scratch/huangw53/dgrp-rna/DGRP_BaseLine
ls *fastq.gz | sed 's/\.fastq\.gz//g' > /mnt/gs18/scratch/users/tansuxu/dgrp/chechID/DGRP_BaseLine.txt

cd /mnt/gs18/scratch/users/tansuxu/dgrp/chechID
split -l 400 DGRP_BaseLine.txt
mv xaa DGRP_BaseLine.1.txt
mv xab DGRP_BaseLine.2.txt
mv xac DGRP_BaseLine.3.txt

while read line
do
  sample=`echo $line | awk '{print $1}'`
  sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,fastq=/mnt/ufs18/scratch/huangw53/dgrp-rna/DGRP_BaseLine/$sample.fastq.gz,tmp=/mnt/gs18/scratch/users/tansuxu/dgrp/chechID/tmp,mem=16G,dgrpExonGeno=/mnt/research/qgg/resource/dgrp/all.jgil.exon.common.sorted.tsv,dgrpLineID=/mnt/research/qgg/resource/dgrp/all.jgil.line.id.txt,outdir=/mnt/research/qgg/dgrp-age-exp/checkID,head=100000,vcf=/mnt/research/qgg/resource/dgrp/all.jgil.exon.common.recode.vcf --out=/mnt/research/qgg/dgrp-age-exp/checkID/log/checkID.$sample.out --error=/mnt/research/qgg/dgrp-age-exp/checkID/log/checkID.$sample.err /mnt/research/qgg/dgrp-age-exp/sbatch/checkID.sbatch
done < DGRP_BaseLine.1.txt

while read line
do
  sample=`echo $line | awk '{print $1}'`
  sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,fastq=/mnt/ufs18/scratch/huangw53/dgrp-rna/DGRP_BaseLine/$sample.fastq.gz,tmp=/mnt/gs18/scratch/users/tansuxu/dgrp/chechID/tmp,mem=16G,dgrpExonGeno=/mnt/research/qgg/resource/dgrp/all.jgil.exon.common.sorted.tsv,dgrpLineID=/mnt/research/qgg/resource/dgrp/all.jgil.line.id.txt,outdir=/mnt/research/qgg/dgrp-age-exp/checkID,head=100000,vcf=/mnt/research/qgg/resource/dgrp/all.jgil.exon.common.recode.vcf --out=/mnt/research/qgg/dgrp-age-exp/checkID/log/checkID.$sample.out --error=/mnt/research/qgg/dgrp-age-exp/checkID/log/checkID.$sample.err /mnt/research/qgg/dgrp-age-exp/sbatch/checkID.sbatch
done < DGRP_BaseLine.2.txt

while read line
do
  sample=`echo $line | awk '{print $1}'`
  sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,fastq=/mnt/ufs18/scratch/huangw53/dgrp-rna/DGRP_BaseLine/$sample.fastq.gz,tmp=/mnt/gs18/scratch/users/tansuxu/dgrp/chechID/tmp,mem=16G,dgrpExonGeno=/mnt/research/qgg/resource/dgrp/all.jgil.exon.common.sorted.tsv,dgrpLineID=/mnt/research/qgg/resource/dgrp/all.jgil.line.id.txt,outdir=/mnt/research/qgg/dgrp-age-exp/checkID,head=100000,vcf=/mnt/research/qgg/resource/dgrp/all.jgil.exon.common.recode.vcf --out=/mnt/research/qgg/dgrp-age-exp/checkID/log/checkID.$sample.out --error=/mnt/research/qgg/dgrp-age-exp/checkID/log/checkID.$sample.err /mnt/research/qgg/dgrp-age-exp/sbatch/checkID.sbatch
done < DGRP_BaseLine.3.txt

# check if run successfully for samples in three files

cd /mnt/gs18/scratch/users/tansuxu/dgrp/chechID

while read line
do
  sample=`echo $line | awk '{print $1}'`
  err=`wc -l /mnt/research/qgg/dgrp-age-exp/checkID/log/checkID.$sample.err | awk '{print $1}'`
  if [ $err -ne 0 ]
    then
	echo $sample
  fi
done < <(cat DGRP_BaseLine.1.txt DGRP_BaseLine.2.txt DGRP_BaseLine.3.txt)

# calculate distance for each 3WK fastq
# ==========================================================================

cd /mnt/ufs18/scratch/huangw53/DGRP_3WK
ls *fastq.gz | sed 's/\.fastq\.gz//g' > /mnt/gs18/scratch/users/tansuxu/dgrp/chechID/DGRP_3WK.txt

cd /mnt/gs18/scratch/users/tansuxu/dgrp/chechID
split -l 450 DGRP_3WK.txt
mv xaa DGRP_3WK.1.txt
mv xab DGRP_3WK.2.txt

while read line
do
  sample=`echo $line | awk '{print $1}'`
  sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,fastq=/mnt/ufs18/scratch/huangw53/DGRP_3WK/$sample.fastq.gz,tmp=/mnt/gs18/scratch/users/tansuxu/dgrp/chechID/tmp,mem=16G,dgrpExonGeno=/mnt/research/qgg/resource/dgrp/all.jgil.exon.common.sorted.tsv,dgrpLineID=/mnt/research/qgg/resource/dgrp/all.jgil.line.id.txt,outdir=/mnt/research/qgg/dgrp-age-exp/checkID,head=100000,vcf=/mnt/research/qgg/resource/dgrp/all.jgil.exon.common.recode.vcf --out=/mnt/research/qgg/dgrp-age-exp/checkID/log/checkID.$sample.out --error=/mnt/research/qgg/dgrp-age-exp/checkID/log/checkID.$sample.err /mnt/research/qgg/dgrp-age-exp/sbatch/checkID.sbatch
done < DGRP_3WK.1.txt

while read line
do
  sample=`echo $line | awk '{print $1}'`
  sbatch --export=env=/mnt/research/qgg/dgrp-age-exp/resource.env,fastq=/mnt/ufs18/scratch/huangw53/DGRP_3WK/$sample.fastq.gz,tmp=/mnt/gs18/scratch/users/tansuxu/dgrp/chechID/tmp,mem=16G,dgrpExonGeno=/mnt/research/qgg/resource/dgrp/all.jgil.exon.common.sorted.tsv,dgrpLineID=/mnt/research/qgg/resource/dgrp/all.jgil.line.id.txt,outdir=/mnt/research/qgg/dgrp-age-exp/checkID,head=100000,vcf=/mnt/research/qgg/resource/dgrp/all.jgil.exon.common.recode.vcf --out=/mnt/research/qgg/dgrp-age-exp/checkID/log/checkID.$sample.out --error=/mnt/research/qgg/dgrp-age-exp/checkID/log/checkID.$sample.err /mnt/research/qgg/dgrp-age-exp/sbatch/checkID.sbatch
done < DGRP_3WK.2.txt

# check if run successfully for samples in two files

cd /mnt/gs18/scratch/users/tansuxu/dgrp/chechID

while read line
do
  sample=`echo $line | awk '{print $1}'`
  err=`wc -l /mnt/research/qgg/dgrp-age-exp/checkID/log/checkID.$sample.err | awk '{print $1}'`
  if [ $err -ne 0 ]
    then
	echo $sample
  fi
done < <(cat DGRP_3WK.1.txt DGRP_3WK.2.txt DGRP_3WK.3.txt)
