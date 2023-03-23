#!/bin/bash

# Job name:
#SBATCH --job-name=checkID
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=2
#
# Memory per node:
#SBATCH --mem=16G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "fastq.gzdays-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=0:30:00
# ============================================================

source "$env"
sample=`echo $fastq | awk -F "/" '{print $NF}' | sed 's/\.fastq\.gz//'`
mkdir $tmp/$sample

# sample reads
# ============================================================

gunzip -c $fastq > $tmp/$sample/tmp.fastq
$FASTQSAMPLE -n $head $tmp/$sample/tmp.fastq -o $tmp/$sample/tmp.sample

# hisat2 alignment
# ============================================================

$HISAT2 -p 2 -k 3 --rna-strandness R -x $HISAT2IDX -U $tmp/$sample/tmp.sample.fastq 2> $tmp/$sample/hisat2.map.log | $SAMTOOLS view -bS - 2> $tmp/$sample/sam.view.log | $SAMTOOLS sort - -o $tmp/$sample/hisat2.map.bam 2> $tmp/$sample/sam.sort.log

# add read group
# ============================================================

java -jar $PICARD \
  AddOrReplaceReadGroups \
  I=$tmp/$sample/hisat2.map.bam \
  O=$tmp/$sample/hisat2.map.rg.bam \
  RGID=sample RGLB=lib RGPL=ILLUMINA RGPU=x RGSM=sample > $tmp/$sample/addrg.log 2>&1

# GATK split N cigar
# ============================================================

java -Xmx"$mem" -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$tmp \
  -jar $GATK SplitNCigarReads \
  -R $REFFA \
  -I $tmp/$sample/hisat2.map.rg.bam \
  -O $tmp/$sample/hisat2.map.splitn.bam > $tmp/$sample/hisat2.map.splitn.log 2>&1

# count ASE
# ============================================================

java -Xmx"$mem" -XX:ParallelGCThreads=1 -Djava.io.tmpdir=$tmp \
  -jar $GATK ASEReadCounter \
  -R $REFFA \
  -I $tmp/$sample/hisat2.map.splitn.bam \
  -O $tmp/$sample/ase.tsv \
  -V $vcf > $tmp/$sample/ase.log 2>&1

# check against genotypes
# calculate the allele frequency (f) of alt/total in the RNA-Seq sample
# take the absolute value of f-geno against each line
# sum these differences up and divide by the total number of sites
# this can be used to measure fraction of difference in genotyps
# do this for each fastq file
# ============================================================

tail -n+2 $tmp/$sample/ase.tsv | awk '$8 > 0 {print $3"\t"$7/$8}' | sort -k1,1 | join -t $'\t' - $dgrpExonGeno | perl -we '$count = 0; $firstline = <>; chomp $firstline; @line = split /\t/, $firstline; @dist = (0) x ($#line - 1); for (my $i = 0; $i <= $#dist; $i++) { $dist[$i] += abs($line[$i+2] - $line[1]); } $count++; while (<>) { chomp $_; @line = split /\t/, $_; for (my $i = 0; $i <= $#dist; $i++) { $dist[$i] += abs($line[$i+2] - $line[1]); }; $count++; }; print join("\n", map { $_/$count } @dist), "\n";' | paste $dgrpLineID - > $outdir/$sample.dist.txt

rm -r $tmp/$sample

echo $(date +"%m-%d-%Y-%T"):info: done.main.process
