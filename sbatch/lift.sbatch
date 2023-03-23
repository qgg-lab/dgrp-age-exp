#!/bin/bash

# Job name:
#SBATCH --job-name=lift
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
#SBATCH --mem=4G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "fastq.gzdays-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00
# ============================================================

source "$env"

mkdir $line
cd $line

# prepare bed file
perl $PREPLINESEQBED --vcf $vcf --dgrp $line > $line.bed 2> prepLineSeqBED.log

# lift
perl $LIFTGTF --gtf $gtf --bed $line.bed --fasta $REFFA --gtfout lift.gtf --varout lift.bed --fastaout lift.fa 2> lift.log

echo $(date) done.main.process