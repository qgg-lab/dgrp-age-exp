#!/bin/bash

# Job name:
#SBATCH --job-name=stringtie
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=9
#
# Memory per node:
#SBATCH --mem=16G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "fastq.gzdays-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=04:00:00
# ============================================================

source "$env"

awk -F "\"" '{print $2"\t"$0}' $gtf | sort -k1,1 | join -v 1 -t $'\t' - /mnt/research/qgg/dgrp-age-exp/gtf/rRNA.7SLRNA.gene.id | cut -f 2- | sort -k1,1 -k2,2n > $dir/lift.no.rRNA.no.7SLRNA.gtf

$STRINGTIE --rf -p 8 -e -B -G $dir/lift.no.rRNA.no.7SLRNA.gtf -o $dir/out.gtf -A $dir/gene_abund.tab $bam 2> $dir/stringtie.eB.log

echo $(date +"%m-%d-%Y-%T"):info: done.main.process
