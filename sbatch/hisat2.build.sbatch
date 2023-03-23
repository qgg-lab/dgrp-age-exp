#!/bin/bash

# Job name:
#SBATCH --job-name=hisat2.build
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=8
#
# Memory per node:
#SBATCH --mem=160G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "fastq.gzdays-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00
# ============================================================

source "$env"

$HISAT2EXON $dir/lift.gtf > $dir/lift.exons
$HISAT2SS $dir/lift.gtf > $dir/lift.splice_sites

$HISAT2BD -p 8 $dir/lift.fa --ss $dir/lift.splice_sites --exon $dir/lift.exons $dir/lift.genome

echo $(date +"%m-%d-%Y-%T"):info: done.main.process