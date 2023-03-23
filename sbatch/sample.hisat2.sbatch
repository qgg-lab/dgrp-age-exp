#!/bin/bash

# Job name:
#SBATCH --job-name=sample.hisat2
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=18
#
# Memory per node:
#SBATCH --mem=32G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "fastq.gzdays-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=04:00:00
# ============================================================

source "$env"

$HISAT2 -p 16 -k 3 --rna-strandness R -x $HISAT2IDX -U $fastq 2> $dir/log/hisat2.map.log | $SAMTOOLS view -bS - 2> $dir/log/sam.view.log | $SAMTOOLS sort - -o $dir/hisat2.map.bam 2> $dir/log/sam.sort.log

echo $(date +"%m-%d-%Y-%T"):info: done.main.process
