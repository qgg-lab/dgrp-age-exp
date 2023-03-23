#!/bin/bash

# Job name:
#SBATCH --job-name=sample.stringtie
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

$STRINGTIE --rf -p 16 -G $gtf -o $dir/pooled.sample.gtf $dir/hisat2.map.bam


echo $(date +"%m-%d-%Y-%T"):info: done.main.process
