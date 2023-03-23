#!/bin/bash

# Job name:
#SBATCH --job-name=mergeBamInd
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
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00

# ============================================================

source "$env"

echo $bamInd | tr "," "\n" | awk '{print "'$tmp'/"$1"/map.bam"}' | xargs $SAMTOOLS merge $dir/sort.bam 2> $dir/log/merge.log

$SAMTOOLS index $dir/sort.bam

echo $(date) done.main.process