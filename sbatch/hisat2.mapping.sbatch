#!/bin/bash

# Job name:
#SBATCH --job-name=hisat2.mapping
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
#SBATCH --mem=8G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "fastq.gzdays-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=4:00:00
# ============================================================

source "$env"

$HISAT2 --rg ID:$fc.$lane --rg PL:illumina --rg PU:$fc.$lane --rg SM:$individual -p 8 --rna-strandness R --summary-file $dir/log/map.summary -x "$index" -U "$fastq" 2> $dir/log/map.log | $SAMTOOLS view -bS - 2> $dir/log/samview.log | $SAMTOOLS sort - -o $tmp/map.bam 2> $dir/log/samsort.log

echo $(date) done.main.process