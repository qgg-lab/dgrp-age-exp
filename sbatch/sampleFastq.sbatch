#!/bin/bash

# Job name:
#SBATCH --job-name=sampleFastq
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
#SBATCH --time=00:30:00
# ============================================================

source "$env"

# sample reads
# ============================================================

gunzip -c $fastq > $tmp/$sample.tmp.fastq
$FASTQSAMPLE -n $head $tmp/$sample.tmp.fastq -o $tmp/$sample.tmp.sample
rm $tmp/$sample.tmp.fastq

echo $(date +"%m-%d-%Y-%T"):info: done.main.process
