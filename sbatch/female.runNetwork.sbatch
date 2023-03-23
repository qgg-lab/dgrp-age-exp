#!/bin/bash

# Job name:
#SBATCH --job-name=runNetwork
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=20
#
# Memory per node:
#SBATCH --mem=64G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "fastq.gzdays-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=00:30:00
# ============================================================

module load GCC/11.2.0  OpenMPI/4.1.1 R/4.1.2

Rscript R/runNetwork.R qg/female.sig.gene.exp.RData figureData/female.wgcna.RData > log/female.wgcna.Rout 2>&1

echo $(date) done.main.process