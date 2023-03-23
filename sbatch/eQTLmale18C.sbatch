#!/bin/bash

# Job name:
#SBATCH --job-name=eQTL.male.18c
#
# Number of nodes needed for use case:
#SBATCH --nodes=1
#
# Tasks per node based on number of cores per node:
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=1
#
# Memory per node:
#SBATCH --mem=4G
#
# Wall clock limit (one of "minutes", "minutes:seconds", "hours:minutes:seconds", "fastq.gzdays-hours", "days-hours:minutes" and "days-hours:minutes:seconds"):
#SBATCH --time=10:00:00
# ============================================================

/mnt/research/qgg/software/plink-v1.90b6.18/plink --silent --bfile /mnt/research/qgg/dgrp-age-exp/eqtl/dgrp.common --pheno eqtl/male.18c.pheno --all-pheno --allow-extra-chr --allow-no-sex --assoc --pfilter 0.00001 --out eqtl/male.18c/eqtl

echo $(date) done.main.process