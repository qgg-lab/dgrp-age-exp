#! /bin/bash
# ============================================================

# job array to simulate expression and correlation
# ============================================================

# request resource and sbatch options
# ============================================================

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

#SBATCH --job-name=male.sim.corr.18c
#SBATCH --array=1-1000
#SBATCH --cpus-per-task=1
#SBATCH --output=/mnt/research/qgg/dgrp-age-exp/log/male.sim.corr.18c.JobArray%a.log
#SBATCH --error=/mnt/research/qgg/dgrp-age-exp/log/male.sim.corr.18c.JobArray%a.log

cd /mnt/research/qgg/dgrp-age-exp
Rscript R/simCorr18C.R qg/male.sig.gene.exp.RData "${SLURM_ARRAY_TASK_ID}" simCorr/male.sim.corr.18c.perm"${SLURM_ARRAY_TASK_ID}".RData > log/male.sim.corr.18c.JobArray"${SLURM_ARRAY_TASK_ID}".log 2>&1

echo $(date) done.main.process
