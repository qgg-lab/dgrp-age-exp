#! /bin/bash

# job array to run eQTL mapping
# ============================================================

# request resource and sbatch options
# ============================================================

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --time=10:00:00

#SBATCH --job-name=eqtl.female.25c
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --output=/mnt/research/qgg/dgrp-age-exp/log/eqtl.female.25c.JobArray%a.log
#SBATCH --error=/mnt/research/qgg/dgrp-age-exp/log/eqtl.female.25c.JobArray%a.err

/mnt/research/qgg/software/plink-v1.90b6.18/plink --silent --bfile /mnt/research/qgg/dgrp-age-exp/eqtl/dgrp.common --pheno /mnt/research/qgg/dgrp-age-exp/eqtl/permPheno/female.25c.perm"${SLURM_ARRAY_TASK_ID}".pheno --all-pheno --allow-extra-chr --allow-no-sex --assoc --pfilter 0.00001 --out /mnt/research/qgg/dgrp-age-exp/eqtl/perm/female.25c.perm"${SLURM_ARRAY_TASK_ID}"/eqtl.perm

echo $(date) done.main.process
