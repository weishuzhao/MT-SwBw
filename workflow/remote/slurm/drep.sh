#!/bin/bash
#SBATCH -J tpm
#SBATCH -p cpu
#SBATCH --ntasks-per-node=40
#SBATCH -o log/drep/%x_%j.out
#SBATCH -e log/drep/%x_%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=Hwrn.aou@sjtu.edu.cn
#SBATCH -N 1
:<<!EOF!
 * @Date: 2022-05-17 15:47:33
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-09-05 11:09:19
 * @FilePath: /2021_09-MT10kSW/workflow/remote/slurm/drep.sh
 * @Description:
!EOF!
set -e && echo "$0 $*" >&2 && source ~/.conda_init
module unload xalt

set -vx
conda_ activate smk


snakemake -c $SLURM_NTASKS \
    --shadow-prefix $SCRATCH/`relapwd . ~` \
    results/MAGs/site_module.csv \
    results/MAGs/gmodule.csv \
    results/MAGs/gene_ko_tpm.csv \
    results/MAGs/Wtdb.relative_abundance.tsv \
    04_bin/04_tpm/Stdb.raw.relative_abundance.tsv \
    04_bin/04_tpm/Wtdb.raw.relative_abundance.tsv \
    04_bin/04_tpm/Wtdb.count -krp
