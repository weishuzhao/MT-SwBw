#!/bin/bash
#SBATCH -J binning
#SBATCH -p cpu
#SBATCH --ntasks-per-node=40
#SBATCH -a 0-16
#SBATCH -o log/binning/%A_%a.out
#SBATCH -e log/binning/%A_%a.err
#SBATCH --mail-type=end
#SBATCH --mail-user=Hwrn.aou@sjtu.edu.cn
#SBATCH -N 1
:<<!EOF!
 * @Date: 2022-05-17 15:47:33
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-08-22 23:34:03
 * @FilePath: /2021_09-MT10kSW/workflow/remote/slurm/binning.sh
 * @Description:
!EOF!
set -e && echo "$0 $*" >&2 && source ~/.conda_init
module unload xalt

set -vx
conda_ activate smk


snakemake -c $SLURM_NTASKS \
    --shadow-prefix $SCRATCH/`relapwd . ~` \
    DASTool_one_site -krp
