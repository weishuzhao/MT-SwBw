#!/bin/bash
#SBATCH -J map_coverm
#SBATCH -p cpu
#SBATCH --ntasks-per-node=40
#SBATCH -a 0-22
#SBATCH -o log/map_coverm/%A_%a.out
#SBATCH -e log/map_coverm/%A_%a.err
#SBATCH --mail-type=end
#SBATCH --mail-user=Hwrn.aou@sjtu.edu.cn
#SBATCH -N 1
:<<!EOF!
 * @Date: 2022-05-17 15:47:33
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-06-29 14:29:57
 * @FilePath: /2021_09-MT10kSW/workflow/remote/slurm/assem.sh
 * @Description:
!EOF!
set -e && echo "$0 $*" >&2 && source ~/.conda_init
module unload xalt

set -vx
conda_ activate smk


#snakemake -c $SLURM_NTASKS \
#    assem_depth_one_site


snakemake -c $SLURM_NTASKS \
    --shadow-prefix $SCRATCH/`relapwd . ~` \
    map_coverm_one_site
