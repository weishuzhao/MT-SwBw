#!/bin/bash
#SBATCH -J tree
#SBATCH -p cpu
#SBATCH --ntasks-per-node=40
#SBATCH -o log/drep/%x_%j.out
#SBATCH -e log/drep/%x_%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=Hwrn.aou@sjtu.edu.cn
#SBATCH -N 1
:<<!EOF!
 * @Date: 2022-07-01 09:36:42
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-09-08 13:38:31
 * @FilePath: /2021_09-MT10kSW/workflow/remote/slurm/tree.sh
 * @Description:
!EOF!
set -e && echo "$0 $*" >&2 && source ~/.conda_init
module unload xalt

set -vx
conda_ activate smk


snakemake -c $SLURM_NTASKS \
    --shadow-prefix $SCRATCH/`relapwd . ~` \
    `# trees_itol_annots` \
    04_bin/03_reference/genome_info.csv \
    fetchMG_to_trees -p
