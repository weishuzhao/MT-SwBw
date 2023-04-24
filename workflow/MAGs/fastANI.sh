#!/bin/bash
#SBATCH -J fastANI
#SBATCH -p cpu
#SBATCH --ntasks-per-node=40
#SBATCH -o fastANI-j-%j.out
#SBATCH -e fastANI-j-%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=Hwrn.aou@sjtu.edu.cn
#SBATCH -N 1
:<<!EOF!
 * @Date: 2022-04-17 00:27:10
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-06-25 12:09:14
 * @FilePath: /2021_09-MT10kSW/workflow/MAGs/fastANI.sh
 * @Description:
    Stdb[(Stdb.Completeness >= 50) & (Stdb.Contamination <= 10)][["genome"]].merge(Bdb)["location"].to_csv("all_path.txt", index=False, header=False)
!EOF!


echo "" > workflow/MAGs/MT_MAGs.tsv
for i in `awk -v FS=',' '{if (NR==1) next}; {print $1}' Stdb.csv| grep -v "^1"`
do
    find ../ -name $i >> workflow/MAGs/MT_MAGs.tsv
done
echo "" > workflow/MAGs/ME_MAGs.tsv
for i in `awk -v FS=',' '{if (NR==1) next}; {print $1}' Stdb.csv| grep "^1"`
do
    find ../ -name $i >> workflow/MAGs/ME_MAGs.tsv
done


    echo $(date +%F%n%T)
    fastANI \
        --ql workflow/MAGs/MT_MAGs.tsv --rl workflow/MAGs/ME_MAGs.tsv \
        -o workflow/MAGs/MT2ME_ANI.out --matrix \
        -t $SLURM_NTASKS --fragLen 1000

    echo $(date +%F%n%T)
    fastANI \
        --ql workflow/MAGs/ME_MAGs.tsv --rl workflow/MAGs/MT_MAGs.tsv \
        -o workflow/MAGs/ME2MT_ANI.out --matrix \
        -t $SLURM_NTASKS --fragLen 1000
