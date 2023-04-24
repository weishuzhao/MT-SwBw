#!/bin/bash
#SBATCH -J test_clust
#SBATCH -p cpu
#SBATCH --ntasks-per-node=40
#SBATCH -o Oerr/%x-%j.out
#SBATCH -e Oerr/%x-%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=Hwrn.aou@sjtu.edu.cn
#SBATCH -N 1
:<<!EOF!
 * @Date: 2022-05-29 20:06:24
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-06-29 14:22:14
 * @FilePath: /2021_09-MT10kSW/workflow/gene_annot/gene_clust.sh
 * @Description:
!EOF!


###############################################################################
##                                                                           ##
##  Base function of SJTU pi                                                 ##
##                                                                           ##
###############################################################################

    function scratch() {
        echo $SCRATCH/`relapwd $1 ~`
    }


    function clean2old() {
        declare dir=""

        for dir in $@
        do
            if [ -d ${dir} ]
            then
                clean2old  ${dir}_old
                #/bin/rm -r ${dir}_old
                mv  ${dir} ${dir}_old
            fi
            #mkdir   ${dir}
        done
    }


###############################################################################
##                                                                           ##
##  classic mmseq methods                                                    ##
##                                                                           ##
###############################################################################


function mmseq_() {
    conda_ activate python36

    declare protein=$1 \
            THREAD=${2:-1}

    declare  clust_dir=`dirname  ${protein}`
    declare clust_base=`basename ${protein}`
            clust_base="${clust_base%%.faa}"

    declare clust_out=${clust_dir}/${clust_base}
    declare clust_dir=${clust_dir}/clust
    /bin/rm -rf `scratch ${clust_dir}`
    mkdir -p ${clust_dir} `scratch ${clust_dir}`
    declare DB=${clust_dir}/${clust_base}

    mmseqs createdb ${protein} ${DB}

    mmseqs cluster ${DB} ${DB}_clu_100 `scratch ${DB}` --cov-mode 1 -c 1 --min-seq-id 1 `#-k 10` --threads ${THREAD}
    mmseqs createtsv ${DB} ${DB} ${DB}_clu_100 ${clust_out}-clu_100.tsv --threads ${THREAD}
    mmseqs createsubdb ${DB}_clu_100 ${DB} ${DB}_100

    mmseqs cluster ${DB}_100 ${DB}_clu `scratch ${DB}` `#--cov-mode 1` -c 0.9 --min-seq-id 0.95 `#-k 10` --threads ${THREAD}
    mmseqs createtsv ${DB}_100 ${DB}_100 ${DB}_clu ${clust_out}-clu.tsv --threads ${THREAD}
    mmseqs createsubdb ${DB}_clu ${DB}_100 ${DB}_clu_rep
    mmseqs convert2fasta ${DB}_clu_rep ${clust_out}-clu_rep.faa

    mmseqs createseqfiledb ${DB} ${DB}_clu ${DB}_clu_seq --threads ${THREAD}
    mmseqs result2flat ${DB} ${DB} ${DB}_clu_seq ${clust_out}-clu_seq.faa.clu

}


source ~/.conda_init
