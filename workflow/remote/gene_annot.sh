#!/bin/bash
#SBATCH -J test_annot
#SBATCH -p small
#SBATCH --ntasks-per-node=20
#SBATCH -o Oerr/%x-%j.out
#SBATCH -e Oerr/%x-%j.err
#SBATCH --mail-type=end
#SBATCH --mail-user=Hwrn.aou@sjtu.edu.cn
#SBATCH -N 1
:<<!EOF!
 * @Date: 2022-01-24 10:46:24
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-06-23 10:08:05
 * @FilePath: /2021_09-MT10kSW/workflow/gene_annot/gene_annot.sh
 * @Description:

###############################################################################
##                                                                           ##
##  Function definition of this class                                        ##
##                                                                           ##
###############################################################################

def {program}_(protein, out_dir, out_base, THREAD):
    declare annot={program}

    declare    protein=\$1 \
             annot_dir=\$2 \
            annot_base=\$3 \
                THREAD=\$4

    asserts:
        protein == os.path.join(out_dir, out_base + ".faa")
        temp_dir.startswith("/scratch")

def annot(protein, THREAD):
    outline:
        1.  out_dir, out_base <- protein
        2.  {{program}_() for {program} in programs}

!EOF!


###############################################################################
##                                                                           ##
##  Base function of SJTU pi                                                 ##
##                                                                           ##
###############################################################################

    function scratch() {
        echo $SCRATCH/`relapwd $1 ~`
    }


    function py_() {
        echo `python -c "print($*)"`
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
##  Various wrapper function for annotation methods                          ##
##                                                                           ##
###############################################################################

function kofam_() {
    conda_ activate python36
    declare annot=kofam

    declare    protein=$1 \
             annot_dir=$2 \
            annot_base=$3 \
                THREAD=$4

    declare annot_out=${annot_dir}/${annot_base}-${annot}
    declare tmp_dir=`scratch ${annot_out}`-tmp
    clean2old ${tmp_dir}; mkdir -p ${tmp_dir}

    echo $(date +%F%n%T)
    exec_annotation \
        --profile=$KOFAM_DB/db/profiles/prokaryote.hal \
        --ko-list=$KOFAM_DB/db/ko_list \
        --tmp-dir ${tmp_dir} \
        -f detail-tsv -E 1e-5 \
        --cpu ${THREAD} \
        -o ${annot_out}.txt \
        ${protein} \
    && /bin/rm -r ${tmp_dir}

    cat ${annot_out}.txt \
    |grep '\*' > ${annot_out}.tsv

}


function eggnog_() {
    conda_ activate python39
    declare annot=eggnog

    declare    protein=$1 \
             annot_dir=$2 \
            annot_base=$3 \
                THREAD=$4

    declare annot_out=${annot_dir}/${annot_base}-${annot}
    declare tmp_dir=`scratch ${annot_out}`-tmp
    clean2old ${tmp_dir}; mkdir -p ${tmp_dir}

    echo $(date +%F%n%T)
    emapper.py \
        --data_dir $EGGNOG_DB \
        --output_dir ${annot_dir}  `# Where output files should be written` \
        --scratch_dir ${tmp_dir}   `# Write output files in a temporary scratch dir, move them to the final output dir when finished.` \
        --temp_dir ${tmp_dir}      `# Where temporary files are created.` \
        -m diamond \
        --cpu ${THREAD} \
        -o ${annot_base}-${annot} \
        -i ${protein}

    ln -sb ${annot_base}-${annot}.emapper.annotations \
        ${annot_out}.tsv

}


function tcdb_() {
    # TODO: use diamond blastp instead
    conda_ activate python39
    declare annot=tcdb

    declare    protein=$1 \
             annot_dir=$2 \
            annot_base=$3 \
                THREAD=$4

    declare annot_out=${annot_dir}/${annot_base}-${annot}

    echo $(date +%F%n%T)
    diamond makedb \
        --in $TCDB_FAA \
        -d ${annot_out}..tcdb \
        -p ${THREAD} --quiet
    diamond blastp \
        -d ${annot_out}..tcdb \
        -q ${protein} \
        -o ${annot_out}.tsv \
        --outfmt 6 \
        -e 1e-20 --id 30 --query-cover 70 --subject-cover 70 \
        --dbsize 1000000000 -p ${THREAD} \
        --sensitive
    rm ${annot_out}..tcdb*

}


function annot() {
    declare protein=$1  `# {sample.faa}` \
            THREAD=$2 \
            programs=${3:-kofam eggnog tcdb}

    declare  annot_dir=`dirname  ${protein}`
    declare annot_base=`basename ${protein}`
            annot_base="${annot_base%%.faa}"

    echo $(date +%F%n%T)
    for program in ${programs}
    do
        ${program}_  ${protein} ${annot_dir} ${annot_base} ${THREAD}
    done
    echo $(date +%F%n%T)

}


source ~/.conda_init
export KOFAM_DB=/lustre/home/acct-clsxx/clsxx/Data/Database2/KofamKOALA
export EGGNOG_DB=/lustre/home/acct-clsxx/clsxx/Data/Database2/eggNOG/emapperdb-5.0.2
export TCDB_FAA=/lustre/home/acct-clsxx/clsxx/Data/Database2/TCDB/tcdb-21_05_30.faa
#annot ${protein} $SLURM_NTASKS
