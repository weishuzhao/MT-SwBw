#!/bin/bash
:<<!EOF!
 * @Date: 2022-03-19 19:34:44
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-05-20 09:02:51
 * @FilePath: /2022_05-ZFMG-release/workflow/gene_annot/prot_id.sh
 * @Description:
!EOF!

function blastp_() {
    conda activate python39

    local dbfaa=$1 query=$2
    local threads=${3:-1}

    local db=${dbfaa}_db
    local prot_out=${query}.blastp.tsv

    diamond makedb \
        --in ${dbfaa} -d ${db} \
        -p ${threads} --quiet

    diamond blastp \
        -d ${db} \
        -q ${query} \
        -o ${prot_out} \
        -k 5 -e 0.00001 \
        --id 30 \
        --query-cover 50 --subject-cover 50 \
        -b 4 --dbsize 1000000000 \
        -p ${threads} --quiet

}

function blastn_() {
    conda activate python36

    local dbfa=$1 query=$2
    local threads=${3:-1}

    local db=${dbfa}_db
    local prot_out=${query}.blastn.tsv

    makeblastdb \
        -in ${dbfa} -out ${db} \
        -dbtype nucl

    blastn \
        -task megablast \
        -db ${db} \
        -query ${query} \
        -strand plus \
        -out ${prot_out} \
        -outfmt 6 \
        -num_alignments 5 \
        -evalue 0.00001 -perc_identity 70 \
        -qcov_hsp_perc 50 \
        -num_threads ${threads}

    echo "finish"
}


for format in faa fna; do

mkdir -p workflow/gene_annot/$format
    declare sequence=workflow/gene_annot/$format/ME.$format
    echo "" > $sequence
    samples=`awk -v FS='\t' '! /^#/ {if($3=="ZF"){print $1}}' 00_data/sample_meta.tsv`
    for sample in ${samples[@]}
    do
        annot_fa=$sample-megahit/03_annot/03_annot..$sample-megahit.$format
        awk -v FS=">k" \
            '{if($1==""){print ">'"${sample}|"'k"$2}else{print $0}}' \
            $annot_fa \
        >> $sequence
    done

    declare sequence=workflow/gene_annot/$format/MT.$format
    echo "" > $sequence
    samples=`awk -v FS='\t' '! /^#/ {if($3!="ZF"){print $1}}' 00_data/sample_meta.tsv`
    for sample in ${samples[@]}
    do
        annot_fa=$sample-megahit/03_annot/03_annot..$sample-megahit.$format
        awk -v FS=">k" \
            '{if($1==""){print ">'"${sample}|"'k"$2}else{print $0}}' \
            $annot_fa \
        >> $sequence
    done

done

declare MT=workflow/gene_annot/faa/MT.faa ME=workflow/gene_annot/faa/ME.faa
blastp_ $ME $MT $SLURM_NTASKS
blastp_ $MT $ME $SLURM_NTASKS

declare MT=workflow/gene_annot/fna/MT.fna ME=workflow/gene_annot/fna/ME.fna
blastn_ $ME $MT $SLURM_NTASKS
blastn_ $MT $ME $SLURM_NTASKS
