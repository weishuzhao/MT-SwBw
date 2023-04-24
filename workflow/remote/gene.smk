contig = str(file_path.contig("{site}", "_cut"))
faa = file_path.gene("{site}", suffix=".faa")
fna = file_path.gene("{site}", suffix=".fna")
gff = file_path.gene("{site}", suffix=".gff")


rule contig_gene_prodigal:
    input:
        contig = contig,
    output:
        faa = faa,
        fna = fna,
        gff = gff,
    params:
        gene_prefix = str(file_path.gene("{site}", suffix=""))
    threads: 1
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python39

        python -m workflow.utils.gene_predict_cut \
            --genome {input.contig} \
            --gene-prefix {params.gene_prefix} \
            --min-faa-len 33 \
            --filter \
            --meta
        """


rule contig_gene_abundance:
    input:
        bam = file_path.bam("{site}", "{layer}"),
        gff = gff,
    output:
        gcount        = file_path.gene("{site}", suffix=".count.txt"),
        count_summary = file_path.gene("{site}", suffix=".count.txt.summary"),
    threads: THREADS-1
    shell:
        """
        source ~/.conda_init
        conda_ activate python36

        featureCounts \
            -a {input.gff} \
            -o {output.gcount} \
            -t CDS -g ID \
            -T {threads} \
            -p {input.bam}
        """


rule gene_clu:
    input:
        all_faa     = "{any}.faa",
    output:
        all_100     = "{any}-clu_100.tsv",
        all_clu_faa = "{any}-clu_rep.faa",
        all_clu     = "{any}-clu.tsv",
        clu_flat    = "{any}-clu_seq.faa.clu",
    threads: THREADS
    shell:
        """
        source workflow/remote/gene_clust.sh
        mmseq_ {input.all_faa} {threads}
        """


rule gene_annot_ko:
    input:
        all_clu_faa  = "{any}.faa",
    output:
        annot  = "{any}-{method}.tsv",
    threads: THREADS
    params:
        method = "{method}",
    wildcard_constraints:
        method = "kofam|eggnog|tcdb"
    shell:
        """
        source workflow/remote/gene_annot.sh
        annot {input.all_clu_faa} {threads} {params.method}
        """


rule gene_annot_mantis:
    input:
        all_clu_faa  = "{any}.faa",
    output:
        annot  = "{any}-{method}.tsv",
    threads: THREADS
    wildcard_constraints:
        method = "mantis"
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate test

        mantis \
            run \
            -mc ~/Data/Database2/MANTIS/MANTIS.cfg \
            -c {threads} \
            -i {input.all_clu_faa} \
            -o {input.all_clu_faa}-mantis

        mv {input.all_clu_faa}-mantis/consensus_annotation.tsv \
           {output.annot}
        """
