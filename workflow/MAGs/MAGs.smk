"""
 * @Date: 2022-05-29 16:19:51
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-04-24 21:25:59
 * @FilePath: /2021_09-MT10kSW/workflow/MAGs/MAGs.smk
 * @Description:
"""
Stdb = (file_path.results("Stdb"),)
Wtdb = (file_path.results("Wtdb"),)
fastani_reference = (file_path.all_bins("fastani_reference.csv"),)

collect_genome = (file_path.all_bins("collect_genome"),)
collect_gene = (file_path.all_bins("collect_gene"),)

genome_abds = file_path.results("Wtdb.relative_abundance.tsv")


rule recover_drep:
    input:
        drep_out=file_path.all_bins("drep"),
        taxon_out=file_path.all_bins("gtdbtk"),
    output:
        Stdb=Stdb,
        Wtdb=Wtdb,
        fastani_reference=fastani_reference,
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python39

        python -m workflow.utils.collect_drep_taxon \
            --drep {input.drep_out} \
            --taxon {input.taxon_out} \
            --Stdb {output.Stdb} \
            --Wtdb {output.Wtdb} \
            --fastani-reference {output.fastani_reference} \
            --depth
        """


# region genome gene annotation
rule genome_prodigal:
    input:
        Stdb=Stdb,
    output:
        collect_genome=directory(collect_genome),
        collect_gene=directory(collect_gene),
        gene_map=file_path.results("gene_map.csv"),
    threads: THREADS
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python39

        python -m workflow.MAGs.prodigal \
            --Stdb {input.Stdb} \
            --genome-dir {output.collect_genome} \
            --gene-dir   {output.collect_gene} \
            --gene-map   {output.gene_map} \
            --threads {threads}
        """


rule genome_all_faa:
    input:
        collect_gene=collect_gene,
    output:
        all_faa=file_path.all_bins("collect_annot") / "all.faa",
    threads: THREADS
    shell:
        """
        cat {input.collect_gene}/*.faa > {output.all_faa}
        """


rule genome_2_ko:
    input:
        all_100=file_path.all_bins("collect_annot") / "all-clu_100.tsv",
        all_clu=file_path.all_bins("collect_annot") / "all-clu.tsv",
        annot_kofam=file_path.all_bins("collect_annot") / "all-clu_rep-kofam.tsv",
        annot_eggnog=file_path.all_bins("collect_annot") / "all-clu_rep-eggnog.tsv",
        annot_tcdb=file_path.all_bins("collect_annot") / "all-clu_rep-tcdb.tsv",
        annot_mantis=file_path.all_bins("collect_annot") / "all-clu_rep-mantis.tsv",
    output:
        gene_annots=file_path.results("all_gene_annots.csv"),
    params:
        annot_prefix=file_path.all_bins("collect_annot") / "all-clu_rep-.*.tsv",
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python39

        python -m workflow.remote.gene_annot \
            --all-100      {input.all_100} \
            --all-clu      {input.all_clu} \
            --annot-prefix {params.annot_prefix} \
            --gene-annots  {output.gene_annots} \
        """


rule get_gmodule:
    input:
        gene_annots=file_path.results("all_gene_annots.csv"),
        gene_map=file_path.results("gene_map.csv"),
    output:
        gmodule=file_path.results("gmodule.csv"),
    run:
        import pandas as pd
        from remote.gene_annot import get_gmodule

        all_gene_annots = pd.read_csv(input.gene_annots, index_col=0)
        gene_map = pd.read_csv(input.gene_map, index_col=1)["genome"]
        genomeko = (
            all_gene_annots.merge(gene_map, left_index=True, right_index=True)
            .dropna()
            .pivot_table(
                values="ko", index="ko", columns="genome", aggfunc=len, fill_value=0
            )
        )
        gmodule = get_gmodule(genomeko)
        gmodule.to_csv(output.gmodule)


rule prodigal2gff:
    input:
        genome_tdb=file_path.results("{genome_tdb}.csv"),
        collect_genome=collect_genome,
        collect_gene=collect_gene,
    output:
        all_fa=file_path.bins_tpm("{genome_tdb}.fa"),
        all_gff=file_path.bins_tpm("{genome_tdb}.gff"),
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python39

        python -m workflow.MAGs.prodigal2gff \
            --genome-tdb     {input.genome_tdb} \
            --collect-genome {input.collect_genome} \
            --collect-gene   {input.collect_gene} \
            --all-fa         {output.all_fa} \
            --all-gff        {output.all_gff}
        """


# endregion genome gene annotation


rule genomeko_class_view:
    input:
        Stdb="results/MAGs/Stdb.csv",
        Wtdb="results/MAGs/Wtdb.csv",
        genomeko="results/MAGs/genomeko.csv",
    output:
        classko_annot="results/MAGs/classko_annot.csv",
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python39

        python workflow/MAGs/genomeko_class_view.py \
            {input.Stdb} {input.Wtdb} \
            {input.genomeko} \
            {output.classko_annot}
        """


# region figures
rule ko_class_heatmap_N:
    input:
        genome_abds=genome_abds,
        Stdb=Stdb,
        Wtdb=Wtdb,
    output:
        fig_ko_class_heatmap_N=file_path.figs("ko_class_heatmap_N"),
    message:
        "percentages of N metabolic genes across class level"
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript workflow/MAGs/ko_class_heatmap_N.r \
            {input.genome_abds} \
            {output.fig_ko_class_heatmap_N}
        """


rule ko_relationship_N:
    input:
        genome_abds=genome_abds,
        Stdb=Stdb,
        Wtdb=Wtdb,
    output:
        fig_ko_relationship_N=file_path.figs("ko_relationship_N_{domain_spec}"),
    params:
        domain_spec="{domain_spec}",
    wildcard_constraints:
        domain_spec="All|Archaea|Bacteria",
    message:
        "co-occurrence network of selected metabolic genes"
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript workflow/MAGs/ko_relationship_N.r \
            {input.genome_abds} \
            {output.fig_ko_relationship_N} \
            {params.domain_spec}
        """


# endregion figures
