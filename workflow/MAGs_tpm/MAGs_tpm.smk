"""
 * @Date: 2022-05-28 23:04:47
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-04-24 21:29:25
 * @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/MAGs_tpm.smk
 * @Description:
"""
# include:
#    "../MAGs/MAGs.smk"


# region remote generate coverm bam
## region remote
rule bwa_index:
    input:
        all_fa="{any}.fa",
    output:
        genome_amb="{any}.fa.amb",
        genome_ann="{any}.fa.ann",
        genome_bwt="{any}.fa.bwt",
        genome_pac="{any}.fa.pac",
        genome_sa="{any}.fa.sa",
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate coverM

        bwa index \
            {input.all_fa}
        """


rule map_as_coverm_bwa:
    input:
        r1=str(file_path.trimmed_reads("{site}", "{layer}", 1)),
        r2=str(file_path.trimmed_reads("{site}", "{layer}", 2)),
        genome_amb=file_path.bins_tpm("{genome_tdb}.fa.amb"),
        genome_ann=file_path.bins_tpm("{genome_tdb}.fa.ann"),
        genome_bwt=file_path.bins_tpm("{genome_tdb}.fa.bwt"),
        genome_pac=file_path.bins_tpm("{genome_tdb}.fa.pac"),
        genome_sa=file_path.bins_tpm("{genome_tdb}.fa.sa"),
    output:
        bam=file_path.bins_tpm("{genome_tdb}", "{site}", "{layer}"),
    threads: THREADS
    params:
        genome=str(file_path.bins_tpm("{genome_tdb}.fa")),
    wildcard_constraints:
        genome_tdb="Stdb|Wtdb",
    shadow:
        "shallow"
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate coverM

        bwa mem \
            -t {threads} \
            {params.genome} \
            {input.r1} {input.r2} \
        > tmp
        cat tmp \
        | samtools sort \
            -T "bwa" \
            -l0 -@ {threads} \
        | samtools view \
            -@ {threads} -b \
            -o {output.bam}
        """


rule coverm_genome_raw:
    input:
        bams=[
            file_path.bins_tpm("{genome_tdb}", f"{site}", f"{layer}")
            for site in site_layer_dict
            for layer in site_layer_dict[site]
        ],
        collect_genome=collect_genome,
    output:
        rela_abd=file_path.bins_tpm("{genome_tdb}.raw.relative_abundance.tsv"),
    log:
        file_path.bins_tpm("{genome_tdb}.raw.relative_abundance.log"),
    threads: THREADS
    wildcard_constraints:
        genome_tdb="Stdb|Wtdb",
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate coverM

        coverm genome \
            -b {input.bams} \
            -d {input.collect_genome} -x .fa \
            -m relative_abundance \
            -t {threads} \
            -o {output.rela_abd} \
            2>{log}
        """


rule coverm_genome:
    input:
        bams=[
            file_path.bins_tpm("{genome_tdb}", f"{site}", f"{layer}")
            for site in site_layer_dict
            for layer in site_layer_dict[site]
        ],
        collect_genome=collect_genome,
    output:
        rela_abd=file_path.bins_tpm("{genome_tdb}.relative_abundance.tsv"),
    threads: THREADS
    wildcard_constraints:
        genome_tdb="Stdb|Wtdb",
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate coverM

        coverm genome \
            -b {input.bams} \
            -d {input.collect_genome} -x .fa \
            -m relative_abundance \
            --min-read-aligned-length 50 \
            --min-read-percent-identity 0.99 \
            --min-covered-fraction 0.1 \
            --proper-pairs-only \
            -t {threads} \
            -o {output.rela_abd}
        """


rule coverm_genome_collect:
    input:
        Stdb=Stdb,
        rela_abd=file_path.bins_tpm("{genome_tdb}.relative_abundance.tsv"),
    output:
        rela_abd=file_path.results("{genome_tdb}.relative_abundance.tsv"),
    run:
        import pandas as pd

        Stdb = pd.read_csv(str(input.Stdb))
        relative_abundance_raw = pd.read_csv(input.rela_abd, sep="\t", index_col=0)

        relative_abundance = (
            relative_abundance_raw.reset_index()
            .merge(
                Stdb.assign(
                    Genome=lambda df: df["genome"].apply(
                        lambda x: x.rsplit(".fa", 1)[0]
                    )
                )[["Genome", "secondary_cluster"]]
            )
            .rename({"secondary_cluster": "Cluster"}, axis=1)
            .melt(
                id_vars=["Cluster", "Genome"],
                var_name="Bam",
                value_name="Relative abundance",
            )
            .groupby(["Cluster", "Bam"])[["Relative abundance"]]
            .apply("sum")
            .pipe(lambda df: df[df["Relative abundance"] > 0])
            .reset_index()
            .assign(
                Bam=lambda df: df["Bam"].apply(
                    lambda x: x.rsplit(" Relative Abundance (%)")[0]
                )
            )
            .assign(
                Bam=lambda df: df["Bam"].apply(
                    lambda x: x[:-8] if x.endswith(".fq.gz") else x
                )
            )
            .assign(Site=lambda df: df["Bam"].apply(lambda x: x.split("..")[1]))
            .assign(Layer=lambda df: df["Bam"].apply(lambda x: x.split("..")[2]))
            .drop("Bam", axis=1)
        )
        relative_abundance.to_csv(output.rela_abd, sep="\t", index=False)


## endregion remote


## region figs
rule coverm_relative_abundance_2d:
    input:
        genome_abds=genome_abds,
        Stdb=Stdb,
        Wtdb=Wtdb,
    output:
        coverm_relative_abundance_2d=file_path.figs("coverm_relabd_2d"),
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript workflow/MAGs_tpm/relative_abundance_nmds.r \
            {input.genome_abds} \
            {output.coverm_relative_abundance_2d}
        """


rule coverm_relative_abundance_taxon:
    input:
        genome_abds=genome_abds,
        Stdb=Stdb,
        Wtdb=Wtdb,
    output:
        coverm_relative_abundance_taxon=file_path.figs("coverm_relabd_bar-color"),
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript workflow/MAGs_tpm/relative_abundance_taxon.r \
            {input.genome_abds} \
            {output.coverm_relative_abundance_taxon}
        """


rule icamp_boot_summary_groups_pie:
    input:
        boot_summary_groups="results/MAGs/icamp.iCAMP.BootSummary.Groups.RPKM.csv",
    output:
        icamp_boot_summary_groups_pie=file_path.figs("icamp_boot_summary_groups_pie"),
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript workflow/MAGs_tpm/boot_summary_groups.r \
            {output.icamp_boot_summary_groups_pie} \
            {input.boot_summary_groups}
        """


rule icamp_boot_summary_groups_bar:
    input:
        boot_summary_groups="results/MAGs/icamp.iCAMP.BootSummary.Groups.RPKM.csv",
    output:
        icamp_boot_summary_groups_pie=file_path.figs("icamp_boot_summary_groups_bar"),
    conda:
        "R4.1"
    shell:
        """
        Rscript workflow/MAGs_tpm/boot_summary_groups_bar.r \
            {output.icamp_boot_summary_groups_pie} \
            {input.boot_summary_groups}
        """


## endregion figs


## region tabs
rule draw_tab_genome_prevalence:
    input:
        genome_abds=genome_abds,
        Stdb=Stdb,
        Wtdb=Wtdb,
    output:
        tab_relative_abundance=file_path.tabs("tab{any}_relative_abundance"),
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript workflow/MAGs_tpm/relative_abundance_table.r \
            {input.genome_abds} \
            {output.tab_relative_abundance}
        """


## endregion tabs
# endregion use coverm to describe MAG abundance and existence


# region calculate gene abundance
## region remote
rule bam_gff_count:
    input:
        bams=[
            file_path.bins_tpm("{genome_tdb}", f"{site}", f"{layer}")
            for site in site_layer_dict
            for layer in site_layer_dict[site]
        ],
        all_gff=file_path.bins_tpm("{genome_tdb}.gff"),
    output:
        all_count=file_path.bins_tpm("{genome_tdb}{if_multi}count"),
        all_count_summary=file_path.bins_tpm("{genome_tdb}{if_multi}count.summary"),
    params:
        m="{if_multi}",
    wildcard_constraints:
        if_multi=".multi.|\\.",
    threads: THREADS
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python36

        if [[ {params.m} == *multi* ]]
        then declare if_multi='-M --fraction'
        else declare if_multi=''
        fi

        featureCounts \
            -a {input.all_gff} \
            -o {output.all_count} \
            -t CDS -g ID \
            $if_multi \
            -T {threads} \
            -p {input.bams}
        """


rule get_all_gene_tpm:
    input:
        all_count=file_path.bins_tpm("Wtdb.count"),
        gene_annots=file_path.results("all_gene_annots.csv"),
        gene_map=file_path.results("gene_map.csv"),
    output:
        gene_ko_tpm=file_path.results("gene_ko_tpm.csv"),
        site_module=file_path.results("site_module.csv"),
    run:
        import pandas as pd
        from MAGs_tpm.MAG_ko_tpm import load_gene_abd, get_gmodule

        gene_annots = pd.read_csv(input.gene_annots, index_col=0)
        gene_map = pd.read_csv(input.gene_map, index_col=1)["genome"]
        gene_tpm = load_gene_abd(input.all_count, "tpm", melt=True, threshold_0=0)
        all_gene_tpm = (
            gene_annots.merge(gene_map, left_index=True, right_index=True)
            .merge(gene_tpm, left_index=True, right_on="Geneid")
            .assign(
                Layer=lambda df: df["bam"].apply(lambda x: x.split("..", 1)[1][:-4])
            )
        )
        gene_ko_tpm = (
            all_gene_tpm.groupby(["ko", "genome", "Layer"])["tpm"]
            .agg("sum")
            .reset_index()
        )
        gene_ko_tpm.to_csv(output.gene_ko_tpm, index=False)

        bamko_tpm = gene_ko_tpm.pivot_table("tpm", "ko", "Layer", sum, 0)
        bmodule = get_gmodule(bamko_tpm)
        bmodule.to_csv(output.site_module)


## endregion remote


## region figs
rule taxonko_venn:
    input:
        genome_abds=genome_abds,
        Wtdb=Wtdb,
        gene_tpm=file_path.results("gene_ko_tpm.csv"),
        site_module=file_path.results("site_module.csv"),
        r="workflow/MAGs_tpm/taxonko_venn.r",
    output:
        fig_mag_venn=file_path.figs("mag_venn"),
    message:
        "Venn of MAG taxonomy and ko"
    conda:
        "R4.1"
    shell:
        """
        Rscript {input.r} \
            {input.genome_abds} \
            {output.fig_mag_venn}
        """


rule ko_tpm_MAG_N_bar:
    input:
        genome_abds=genome_abds,
        gene_tpm=file_path.results("gene_ko_tpm.csv"),
        r="workflow/MAGs_tpm/ko_tpm_MAG_N_bar.r",
    output:
        fig_ko_tpm_MAG_N_bar=file_path.figs("ko_tpm_MAG_N_bar"),
    shadow:
        "shallow"
    conda:
        "R4.1"
    shell:
        """
        Rscript {input.r} \
            {input.genome_abds} \
            {output.fig_ko_tpm_MAG_N_bar}
        """


rule ko_tpm_MAG_N_pie:
    input:
        gene_tpm=file_path.results("gene_ko_tpm.csv"),
        r="workflow/MAGs_tpm/ko_tpm_MAG_N_pie.r",
    output:
        fig_ko_tpm_MAG_N_pie=file_path.figs("ko_tpm_MAG_N_pie"),
    shadow:
        "shallow"
    conda:
        "R4.1"
    shell:
        """
        Rscript {input.r} \
            {output.fig_ko_tpm_MAG_N_pie}
        """


## endregion figs


## region tabs
rule draw_tab_ko_tpm_genome_N_kayer_ko:
    input:
        genome_abds=genome_abds,
        Stdb=Stdb,
        Wtdb=Wtdb,
    output:
        tab_ko_tpm_genome_N_kayer_ko=file_path.tabs(
            "tab{any}_tab_ko_tpm_genome_N_kayer_ko"
        ),
    conda:
        "R4.1"
    shell:
        """
        Rscript workflow/MAGs_tpm/ko_tpm_MAG_N_table.r \
            {input.genome_abds} \
            {output.tab_ko_tpm_genome_N_kayer_ko}
        """


## endregion tabs
# endregion calculate gene abundance
