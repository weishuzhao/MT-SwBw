"""
 * @Date: 2022-06-28 20:33:42
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-04-23 17:12:47
 * @FilePath: /2021_09-MT10kSW/workflow/others/others.smk
 * @Description:
"""


rule download_gebco_bathymetry_data:
    output:
        asc="data/GEBCO_28_Jun_2022_d37548b61e36/gebco_2021_n12.5_s10.0_w140.5_e143.5.asc",
    shell:
        """
        ls {output.asc}
        """


rule draw_map:
    input:
        asc="data/GEBCO_28_Jun_2022_d37548b61e36/gebco_2021_n12.5_s10.0_w140.5_e143.5.asc",
    output:
        fig_sample_site_map=file_path.figs("fig1_sample_site_map"),
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript workflow/others/draw_map.r \
            {input.asc} \
            {output.fig_sample_site_map}
        """


rule draw_fig2:
    input:
        genome_abds=genome_abds,
        Stdb=Stdb,
        Wtdb=Wtdb,
    output:
        fig_relative_nmds_class=file_path.figs("fig2_relative_nmds_class"),
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript workflow/others/draw_fig2.r \
            {input.genome_abds} \
            {output.fig_relative_nmds_class}
        """


rule draw_fig4:
    input:
        gene_tpm=file_path.results("gene_ko_tpm.csv"),
    output:
        fig_ko_tpm_site_signif_N=file_path.figs("fig4_ko_tpm_site_signif_N"),
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript workflow/others/draw_fig4.r \
            {output.fig_ko_tpm_site_signif_N}
        """


rule draw_supp_fig1:
    input:
        div_raw=file_path.otus("phyloFlash_raw"),
    output:
        fig_phylpflash_unannot=file_path.figs("supp.fig1_phylpflash_unannot"),
    message:
        "show mOTU annotated level and rarefy it"
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript workflow/others/draw_supp_fig1.r \
            {output.fig_phylpflash_unannot}
        """


rule draw_supp_fig2:
    input:
        genome_abds=genome_abds,
        Stdb=Stdb,
        Wtdb=Wtdb,
    output:
        fig_coverm_relabd_signif_N=file_path.figs("supp.fig2_coverm_relabd_signif_N"),
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript workflow/others/draw_supp_fig2.r \
            {input.genome_abds} \
            {output.fig_coverm_relabd_signif_N}
        """


rule draw_supp_fig4:
    input:
        genome_abds=genome_abds,
        gene_tpm=file_path.results("gene_ko_tpm.csv"),
        r="workflow/others/draw_supp_fig4.r",
    output:
        fig_ko_tpm_signif_site=file_path.figs("supp.fig4_ko_tpm_signif_site"),
    shadow:
        "shallow"
    conda:
        "R4.1"
    shell:
        """
        Rscript {input.r} \
            {output.fig_ko_tpm_signif_site}
        """
