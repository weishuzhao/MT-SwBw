"""
 * @Date: 2022-06-28 20:33:42
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-09-14 21:26:24
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
        r="workflow/others/draw_fig2.r",
    output:
        fig_relative_nmds_class=file_path.figs("fig2_relative_nmds_class"),
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript {input.r} \
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
        div_raw=file_path.otus("abundance.csv"),
        r="workflow/others/draw_supp_fig1.r",
    output:
        fig_class_16s=file_path.figs("supp.fig1_class_16s"),
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript {input.r} \
            {input.div_raw} \
            {output.fig_class_16s}
        """


rule draw_supp_fig3:
    input:
        genome_abds=genome_abds,
        Stdb=Stdb,
        Wtdb=Wtdb,
        r="workflow/others/draw_supp_fig3.r",
    output:
        fig_coverm_relabd_signif_N=file_path.figs("supp.fig3_coverm_relabd_signif_N"),
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript {input.r} \
            {input.genome_abds} \
            {output.fig_coverm_relabd_signif_N}
        """


rule draw_supp_fig2:
    input:
        genome_abds=genome_abds,
        Stdb=Stdb,
        Wtdb=Wtdb,
        otu=file_path.otus("otu.tsv"),
        r="workflow/others/draw_supp_fig2.r",
    output:
        fig_share_water_sed=file_path.figs("supp.fig2_share_water_sed"),
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript {input.r} \
            {input.genome_abds} \
            {input.otu} \
            {output.fig_share_water_sed}
        """


rule draw_supp_fig5:
    input:
        genome_abds=genome_abds,
        gene_tpm=file_path.results("gene_ko_tpm.csv"),
        r="workflow/others/draw_supp_fig5.r",
    output:
        fig_ko_tpm_signif_site=file_path.figs("supp.fig5_ko_tpm_signif_site"),
    shadow:
        "shallow"
    conda:
        "R4.1"
    shell:
        """
        Rscript {input.r} \
            {output.fig_ko_tpm_signif_site}
        """
