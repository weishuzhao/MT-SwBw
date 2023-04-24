phyloflash_rarefy = "workflow/reads_diversity/div.phyloflash.rarefy.csv"
ko_reads_zone = ("workflow/gene_annot/ko_reads_zone.csv",)
module_reads_zone = ("workflow/gene_annot/module_reads_zone.csv",)


rule collect_phyloflash:
    input:
        NTUfull_abundances=[
            file_path.stat_reads(f"{site}", f"{layer}", "phyloFlash", ".csv")
            for site in site_layer_dict
            for layer in site_layer_dict[site]
        ],
    output:
        div_raw=file_path.otus("phyloFlash_raw"),
    run:
        import pandas as pd

        pd.concat(
            [
                pd.read_csv(file, names=["OTU", "Reads"]).assign(File=file)
                for file in input.NTUfull_abundances
            ]
        ).to_csv(output.div_raw, index=False)


rule rarefy_phyloFlash:
    input:
        div_raw=file_path.otus("phyloFlash_raw"),
    output:
        fig_phylpflash_unannot=file_path.figs("phylpflash_unannot"),
        tab_phyloflash_rarefy=phyloflash_rarefy,
    message:
        "show mOTU annotated level and rarefy it"
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate R4.1

        Rscript workflow/reads_diversity/rarefy_phyloFlash.r \
            {output.fig_phylpflash_unannot} \
            {output.tab_phyloflash_rarefy}
        """
