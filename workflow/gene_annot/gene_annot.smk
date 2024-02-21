"""
 * @Date: 2022-05-28 23:04:47
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-02-19 13:20:19
 * @FilePath: /2021_09-MT10kSW/workflow/gene_annot/gene_annot.smk
 * @Description:
"""


rule contig_all_faa:
    input:
        collect_gene=[
            file_path.gene(f"{site}", suffix=".faa") for site in site_layer_dict
        ],
    output:
        all_faa=file_path.env_genes("all.faa"),
    threads: THREADS
    shell:
        """
        cat {input.collect_gene} > {output.all_faa}
        """


rule get_contig_rep_mantis:
    input:
        annot_msntis=file_path.env_genes("all-clu_rep-mantis.tsv"),
