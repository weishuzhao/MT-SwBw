
rule md5sum:
    input:
        file = "{any}"
    output:
        file = "{any}.md5sum"
    threads: 1
    shell:
        """
        md5sum {input.file} > {output.file}
        """

include:
    "trim.smk"
include:
    "assembly.smk"
include:
    "binning.smk"
include:
    "gene.smk"
include:
    "drep.smk"


if SLURM_ARRAY_TASK_ID >= 0:
    site = sorted(site_layer_dict)[SLURM_ARRAY_TASK_ID]
    print(site)
    rule stat_trim_one_site:
        input:
            npo = [file_path.stat_reads(f"{site}", f"{layer}", "nonpareil", ".npo") for layer in site_layer_dict[site]],
            NTUfull_abundance = [file_path.stat_reads(f"{site}", f"{layer}", "phyloFlash", ".csv") for layer in site_layer_dict[site]],
        shell:
            """
            ls {input}
            """

    rule trim_one_site:
        input:
            r1 = [file_path.trimmed_reads(site, layer, 1) for layer in site_layer_dict[site]],
            r2 = [file_path.trimmed_reads(site, layer, 2) for layer in site_layer_dict[site]],
        shell:
            """
            ls {input.r1} {input.r2}
            """

    rule assem_depth_one_site:
        input:
            contig = str(file_path.contig(f"{site}", "_cut")),
            jgi = file_path.contig(f"{site}", "_cut", "-jgi.depth")
        run:
            """
            ls -l {input.contig} {input.jgi}
            """

    rule DASTool_one_site:
        input:
            scf2bin = file_path.dastool(f"{site}", "_scaffolds2bin.tsv"),
            checkm = file_path.checkm_of(file_path.dastool_bins(f"{site}")),
        run:
            """
            ls -l {input.scf2bin} {input.checkm}
            """

else:
    rule assem_all_site:
        input:
            contigs = [str(file_path.contig(f"{site}", "_cut")) for site in site_layer_dict],
        run:
            """
            ls {input.contigs}
            """
