#include:
#    "MAGs.smk"
#include:
#    "../MAGs_tpm/MAGs_tpm.smk"


rule get_ref_genomes:
    input:
        Wtdb = Wtdb,
    output:
        reference_genome_info = file_path.all_bins("ref_genome_info.csv"),
        reference_genome_dir  = directory(file_path.all_bins("ref_genome")),
        reference_gene_dir    = directory(file_path.all_bins("ref_gene")),
        reference_gene_map    = file_path.all_bins("ref_gene_map.csv"),
    threads: THREADS
    message: "check and DELETE the STUPID GCA_013288705.1.fna"
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python39

        python -m workflow.MAGs.get_ref_genomes \
            --Wtdb {input.Wtdb} \
            --ref-genome-info {output.reference_genome_info} \
            --ref-genome      {output.reference_genome_dir} \
            --ref-gene        {output.reference_gene_dir} \
            --ref-gene-map    {output.reference_gene_map} \
            --threads {threads}
        """


tree_prefix_expr = str(file_path.bins_tree("{bottom_or_slope}_mid_{min_exist_marker}"))
rule fetchMG_to_one_tree:
    input:
        Wtdb = Wtdb,
        collect_gene = collect_gene,
        reference_gene_dir = file_path.all_bins("ref_gene"),
        genome_abds = genome_abds,
    output:
        afa    = f"{tree_prefix_expr}.afa",
        trimal = f"{tree_prefix_expr}.afa.trimal",
        tree   = f"{tree_prefix_expr}.afa.trimal.fasttree",
    params:
        tree   = f"{tree_prefix_expr}",
        bottom_or_slope  = "{bottom_or_slope}",
        min_exist_marker = "{min_exist_marker}",
    threads: THREADS
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python39

        python -m workflow.MAGs.fetchMGtree \
            --Wtdb               {input.Wtdb} \
            --gene               {input.collect_gene} \
            --ref-gene           {input.reference_gene_dir} \
            --relative-abundance {input.genome_abds} \
            --tree               {params.tree} \
            --min-exist-marker   {params.min_exist_marker} \
            -t {threads} \
            --{params.bottom_or_slope}
        """


rule fetchMG_to_trees:
    input:
        trees = [
            f"{tree_prefix_expr}.afa.trimal.fasttree".format(
                bottom_or_slope=bottom_or_slope, min_exist_marker=min_exist_marker
            )
            for bottom_or_slope in ["all", "slope", "bottom"]
            for min_exist_marker in [20]
        ],
    shell:
        """
        ls {input}
        """


rule trees_itol_annots:
    input:
        Wtdb_abds = file_path.results(f"Wtdb.relative_abundance.tsv"),
    output:
        table2itol = directory(file_path.bins_tree("table2itol")),
    message:
        "Warning: may modify file '04_bin/04_tree/table2itol.csv' manually"
    shell:
        """
            set +u
            source workflow/utils/.conda_init
            conda activate python39

        python -c "from workflow.MAGs.fetchMGtree import collect_itol_table, file_path; collect_itol_table('{input.Wtdb_abds}')"

        bash workflow/MAGs/table2itol.sh
        """
