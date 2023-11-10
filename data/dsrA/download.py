# -*- coding: utf-8 -*-
"""
 * @Date: 2023-11-09 21:54:56
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-11-10 14:44:52
 * @FilePath: /2021_09-MT10kSW/data/dsrA/download.py
 * @Description:
"""
# """

from typing import TextIO

from Bio import Entrez, Phylo, SeqIO

Entrez.email = "hwrn.aou@sjtu.edu.cn"
RETMAX = 100


def download_ncbi(accs: list[str], fo: TextIO):
    GCx_accs: list[str] = []
    for acc in accs:
        if acc.startswith("GCA_") or acc.startswith("GCF_"):
            GCx_accs.append(acc)
        else:
            with Entrez.efetch(db="Protein", id=acc, rettype="fasta") as handle:
                for line in handle:
                    print(line, end="", file=fo)
    return GCx_accs


def extract_genome_gene_ko_faa(ko2gene: dict[str, str]):
    from workflow.utils.file_path import file_path

    all_clu_faa = file_path.all_bins("collect_annot") / "all-clu_rep.faa"

    gene_annots = pd.read_csv(file_path.results("all_gene_annots.csv"))
    putate_dsrA_gene_annots: pd.DataFrame = gene_annots.pipe(
        lambda df: df[df["ko"].apply(lambda x: x in ko2gene)]
    )

    ko2faa = {
        ko: file_path.data / "dsrA" / f"putate_{gene}_gene.faa"
        for ko, gene in ko2gene.items()
    }

    for i, df in putate_dsrA_gene_annots.groupby("ko"):
        SeqIO.write(
            [
                record
                for record in SeqIO.parse(all_clu_faa, "fasta")
                if record.id in df["all"].values
            ],
            ko2faa[str(i)],
            "fasta-2line",
        )
    return ko2faa


if __name__ == "__main__":
    import pandas as pd

    with Entrez.einfo() as handle:
        record: dict = Entrez.read(handle)  # type: ignore

    genes = pd.read_csv("data/dsrA/dsrA.tsv", sep="\t").assign(
        acc=lambda df: df["Accession"].apply(lambda i: i.split(".")[0])
    )
    label_genes = genes.groupby("Clade")["acc"].agg(list).to_dict()
    with open("data/dsrA/ref_dsrA_gene.faa", "w") as fo:
        GCx_accs = download_ncbi([i for i in genes["acc"]], fo)
        # GCx_accs = [i for i in genes["acc"] if i.startswith("GC")]

    genome_dsrA_faa_file = extract_genome_gene_ko_faa({"K11180": "dsrA"})["K11180"]
    """{bash}
    cat data/dsrA/ref_dsrA_gene.faa data/dsrA/putate_dsrA_gene.faa > data/dsrA/dsrA.faa
    # ~/Templates/meta-snakemake-minimal$ python -m snakemake \
    #     /home/hwrn/Work/2021_09-MT10kSW/data/dsrA/dsrA.afa.trimal.iqtree \
    #     --profile profile/local/ -rp -c 10
    """

    from Bio import Phylo
    from Bio.Phylo import Newick

    tree: Newick.Tree = Phylo.read("data/dsrA/dsrA.afa.trimal.contree", "newick")
    tree_names: dict[str, str] = {
        i.name.split(".")[0]: i.name for i in tree.get_terminals()
    }
    tree.root_with_outgroup(
        [tree_names[i] for i in label_genes["Moorella"] if i not in GCx_accs]
    )
    assert sorted(
        [
            i.name
            for i in tree.common_ancestor(
                *[tree_names[i] for i in label_genes["SRB"] if i not in GCx_accs]
            ).get_terminals()
        ]
    ) == sorted([tree_names[i] for i in label_genes["SRB"] if i not in GCx_accs])
    sob_ca: Newick.Clade = tree.common_ancestor(
        [tree_names[i] for i in label_genes["SOB"] if i not in GCx_accs]
    )
    sob_gene = [
        i.name
        for i in sob_ca.get_terminals()
        if i.name.split(".")[0] not in label_genes["SOB"]
    ]
    sra_gene = [
        i.name
        for i in tree.common_ancestor(
            [tree_names[i] for i in label_genes["SRA"] if i not in GCx_accs]
        ).get_terminals()
        if i.name.split(".")[0] not in label_genes["SRA"]
    ]
    srb_gene = [
        i.name
        for i in tree.common_ancestor(
            [tree_names[i] for i in label_genes["SRB"] if i not in GCx_accs]
        ).get_terminals()
        if i.name.split(".")[0] not in label_genes["SRB"]
    ]
