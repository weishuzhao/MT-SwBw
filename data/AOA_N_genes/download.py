# -*- coding: utf-8 -*-
"""
 * @Date: 2024-02-27 17:52:35
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2024-02-29 14:10:27
 * @FilePath: /2021_09-MT10kSW/data/AOA_N_genes/download.py
 * @Description:
"""


import pandas as pd
from Bio import Phylo, SeqIO
from Bio.Phylo import Newick
from workflow.utils.file_path import file_path


def extract_genome_gene_ko_faa(
    ko2gene: dict[str, str],
    faa_file=file_path.all_bins("collect_annot") / "all.faa",
    ko_file=lambda gene: file_path.data / "AOA_N_genes" / f"{gene}_all.faa",
):
    gene_annots = pd.read_csv(file_path.results("all_gene_annots.csv"))
    putate_dsrA_gene_annots: pd.DataFrame = gene_annots.pipe(
        lambda df: df[df["ko"].apply(lambda x: x in ko2gene)]
    )

    ko2faa = {ko: ko_file(gene) for ko, gene in ko2gene.items()}

    for i, df in putate_dsrA_gene_annots.groupby("ko"):
        SeqIO.write(
            [
                record
                for record in SeqIO.parse(faa_file, "fasta")
                if record.id in df["all"].values
            ],
            ko2faa[str(i)],
            "fasta-2line",
        )
    return ko2faa


if __name__ == "__main__":

    genome_amoB_faa_file = extract_genome_gene_ko_faa(
        {
            "K10944": "amoA",
            "K10945": "amoB",
            "K00261": "gdhA",
            "K01428": "ureC",
            "K01905": "glnA",
            "K03320": "amt",
        },
        file_path.all_bins("collect_annot") / "all-clu_rep.faa",
        ko_file=lambda gene: file_path.data / "AOA_N_genes" / f"{gene}.faa",
    )
    genome_amoB_faa_file = extract_genome_gene_ko_faa(
        {
            "K00368": "nirK",
        },
        file_path.all_bins("collect_annot") / "all-clu_rep.faa",
        ko_file=lambda gene: file_path.data / "AOA_N_genes" / f"{gene}_1.faa",
    )
    SeqIO.write(
        (
            i
            for i in SeqIO.parse("data/AOA_N_genes/nirK_1.faa", "fasta")
            if len(i) >= 90
        ),
        "data/AOA_N_genes/nirK_ge90.faa",
        "fasta-2line",
    )
    faa_file = extract_genome_gene_ko_faa(
        {
            "K10944": "amoA",
            "K10945": "amoB",
            "K00261": "gdhA",
            "K00368": "nirK",
            "K01428": "ureC",
            "K01905": "glnA",
            "K03320": "amt",
        },
        ko_file=lambda gene: file_path.data / "AOA_N_genes" / f"{gene}_all.faa",
    )
    """{bash ~/Templates/meta-snakemake-minimal$}
    cat data/amoB/putate_amoB_gene.faa > data/AOA_N_genes/amoB.faa
    python -m snakemake \
        /home/hwrn/Work/2021_09-MT10kSW/data/AOA_N_genes/{amoA,amoB,gdhA,nirK,ureC,glnA,amt}.afa.trimal.iqtree \
        --profile profile/local/ -c 10
    """

    tree: Newick.Tree = Phylo.read("data/amoB/amoB.afa.trimal.contree", "newick")
    tree_names: dict[str, str] = {
        i.name.split("|")[0]: i.name for i in tree.get_terminals()
    }
    Phylo.draw(tree)
