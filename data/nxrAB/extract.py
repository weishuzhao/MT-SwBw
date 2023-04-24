# -*- coding: utf-8 -*-
"""
 * @Date: 2022-09-21 21:03:00
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-09-22 10:16:32
 * @FilePath: /2021_09-MT10kSW/data/nxrAB/extract.py
 * @Description:
    files:
        - ref genes from TrEMBL of uniprot have been download: data/nxrAB/{gene}.faa
            - narH: (gene:narH) NOT nitrite NOT "uncultured bacterium"
            - nxrB: (gene:nxrB) NOT nitrate NOT "uncultured bacterium"
            - narG: (gene:narG) NOT nitrite NOT "uncultured bacterium"
            - nxrA: (gene:nxrA) NOT nitrate NOT "uncultured bacterium"
        - genes of MAGs: 04_bin/03_collection/annot/all-clu_rep.faa
        - annotations of genes: results/MAGs/all_gene_annots.csv

    workflow:
        1. collect potential genes annotated as K00370 (narG/nxrA) and K00371 (narH/nxrB)
            - putate_nxrAB_gene_file = file_path.data / "nxrAB" / "putate_{gene}_gene.faa"
        2. merge these genes with rename
            - file_path.data / "nxrAB" / f"ref_{gene}.faa"
        3. assign MAG genes to reference genes (cross-blast? build tree?)
"""

from typing import Union
from Bio import SeqIO
import pandas as pd
from workflow.utils.file_path import file_path


def check_diamond_out(diamond_out, head_n=1):
    """
    >>> diamond_1to1(
    ...     f"putate_{gene}_gene.faa",
    ...     f"ref_{gene}.faa",
    ...     diamond_tsv,
    ... )
    >>> assign_gene_info: pd.DataFrame = (
    ...     check_diamond_out(diamond_tsv)
    ...     .pipe(lambda df: df[df["sseqid"].apply(lambda x: x.startswith("ref_"))])
    ...     .assign(gene=lambda df: df["sseqid"].apply(lambda x: x[4:].split("|")[0]))
    ...     .merge(ref_gene_infos, left_on="sseqid", right_on="id")
    ...     .drop("id", axis=1)
    ... )
    """
    if not isinstance(diamond_out, pd.DataFrame):
        diamond_out = pd.read_csv(
            diamond_out,
            sep="\t",
            names="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split(),
        )

    diamond_out_index = (
        diamond_out.sort_values("pident", ascending=False)
        .groupby("qseqid")
        .apply(lambda df: df.reset_index().drop("index", axis=1).reset_index())
        .reset_index(drop=True)
    )
    diamond_out_merge_filter = (
        diamond_out_index.iloc[:, :4]
        .merge(
            diamond_out_index.iloc[:, :4],
            left_on=["qseqid", "sseqid"],
            right_on=["sseqid", "qseqid"],
            suffixes=("_q", "_s"),
        )
        .drop(["qseqid_s", "sseqid_s"], axis=1)
        .rename({"qseqid_q": "qseqid", "sseqid_q": "sseqid"}, axis=1)
        .groupby("qseqid")
        .apply(lambda df: df)
    )
    return diamond_out_merge_filter.groupby("qseqid").head(head_n)


def assign_unipprt_fa_features(
    headline=">ref_narH|tr|X5L138|X5L138_9MYCO Nitrate reductase subunit beta OS=Mycolicibacterium mageritense DSM 44476 = CIP 104973 OX=1209984 GN=narH PE=4 SV=1",
):
    features: dict[str, str] = {}
    if headline.startswith(">"):
        headline = headline[1:]
    headline = headline.strip()
    seqid, remain = headline.split(maxsplit=1)
    features["id"] = seqid
    key = "desc"
    while "=" in remain:
        value_, remain = remain.split("=", 1)
        while value_[-1] == " " or remain[0] == " ":
            remain_, remain = remain.split("=", 1)
            value_ = value_ + "=" + remain_
        features[key], key = value_.rsplit(maxsplit=1)
    features[key] = remain
    return features


# region 1 collect
all_clu_faa = file_path.all_bins("collect_annot") / "all-clu_rep.faa"

gene_annots = pd.read_csv(file_path.results("all_gene_annots.csv"))
putate_nxrAB_gene_annots: pd.DataFrame = gene_annots.pipe(
    lambda df: df[df["ko"].apply(lambda x: x in ["K00370", "K00371"])]
).assign(gene=lambda df: df["ko"].apply(lambda x: "nxrA" if x == "K00370" else "nxrB"))

putate_nxrAB_gene_file = file_path.data / "nxrAB" / "putate_{gene}_gene.faa"

for i, df in putate_nxrAB_gene_annots.groupby("gene"):
    SeqIO.write(
        [
            record
            for record in SeqIO.parse(all_clu_faa, "fasta")
            if record.id in df["all"].values
        ],
        str(putate_nxrAB_gene_file).format(gene=i),
        "fasta-2line",
    )
# endregion 1 collect

# region 2 merge and rename
ref_gene_infos_dict: dict = {}
for homogenes in (("nxrA", "narG"), ("nxrB", "narH")):
    with open(file_path.data / "nxrAB" / f"ref_{homogenes[0]}.faa", "w") as fo:
        for gene in homogenes:
            with open(file_path.data / "nxrAB" / f"{gene}.faa") as fi:
                for line in fi:
                    if line.startswith(">"):
                        line = f">ref_{gene}|" + line[1:]
                        features = assign_unipprt_fa_features(line)
                        ref_gene_infos_dict[features["id"]] = features
                    fo.write(line)
ref_gene_infos = pd.DataFrame(ref_gene_infos_dict).T
# endregion 2 merge and rename

# region 3 crossblast
from PyLib.biotool.diamond import diamond_1to1

for gene in ("nxrA", "nxrB"):
    diamond_tsv = file_path.data / "nxrAB" / f"diamond_{gene}.tsv"
    diamond_1to1(
        file_path.data / "nxrAB" / f"putate_{gene}_gene.faa",
        file_path.data / "nxrAB" / f"ref_{gene}.faa",
        diamond_tsv,
    )

    assign_gene_info: pd.DataFrame = (
        check_diamond_out(diamond_tsv)
        .pipe(lambda df: df[df["sseqid"].apply(lambda x: x.startswith("ref_"))])
        .assign(gene=lambda df: df["sseqid"].apply(lambda x: x[4:].split("|")[0]))
        .merge(ref_gene_infos, left_on="sseqid", right_on="id")
        .drop("id", axis=1)
    )
    assign_gene_info.to_csv(
        file_path.data / "nxrAB" / f"checked_{gene}.csv", index=False
    )

# endregion 3 crossblast


# region 4 nxrAB in genome contents
# de-novo init
import pandas as pd
from workflow.utils.file_path import file_path
from workflow.remote.gene_annot import load_rep2all

assign_gene_info = pd.concat(
    [
        pd.read_csv(file_path.data / "nxrAB" / f"checked_{gene}.csv")
        for gene in ("nxrA", "nxrB")
    ]
)
gene_info_nxrAB = assign_gene_info.pipe(
    lambda df: df[df["gene"].apply(lambda x: x in ("nxrA", "nxrB"))]
)
Wtdb = pd.read_csv(file_path.results("Wtdb"))
all_100 = file_path.all_bins("collect_annot") / "all-clu_100.tsv"
all_clu = file_path.all_bins("collect_annot") / "all-clu.tsv"
rep2all = load_rep2all(all_100, all_clu)
all_gene_annots = pd.read_csv(file_path.results("all_gene_annots.csv"), index_col=0)

gene_map: pd.Series = pd.read_csv(file_path.results("gene_map.csv"), index_col=1)[
    "genome"
]

NOB_genomes = (
    rep2all.merge(gene_info_nxrAB, left_on="rep", right_on="qseqid")
    .merge(gene_map, left_on="all", right_index=True)[
        ["genome", "all", "gene", "desc", "OS", "GN", "index_q", "pident_q"]
    ]
    .rename({"all": "gene", "gene": "GENE"}, axis=1)
)

NOB_genomes_putate_narGH: pd.DataFrame = (
    gene_map.reset_index()
    .merge(NOB_genomes["genome"])
    .merge(all_gene_annots, left_on="gene", right_index=True)
    .merge(pd.Series(["K00370", "K00371"], name="ko"))
    .pipe(lambda df: df[1 ^ df.duplicated()])
)
NOB_genomes_putate_narGH.merge(NOB_genomes, how="left").sort_values("genome").merge(
    Wtdb
).to_csv(file_path.data / "nxrAB" / f"final_nxrAB.csv", index=False)

# endregion 4 nxrAB in genome contents
