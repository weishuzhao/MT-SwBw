# -*- coding: utf-8 -*-
"""
 * @Date: 2022-07-15 19:30:40
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-07-17 09:46:10
 * @FilePath: /2021_09-MT10kSW/workflow/others/check_gene_genome_tpm.py
 * @Description:
    check:
        \\any i(gene) from a(genome), \\any x(sample) contains a: assert i in x
        \\exists a(genome) contains i(gene), \\any x(sample) contains a: assert a in x
"""
from pathlib import Path
from typing import Literal

import click
import pandas as pd
from PyLib.PyLibTool.file_info import verbose_import
from workflow.utils.file_path import FilePath

from workflow.MAGs_tpm.MAG_ko_tpm import load_gene_abd, get_all_gene_tpm
from workflow.remote.gene_annot import load_rep2all


logger = verbose_import(__name__, __doc__)
file_path = FilePath()


def load_rep_sample():
    gene_count = load_gene_abd(
        file_path.bins_tpm("Wtdb.gff.count"),
        "count",
        melt=True,
        threshold_0=0,
    ).assign(
        layer=lambda df: df["bam"].apply(lambda x: Path(x).name[:-4].split("..", 1)[1])
    )
    rep2all = load_rep2all(
        file_path.all_bins("collect_annot") / "all-clu_100.tsv",
        file_path.all_bins("collect_annot") / "all-clu.tsv",
    )
    rep_count: pd.DataFrame = (
        pd.merge(rep2all, gene_count, how="left", left_on="all", right_on="Geneid")
        .pivot_table(values="count", index="rep", columns="layer", aggfunc=sum)
        .melt(value_name="count", ignore_index=False)
        .reset_index()
        .dropna(subset=["count"])
    )
    Wtdb_abds: pd.DataFrame = (
        pd.read_csv(file_path.results("Wtdb.relative_abundance.tsv"), sep="\t")
        .melt(id_vars="Genome", var_name="site", value_name="relative_abundance")
        .assign(
            relative_abundance=lambda df: df["relative_abundance"].apply(
                lambda x: x or float("NaN")
            )
        )
        .dropna()
        .assign(genome=lambda df: df["Genome"].apply(lambda x: f"{x}.fa"))
        .assign(layer=lambda df: df["site"].apply(lambda x: x.split("..", 1)[1][:-23]))
        .pipe(lambda df: df[df["Genome"] != "unmapped"])
    )[["genome", "layer", "relative_abundance"]]
    gene_map = pd.read_csv(file_path.results("gene_map.csv"))
    rep_genome: pd.DataFrame = (
        pd.merge(rep2all, gene_map, left_on="all", right_on="gene")
        .pipe(lambda df: df[df["rep"] == df["all"]])
        .assign(len=lambda df: (df["end"] - df["start"] + 1))
    )[["rep", "genome", "partial", "len"]]

    rep_sample: pd.DataFrame = (
        pd.merge(rep_genome, Wtdb_abds, on="genome")[["rep", "layer", "len"]]
        .drop_duplicates()
        .merge(rep_count, how="outer", on=["rep", "layer"])
        .assign(miss=lambda df: df["len"].isna() * 2 + df["count"].isna())
    )
    return rep_sample


if __name__ == "__main__":
    rep_sample = load_rep_sample()
    rep_sample.pivot("rep", "layer", "miss")
    rep_sample.assign(
        wrong=lambda df: df["miss"].apply(
            lambda x: "" if x == 0 else "genome" if x == 1 else "ko"
        )
    )["wrong"].value_counts()
    # ko          9938986
    # ""          6247040
    # genome      1602312
    rep_len = (
        rep_sample[["rep", "len"]]
        .dropna(subset=["len"])
        .drop_duplicates()
        .rename({"len": "_len"}, axis=1)
    )
    (
        rep_sample.merge(rep_len)
        .assign(
            wrong=lambda df: df["miss"].apply(
                lambda x: "" if x == 0 else "genome" if x == 1 else "ko"
            )
        )
        .assign(rpb=lambda df: df["count"] / df["_len"])
        .groupby("wrong")["rpb"]
        .describe()
    )
    rep_sample_tpm: pd.DataFrame = (
        rep_sample.merge(rep_len)
        .assign(
            wrong=lambda df: df["miss"].apply(
                lambda x: "" if x == 0 else "genome" if x == 1 else "ko"
            )
        )
        .assign(rpb=lambda df: df["count"] / df["_len"])
        .pipe(
            lambda df: df.merge(
                df.groupby("layer")["rpb"].sum().rename("rpb_sample"),
                left_on="layer",
                right_index=True,
            )
        )
        .assign(tpm=lambda df: df["rpb"] / df["rpb_sample"] * 1e6)
    )
    rep_sample_tpm[["rep", "layer", "count", "rpb", "tpm", "miss", "wrong"]].to_csv(
        file_path.cache("rep_sample_tpm.csv"), index=False
    )
    # wrong      count      mean       std       min       25%       50%       75%         max
    #         6247040.0  0.041436  0.311519  0.000095  0.003481  0.008292  0.022084  180.339286
    # genome        0.0       NaN       NaN       NaN       NaN       NaN       NaN         NaN
    # ko      5797715.0  0.007514  0.099441  0.000039  0.001182  0.002347  0.005291   77.304598
