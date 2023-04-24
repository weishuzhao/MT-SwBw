# -*- coding: utf-8 -*-
"""
 * @Date: 2022-05-04 18:58:25
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-07-15 19:45:45
 * @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/MAG_ko_tpm.py
 * @Description:
"""

from pathlib import Path
from typing import Literal

import click
import pandas as pd
from PyLib.PyLibTool.file_info import verbose_import
from workflow.utils.file_path import FilePath

from workflow.remote.gene_annot import get_gmodule


logger = verbose_import(__name__, __doc__)
file_path = FilePath()


GENE_ABD_METHOD = Literal["count", "rpb", "tpm"]


def load_gene_abd(
    gene_count, method: GENE_ABD_METHOD = "rpb", melt=False, threshold_0=0.001
) -> pd.DataFrame:
    if isinstance(gene_count, pd.DataFrame):
        gene_count_raw = gene_count
    else:
        gene_count_raw = pd.read_csv(gene_count, index_col=0, sep="\t", comment="#")

    if method == "count":
        gene_abd = gene_count_raw.iloc[:, 5:]
    elif method == "rpb":
        gene_abd = gene_count_raw.pipe(lambda df: df.iloc[:, 5:].T / df.iloc[:, 4]).T
    elif method == "tpm":
        gene_abd = load_gene_abd(gene_count_raw, "rpb", threshold_0=threshold_0).pipe(
            lambda df: df * 1e6 / df.agg("sum")
        )
    else:
        raise ValueError("unknown method")

    gene_abd_cut: pd.DataFrame = gene_abd.pipe(lambda df: df * (df > threshold_0))

    if melt:
        return (
            gene_abd_cut.reset_index()
            .pipe(
                lambda df: df.melt(
                    id_vars=["Geneid"],
                    value_vars=df.columns,
                    var_name="bam",
                    value_name=method,
                )
            )
            .pipe(lambda df: df[df[method] > 0])
        )
    else:
        return gene_abd


def get_all_gene_tpm(
    all_gene_ko: pd.DataFrame,
    genome_map: pd.DataFrame,
    gene_abd: pd.DataFrame,
):
    all_gene_abd = all_gene_ko.merge(
        genome_map, left_index=True, right_index=True
    ).merge(gene_abd, left_index=True, right_on="Geneid")
    return all_gene_abd


def main(
    gene_annots_path: Path,
    gene_map_path: Path,
    gene_count_path: Path,
    gene_tpm_path: Path,
    genomeko_path: Path,
    gmodule_path: Path,
):
    gene_annots = pd.read_csv(gene_annots_path, index_col=0)
    gene_map = pd.read_csv(gene_map_path, index_col=1)["genome"]
    # considering 1 read maximum mapping to 150 reads (rpb = 0.00667),
    #  a rpb < 0.001 means reads mapped to a gene (1. longer than ~900 bp) and (2. at most 1/6 of gene are recognized)
    gene_tpm = load_gene_abd(gene_count_path, "tpm", melt=True, threshold_0=0)  # 0.001

    # gene_rpb.assign(
    #    zone=lambda df: pd.cut(df["rpb"], [i * 0.001 for i in range(20)])
    # ).groupby("zone")["rpb"].agg(len)
    all_gene_tpm = (
        gene_annots.merge(gene_map, left_index=True, right_index=True)
        .merge(gene_tpm, left_index=True, right_on="Geneid")
        .assign(
            layer=lambda df: df["bam"].apply(
                lambda x: Path(x).name[:-4].split("..", 1)[1]
            )
        )
    )

    ko_genome_bam_tpm = (
        all_gene_tpm.groupby(["ko", "genome", "layer"])["tpm"].agg("sum").reset_index()
    )
    ko_genome_bam_tpm.to_csv(gene_tpm_path, index=False)

    bamko_tpm = ko_genome_bam_tpm.pivot_table("tpm", "ko", "layer", sum, 0)
    bamko_tpm.to_csv(genomeko_path)

    bmodule = get_gmodule(bamko_tpm)
    bmodule.to_csv(gmodule_path)


@click.command()
@click.option("--loglevel", default="INFO", type=str, help="set level of logger")
@click.option(
    "--gene-map",
    default=file_path.results("gene_map.csv"),
    type=Path,
    help="genome gene map",
)
@click.option(
    "--gene-annots",
    default=file_path.results("all_gene_annots.csv"),
    type=Path,
    help="",
)
@click.option(
    "--gene-count",
    default=file_path.bins_tpm("Wtdb.gff.count"),
    type=Path,
    help="",
)
@click.option(
    "--gene-tpm",
    default=file_path.bins_tpm("Wtdb_gene_tpm.csv"),
    type=Path,
    help="",
)
@click.option(
    "--genomeko", default=file_path.results("Wtdb_ko_tpm.csv"), type=Path, help=""
)
@click.option(
    "--gmodule", default=file_path.results("Wtdb_module_tpm.csv"), type=Path, help=""
)
def run(
    loglevel: str,
    gene_map: Path,
    gene_annots: Path,
    gene_count: Path,
    gene_tpm: Path,
    genomeko: Path,
    gmodule: Path,
):
    logger.setLevel(level=loglevel.upper())  # info

    gene_tpm.parent.mkdir(parents=True, exist_ok=True)
    genomeko.parent.mkdir(parents=True, exist_ok=True)
    gmodule.parent.mkdir(parents=True, exist_ok=True)

    logger.warning(">>> job start")
    state = main(
        gene_annots,
        gene_map,
        gene_count,
        gene_tpm,
        genomeko,
        gmodule,
    )
    logger.warning(">>> job finish")
    if state == 0:
        logger.info("success!")


if __name__ == "__main__":
    run()
