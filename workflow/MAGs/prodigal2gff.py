# -*- coding: utf-8 -*-
"""
 * @Date: 2022-05-04 18:58:25
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-06-24 19:48:42
 * @FilePath: /2021_09-MT10kSW/workflow/MAGs_tpm/prodigal2gff.py
 * @Description:
"""

from pathlib import Path

import click
import pandas as pd
from PyLib.PyLibTool.file_info import verbose_import
from workflow.utils.file_path import FilePath


logger = verbose_import(__name__, __doc__)
file_path = FilePath()


def prodigal2gff(prodigal_faa_file: Path):
    with prodigal_faa_file.open() as pfain:
        for line in pfain:
            if line.startswith(">"):
                head, start, end, strand, other_ = line.strip().split(" # ")
                contigId = head[1:].rsplit("_", 1)[0]
                other = "ID=" + head[1:] + ";" + other_.split(";", 1)[1]
                yield (
                    contigId,
                    "Prodigal_v2.6.3-modify",
                    "CDS",
                    start,
                    end,
                    ".",
                    strand[:-1] or "+",
                    0,
                    other,
                )


def main(
    genomes: list[str],
    collect_genome: Path,
    all_fa: Path,
    collect_gene: Path,
    all_gff: Path,
):
    genome_paths = [collect_genome / genome for genome in genomes]
    assert all([path.is_file for path in genome_paths])
    with all_fa.open("w") as fo:
        for path in genome_paths:
            with path.open() as fi:
                fo.write(fi.read())
    #
    gene_paths = [collect_gene / f"{genome}.faa" for genome in genomes]
    assert all([path.is_file for path in gene_paths])
    with all_gff.open("w") as fo:
        for path in gene_paths:
            for item in prodigal2gff(path):
                print(*item, sep="\t", file=fo)
    #
    return 0


@click.command()
@click.option("--loglevel", default="INFO", type=str, help="set level of logger")
@click.option("--genome-tdb", default=file_path.results("Stdb.csv"), type=Path, help="")
@click.option(
    "--collect-genome", default=file_path.all_bins("collect_genome"), type=Path, help=""
)
@click.option(
    "--collect-gene", default=file_path.all_bins("collect_gene"), type=Path, help=""
)
@click.option("--all-fa", default=file_path.bins_tpm("Stdb.fa"), type=Path, help="")
@click.option("--all-gff", default=file_path.bins_tpm("Stdb.gff"), type=Path, help="")
def run(
    loglevel: str,
    genome_tdb: Path,
    collect_genome: Path,
    collect_gene: Path,
    all_fa: Path,
    all_gff: Path,
):
    logger.setLevel(level=loglevel.upper())  # info
    # loglevel = "INFO"
    # genome_tdb = Path("results/MAGs/Stdb.csv")
    # collect_genome = Path("04_bin/03_collection/genome")
    # collect_gene = Path("04_bin/03_collection/gene")
    # all_fa = Path("04_bin/04_tpm/Stdb.fa")
    # all_gff = Path("04_bin/04_tpm/Stdb.gff")

    genomes = pd.read_csv(genome_tdb, usecols=[0])["genome"]

    all_fa.parent.mkdir(parents=True, exist_ok=True)
    all_gff.parent.mkdir(parents=True, exist_ok=True)

    logger.warning(">>> job start")
    state = main(genomes, collect_genome, all_fa, collect_gene, all_gff)
    logger.warning(">>> job finish")
    if state == 0:
        logger.info("success!")


if __name__ == "__main__":
    run()
