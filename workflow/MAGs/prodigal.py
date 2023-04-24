# -*- coding: utf-8 -*-
"""
 * @Date: 2022-04-10 21:02:42
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-08-19 22:41:02
 * @FilePath: /2021_09-MT10kSW/workflow/MAGs/prodigal.py
 * @Description:
"""

# first get high quality MAGs and reference MAGs as well

from pathlib import Path

import click
import pandas as pd
from Bio import SeqIO
from joblib import Parallel, delayed
from PyLib.PyLibTool.file_info import verbose_import
from workflow.utils.gene_predict_cut import prodigal
from workflow.utils.file_path import FilePath, SnakeConfigGet


logger = verbose_import(__name__, __doc__)
snake_config = SnakeConfigGet({})
file_path = FilePath()


def prodigal_faa_to_genome(prodigal_prefix: Path):
    genome = prodigal_prefix.name
    gene_map: list[tuple] = []
    with open(f"{prodigal_prefix}.faa") as fi:
        for line in fi:
            if line.startswith(">"):
                values = line[1:].split(" # ")
                gene, start, end, strand_, values_ = values
                strand = "+" if strand_ == "1" else "-"
                partial = dict(i.split("=") for i in values_.split(";") if i).get(
                    "partial", ""
                )
                gene_map.append((genome, gene, start, end, strand, partial))
    return pd.DataFrame(
        gene_map, columns=["genome", "gene", "start", "end", "strand", "partial"]
    )


def main(dir2MAGs, genome_dir, gene_dir, gene_map_path, threads):
    def main_thread1(genome_path: Path):
        genomeid = genome_path.name
        genome_prefix = genome_dir / genomeid

        genome = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
        SeqIO.write(genome.values(), genome_prefix, "fasta-2line")

        annot_prefix = prodigal(genome_prefix, gene_dir)
        return annot_prefix

    annot_prefixs: list[Path] = Parallel(threads, verbose=20)(
        delayed(main_thread1)(Path(genome)) for genome in (dir2MAGs)
    )

    # annot_prefixs = [
    #     Path("04_bin/03_collection/gene") / i.name
    #     for i in Path("04_bin/03_collection/genome").glob("*.fa")
    # ]
    gene_map = pd.concat(
        [prodigal_faa_to_genome(annot_prefix) for annot_prefix in annot_prefixs]
    )
    gene_map.to_csv(gene_map_path, index=False)

    return 0


def collect_dir2MAGs(Stdb: pd.DataFrame):
    return Stdb["genome_path"]


@click.command()
@click.option("--loglevel", default="INFO", type=str, help="set level of logger")
@click.option(
    "--Stdb", default=file_path.results("Stdb"), type=Path, help="output Stdb"
)
@click.option(
    "--genome-dir",
    default=file_path.all_bins("collect_genome"),
    type=Path,
    help="collect genomes",
)
@click.option(
    "--gene-dir",
    default=file_path.all_bins("collect_gene"),
    type=Path,
    help="collect genes",
)
@click.option(
    "--gene-map",
    default=file_path.results("gene_map.csv"),
    type=Path,
    help="genome gene map",
)
@click.option("-t", "--threads", default=snake_config.THREADS)
def run(
    loglevel: str,
    stdb: Path,
    genome_dir: Path,
    gene_dir: Path,
    gene_map: Path,
    threads: int,
):
    logger.setLevel(level=loglevel.upper())  # info

    Stdb = pd.read_csv(stdb)
    dir2MAGs = collect_dir2MAGs(Stdb)
    genome_dir.mkdir(parents=True, exist_ok=True)
    gene_dir.mkdir(parents=True, exist_ok=True)

    gene_map.parent.mkdir(parents=True, exist_ok=True)

    logger.warning(">>> job start")
    state = main(dir2MAGs, genome_dir, gene_dir, gene_map, int(threads))
    logger.warning(">>> job finish")
    if state == 0:
        logger.info("success!")


if __name__ == "__main__":
    run()
