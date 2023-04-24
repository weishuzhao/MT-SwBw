# -*- coding: utf-8 -*-
"""
 * @Date: 2022-06-26 21:55:54
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-06-30 10:27:06
 * @FilePath: /2021_09-MT10kSW/workflow/MAGs/get_ref_genomes.py
 * @Description:
"""

from pathlib import Path

import click
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from PyLib.biotool.download import download_fna
from PyLib.PyLibTool.file_info import verbose_import
from workflow.MAGs.prodigal import prodigal, prodigal_faa_to_genome
from workflow.utils.file_path import FilePath, SnakeConfigGet
from PyLib.biotool.fna_msg import statistic_fna
from PyLib.reader.read_outputs import fasta

logger = verbose_import(__name__, __doc__)
snake_config = SnakeConfigGet({})
file_path = FilePath()


GTDBTK_DATA_PATH = Path("~/Data/Database2/GTDB/release202/").expanduser()
gtdb_taxonomy = pd.read_csv(
    GTDBTK_DATA_PATH / "taxonomy" / "gtdb_taxonomy.tsv",
    sep="\t",
    names=["reference", "taxonomy"],
)


def infer_references(classification):
    my_phylum = set(taxon.rsplit(";", 7 - 2)[0] for taxon in classification)
    rs = np.random.RandomState(804)
    selected = (
        gtdb_taxonomy.assign(
            p__=lambda df: df["taxonomy"].apply(
                lambda taxon: taxon.rsplit(";", 7 - 2)[0]
            )
        )
        .merge(pd.DataFrame({"p__": list(my_phylum)}))
        .assign(
            c__=lambda df: df["taxonomy"].apply(
                lambda taxon: taxon.rsplit(";", 7 - 3)[0]
            )
        )
        .groupby("c__")
        .apply(lambda df: df.sample(1, random_state=rs))
        .drop(["p__", "c__"], axis=1)
    )
    return selected


def collect_refs(genome_dir: Path, genomes: list[str]):
    genome_dir.mkdir(parents=True, exist_ok=True)
    ref_genomes: list[str] = Parallel(THREADS, verbose=20)(
        delayed(download_fna)(i, genome_dir, retry=3) for i in genomes
    )
    return [Path(i) for i in ref_genomes]


def annot_genome(genomes_dirs: list, gene_dir: Path, gene_map_path: Path):
    annot_prefixs: list[Path] = Parallel(THREADS, verbose=20)(
        delayed(prodigal)(Path(genome), gene_dir) for genome in (genomes_dirs)
    )
    gene_map = pd.concat(
        [prodigal_faa_to_genome(annot_prefix) for annot_prefix in annot_prefixs]
    )
    gene_map.to_csv(gene_map_path, index=False)


def main(
    Wtdb: pd.DataFrame,
    reference_genome_dir: Path,
    reference_genome_info_dir: Path,
    reference_gene_dir: Path,
    reference_gene_map: Path,
):
    Wtdb = pd.read_csv(file_path.results("Wtdb"))

    references_table = infer_references(Wtdb["classification"])
    references = references_table["reference"].apply(
        lambda GB_GCX: GB_GCX.split("_", 1)[1]
    )

    reference_genomes = collect_refs(reference_genome_dir, references)
    ref_idb = (
        pd.Series(reference_genomes)
        .apply(
            lambda genome_path: pd.Series(
                data=[genome_path, *statistic_fna(fasta(genome_path))],
                index=[
                    "genome_path",
                    "SeqNumbers",
                    "MaxLength",
                    "GenomeSize",
                    "GC",
                    "N50",
                    "L50",
                ],
            )
        )
        .assign(genome=lambda df: df["genome_path"].apply(lambda x: Path(x).name))
    )

    reference_genome_info = references_table.assign(
        genome=lambda df: df["reference"].apply(lambda x: x.split("_", 1)[1] + ".fna")
    ).merge(ref_idb)
    reference_genome_info.to_csv(reference_genome_info_dir, index=False)
    annot_genome(reference_genomes, reference_gene_dir, reference_gene_map)


THREADS = 10


@click.command()
@click.option("--loglevel", default="INFO", type=str, help="set level of logger")
@click.option(
    "--Wtdb", default=file_path.results("Wtdb"), type=Path, help="output Wtdb"
)
@click.option(
    "--ref-genome-info",
    default=file_path.all_bins("ref_genome_info.csv"),
    type=Path,
    help="infered reference genomes",
)
@click.option(
    "--ref-genome",
    default=file_path.all_bins("ref_genome"),
    type=Path,
    help="collect genomes",
)
@click.option(
    "--ref-gene",
    default=file_path.all_bins("ref_gene"),
    type=Path,
    help="collect genes",
)
@click.option(
    "--ref-gene-map",
    default=file_path.all_bins("ref_gene_map.csv"),
    type=Path,
    help="output ref genome gene map",
)
@click.option("-t", "--threads", default=snake_config.THREADS)
def run(
    loglevel: str,
    wtdb: Path,
    ref_genome: Path,
    ref_genome_info: Path,
    ref_gene: Path,
    ref_gene_map: Path,
    threads: int,
):
    logger.setLevel(level=loglevel.upper())  # info

    Wtdb = pd.read_csv(wtdb)

    ref_genome.mkdir(parents=True, exist_ok=True)
    ref_gene.mkdir(parents=True, exist_ok=True)

    ref_genome_info.parent.mkdir(parents=True, exist_ok=True)
    ref_gene_map.parent.mkdir(parents=True, exist_ok=True)

    global THREADS
    THREADS = threads

    logger.warning(">>> job start")
    state = main(
        Wtdb,
        ref_genome,
        ref_genome_info,
        ref_gene,
        ref_gene_map,
    )
    logger.warning(">>> job finish")
    if state == 0:
        logger.info("success!")


if __name__ == "__main__":
    run()
