# -*- coding: utf-8 -*-
"""
 * @Date: 2021-09-12 21:45:24
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-06-29 10:10:51
 * @FilePath: /2021_09-MT10kSW/workflow/utils/collect_drep_taxon.py
 * @Description:
"""

from pathlib import Path

import click
import pandas as pd
from drep.WorkDirectory import WorkDirectory

from .file_path import FilePath
from PyLib.PyLibTool.file_info import verbose_import

logger = verbose_import(__name__, __doc__)
file_path = FilePath()


def join_Widb_gtdb(wd: WorkDirectory, taxon_dir: Path):
    """Try to get Wtdb at primary_cluster."""
    Widb = wd.get_db("Widb")
    Widb = Widb[[i for i in Widb.columns if not Widb[i].isna().all()]]
    #
    gtdbtk_summary = pd.concat(
        [
            pd.read_csv(
                taxon_dir / f"gtdbtk.{marker_gene}.summary.tsv",
                sep="\t",
                na_values="N/A",
            )
            for marker_gene in ["ar122", "bac120"]
        ]
    )
    gtdbtk_summary["genome"] = gtdbtk_summary["user_genome"].apply(lambda x: f"{x}.fa")
    fastani_reference = gtdbtk_summary[1 ^ gtdbtk_summary["fastani_reference"].isna()][
        [
            "fastani_reference",
            "fastani_reference_radius",
            "fastani_taxonomy",
            "fastani_ani",
        ]
    ]
    Wtdb = Widb.merge(gtdbtk_summary[["genome", "classification"]])
    return Wtdb, fastani_reference


def join_SChdb(wd: WorkDirectory):
    Sdb = wd.get_db("Sdb")
    Cdb = wd.get_db("Cdb")
    SCdb = Sdb.merge(Cdb)
    genomeInfo = wd.get_db("genomeInfo").rename(
        {
            "genome": "Bin Id",
            "completeness": "Completeness",
            "contamination": "Contamination",
            "strain_heterogeneity": "Strain heterogeneity",
        },
        axis=1,
    )
    Stdb = SCdb.merge(genomeInfo, left_on="genome", right_on="Bin Id")[
        [
            "genome",
            "score",
            "secondary_cluster",
            "Completeness",
            "Contamination",
            "Strain heterogeneity",
        ]
    ]
    Stdb["site"] = Stdb["genome"].apply(file_path.infer_site_of)
    return Stdb


def infer_contig_from_genome_name(genome_path: Path, __file_depths={}):
    from PyLib.biotool.fna_msg import seq_total_depth, statistic_fna
    from PyLib.reader.read_outputs import fasta, jgi_depths

    #
    site = file_path.infer_site_of(genome_path)
    depth_path = file_path.contig(f"{site}", "_cut", "-jgi.depth")
    #
    if depth_path not in __file_depths:
        with open(depth_path) as fi:
            sample_list, ctg_depth = jgi_depths(fi)
            __file_depths[depth_path] = ctg_depth
            logger.info(depth_path, "load to memory")
    ctg_depth = __file_depths[depth_path]
    #
    totalAvgDepth, depths = seq_total_depth(ctg_depth, fasta(genome_path))
    return pd.Series(
        data=[genome_path, *statistic_fna(fasta(genome_path)), totalAvgDepth],
        index=[
            "genome_path",
            "SeqNumbers",
            "MaxLength",
            "GenomeSize",
            "GC",
            "N50",
            "L50",
            "totalAvgDepth",
        ],
    )


def recover_Stdb_geonme(wd: WorkDirectory, Stdb: pd.DataFrame):
    MAG_dir = wd.get_db("Bdb")[["genome", "location"]]
    dir2MAGs = MAG_dir.set_index("genome").loc[Stdb["genome"], "location"]
    return dir2MAGs


def main(wd, taxon_dir, Stdb_dir, Wtdb_dir, fastani_reference_dir, depth):
    Stdb = join_SChdb(wd)
    if depth:
        dir2MAGs = recover_Stdb_geonme(wd, Stdb)
        Sidb = dir2MAGs.apply(Path).apply(infer_contig_from_genome_name)
        Stdb = Stdb.merge(Sidb.reset_index())
    Stdb.to_csv(Stdb_dir, index=False)

    Wtdb, fastani_reference = join_Widb_gtdb(wd, taxon_dir)
    Wtdb.to_csv(Wtdb_dir, index=False)
    fastani_reference.to_csv(fastani_reference_dir, index=False)

    return 0


@click.command()
@click.option("--loglevel", default="INFO", type=str, help="set level of logger")
@click.option(
    "--drep", default=file_path.all_bins("drep"), type=Path, help="input drep"
)
@click.option(
    "--taxon", default=file_path.all_bins("gtdbtk"), type=Path, help="input taxon"
)
@click.option(
    "--Stdb", default=file_path.results("Stdb"), type=Path, help="output Stdb"
)
@click.option(
    "--Wtdb", default=file_path.results("Wtdb"), type=Path, help="output Wtdb"
)
@click.option(
    "--fastani-reference",
    default=file_path.results("fastani_reference"),
    type=Path,
    help="output fastani_reference",
)
@click.option(
    "--depth",
    is_flag=True,
    help="collect depth",
)
def run(
    loglevel: str,
    drep,
    taxon,
    stdb: Path,
    wtdb: Path,
    fastani_reference: Path,
    depth=False,
):
    logger.setLevel(level=loglevel.upper())  # info

    wd = WorkDirectory(drep)
    stdb.parent.mkdir(parents=True, exist_ok=True)
    wtdb.parent.mkdir(parents=True, exist_ok=True)
    fastani_reference.parent.mkdir(parents=True, exist_ok=True)

    logger.warning(">>> job start")
    state = main(wd, taxon, stdb, wtdb, fastani_reference, depth)
    logger.warning(">>> job finish")
    if state == 0:
        logger.info("success!")


if __name__ == "__main__":
    run()
