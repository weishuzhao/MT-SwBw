# -*- coding: utf-8 -*-
"""
 * @Date: 2022-06-19 14:59:53
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-16 14:12:06
 * @FilePath: /2021_09-MT10kSW/workflow/utils/file_path.py
 * @Description:
"""

import os
from pathlib import Path
import pandas as pd


def load_sample_meta():

    tmp = (
        pd.read_csv("data/sample_meta.tsv", sep="\t", na_values="——")
        .assign(
            Latitude=lambda df: df.iloc[:, 2]
            .apply(lambda x: x.split(" ")[0].split("'N", 1)[0].split("°", 1))
            .apply(lambda x: (int(x[0]) + float(x[1]) / 60))
        )
        .assign(
            Longitude=lambda df: df.iloc[:, 2]
            .apply(lambda x: x.split(" ")[1].split("'E", 1)[0].split("°", 1))
            .apply(lambda x: (int(x[0]) + float(x[1]) / 60))
        )
    )
    tmp1 = pd.read_csv("data/sample_meta_sediment.tsv", sep="\t")

    sample_meta = (
        pd.concat(
            [
                pd.DataFrame(
                    {
                        "Site": tmp.iloc[:, 0],
                        "Layers": "water",
                        "Depth": tmp.iloc[:, 3],
                        "Latitude": tmp["Latitude"],
                        "Longitude": tmp["Longitude"],
                    }
                ),
                pd.DataFrame(
                    {
                        "Site": tmp1.iloc[:, 3],
                        "Layers": tmp1.iloc[:, 1].apply(lambda x: x[:-5]),
                        "Depth": tmp1.iloc[:, 7].apply(
                            lambda x: int(x.replace(",", ""))
                        ),
                        "Latitude": tmp1.iloc[:, 5],
                        "Longitude": tmp1.iloc[:, 6],
                    }
                ),
            ]
        )
        .assign(
            Location=lambda df: df["Depth"].apply(
                lambda x: "Bottom" if 10000 < x else "Slope",
            ),
        )
        .assign(
            Group=lambda df: df.apply(
                lambda s: ("B" if s["Location"] == "Bottom" else "S")
                + ("w" if s["Layers"] == "water" else "s"),
                axis=1,
            ),
        )
    )
    return sample_meta


def exec_site_layer_dict():
    """
    >>> from pprint import pprint
    >>> pprint(exec_site_layer_dict())
    {'D1T1': ['0-2', '4-6'],
     'D1T2': ['4-6', '24-26'],
     'MC02': ['8-10', '28-30'],
     'T1B10': ['0-2', '36-38', '44-46'],
     'T1B11': ['0-3'],
     'T1B3': ['R1.0-3', 'R2.0-3'],
     'T1B5': ['0-2', '8-10', '38-40', '28-30'],
     'T1B8': ['2-4', '16-18'],
     'T1L10': ['0-3', '6-9', '12-15', '18-21'],
     'T1L6': ['R1.0-3', 'R2.0-3'],
     'T3L11': ['0-3', '6-9', '12-15', '18-21'],
     'T3L14': ['0-2', '4-6', '6-8', '12-14', '18-20'],
     'T3L8': ['0-3', '6-9', '12-15', '18-21'],
     'TY.040': ['water'],
     'TY.041': ['water'],
     'TY.044': ['water'],
     'WQ.018': ['water'],
     'WQ.021': ['water'],
     'WQ.022': ['water'],
     'WQ.023': ['water'],
     'WQ.024': ['water'],
     'YW.019': ['water'],
     'YW.020': ['water'],
     'YW.021': ['water'],
     'YW.023': ['water']}
    """
    return sample_meta.groupby("Site").agg({"Layers": list}).to_dict()["Layers"]


class SnakeConfigGet:
    def __init__(self, config):
        """
        >>> snake_config = SnakeConfigGet({})
        >>> snake_config.THREADS
        4
        >>> snake_config.SLURM_ARRAY_TASK_ID
        22
        """
        self.config = config

    @property
    def THREADS(self):
        return int(
            self.config.get(
                "THREADS",
                self.config.get("SLURM_NTASKS", os.environ.get("SLURM_NTASKS", 4)),
            )
        )

    @property
    def SLURM_ARRAY_TASK_ID(self):
        return int(
            self.config.get(
                "SLURM_ARRAY_TASK_ID",
                os.environ.get("SLURM_ARRAY_TASK_ID", -1),
            )
        )


class FilePath:
    def __init__(self, WORK_DIR=Path("./")):
        self.WORK_DIR: Path = WORK_DIR
        self.data = WORK_DIR / "data"
        self.pipe = WORK_DIR / "pipe"

    def log(self, name, site, layer=None) -> Path:
        if layer:
            return self.pipe / site / "log" / f"site..{layer}" / f"{name}.log"
        return self.pipe / site / "log" / f"{name}.log"

    def trimmed_reads(self, site, layer, i) -> Path:
        return self.pipe / site / f"01_trim..{site}..{layer}_{i}.fq.gz"

    def stat_reads(self, site, layer, work, suffix) -> Path:
        return self.pipe / site / "01_alpha" / f"{work}..{site}..{layer}{suffix}"

    def contig(self, site, label="", suffix=".fa") -> Path:
        return self.pipe / site / f"02_assem..{site}{label}{suffix}"

    def bam(self, site, layer, suffix=".bam") -> Path:
        return self.pipe / site / f"02_assem..{site}..{layer}{suffix}"

    def gene(self, site, label="", suffix=".faa") -> Path:
        return self.pipe / site / f"03_annot..{site}{label}{suffix}"

    def ctg2mag(self, site, method="") -> Path:
        return self.pipe / site / "04_bin" / "ctg2mag" / f"{method}.tsv"

    def dastool(self, site, suffix="") -> Path:
        return self.pipe / site / "04_bin" / "ctg2mag" / f"DASTool{suffix}"

    def dastool_bins(self, site) -> Path:
        return self.pipe / site / "04_bin" / "bins"

    def env_genes(self, work: str) -> Path:
        return self.WORK_DIR / "results" / "env_genes" / work

    def all_bins(self, work: str) -> Path:
        if work == "drep":
            return self.WORK_DIR / "04_bin" / "02_drep"
        if work == "gtdbtk":
            return self.WORK_DIR / "04_bin" / "02_taxon"
        if work.startswith("collect"):
            work = work.split("_", 1)[1]
            if work == "genome":
                return self.WORK_DIR / "04_bin" / "03_collection" / "genome"
            if work == "gene":
                return self.WORK_DIR / "04_bin" / "03_collection" / "gene"
            return self.WORK_DIR / "04_bin" / "03_collection" / work
        if work.startswith("ref"):
            work = work.split("_", 1)[1]
            if work == "genome":
                return self.WORK_DIR / "04_bin" / "03_reference" / "genome"
            if work == "gene":
                return self.WORK_DIR / "04_bin" / "03_reference" / "gene"
            return self.WORK_DIR / "04_bin" / "03_reference" / work
        return self.WORK_DIR / "results" / "MAGs" / work

    def bins_tree(self, work: str) -> Path:
        return self.WORK_DIR / "04_bin" / "04_tree" / work

    def bins_tpm(self, work: str, site=None, layer=None, suffix=".bam") -> Path:
        if site and layer:  # bam file
            return (
                self.WORK_DIR / "04_bin" / "04_tpm" / f"{work}..{site}..{layer}{suffix}"
            )
        return self.WORK_DIR / "04_bin" / "04_tpm" / work

    def bins_salmon(self, work: str, site=None, layer=None, suffix=".quant") -> Path:
        if site and layer:  # bam file
            return (
                self.WORK_DIR
                / "04_bin"
                / "04_salmon"
                / f"{work}..{site}..{layer}{suffix}"
            )
        return self.WORK_DIR / "04_bin" / "04_salmon" / work

    def otus(self, work: str) -> Path:
        if work.startswith("phyloFlash"):
            method, work = work.split("_", 1)
            return (
                self.WORK_DIR
                / "results"
                / "reads_diversity"
                / f"div.{method}.{work}.csv"
            )
        return self.WORK_DIR / "results" / "reads_diversity" / work

    def results(self, work: str) -> Path:
        if work == "Stdb":
            return self.WORK_DIR / "results" / "MAGs" / "Stdb.csv"
        if work == "Wtdb":
            return self.WORK_DIR / "results" / "MAGs" / "Wtdb.csv"
        if work == "fastani_reference":
            return self.WORK_DIR / "results" / "MAGs" / "fastani_reference.csv"
        return self.WORK_DIR / "results" / "MAGs" / work

    def cache(self, work: str) -> Path:
        return self.WORK_DIR / "results" / "cache" / work

    def figs(self, work: str, suffix=".svg") -> Path:
        return self.WORK_DIR / "results" / "figs" / f"{work}{suffix}"

    def tabs(self, work: str, suffix=".tsv") -> Path:
        return self.WORK_DIR / "results" / "tabs" / f"{work}{suffix}"

    @staticmethod
    def infer_site_of(genome):
        genome = Path(genome).name
        site = genome.rsplit("-", 1)[0]
        return site

    @staticmethod
    def annot_tsv_of(gene: Path, annot=""):
        path, name_ = gene.parent, gene.name
        name = name_.rsplit(".faa", 1)[0] + (f"-{annot}" if annot else "") + ".tsv"
        return path / name

    @staticmethod
    def checkm_of(bins: Path):
        path, name_ = bins.parent, bins.name
        name = name_ + "-checkm.tsv"
        return path / name

    @staticmethod
    def scratch(path, here="~/", there="/scratch/home/acct-clsxx/clsxx"):
        full_path = Path(path).expanduser().absolute()
        full_here = Path(here).expanduser().absolute()
        full_there = Path(there).expanduser().absolute()
        here_common = os.path.commonprefix([full_path, full_here])
        path_relative = str(full_path)[len(here_common) :].strip("/")
        scratch_path = full_there / path_relative
        return scratch_path


sample_meta = load_sample_meta()
site_layer_dict = exec_site_layer_dict()
file_path = FilePath()
