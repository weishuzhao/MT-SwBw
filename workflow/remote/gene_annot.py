# -*- coding: utf-8 -*-
"""
 * @Date: 2022-04-15 13:56:44
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-09-12 16:28:37
 * @FilePath: /2021_09-MT10kSW/workflow/remote/gene_annot.py
 * @Description:
"""

import os
import re
from pathlib import Path
from typing import Generator, Union

import click
import pandas as pd
from PyLib.biotool.kegg import load_ko00002
from PyLib.PyLibTool.file_info import verbose_import
from workflow.utils.file_path import FilePath


logger = verbose_import(__name__, __doc__)
file_path = FilePath()


def read_table(text, sep="\t", annot="#", title: list = None, openit=False):
    if openit:
        text = open(text)  # type: ignore
    for line in text:
        if line.startswith(annot):
            if title is not None:
                title.clear()
                title.extend(line[len(annot) :].rstrip().split(sep))
            continue
        values = line.strip().split(sep)
        if values:
            yield values
    if openit:
        text.close()  # type: ignore


class gene2KO:
    class gene_ko_iter:
        def __init__(self, filename: Path):
            self.filename = filename

        def __call__(self) -> Generator[tuple[str, str], None, None]:
            raise NotImplementedError

    class ghost(gene_ko_iter):
        def __call__(self) -> Generator[tuple[str, str], None, None]:
            with open(self.filename) as text:
                for values in read_table(text):
                    if len(values) > 1:
                        gene, ko = values[0:2]
                        yield gene, ko

    class kofam(gene_ko_iter):
        def __call__(self) -> Generator[tuple[str, str], None, None]:
            with open(self.filename) as text:
                for values in read_table(text):
                    if len(values) > 1:
                        gene, ko = values[1:3]
                        yield gene, ko

    class eggnog(gene_ko_iter):
        def __call__(self) -> Generator[tuple[str, str], None, None]:
            i_KEGG_ko = 11
            with open(self.filename) as text:
                for values in read_table(text):
                    kos = values[i_KEGG_ko]
                    if len(kos) > 1:
                        gene = values[0]
                        for ko in kos.split(","):
                            yield gene, ko[3:]

    class mantis(gene_ko_iter):
        KO_PATTERN = re.compile("(^|;)(K\\d{5})(;|$)")

        def __call__(self) -> Generator[tuple[str, str], None, None]:
            with open(self.filename) as text:
                for values in read_table(text):
                    Ref_Hits = values[2]
                    kos: list[tuple[str, str, str]] = re.findall(
                        self.KO_PATTERN, Ref_Hits
                    )
                    if len(kos) > 0:
                        gene = values[0]
                        for ko in kos:
                            yield gene, ko[1]

    ## collect gene KO
    annoters = [ghost, kofam, eggnog, mantis]

    def get_gene_KOs(self) -> dict[str, str]:
        """Only keep the first match:
        >>> gene_KOs.setdefault(gene, ko)"""
        gene_KOs: dict[str, str] = {}
        for annoter, file in zip(self.annoters, self.ann_files):
            if not os.path.isfile(file):
                continue

            gene_KOs_: dict[str, list[str]] = {}
            for gene, ko in annoter(file)():
                gene_KOs_.setdefault(gene, []).append(ko)
            for gene in gene_KOs_:
                if gene not in gene_KOs:
                    gene_KOs[gene] = ":".join(gene_KOs_[gene])

        return gene_KOs

    def get_gene_annots(self):
        return pd.Series(self.get_gene_KOs(), name="ko")

    def __init__(self, pattern: Union[Path, list[Path]]):

        if not isinstance(pattern, list):
            pattern = [pattern]
        ann_files_ = self._infer_ann_files(pattern)
        ann_files = [ann_files_.get(i, Path()) for i, _ in enumerate(self.annoters)]

        if not any(ann_files):
            raise FileNotFoundError(f"pattren(s) '{pattern}' donot match any file!")

        self.ann_files = ann_files

    def _infer_ann_files(self, patterns: list[Path]):
        ann_files = {i: Path() for i, _ in enumerate(self.annoters)}

        for pattern in reversed(patterns):
            pattern_re = re.compile(pattern.name)
            for file in pattern.parent.iterdir():
                if pattern_re.search(file.name):
                    for i, source in enumerate(self.annoters):
                        if source.__name__ in file.name.lower():
                            ann_files[i] = file
        return ann_files


def load_rep2all(all_100_: Path, all_clu_: Path):
    all_100 = pd.read_csv(all_100_, sep="\t", header=None, names=["rep100", "all"])
    all_clu = pd.read_csv(all_clu_, sep="\t", header=None, names=["rep", "rep100"])
    rep2all = all_100.merge(all_clu)[["rep", "all"]]
    return rep2all


def get_all_gene_annots(gene_annots: pd.Series, rep2all: pd.DataFrame):
    ko_exploded = (
        gene_annots.apply(lambda x: x.split(":") if x.startswith("K") else [])
        .explode()
        .dropna()
        .reset_index()
    )
    all_gene_annots = ko_exploded.merge(rep2all, left_on="index", right_on="rep")[
        ["all", "ko"]
    ].set_index("all")

    return all_gene_annots


def get_gmodule(genomeko: pd.DataFrame):
    KEGG_DIR = "~/Data/Database2/KEGG"
    _, modules = load_ko00002(KEGG_DIR)
    gmodule_ = (
        genomeko.apply(lambda x: x[x > 0].index, axis=0)
        .apply(
            lambda x: {
                mname: module.completeness(x) for mname, module in modules.items()
            }
        )
        .apply(lambda x: pd.Series(x))
    )
    gmodule = gmodule_.T[gmodule_.apply(sum, 0) > 0]

    return gmodule


def main(
    annot_prefix,
    all_100_path: Path,
    all_clu_path: Path,
    all_gene_annots_path: Path,
):
    gene_annots = gene2KO(annot_prefix).get_gene_annots()

    # gene_annots = pd.read_csv(gene_annots_path, index_col=0)["ko"]
    rep2all = load_rep2all(all_100_path, all_clu_path)

    all_gene_annots = get_all_gene_annots(gene_annots, rep2all)
    all_gene_annots.to_csv(all_gene_annots_path)


@click.command()
@click.option("--loglevel", default="INFO", type=str, help="set level of logger")
@click.option(
    "--all-100",
    default=file_path.all_bins("collect_annot") / "all-clu_100.tsv",
    type=Path,
    help="",
)
@click.option(
    "--all-clu",
    default=file_path.all_bins("collect_annot") / "all-clu.tsv",
    type=Path,
    help="",
)
@click.option(
    "--annot-prefix",
    default=file_path.all_bins("collect_annot") / "all-clu_rep-kofam.tsv",
    type=Path,
    help="",
)
@click.option(
    "--gene-annots",
    default=file_path.results("all_gene_annots.csv"),
    type=Path,
    help="",
)
def run(
    loglevel: str,
    all_100: Path,
    all_clu: Path,
    annot_prefix: Path,
    gene_annots: Path,
):
    logger.setLevel(level=loglevel.upper())  # info

    gene_annots.parent.mkdir(parents=True, exist_ok=True)

    logger.warning(">>> job start")
    state = main(annot_prefix, all_100, all_clu, gene_annots)
    logger.warning(">>> job finish")
    if state == 0:
        logger.info("success!")


if __name__ == "__main__":
    run()
