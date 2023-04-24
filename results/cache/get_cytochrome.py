# -*- coding: utf-8 -*-
"""
 * @Date: 2022-09-12 09:37:49
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-09-13 10:31:12
 * @FilePath: /2021_09-MT10kSW/results/cache/get_cytochrome.py
 * @Description:
    recognize cytochrome-related genomes
"""

from typing import Union, TextIO, Optional
import re
from workflow.utils.file_path import file_path
from workflow.remote.gene_annot import load_rep2all
import pandas as pd


def read_table(
    text: Union[list[str], str, TextIO],
    sep="\t",
    annot="#",
    title: Optional[list[str]] = None,
    openit=False,
):
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


class Mantis:
    KO_PATTERN = re.compile("(^|;)(K\\d{5})(;|$)")

    @classmethod
    def _Ref_Hits_to_KO(cls, ref_hits_str) -> list[str]:
        kos: list[tuple[str, str, str]] = re.findall(cls.KO_PATTERN, ref_hits_str)
        if len(kos) > 0:
            return [ko[1] for ko in kos]
        return []


annot_mantis = file_path.all_bins("collect_annot") / "all-clu_rep-mantis.tsv"
all_100 = file_path.all_bins("collect_annot") / "all-clu_100.tsv"
all_clu = file_path.all_bins("collect_annot") / "all-clu.tsv"
gene_map = pd.read_csv(file_path.results("gene_map.csv"))
Wtdb = pd.read_csv(file_path.results("Wtdb"))

exclude_kos = {
    "K00101": "lldD; L-lactate dehydrogenase (cytochrome) [EC:1.1.2.3]",
    "K00102": "LDHD, dld; D-lactate dehydrogenase (cytochrome) [EC:1.1.2.4]",
    "K00114": "exaA; alcohol dehydrogenase (cytochrome c) [EC:1.1.2.8]",
    "K00142": "AASDH; acyl-CoA synthetase [EC:6.2.1.-]",
    "K00266": "gltD; glutamate synthase (NADPH) small chain [EC:1.4.1.13]",
    "K00380": "cysJ; sulfite reductase (NADPH) flavoprotein alpha-component [EC:1.8.1.2]",
    "K00386": "sorB; sulfite dehydrogenase (cytochrome) subunit B [EC:1.8.2.1]",
}

mantis_annot_ = {}
with annot_mantis.open() as ain:
    title = next(read_table(ain))
    for values in read_table(ain):
        if "cytochrome" in "\t".join(values):
            mantis_annot_[values[0]] = (*values[1:5], "\t".join(values[6:]).split("|"))


mantis_annot = pd.DataFrame(mantis_annot_, index=[*title[1:5], title[6]]).T
mantis_annot["KO"] = mantis_annot["Ref_Hits"].apply(lambda x: Mantis._Ref_Hits_to_KO(x))

rep2all = load_rep2all(all_100, all_clu)
genome_mantis = (
    mantis_annot.merge(rep2all, left_index=True, right_on="rep")
    .merge(gene_map, left_on="all", right_on="gene")
    .drop(["rep", "all"], axis=1)
)

Wtdb_mantis = Wtdb[["genome", "classification"]].merge(genome_mantis)
Wtdb_mantis[["genome", "classification", "partial", "Ref_Hits", "Links"]].to_csv(
    "results/cache/Wtdb_mantis_cytochrome.csv"
)
Wtdb_mantis.explode("KO").pivot_table("Links", "genome", "KO", len, 0).to_csv(
    "results/cache/genomeko_cytochrome.csv"
)
