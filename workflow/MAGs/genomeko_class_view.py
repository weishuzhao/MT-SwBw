# -*- coding: utf-8 -*-
"""
 * @Date: 2022-06-14 17:10:52
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-06-14 18:27:48
 * @FilePath: /2021_09-MT10kSW/workflow/MAGs/genomeko_class_view.py
 * @Description:
"""

import sys
from pathlib import Path
from time import sleep

import pandas as pd
from Bio.KEGG import REST


def get_ko_annots(kos, chunk_size=100):
    unannot_ko = set(kos)
    ko_annots = pd.DataFrame(columns=["ko", "desc"])
    while unannot_ko:
        print(len(unannot_ko))
        kos = list(unannot_ko)[:chunk_size]
        try:
            rest = REST.kegg_list("+".join(kos))
        except Exception:
            sleep(5)
            continue
        else:
            a = pd.read_csv(rest, names=["ko", "desc"], sep="\t")
            ko_annots = pd.concat([ko_annots, a])
            unannot_ko = unannot_ko - set(kos)

    return ko_annots.assign(ko=lambda df: df["ko"].apply(lambda x: x[3:]))


if __name__ == "__main__":
    Stdb_dir = Path(sys.argv[1])
    Wtdb_dir = Path(sys.argv[2])
    genomeko_dir = Path(sys.argv[3])
    classkoannot_dir = Path(sys.argv[4])

    Stdb = pd.read_csv(Stdb_dir)
    Wtdb = pd.read_csv(Wtdb_dir)
    genomeko = pd.read_csv(genomeko_dir)

    genomeclass = (
        Wtdb[["cluster", "classification"]]
        .merge(Stdb, left_on="cluster", right_on="secondary_cluster")
        .assign(
            c__=lambda df: df["classification"].apply(lambda x: x.split(";")[2][3:])
        )
    )[["genome", "c__"]]
    classko = (
        genomeko.set_index("ko")
        .T.merge(genomeclass, left_index=True, right_on="genome")
        .assign(genome=lambda df: 1)
        .groupby("c__")
        .agg((lambda x: sum(x > 0)))
        .T
    )

    ko_annots = get_ko_annots(classko.index)
    classkoannot = classko.merge(
        ko_annots,
        left_index=True,
        right_on="ko",
        how="left",
    ).set_index("ko")
    classkoannot.to_csv(classkoannot_dir)
