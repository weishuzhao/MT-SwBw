# -*- coding: utf-8 -*-
"""
 * @Date: 2021-11-16 13:49:03
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-11-18 15:54:09
 * @FilePath: /2021_09-MT10kSW/Analyze/binannot/map_bin_gene.py
 * @Description:
"""


import os
import pandas as pd
from Bio import SeqIO

SAMPLES = []


def load_samples():
    sample_meta = pd.read_csv("00_data/sample_meta.tsv", sep="\t")
    global SAMPLES
    SAMPLES = type(SAMPLES)(sample_meta["Station"])


load_samples()


def load_bin_names(sample, _Stdb={}):
    if not _Stdb:
        _Stdb = {"sample": pd.read_csv("Analyze/drep/Stdb.csv")}
    Stdb = _Stdb["sample"]
    bin_names = Stdb[Stdb["sample"] == sample]["genome"]
    return list(bin_names)


def load_bin_contigs(sample, bin_name):
    method = f"{sample}-megahit"
    bin_file = os.path.join("Pipe", method,
                            "04_bin", "DASTool_bins",
                            bin_name)
    with open(bin_file) as fa_in:
        contigs = [contigId.split()[0] for (contigId, seq) in
                   SeqIO.FastaIO.SimpleFastaParser(fa_in)]
    return contigs


def load_RPb_bin(sample, bin_name, _RPb={}):
    if sample not in _RPb:
        method = f"{sample}-megahit"
        RPb_file = os.path.join("Pipe", method,
                                "03_annot",
                                f"03_annot..{method}.tpm")
        _RPb = {sample: pd.read_csv(RPb_file, sep="\t",
                                    names=["gene", "KO", "contig", "RPb"],
                                    header=0)}
    RPb = _RPb[sample]
    bin_contigs = load_bin_contigs(sample, bin_name)
    bin_RPb = RPb[RPb["contig"].apply(lambda x: x in bin_contigs)]
    bin_RPb.loc[:, "binId"] = bin_name
    return bin_RPb


def load_mrk_RPb(sample, bin_RPb: pd.DataFrame, _mrk = {}):
    if sample not in _mrk:
        method = f"{sample}-megahit"
        checkm_marker_raw: pd.DataFrame = pd.read_csv(
            os.path.join("Pipe", method,
                         "03_annot", f"03_annot..{method}_cut-checkm_marker.tsv"),
            sep="\t", names=["sample", "marker", "gene"])
        checkm_marker_raw["dup"] = ["&&" in gene for gene in checkm_marker_raw["gene"]]
        checkm_marker = checkm_marker_raw. \
            join(checkm_marker_raw['gene']. \
                     str.split('&&', expand=True). \
                     stack(). \
                     reset_index(level=1, drop=True). \
                     rename('gene_single'))
        _mrk = {sample: checkm_marker_raw. \
                join(checkm_marker_raw['gene']. \
                str.split('&&', expand=True). \
                stack(). \
                reset_index(level=1, drop=True). \
                rename('gene_single'))}
    checkm_marker = _mrk[sample]
    #
    bin_markers = checkm_marker.merge(bin_RPb, how="left",
                                      left_on="gene_single",
                                      right_on="gene")
    mrk_bin = bin_markers.pivot_table(index="marker", columns="binId",
                                      values="RPb", aggfunc=sum)
    #
    return mrk_bin


def load_PF_bin(sample, bin_name):
    bin_RPb: pd.DataFrame = load_RPb_bin(sample, bin_name)
    mrk_bin = load_mrk_RPb(sample, bin_RPb)
    KO_bin = bin_RPb.pivot_table(index="KO", columns="binId",
                                 values="RPb", aggfunc=sum)
    return mrk_bin.append(KO_bin)


def load_PF_sample(sample):
    bin_names = load_bin_names(sample)
    PF_sample = pd.DataFrame()
    for bin_name in bin_names:
        PF_sample = pd.merge(PF_sample, load_PF_bin(sample, bin_name),
                             how="outer", left_index=True, right_index=True)
    return PF_sample


if __name__ == "__main__":
    for sample in SAMPLES:
        PF_sample = load_PF_sample(sample)
        PF_sample.to_csv(
            "Analyze/binannot/" + sample + "-PF_bin.tsv",
            sep="\t")
