# -*- coding: utf-8 -*-
"""
 * @Date: 2022-06-18 16:50:39
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-06-19 16:10:34
 * @FilePath: /2021_09-MT10kSW/workflow/utils/ncbi_api.py
 * @Description:
"""

from typing import Dict
import pandas as pd
from Bio import Entrez
import xml.dom.minidom

from workflow.utils.PyLib.PyLibTool.file_info import verbose_import
from pprint import pprint


logger = verbose_import(__name__, __doc__)
logger.setLevel("INFO")


Entrez.email = "hwrn.aou@sjtu.edu.cn"
RETMAX = 100

warning_samples_cache: Dict = {}


def get_remote_sra():
    SraList_file = "SraAccList.txt"
    SraList = open(SraList_file).read().split()
    from joblib import Parallel, delayed

    def get_remote_sra(srr):
        import os

        os.system(f"prefetch {srr}; cd {srr}; fastq-dump --split-3 {srr}.sra")

    def gzip_remote_fq(srr, i):
        import os

        os.system(f"gzip {srr}/{srr}_{i}.fastq")

    Parallel(40, verbose=20)(delayed(get_remote_sra)(srr) for srr in SraList)
    Parallel(40, verbose=20)(
        delayed(gzip_remote_fq)(srr, i) for srr in SraList for i in (1, 2)
    )


def parse_ExpXml(ExpXml) -> tuple[str, str]:
    DOMTree: xml.dom.minidom.Document = xml.dom.minidom.parseString(
        "<record>" + ExpXml + "</record>"
    )
    srx_acc = DOMTree.documentElement.childNodes[2].attributes["acc"].value
    srx_nane = DOMTree.documentElement.childNodes[2].attributes["name"].value
    return (srx_acc, srx_nane)


def SRX2SRR(srx):
    with Entrez.esearch(db="sra", term=srx) as handle:
        record = Entrez.read(handle)
    Id = record["IdList"][0]
    with Entrez.esummary(db="sra", id=Id) as handle:
        record = Entrez.read(handle)
    if isinstance(record, Entrez.Parser.ListElement):
        if len(record) == 1:
            record = record[0]
            Run = record["Runs"]
        else:
            warning_samples_cache[Id] = "multiple records of srx", record
    else:
        Run = record["Runs"]
    srr = Run.split('acc="')[1].split('"')[0]
    return srr


if __name__ == "__main__":
    with Entrez.einfo() as handle:
        record = Entrez.read(handle)

    # esearch.fcgi?db=assembly&term=SAMEA2619802
    prj = "PRJNA635214"
    with Entrez.esearch(db="sra", term=prj, retmax=RETMAX) as handle:
        record = Entrez.read(handle)
    assert len(record["IdList"]) < RETMAX

    records = {}
    for Id in record["IdList"]:
        logger.info("collect info by key '%s'", Id)
        with Entrez.esummary(db="sra", id=Id) as handle:
            record = Entrez.read(handle)
            if isinstance(record, Entrez.Parser.ListElement):
                if len(record) == 1:
                    record = record[0]
                    records[Id] = record
                else:
                    warning_samples_cache[Id] = "multiple records", record
            else:
                records[Id] = record

    logger.info("total items number: %d", len(records))
    srx_message = dict(parse_ExpXml(record["ExpXml"]) for record in records.values())
    srr_message = {SRX2SRR(srx): message for srx, message in srx_message.items()}
    len(srr_message)
    tmp = {
        k: w.split("sediment ")[1]
        for k, w in srr_message.items()
        if "metatranscriptomic" not in w.lower()
    }
    len(tmp)

    site_message: dict[str, dict[str, str]] = {}
    for srr, message in srr_message.items():
        if "metatranscriptomic" in message:
            continue
        message = message.split("sediment ")[1]
        if message.startswith("dive "):
            message = message[5:]
        site, subsample = message.split(" ", 1)
        site_message.setdefault(site, {})[subsample] = srr
    pprint(site_message)
    print(*((srr, tmp[srr]) for srr in sorted(tmp)), sep="\n")

    # now srr_message gone to data/sample_meta_sediment.tsv

    sample_meta_sediment = pd.read_csv("data/sample_meta_sediment.tsv", sep="\t")
    sample_meta_sediment.apply(
        lambda x: (x.iloc[0], x.iloc[1], x.iloc[3]), axis=1
    ).apply(
        lambda x: [
            f"ln -s "
            f"~/Data/Database2/metagenome/PRJNA635214/{x[0]}/{x[0]}_{i}.fastq.gz "
            f"data/{x[2]}..{x[1]}_{i}.fq.gz"
            for i in (1, 2)
        ]
    ).apply(
        lambda x: print(*x, sep="\n")
    )
