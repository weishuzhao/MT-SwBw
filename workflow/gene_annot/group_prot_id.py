# -*- coding: utf-8 -*-
"""
 * @Date: 2022-03-28 15:43:25
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-05-20 19:22:25
 * @FilePath: /2022_05-ZFMG-release/workflow/gene_annot/group_prot_id.py
 * @Description:
"""

import os
import re
from pathlib import Path
from typing import Dict, Generator, List, Tuple

import pandas as pd


def read_table(text, sep="\t", annot="#", title=None, openit=False):
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
        def __init__(self, filename):
            self.filename = filename

        def __call__(self) -> Generator[Tuple[str, str], None, None]:
            raise NotImplementedError

    class ghost(gene_ko_iter):
        def __call__(self) -> Generator[Tuple[str, str], None, None]:
            with open(self.filename) as text:
                for values in read_table(text):
                    if len(values) > 1:
                        gene, ko = values[0:2]
                        yield gene, ko

    class kofam(gene_ko_iter):
        def __call__(self) -> Generator[Tuple[str, str], None, None]:
            with open(self.filename) as text:
                for values in read_table(text):
                    if len(values) > 1:
                        gene, ko = values[1:3]
                        yield gene, ko

    class eggnog(gene_ko_iter):
        def __call__(self) -> Generator[Tuple[str, str], None, None]:
            i_KEGG_ko = 11
            with open(self.filename) as text:
                for values in read_table(text):
                    kos = values[i_KEGG_ko]
                    if len(kos) > 1:
                        gene = values[0]
                        for ko in kos.split(","):
                            yield gene, ko[3:]

    ## collect gene KO
    annoters = [ghost, kofam, eggnog]

    def get_gene_KOs(self) -> Dict[str, str]:
        """Only keep the first match:
        >>> gene_KOs.setdefault(gene, ko)"""
        gene_KOs: Dict[str, str] = {}
        for format, file in zip(self.annoters, self.ann_files):
            if not os.path.exists(file):
                continue

            gene_KOs_: Dict[str, List[str]] = {}
            for gene, ko in format(file)():
                gene_KOs_.setdefault(gene, []).append(ko)
            for gene in gene_KOs_:
                if gene not in gene_KOs:
                    gene_KOs[gene] = ":".join(gene_KOs_[gene])

        return gene_KOs

    def __init__(self, pattern: Path):
        ann_files = [Path() for i in self.annoters]

        pattern_re = re.compile(pattern.name)
        for file in pattern.parent.iterdir():
            if pattern_re.search(file.name):
                for i, source in enumerate(self.annoters):
                    if source.__name__ in file.name.lower():
                        ann_files[i] = file

        if not any(ann_files):
            raise FileNotFoundError(f"pattren '{pattern}' donot match any file!")

        self.ann_files = ann_files


entry2ko = pd.read_csv(
    "workflow/gene_annot/entry2ko.csv", index_col=0, names=["module", "ko"], header=0
)


def group_prot_id_iter(filename, annot_out_f, bins=range(102)):
    MT2ME: pd.DataFrame = pd.read_csv(
        filename,
        sep="\t",
        names="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split(),
        usecols=["qseqid", "pident"],
    )
    MT2ME["sample"] = MT2ME["qseqid"].apply(lambda x: x.split("|")[0])
    for sample, MT2ME_i in MT2ME.groupby("sample"):
        annot_out = annot_out_f.format(sample=sample)
        count = pd.read_csv(f"{annot_out}-count.txt", sep="\t", skiprows=1)
        zone_gene = (
            pd.DataFrame(
                {
                    "gene": count.apply(
                        lambda x: f'{x["Chr"]}_{x["Geneid"].split("_")[1]}', axis=1
                    ),
                    "qseqid": count.apply(
                        lambda x: f'{sample}|{x["Chr"]}_{x["Geneid"].split("_")[1]}',
                        axis=1,
                    ),
                    "reads": count[count.columns[6]],
                    "tpm": count[count.columns[6]] / count[count.columns[5]],
                    "count": 1,
                }
            )
            .merge(MT2ME_i, on="qseqid", how="left")
            .fillna(value={"pident": 0})
        )
        zone_gene["tpm"] = zone_gene["tpm"] / sum(zone_gene["tpm"]) * 1e6
        zone_gene["zone"] = pd.cut(
            zone_gene["pident"],
            bins=bins,
            precision=0,
            include_lowest=True,
            right=False,
        )

        gene_KOs = (
            pd.Series(gene2KO(f"{annot_out}-.*.tsv").get_gene_KOs(), name="ko")
            .apply(lambda x: x.split(":"))
            .explode()
            .reset_index()
            .rename({"index": "gene"}, axis=1)
        )
        koreadszone = (
            zone_gene.merge(gene_KOs, how="left")
            .fillna({"ko": ""})
            .groupby(["ko", "zone"])[["reads", "tpm", "count"]]
            .apply(sum)
            .dropna()
            .reset_index()
        )
        koreadszone["sample"] = sample
        modulereadszone = (
            zone_gene.merge(gene_KOs.merge(entry2ko)[["gene", "module"]], how="left")
            .fillna({"module": ""})
            .groupby(["module", "zone"])[["reads", "tpm", "count"]]
            .apply(sum)
            .dropna()
            .reset_index()
        )
        modulereadszone["sample"] = sample

        yield koreadszone, modulereadszone


def get_reads_zone(filename):
    annot_out_f = "{sample}-megahit/03_annot/03_annot..{sample}-megahit"
    koreadszones = []
    modulereadszones = []
    for koreadszone, modulereadszone in group_prot_id_iter(filename, annot_out_f):
        koreadszones.append(koreadszone)
        modulereadszones.append(modulereadszone)

    return pd.concat(koreadszones), pd.concat(modulereadszones)


if os.path.isfile("workflow/gene_annot/ko_reads_zone.csv"):
    ko_reads_zone = pd.read_csv("workflow/gene_annot/ko_reads_zone.csv", index_col=0)
else:
    filename = "workflow/gene_annot/faa/MT.faa.blastp.tsv"
    ko_reads_zone_MT, module_reads_zone_MT = get_reads_zone(filename)
    filename = "workflow/gene_annot/faa/ME.faa.blastp.tsv"
    ko_reads_zone_ME, module_reads_zone_ME = get_reads_zone(filename)
    pd.concat([ko_reads_zone_MT, ko_reads_zone_ME]).to_csv(
        "workflow/gene_annot/ko_reads_zone.csv"
    )
    pd.concat([module_reads_zone_MT, module_reads_zone_ME]).to_csv(
        "workflow/gene_annot/module_reads_zone.csv"
    )

# part 2: best match genes
if os.path.isfile("workflow/gene_annot/species_ident_gene.csv"):
    blast_geneko = pd.read_csv("workflow/gene_annot/species_ident_gene.csv")
else:
    blastn_ident_ = pd.concat(
        [
            pd.read_csv(
                f"workflow/gene_annot/fna/{location}.fna.blastn.tsv",
                sep="\t",
                names="qaccver saccver pident length mismatch gapopen qstart qend sstart send\n   evalue bitscore".split(),
                usecols=["qaccver", "pident"],
            )
            for location in ("ME", "MT")
        ]
    )
    blastn_ident = (
        blastn_ident_[blastn_ident_["pident"] >= 95]["qaccver"]
        .drop_duplicates()
        .apply(lambda x: pd.Series(x.split("|")))
        .rename({0: "sample", 1: "gene"}, axis=1)
    )

    annot_out_f = "{sample}-megahit/03_annot/03_annot..{sample}-megahit"
    blast_geneko = pd.DataFrame()
    for sample, genei in blastn_ident.groupby("sample"):
        annot_out = annot_out_f.format(sample=sample)
        gene_KOs = pd.Series(
            gene2KO(Path(f"{annot_out}-.*.tsv")).get_gene_KOs(), name="ko"
        )
        blast_geneko = pd.concat(
            [
                blast_geneko,
                genei.merge(gene_KOs, left_on="gene", right_index=True),
            ]
        )

    blast_geneko.to_csv("workflow/gene_annot/species_ident_gene.csv", index=False)
