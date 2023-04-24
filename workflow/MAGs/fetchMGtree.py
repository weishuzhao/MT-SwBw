# -*- coding: utf-8 -*-
"""
 * @Date: 2022-04-10 21:02:42
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-09-08 11:35:19
 * @FilePath: /2021_09-MT10kSW/workflow/MAGs/fetchMGtree.py
 * @Description:
"""

import os
import tempfile
from pathlib import Path
from typing import Iterable, Literal

import click
import pandas as pd
import numpy as np
from Bio import AlignIO, SeqIO
from joblib import Parallel, delayed
from PyLib.PyLibTool.file_info import verbose_import
from workflow.utils.file_path import file_path, SnakeConfigGet, sample_meta

logger = verbose_import(__name__, __doc__)
snake_config = SnakeConfigGet({})


def fetchMGs(annot_prefix: Path, out_dir: Path = None, threads=1, rerun=False):
    if out_dir is None:
        out_dir = annot_prefix.parent
    out_prefix = out_dir / annot_prefix.name
    if not rerun and list(out_prefix.glob("*.marker_genes_scores.table")):
        return out_prefix
    try:
        out_prefix.rmdir()
    except OSError:
        os.system(f"rm -rf {out_prefix}")
    os.system(
        f"~/software/fetchMGs/fetchMGs.pl "
        f"    -m extraction "
        f"    {annot_prefix}.faa "
        f"    -o {out_prefix} "
        f"    -d {annot_prefix}.fna "
        f"    -v -i "
        f"    -t {threads}"
    )
    if not list(out_prefix.glob("*.marker_genes_scores.table")):
        raise Exception(f"fetchMGs failed for {annot_prefix}")
    return out_prefix


def collect_align(COG: str, fetchMG_dirs: list[Path], output_dir: Path):
    clust_prefix = output_dir / f"{COG}"
    seqs = []
    for fetchMG_dir in fetchMG_dirs:
        try:
            record = SeqIO.read(fetchMG_dir / f"{COG}.faa", "fasta")
        except ValueError:
            continue
        record.id = fetchMG_dir.name
        seqs.append(record)
    SeqIO.write(seqs, f"{clust_prefix}.faa", "fasta")
    return Path(f"{clust_prefix}.faa")


def align_trim(
    handles=tuple(),
    seqs=tuple(),
    name=Path("align_trim"),
    tempdir: Path = None,
):
    prefix = tempdir / name.name
    SeqIO.write(
        (
            *(seq for f in handles for seq in SeqIO.parse(f, "fasta")),
            *(seq for seq in seqs),
        ),
        f"{prefix}.faa",
        "fasta",
    )
    os.system(
        f"mafft "
        f"    --maxiterate 1000 --localpair "
        f"    {prefix}.faa "
        f">   {prefix}.afa"
    )
    return f"{prefix}.afa"


def marker2tree(genomes, align_dir: Path, output_prefix: Path):
    cat_seqs = {genome: "" for genome in genomes}

    with tempfile.TemporaryDirectory() as tempdir_:
        tempdir = Path(tempdir_)
        for out in Parallel(n_jobs=THREADS, verbose=20)(
            delayed(align_trim)(
                seqs=[i for i in SeqIO.parse(faa, "fasta") if i.id in genomes],
                name=faa,
                tempdir=tempdir,
            )
            for faa in align_dir.glob("*.faa")
        ):
            align = AlignIO.read(out, "fasta")
            tmp = {k: "" for k in cat_seqs.keys()}
            for record in align:
                cat_seqs[record.id] += record.seq
                tmp.pop(record.id)
            for id in tmp:
                cat_seqs[id] += "-" * len(record.seq)

    with open(f"{output_prefix}.afa", "w") as fout:
        for genome, seq in cat_seqs.items():
            print(f">{genome}\n{seq}", file=fout)
    os.system(
        f"trimal -in {output_prefix}.afa "
        f"       -out {output_prefix}.afa.trimal -automated1"
    )
    tre_file = f"{output_prefix}.afa.trimal.fasttree"
    os.system(f"FastTree -gamma {output_prefix}.afa.trimal > {tre_file}")
    return tre_file


def fetchMG_2_alignment_dir(
    annot_prefixs: Iterable[Path], fetchMG_dir: Path, alignments_dir: Path
):
    fetchMG_dirs: list[Path] = Parallel(THREADS, verbose=20)(
        delayed(fetchMGs)(annot_prefix, fetchMG_dir) for annot_prefix in annot_prefixs
    )

    COGs = [i.name[:-4] for i in fetchMG_dirs[0].glob("*.faa")]
    COG_files: list[Path] = Parallel(THREADS, verbose=20)(
        delayed(collect_align)(COG, fetchMG_dirs, alignments_dir) for COG in COGs
    )
    markers_: dict[str, dict[str, int]] = {}
    for file in COG_files:
        markers_[file.name] = {}
        for record in SeqIO.parse(file, "fasta"):
            markers_[file.name][record.id] = 1
    markers = pd.DataFrame(markers_).fillna(0)
    return markers


def get_cross(wtdb: Path, relative_abundance_path: Path):
    Wtdb = pd.read_csv(wtdb)
    relative_abundance = pd.read_csv(relative_abundance_path, sep="\t").merge(
        Wtdb[["genome", "cluster"]].rename(
            {"genome": "Genome", "cluster": "Cluster"}, axis=1
        )
    )

    Wtdb_cross = (
        relative_abundance.merge(sample_meta)
        .groupby("Genome")["Group"]
        .apply(lambda s: s.drop_duplicates())
        .reset_index("Genome")
        .pivot_table("Group", "Genome", "Group", lambda x: True, False)
    )
    return Wtdb_cross


def assign_table2itol_color(table2itol: pd.DataFrame):
    from collections import Counter

    table2itol_pc = table2itol.assign(
        p__=lambda df: df["Taxonomy"].apply(lambda taxon: taxon.rsplit(";", 7 - 2)[0]),
        c__=lambda df: df["Taxonomy"].apply(lambda taxon: taxon.rsplit(";", 7 - 3)[0]),
    )

    phylum_freq = table2itol_pc.groupby("p__")["p__"].agg(len).sort_values()
    show_phylum = {
        *phylum_freq.pipe(lambda s: s[s >= 5]).index,
        "d__Bacteria;p__SAR324",
        "d__Bacteria;p__Desulfobacterota_D",
        "Bacteria;Gemmatimonadota",
        "d__Bacteria;p__Zixibacteria",
    }
    phylum_freq.pipe(lambda s: s[s < 5])

    class_freq = table2itol_pc.groupby("c__")["p__"].agg(len).sort_values()
    Counter(class_freq)
    freq_class = [
        *class_freq.pipe(lambda s: s[s > 10]).index,
    ]
    freq_class_phylum = [c__.rsplit(";", 1)[0] for c__ in freq_class]

    # replace phylum to its only class
    single_class_phylums = (
        class_freq.reset_index()
        .assign(
            p__=lambda df: df["c__"].apply(lambda taxon: taxon.rsplit(";", 3 - 2)[0]),
        )
        .groupby("p__")
        .agg(len)["c__"]
        .pipe(lambda df: df[df == 1])
        .index
    )

    def assign_color(p__, c__):
        color_list = [i[3:] for i in c__.split(";")[:3]]
        if p__ not in show_phylum:
            return "others"
        if p__ in single_class_phylums:
            return ";".join(color_list)
        if p__ in freq_class_phylum:
            if c__ in freq_class:
                return ";".join(color_list)
            # return ";".join([*color_list[:2], "others"])
        return ";".join(color_list[:2])

    color = table2itol_pc.apply(
        lambda x: assign_color(x["p__"], x["c__"]),
        axis=1,
    )
    print(len(Counter(color)))
    print(pd.Series(Counter(color)).sort_index())
    return color


def collect_itol_table(relative_abundance: Path):
    """
    columns:
        genome, class, GC, genome size, compl.., conta.., cross
    """
    Wtdb = pd.read_csv(file_path.results("Wtdb"))
    Stdb = pd.read_csv(file_path.results("Stdb"))
    reference_genome_info = pd.read_csv(file_path.all_bins("ref_genome_info.csv"))
    Wtdb_cross = get_cross(file_path.results("Wtdb"), relative_abundance)

    Wtdb_cross_1 = (
        Wtdb["genome"]
        .pipe(lambda x: Wtdb_cross.loc[x, :])
        .assign(
            Bottom=lambda df: df.apply(
                lambda x: "cross"
                if (x["Bw"] and x["Bs"])
                else ("water" if x["Bw"] else ("sediment" if x["Bs"] else "")),
                axis=1,
            ),
            Slope=lambda df: df.apply(
                lambda x: "cross"
                if (x["Sw"] and x["Ss"])
                else ("water" if x["Sw"] else ("sediment" if x["Ss"] else "")),
                axis=1,
            ),
        )
        .reset_index()[["Genome", "Bottom", "Slope"]]
    )

    genome_info = Wtdb_cross_1.merge(
        Wtdb[["genome", "classification"]]
        .merge(
            Stdb[
                [
                    "genome",
                    "Completeness",
                    "Contamination",
                    "GenomeSize",
                    "GC",
                ]
            ],
        )
        .rename({"classification": "Taxonomy", "genome": "Genome"}, axis=1)
    )
    table2itol = (
        pd.concat(
            [
                genome_info,
                reference_genome_info.rename(
                    {"taxonomy": "Taxonomy", "genome": "Genome"}, axis=1
                )[["Genome", "Taxonomy", "GenomeSize", "GC"]],
            ]
        )
        .assign(**{"log10(GenomeSize)":lambda df: np.log10(df["GenomeSize"])})
        .assign(Color=lambda df: assign_table2itol_color(df))
        .assign(
            Domain=lambda df: df.apply(
                lambda x: x["Taxonomy"].split(";")[0][3:],
                axis=1,
            ),
            Label=lambda df: df.apply(
                lambda x: [i for i in x["Taxonomy"].split(";") if i[3:]][-1],
                axis=1,
            ),
        )
    )

    table2itol.to_csv(file_path.bins_tree("table2itol.csv"), index=False)
    return table2itol


THREADS = 14


def main(
    gene_prefix: pd.Series,
    reference_gene_prefix: pd.Series,
    fetchMG_dir: Path,
    alignments_dir: Path,
    Wtdb_cross: pd.DataFrame,
    tree: Path,
    bottom_or_slope: Literal["S", "B", "all"],
    min_exist_marker: int,
):

    markers = fetchMG_2_alignment_dir(
        (*gene_prefix, *reference_gene_prefix), fetchMG_dir, alignments_dir
    )

    def passed_genomes(min_markers):
        _ = markers.isna().apply(lambda x: sum(1 ^ x) < min_markers, axis=1)
        return set(_.index[1 ^ _])

    if bottom_or_slope == "all":
        genomes = passed_genomes(min_exist_marker)
    else:
        genomes_slope = [
            i
            for i in Wtdb_cross.pipe(
                lambda df: df.index[
                    df[f"{bottom_or_slope}s"] | df[f"{bottom_or_slope}w"]
                ]
            )
        ]
        genomes = [
            i
            for i in passed_genomes(min_exist_marker)
            if (i.endswith(".fna") or i in genomes_slope)
        ]

    marker2tree(genomes, alignments_dir, tree)
    return 0


@click.command()
@click.option("--loglevel", default="INFO", type=str, help="set level of logger")
@click.option(
    "--Wtdb", default=file_path.results("Wtdb"), type=Path, help="output Wtdb"
)
@click.option(
    "--gene",
    default=file_path.all_bins("collect_gene"),
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
    "--fetch-mg",
    default=file_path.bins_tree("fetchMGs"),
    type=Path,
    help="infered reference genomes",
)
@click.option(
    "--alignments",
    default=file_path.bins_tree("alignments"),
    type=Path,
    help="output ref genome gene map",
)
@click.option(
    "--relative-abundance",
    default=file_path.results("Wtdb.relative_abundance.tsv"),
    type=Path,
    help="output ref genome gene map",
)
@click.option(
    "--tree",
    required=True,
    # default=file_path.bins_tree("bottom_mid_20"),
    type=Path,
    help="output ref genome gene map",
)
@click.option(
    "--all",
    "bottom_or_slope",
    flag_value="all",
)
@click.option(
    "--bottom",
    "bottom_or_slope",
    flag_value="B",
)
@click.option(
    "--slope",
    "bottom_or_slope",
    flag_value="S",
)
@click.option("--min-exist-marker", default=20, type=click.IntRange(1, 40, clamp=True))
@click.option("-t", "--threads", default=snake_config.THREADS)
def run(
    loglevel: str,
    wtdb: Path,
    gene: Path,
    ref_gene: Path,
    fetch_mg: Path,
    alignments: Path,
    relative_abundance: Path,
    tree: Path,
    bottom_or_slope: Literal["S", "B", "all"],
    min_exist_marker: int,
    threads: int,
):
    logger.setLevel(level=loglevel.upper())  # info

    Wtdb = pd.read_csv(wtdb)
    gene_prefix = Wtdb["genome"].apply(lambda x: gene / x)
    reference_gene_prefix = [i.parent / i.name[:-4] for i in ref_gene.glob("*.faa")]

    Wtdb_cross = get_cross(wtdb, relative_abundance)

    fetch_mg.mkdir(parents=True, exist_ok=True)
    alignments.mkdir(parents=True, exist_ok=True)

    global THREADS
    THREADS = threads

    logger.warning(">>> job start")
    state = main(
        gene_prefix,
        reference_gene_prefix,
        fetch_mg,
        alignments,
        Wtdb_cross,
        tree,
        bottom_or_slope,
        min_exist_marker,
    )
    logger.warning(">>> job finish")
    if state == 0:
        logger.info("success!")


if __name__ == "__main__":
    run()
