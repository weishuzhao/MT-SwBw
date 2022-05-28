# -*- coding: utf-8 -*-
"""
 * @Date: 2021-09-12 21:45:24
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-12-21 20:20:52
 * @FilePath: /2021_09-MT10kSW/Analyze/drep/recover.py
 * @Description:
"""

import os
from typing import Dict

import drep.d_cluster.utils
import drep.d_evaluate
import drep.d_choose
from drep.WorkDirectory import WorkDirectory
import pandas as pd

from PyLib.biotool.fna_msg import seq_total_depth, statistic_fna
from PyLib.reader.read_outputs import fasta, jgi_depths
from PyLib.reader.iters import gtdbtk_iter, read_table


def get_CW_by_Sani(wd: WorkDirectory, S_ani=0.95):
    Ndb = wd.get_db("Ndb")
    Gdb = wd.get_db("genomeInformation")
    Cdb, c2ret = drep.d_cluster.utils._cluster_Ndb(Ndb, comp_method="ANImf", S_ani=0.95)
    #
    Sdb, Wdb = drep.d_choose.choose_winners(Cdb, Gdb)
    return Cdb, Wdb


def store_db_wrapper(wd: WorkDirectory, modified_dir, **kwargs):
    if not os.path.exists(modified_dir):
        os.makedirs(modified_dir)
    for name, db in kwargs.items():
        wd.data_tables[name] = os.path.join(modified_dir, name + ".csv")
        db.to_csv(wd.data_tables[name], index=False)


def evaluate_winners_wrapper(wd: WorkDirectory, S_ani=0.95):
    modified_dir = f"Analyze/drep/sa_{S_ani}"
    #
    Cdb, Wdb = get_CW_by_Sani(wd, S_ani)
    last_db = wd.data_tables.copy()
    #
    store_db_wrapper(wd, modified_dir, **{"Cdb": Cdb, "Wdb": Wdb})
    # evaluate_winners
    Widb = drep.d_evaluate.evaluate_winners()
    store_db_wrapper(wd, modified_dir, **{"Widb": Widb})
    #
    return Widb, last_db


def join_Widb_gtdb(wd: WorkDirectory):
    Wtdb = wd.get_db("Widb")
    #
    taxonomy = {genome: "" for genome in Wtdb["genome"]}
    #
    for marker_gene in ["ar122", "bac120"]:
        with open(f"04_bin/02_taxon/gtdbtk.{marker_gene}.summary.tsv") as tin:
            for value in gtdbtk_iter(tin):
                genome = value[0] + ".fna"
                classification = value[1][0]
                taxonomy[genome] = classification
    #with open("Analyze/comapre_bins/ZF_Bins_info.tsv") as tin:
    #    for value in read_table(tin):
    #        taxonomy[value[0].rstrip("*")] = value[8]
    #
    Wtdb["taxonomy"] = [taxonomy[genome] for genome in Wtdb["genome"]]
    #
    # join sample
    Cdb = wd.get_db("Cdb")
    Cdb["sample"] = [genome.split("-megahit")[0]
    #                 if genome.endswith(".fna") else genome.split("_")[0]
                     for genome in Cdb["genome"]]
    S_c = Cdb.pivot_table(index="secondary_cluster", columns="sample", values="genome",
                          aggfunc=len, fill_value=0)
    samples = S_c.columns
    Wtdb = Wtdb.join(S_c, on="cluster")
    #
    return Wtdb, samples


def Wtdb_2_krona(Wtdb, samples, filebasename, ifcount = True):
    with open(filebasename +  ".tsv", "w") as fout:
        for i, row in Wtdb.iterrows():
            cluster = row["cluster"]
            taxonomy = row["taxonomy"]
            samples_exist = row[samples]
            print(sum(samples_exist) if ifcount else 1,
                  *[taxon[3:] for taxon in taxonomy.split(";")],
                  cluster,
                  sep="\t", file=fout)
    os.system(f"ktImportText {filebasename}.tsv -o {filebasename}.krona.html")


def join_Widb_gtdb_1(wd: WorkDirectory):
    """ Try to get Wtdb at primary_cluster. """
    Wtdb = wd.get_db("Widb")
    #
    taxonomy = {genome: "" for genome in Wtdb["genome"]}
    #
    for marker_gene in ["ar122", "bac120"]:
        with open(f"04_bin/02_taxon/gtdbtk.{marker_gene}.summary.tsv") as tin:
            for value in gtdbtk_iter(tin):
                genome = value[0] + ".fna"
                classification = value[1][0]
                taxonomy[genome] = classification
    with open("Analyze/comapre_bins/ZF_Bins_info.tsv") as tin:
        for value in read_table(tin):
            taxonomy[value[0].rstrip("*")] = value[8]
    #
    Wtdb["primary_cluster"] = [int(primary_cluster.split("_")[0])
                               for primary_cluster in Wtdb["cluster"]]
    Wtdb = Wtdb.iloc[Wtdb.groupby('primary_cluster')['score'].idxmax()]
    #
    Wtdb["taxonomy"] = [taxonomy[genome] for genome in Wtdb["genome"]]
    #
    # join sample
    Cdb = wd.get_db("Cdb")
    Cdb["sample"] = [genome.split("-megahit")[0]
                     if genome.endswith(".fna") else genome.split("_")[0]
                     for genome in Cdb["genome"]]
    S_c = Cdb.pivot_table(index="primary_cluster", columns="sample", values="genome",
                          aggfunc=len, fill_value=0)
    samples = S_c.columns
    Wtdb = Wtdb.join(S_c, on="primary_cluster")
    #
    return Wtdb, samples


def print_Wtdb_info(Wtdb):
    phylum_dict: Dict[str, int] = {}
    for phylum in [taxon.split(";")[1][3:] for taxon in Wtdb["taxonomy"]]:
        phylum_dict[phylum] = phylum_dict.get(phylum, 0) + 1
    #
    for k, v in sorted(phylum_dict.items()):
        print(k, v, f"{(v/len(Wtdb['genome'])):.3f}", sep="\t")
    #
    taxon_levels = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    for i in range(7):
        number = sum([not taxon.split(";")[i][3:] for taxon in Wtdb["taxonomy"]])
        print(taxon_levels[i], number, f"{(number/len(Wtdb['genome'])):.3f}", sep="\t")


def get_depth(BinId: str, __file_depths = {}):
    BinId = BinId.split(os.path.sep)[-1]
    runname = BinId.rsplit(".", 2)[0]
    bin_filename = f"{PIPE_PATH}/{runname}/04_bin/DASTool_bins/{BinId}"
    depth_file = f"{PIPE_PATH}/{runname}/02_assem..{runname}-jgi.depth"

    if depth_file not in __file_depths:
        with open(depth_file) as fi:
            sample_list, ctg_depth = jgi_depths(fi)
            __file_depths[depth_file] = ctg_depth
            print(depth_file, "load to memory")
    ctg_depth = __file_depths[depth_file]
    #
    #fna_msg = statistic_fna(fasta(bin_filename))
    totalAvgDepth, depths = seq_total_depth(ctg_depth, fasta(bin_filename))
    #return list(fna_msg) + depths
    return depths


def join_SChdb(wd: WorkDirectory):
    Stdb = wd.get_db("Sdb").merge(wd.get_db("Cdb")).merge(wd.get_db("genomeInfo")).merge(
        wd.get_db("Chdb"), left_on="genome", right_on="Bin Id")[[
            "genome", "score", "secondary_cluster",
            "Completeness", "Contamination", "Strain heterogeneity", "length"
        ]]
    Stdb["sample"] = [genome.split("-megahit")[0] for genome in Stdb["genome"]]
    print(len(Stdb))
    Bin_states = [get_depth(BinId) for BinId in Stdb["genome"]]
    print(Bin_states[1])
    for i, key in enumerate(["AvgDepth"]):  #"SeqNumbers", "MaxLength", "GenomeSize", "GC", "N50", "L50", "AvgDepth"]):
        Stdb[key] = [values[i] for values in Bin_states]
        print(key)
    return Stdb


PIPE_PATH = "Pipe/"

modified_dir = f"Analyze/drep"


if __name__ == '__main__':
    wd = WorkDirectory("04_bin/02_drep")

    Wtdb = pd.read_csv("Analyze/drep/Wtdb.csv")
    samples = [
        "TY.040", "TY.041", "TY.044",
        "WQ.018", "WQ.021", "WQ.022", "WQ.023", "WQ.024",
        "YW.019", "YW.020", "YW.021", "YW.023"]
    #Wtdb, samples = join_Widb_gtdb(wd)
    #store_db_wrapper(wd, modified_dir, Wtdb=Wtdb)

    Stdb = join_SChdb(wd)
    store_db_wrapper(wd, modified_dir, Stdb=Stdb)


    #filebasename = os.path.join(modified_dir, "bins")
    #Wtdb_2_krona(Wtdb, samples, filebasename, False)

    #print_Wtdb_info(Wtdb)

    #Wtdb_1, _ = join_Widb_gtdb_1(wd)
    #store_db_wrapper(wd, modified_dir, Wtdb_1=Wtdb_1)
