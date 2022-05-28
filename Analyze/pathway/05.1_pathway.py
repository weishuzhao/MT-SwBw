# -*- coding: utf-8 -*-
"""
 * @Date: 2021-06-15 09:49:37
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-02-10 10:46:57
 * @FilePath: /2022_01-ZFMG/home/hwrn/Work/2021_09-MT10kSW/Analyze/pathway/05.1_pathway.py
 * @Description:
"""

import os
from typing import Iterable, Dict, List, Tuple

import pandas as pd

from PyLib.biotool.kegg import KModule, load_KEGG_module_raw, module_from_brite
from PyLib.reader.iters import DASTool_scaffolds2bin_iter, read_table

# from PyLib.tool.path import appendcwd


KEGG_DIR = "~/Data/Database2/KEGG"

sample_meta = pd.read_csv("00_data/sample_meta.tsv", sep="\t", header=0, na_values="——")
SAMPLES = sample_meta["Station"]
fo_KO_RPb = f"Analyze/pathway/KO_sample_RPb.tsv"


def GENE_LNK_FORMAT(method):
    """
    header like: #gene\tKO\tcontig\t...bam
    """
    return f"Pipe/{method}/03_annot/03_annot..{method}.tpm"


def KO_RPb_sample(
    method: str, contig_subset: Iterable = {}, KO_subset: Iterable = {}
) -> Dict[str, float]:
    # first load KO and genes
    KO_RPb: Dict[str, float] = {}
    with open(GENE_LNK_FORMAT(method)) as tab_in:
        for line in read_table(tab_in):
            gene, KO, contig, RPb = line
            if contig_subset and contig not in contig_subset:
                continue
            if KO_subset and KO not in KO_subset:
                continue
            KO_RPb[KO] = KO_RPb.get(KO, 0.0) + float(RPb)
    return KO_RPb


def pathway_cpl_RPb(KO_RPb: Dict[str, float], modules: List[Tuple[str, KModule]]):
    KO_RPb = {KO: RPb for KO, RPb in dict(KO_RPb).items() if RPb > 0}
    return (
        module.abundance(KO_RPb) * (1 if int(module.completeness(KO_RPb)) else -1)
        for _, module in modules
    )


def main():
    fo_module_name = f"Analyze/pathway/module_name.tsv"
    fo_KO_RPb = f"Analyze/pathway/KO_sample_RPb.tsv"
    fo_abd_cpl = f"Analyze/pathway/abd_cpl.tsv"
    # first load KEGG modules
    module_levels, modules = module_from_brite(
        "br:ko00002",
        os.path.join(KEGG_DIR, "brite", "ko00002.json"),
        os.path.join(KEGG_DIR, "module"),
    )
    with open(fo_module_name, "w") as file_out:
        print(
            "A\tB\tC\tentry\tname",
            *("\t".join(module_level) for module_level in module_levels),
            sep="\n",
            file=file_out,
        )
    KO_sample_RPb = pd.DataFrame()
    for sample in SAMPLES:
        KO_RPb = KO_RPb_sample(f"{sample}-megahit")
        KO_sample_RPb[sample] = pd.Series(KO_RPb)
    KO_sample_RPb.fillna(0, inplace=True)
    KO_sample_RPb.to_csv(fo_KO_RPb, sep="\t", index_label="KO")
    # with open(fo_abd_cpl, 'w') as RPb_out:
    #    print('# sample',
    #          *(entry for entry, _ in modules),
    #          sep='\t', file=RPb_out)
    #    for sample in SAMPLES:
    #        print(sample,
    #              *pathway_cpl_RPb(KO_sample_RPb[sample], modules),
    #              sep='\t', file=RPb_out)


def main3():
    mag_pathway_file = "Analyze/pathway/MAG_KO.tsv"

    module_levels, modules = module_from_brite(
        "br:ko00002",
        os.path.join(KEGG_DIR, "brite", "ko00002.json"),
        os.path.join(KEGG_DIR, "module"),
    )

    # with open(module_name, 'w') as file_out:
    i = 0
    with open(mag_pathway_file, "w") as RPb_out:
        print("# sample", *(entry for entry, _ in modules), sep="\t", file=RPb_out)
        for sample in SAMPLES:
            method = f"{sample}-{TRIM}-{ASSEM}"
            with open(SCAFFOLD2BIN_FORMAT(method)) as s2b_in:
                for MAG, scaffolds in DASTool_scaffolds2bin_iter(s2b_in):
                    MAG = f"{method}|{MAG}"
                    if "log":
                        i += 1
                        print(MAG, i, end="\r")
                    KO_RPb = KO_RPb_sample(method, contig_subset=scaffolds)
                    print(
                        MAG, *pathway_cpl_RPb(KO_RPb, modules), sep="\t", file=RPb_out
                    )


def main2():
    modules: List[Tuple[str, KModule]] = []
    focus_ko_abundance: Dict[str, List[float]] = {}
    for module_name in (
        "M00307",
        "M00580",
        "M00308",
        "M00165",
        "M00166",
        "M00358",
        "M00174",
        "M00346",
        "M00345",
        "M00935",
        "M00528",
        "M00804",
        "M00596",
    ):
        raw_module = load_KEGG_module_raw(module_name, KEGG_DIR)
        modules.append(
            (
                module_name,
                KModule(
                    "".join(raw_module["DEFINITION"]),
                    additional_info="".join(raw_module["NAME"]),
                ),
            )
        )
        for KO in modules[-1][1].list_ko():
            focus_ko_abundance[KO] = [0.0 for _ in SAMPLES]

    for i, sample in enumerate(SAMPLES):
        method = f"{sample}-{TRIM}-{ASSEM}"
        KO_abundance = KO_RPb_sample(method, focus_ko_abundance)
        for KO in focus_ko_abundance:
            focus_ko_abundance[KO][i] = KO_abundance.get(KO, 0.0)

    focus = f"Analyze/pathway/diff_KO_tpm.tsv"
    with open(focus, "w") as file_out:
        print("KO", "module", *SAMPLES, sep="\t", file=file_out)
        for module_name, module in modules:
            for KO in module.list_ko():
                print(KO, module_name, *focus_ko_abundance[KO], sep="\t", file=file_out)


def main1():
    fo_abd_cpl = f"Analyze/pathway/abd_cpl.tsv"
    # first load KEGG modules
    module_levels, modules = module_from_brite(
        "br:ko00002",
        os.path.join(KEGG_DIR, "brite", "ko00002.json"),
        os.path.join(KEGG_DIR, "module"),
    )
    KO_sample_RPb = pd.read_csv(fo_KO_RPb, sep="\t", index_col="KO")
    with open(fo_abd_cpl, "w") as RPb_out:
        print("# sample", *(entry for entry, _ in modules), sep="\t", file=RPb_out)
        for sample in KO_sample_RPb.columns:
            print(
                sample,
                *pathway_cpl_RPb(KO_sample_RPb[sample], modules),
                sep="\t",
                file=RPb_out,
            )


def report_abd_cpl(KO_sample_RPb, module_name):
    def avail_RPb(sample):
        return {
            k: v
            for k, v in dict(KO_sample_RPb[sample]).items()
            if v > 0 and not pd.isna(v)
        }

    from PyLib.biotool.kegg import KModule, load_KEGG_module, map_KO_dict, map_KO_substr

    KEGG_DIR = "~/Data/Database2/KEGG"
    SAMPLES = KO_sample_RPb.columns
    module: KModule = load_KEGG_module(module_name, os.path.join(KEGG_DIR, "module"))
    paths = {
        path
        for path in module.all_paths()
        if path.replace(" ", "")
        in {
            pi.replace(" ", "")
            for sample in SAMPLES
            for pi in module.all_paths(avail_RPb(sample))
        }
    }
    module_def = str(module)
    print(module_def, "|", module_name, "; ".join(module.additional_info["NAME"]))
    for i, path in enumerate(sorted(paths)):
        print(map_KO_substr(module_def, path), i + 1)
    for sample in SAMPLES:
        print(map_KO_dict(module_def, avail_RPb(sample)), "|", sample)


def console():
    import pandas as pd

    fo_KO_RPb = f"Analyze/pathway/KO_sample_RPb.tsv"
    fo_mrk = f"Analyze/annot/marker/mkr_stat.tsv"
    mrk_RPb = pd.read_csv(fo_mrk, sep="\t")["mrk.med"]
    KO_sample_RPb = pd.read_csv(fo_KO_RPb, sep="\t", index_col="KO")
    for sample in KO_sample_RPb.columns:
        KO_sample_RPb[sample] = KO_sample_RPb[sample] / mrk_RPb[sample]
    module_names = """M00028, M00844, M00015
        M00029
        M00040, M00533, M00038
        M00545
        M00551, M00623
        M00144, M00149, M00150, M00151, M00417, M00157, M00159
        M00153, M00156
        M00790
        M00019, M00535, M00432, M00036
        M00168, M00579
        M00165, M00166
        M00010, M00004, M00006, M00008
        M00307, M00308
        M00127, M00896, M00125, M00124, M00916, M00115, M00622, M00914, M00123, M00140, M00846, M00868, M00121, M00924, M00122, M00117, M00116
        M00141, M00926, M00847, M00930
        M00899, M00925
        M00035
        M00627, M00745
        M00642
        M00082, M00087
        M00085
        M00077
        M00026
        M00088, M00091
        M00066, M00100
        M00064
        M00016, M00526, M00527, M00031
        M00378
        M00358
        M00531, M00530
        M00528, M00804
        M00027, M00118
        M00631, M00061, M00632, M00855, M00549, M00909, M00012, M00741
        M00014, M00129, M00373
        M00597
        M00048, M00049, M00050
        M00051, M00052, M00938, M00046, M00939
        M00020, M00033
        M00595
        M00596""".replace(
        ",", " "
    ).split()
    for module_name in module_names:
        report_abd_cpl(KO_sample_RPb, module_name)


if __name__ == "__main__":
    console()  # in remote
    # main1() in local
