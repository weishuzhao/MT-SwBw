# -*- coding: utf-8 -*-
"""
 * @Date: 2021-10-27 21:27:18
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-10-27 21:43:58
 * @FilePath: /2021_09-MT10kSW/Analyze/funct/ex_ab.py
 * @Description:
    Existance and abundance of given pathways
"""

import os
from typing import Iterable, Dict
from PyLib.biotool.kegg import KModule, load_KEGG_module_raw, module_from_brite
from PyLib.reader.iters import read_table


KEGG_DIR = os.path.expanduser('~/Data/Database2/KEGG')

SAMPLES = [
    "TY.040", "TY.041", "TY.044",
    "WQ.018", "WQ.021", "WQ.022", "WQ.023", "WQ.024",
    "YW.019", "YW.020", "YW.021", "YW.023"]
TRIM = 'sickle'
ASSEM = 'megahit'


def KO_abd_sample(method: str,
                  contig_subset: Iterable = {},
                  KO_subset: Iterable = {}
                  ) -> Dict[str, float]:
    # first load KO and genes
    KO_abd: Dict[str, float] = {}
    with open(GENE_LNK_FORMAT(method)) as tab_in:
        for line in read_table(tab_in):
            gene, KO, abd, contig = line
            if contig_subset and contig not in contig_subset:
                continue
            if KO_subset and KO not in KO_subset:
                continue
            KO_abd[KO] = KO_abd.get(KO, 0.0) + float(abd)
    return KO_abd


def main():
    module_name = f'Analyze/funct/module_name.tsv'
    abundance = f'Analyze/funct/abundance_1.tsv'

    # first load KEGG modules
    module_levels, modules = module_from_brite('br:ko00002',
                                               os.path.join(KEGG_DIR, 'brite', 'ko00002.json'),
                                               os.path.join(KEGG_DIR, 'module'))
    with open(module_name, 'w') as file_out:
        print('A\tB\tC\tentry\tname',
              *('\t'.join(module_level) for module_level in module_levels),
              sep='\n', file=file_out)
    with open(abundance, 'w') as abd_out:
        print('# sample',
              *(entry for entry, _ in modules),
              sep='\t', file=abd_out)
        for sample in SAMPLES:
            method = f"{sample}-{TRIM}-{ASSEM}"
            KO_abd = KO_abd_sample(method)
            #compare_local_online(method, modules, KO_abundance)
            print(method,
                  *pathway_cpl_abd(KO_abd, modules),
                  sep='\t', file=abd_out)
