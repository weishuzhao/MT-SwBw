# -*- coding: utf-8 -*-
"""
 * @Date: 2022-06-19 22:29:49
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-06-19 22:45:11
 * @FilePath: /2021_09-MT10kSW/workflow/utils/contigRenameCut.py
 * @Description:
"""

import argparse
from pathlib import Path
import sys
from typing import TextIO, Tuple

from Bio import SeqIO

from PyLib.PyLibTool.file_info import verbose_import

logger = verbose_import(__name__, __doc__)


def rename_cut(fi: TextIO, fo: TextIO, threshold=500, name: str = None):  # , report):
    discard_seqs, discard_bases = 0, 0
    header = f">{name}|" if name else ">"

    for line in SeqIO.parse(fi, "fasta"):
        if len(line.seq) >= threshold:
            fo.write(header + str(line.id) + "\n" + str(line.seq) + "\n")
        else:
            discard_seqs += 1
            discard_bases += len(line.seq)

    logger.warn(
        "    {seqs_n} seqs ({bases_n} bases) are discarded".format(
            seqs_n=discard_seqs, bases_n=discard_bases
        ),
    )

    return discard_seqs, discard_bases


def main(fi: TextIO, fo: TextIO, threshold, name):  # , report):
    rename_cut(fi, fo, threshold, name)

    fi.close()
    fo.close()

    return 0


def get_args() -> Tuple:
    parser = argparse.ArgumentParser(description=__doc__)
    set_args(parser)
    args = parser.parse_args()
    logger.setLevel(level=args.loglevel.upper())  # info

    _input = args.input
    _output = args.output
    threshold = args.threshold
    name = args.name

    input = Path(_input).expanduser().absolute()
    logger.info(f"intput: {input}")
    in_file = open(input)

    if _output == "output":
        out_file = sys.stdout
    else:
        output = Path(_output).expanduser().absolute()
        out_file = open(output, "w")
    logger.info(f"output: {output}")

    logger.info(f"threshold: {threshold}")
    return in_file, out_file, threshold, name


def set_args(parser: argparse.ArgumentParser):
    parser.add_argument(
        "--loglevel", default="INFO", type=str, help="set level of logger"
    )
    parser.add_argument(
        "-i", "--input", type=str, required=True, help="input file in FASTA format"
    )
    parser.add_argument(
        "-o", "--output", type=str, default="stdout", help="output file"
    )
    parser.add_argument(
        "threshold",
        default=500,
        type=int,
        nargs="?",
        help="threshold. " "Any sequence < threshold will be discard",
    )
    parser.add_argument(
        "name",
        default=None,
        type=str,
        nargs="?",
        help="sample name. " "Will be add forward contigId seperated by an '|'.",
    )


def run():
    args = get_args()

    logger.warning(">>> job start at " + now.strftime("%Y-%m-%d %H:%M:%S"))
    state = main(*args)
    logger.warning(">>> job run time: " + str(datetime.now() - now))

    if state == 0:
        logger.info("success!")


if __name__ == "__main__":
    run()
