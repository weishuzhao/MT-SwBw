# -*- coding: utf-8 -*-
"""
 * @Date: 2022-01-12 14:08:09
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-01-12 22:07:00
 * @FilePath: /Analyze/2022/01/12.py
 * @Description:
    Align genes
"""
# %%
# %%
import os
from sys import stdout

from BCBio import GFF
from Bio import Seq, SeqFeature, SeqIO, SeqRecord
from PyLib.biotool.download import (download, download_fna,
                                    retrieve_refseq_url_ls)
from PyLib.tool.fileIO import open_r

# %%
#download_fna("GCA_012928605.1")  # MTA1
#tmp = []
#tmp = tmp or retrieve_refseq_url_ls("GCA_012928605.1")

tmp1 = [
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/012/928/605/GCA_012928605.1_ASM1292860v1/GCA_012928605.1_ASM1292860v1_genomic.gff.gz",
]
for url in tmp1:
    download_file = download(url)
    gff_file_gz = os.path.abspath(download_file)

# %%
download_file = "GCA_012928605.1_ASM1292860v1_genomic.gff.gz"
gff_file_gz = os.path.abspath(download_file)

# %%
genome_file = "GCA_012928605.1.fna.gz"
with open_r(genome_file, mode="rt") as fi:
    genome = SeqIO.to_dict(SeqIO.parse(fi, "fasta"))

# %%
genes = {}

with open_r(gff_file_gz, mode="rt") as in_handle:
    rec: SeqRecord.SeqRecord = None
    for rec in GFF.parse(in_handle, base_dict=genome):
        r: SeqFeature.SeqFeature = None
        for r in rec.features:
            if r.type == "gene":
                cds: SeqFeature.SeqFeature = None
                for cds in r.sub_features:
                    description = "; ".join(cds.qualifiers.get("product", []))
                    genes[cds.id] = SeqRecord.SeqRecord(
                        cds.extract(genome[rec.id]).seq, id = cds.id,
                        description=description, annotations = cds.qualifiers)


# %%
SeqIO.write(genes.values(), "GCA_012928605.1_gene.fna", "fasta-2line")
# and next, run blast on these genes

# %%
#from PyLib.biotool.diamond import diamond_blastp, diamond_check_db
from PyLib.tool.shell import runsh_safe


def blastn_check_db(fastdir, fastafile, outdir):
    fasta_file = os.path.join(fastdir, fastafile)
    blastn_db = os.path.join(outdir, fastafile)
    if not os.access(blastn_db, os.F_OK):
        stdout, stderr = runsh_safe(
            f"makeblastdb -dbtype nucl "
            f"            -in {fasta_file} "
            f"            -out {blastn_db} "
        )
    return blastn_db


def blastn(fastafile, blast_db, blast_out, thread=1, **kwargs):
    """bash
    blastn -query GCA_012928605.1_gene.fna \
           -db    blast/TY.040-megahit.128 \
           -out   blast/TY.040-megahit.128.tsv \
           -outfmt 6 \
           -num_threads 12
    """
    if not os.access(blast_out, os.F_OK):
        stdout, stderr = runsh_safe(
            f"blastn -query {fastafile} "
            f"       -db {blast_db} "
            f"       -out {blast_out} "
            f"       -num_threads {thread} "
            f"       " + (" ".join((f"-{k} {v}" for k, v in kwargs.items())))
        )
    return blast_out

# %%
genome_MAGs = []
genome_dir = "/home/hwrn/Work/2021_12-g__pangenome/results/21-12-26/selected_taxons/140_1"
for file in os.listdir(genome_dir):
    if file.endswith(".fasta"):
        genome_MAGs.append(file)

#for genome in genome_MAGs:
#    blastn_db = blastn_check_db(genome_dir, genome, "blast")
#    blast_out = blastn("GCA_012928605.1_gene.fna", blastn_db, f"blast/{genome}.tsv",
#                       thread=12, outfmt=6)

# %%
matched_genes = {}

for genome in genome_MAGs:
    blast_out = f"blast/{genome}.tsv"
    with open(blast_out) as bo:
        for line in bo.readlines():
            values = line.split()
            gene, pident = values[0], float(values[2])
            matched_genes[gene] = max(matched_genes.get(gene, 0), pident)

# %%
min(matched_genes.values())  # 93.452

{genes[gene].description for gene in matched_genes}

# %%
with open("40168_2020_849_MOESM1_ESM.table.s2.txt") as fs2:
    s2 = fs2.read()
    miss_genes = {gene for gene in matched_genes if genes[gene].description not in s2}

# %%
s2.split("Table S3")[0].split("\n")[1:]

# %%
{genes[gene].description for gene in genes.keys() - matched_genes.keys()}

{gene for gene in genes.keys() if genes[gene].description == "MBL fold metallo-hydrolase"}
# %%
