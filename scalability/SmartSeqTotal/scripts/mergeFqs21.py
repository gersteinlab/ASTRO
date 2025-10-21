#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import csv
import gzip
import glob
from pathlib import Path
from typing import List, Iterable, TextIO


def int_to_base4_fixed(n: int, width: int = 13, alphabet: str = "ACGT") -> str:
    digits = ["A"] * width
    i = width - 1
    while n > 0 and i >= 0:
        n, r = divmod(n, 4)
        digits[i] = alphabet[r]
        i -= 1
    return "".join(digits)

def read_srr_order_from_tsv(tsv_path: Path, srr_col_name: str = "SRR Accession") -> List[str]:
    srr_order: List[str] = []
    with open(tsv_path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            srr = (row.get(srr_col_name) or "").strip()
            if srr:
                srr_order.append(srr)
    return srr_order

def iter_fastq_records(handle: TextIO):

    while True:
        header = handle.readline()
        if not header:
            break  # EOF
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()

        if not (seq and plus and qual):
            raise ValueError("FASTQ stru wrong")

        h = header.rstrip("\r\n")
        s = seq.rstrip("\r\n")
        p = plus.rstrip("\r\n")
        q = qual.rstrip("\r\n")

        yield (h, s, p, q)



def merge_and_rename(
    tsv_path: Path,
    fastq_dir: Path,
    out_fastq_path: Path,
    barcode_len: int = 13,
    alphabet: str = "ACGT",) -> None:

    srr_order = read_srr_order_from_tsv(tsv_path)

    with open(out_fastq_path, "wt") as fout:
        total_files = 0


        for srr in srr_order:
            fq_path = fastq_dir + '/' + srr + '_1.fastq'
            if not os.path.isfile(fq_path):
                print(f"no fastq SRR: {srr}", file=sys.stderr)
                sys.exit()

            total_reads = 0
            with open(fq_path, "r") as fin:
                for _, seq, plus, qual in iter_fastq_records(fin):
                    barcode = int_to_base4_fixed(total_reads, width=barcode_len, alphabet=alphabet)
                    new_header = f"@{total_files}_{total_reads}|:_:|{srr}:{barcode}"

                    fout.write(new_header + "\n")
                    fout.write(seq + "\n")
                    fout.write("+\n")
                    fout.write(qual + "\n")

                    total_reads += 1
            total_files += 1
            



tsv_path  = 'cache/SraRunTable_human_smartseq.tsv'
fastq_dir = 'fastq/'
out_fastq_path = 'result/combine.fq'
merge_and_rename(tsv_path, fastq_dir, out_fastq_path)

