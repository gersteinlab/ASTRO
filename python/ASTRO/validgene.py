#!/usr/bin/env python3
import sys
import subprocess
import re
import math
import os
import multiprocessing as mp
from statsmodels.stats.rates import test_poisson_2indep


def alignment_end(pos, cigar):
    """
    Calculate the end position of an alignment from start position and CIGAR string.

    Parses a CIGAR string to determine how many reference bases are consumed
    by the alignment, then calculates the end position.

    Args:
        pos (int): Start position of alignment (1-based)
        cigar (str): CIGAR string describing alignment operations

    Returns:
        int: End position of alignment on reference (1-based, inclusive)

    Note:
        Only operations that consume reference bases (M, D, N, =, X) are counted.
        Insertions (I), soft clips (S), hard clips (H), and padding (P) are ignored.
    """
    ref_len = 0
    matches = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    for length_str, op in matches:
        length = int(length_str)
        if op in ["M", "D", "N", "=", "X"]:
            ref_len += length
    return pos + ref_len - 1


def count_reads_end_in_region(bam, chrom, start, end, bkstart, bkend, bkstart0, bkend0):
    region = f"{chrom}:{start}-{end}"
    proc = subprocess.Popen(
        ["samtools", "view", bam, region], stdout=subprocess.PIPE, text=True
    )
    countbegin = 0
    countend = 0
    bkcountbegin = 0
    bkcountend = 0
    bkcountbegin0 = 0
    bkcountend0 = 0
    for line in proc.stdout:
        line = line.strip()
        if not line:
            continue
        fields = line.split("\t")
        pos = int(fields[3])
        cigar_str = fields[5]
        aln_begin = pos
        aln_end = alignment_end(pos, cigar_str)

        if start <= aln_begin <= end and start <= aln_end <= end:
            countbegin += 1
            countend += 1
        else:
            if bkstart <= aln_begin < start:
                bkcountbegin += 1
            if bkstart0 <= aln_begin < start:
                bkcountbegin0 += 1
            if end < aln_end <= bkend:
                bkcountend += 1
            if end < aln_end <= bkend0:
                bkcountend0 += 1
        if start <= aln_begin <= end and start <= aln_end <= end:
            countbegin += 1
            countend += 1
        else:
            if bkstart <= aln_begin < start:
                bkcountbegin += 1
            if bkstart0 <= aln_begin < start:
                bkcountbegin0 += 1
            if end < aln_end <= bkend:
                bkcountend += 1
            if end < aln_end <= bkend0:
                bkcountend0 += 1

    proc.stdout.close()
    proc.wait()
    return countbegin, countend, bkcountbegin, bkcountend, bkcountbegin0, bkcountend0


def main():
    if len(sys.argv) < 3:
        print(
            "Usage: python count_bam_end_in_region.py <bedfile> <bamfile>",
            file=sys.stderr,
        )
        sys.exit(1)

    bedfile = sys.argv[1]
    bamfile = sys.argv[2]
    outputprefix = sys.argv[3]
    if len(sys.argv) == 5:
        hangout = int(sys.argv[4])
    else:
        hangout = 5
    outputtrue = outputprefix + ".map"
    outputfalse = outputprefix + ".unmap"
    with open(bedfile, "r") as bf, open(outputtrue, "w") as tm, open(
        outputfalse, "w"
    ) as fm:
        for line in bf:
            line = line.strip()
            parts = line.split()

            chrom = parts[0]
            start = int(parts[1]) + 1 - hangout
            end = int(parts[2]) + hangout
            bedlen = end - start + 1

            bkstart = start - math.floor(bedlen / 2)
            bkend = end + math.ceil(bedlen / 2)
            bkstart0 = start - hangout
            bkend0 = end + hangout
            countbegin, countend, _, _, bkcountbegin0, bkcountend0 = (
                count_reads_end_in_region(
                    bamfile, chrom, start, end, bkstart, bkend, bkstart0, bkend0
                )
            )
            ratiobegin = countbegin / bedlen
            ratiobkbegin0 = bkcountbegin0 / (start - bkstart0)
            # print(line)
            ratioend = countend / bedlen
            ratiobkend0 = bkcountend0 / (bkend0 - end)
            # print((ratiobegin, ratioend, ratiobkbegin0, ratiobkend0, countbegin, countend, bkcountbegin0, bkcountend0))
            if ratiobegin > 2 * ratiobkbegin0:
                # print(line)
                # print((countbegin, bedlen,bkcountbegin0, start - bkstart0))
                # print((countend, bedlen,bkcountbegin0, bkend0 - end))
                res = test_poisson_2indep(
                    count1=countbegin + 1e-7,
                    exposure1=bedlen,
                    count2=bkcountbegin0 + 1e-7,
                    exposure2=start - bkstart0,
                    compare="diff",
                    method="score",
                    alternative="larger",
                )
                p_value = res.pvalue
                # print(p_value)
                if p_value < 0.05:
                    ratioend = countend / bedlen
                    ratiobkend0 = bkcountend0 / (bkend0 - end)
                    if ratioend > 2 * ratiobkend0:
                        res = test_poisson_2indep(
                            count1=countend + 1e-7,
                            exposure1=bedlen,
                            count2=bkcountend0 + 1e-7,
                            exposure2=bkend0 - end,
                            compare="diff",
                            method="score",
                            alternative="larger",
                        )
                        p_value = res.pvalue
                        # print(p_value)
                        if p_value < 0.05:
                            tm.write(line + "\n")
                            continue
            # ratiobkbegin = bkcountbegin / (start - bkstart  + bkend - end)
            # if ratiobegin > 1.5 * ratiobkbegin:
            #    table = np.array([[countbegin, bkcountbegin], [bedlen, start - bkstart]])
            #    _, p_value = fisher_exact(table, alternative='greater')
            #    if p_value < 0.001:
            #        ratioend = countend / bedlen
            #        ratiobkend = bkcountend / bedlen
            #        if ratioend > 1.5 * ratiobkend:
            #            table = np.array([[countend, bkcountend], [bedlen, bkend - end ]])
            #            _, p_value = fisher_exact(table, alternative='greater')
            #            if p_value < 0.001:
            #                tm.write(line + '\n')
            #                continue
            fm.write(line + "\n")


def getvalidedgtf(gtfin, outputfolder, genes2check, hangout=3):
    genes2checkdict = {}
    with open(genes2check, "r") as bf:
        for line in bf:
            line = line.strip()
            genes2checkdict[line] = 1
    bamfile = os.path.join(outputfolder, "STAR/tempfiltered.bam")
    if not os.path.exists(bamfile + ".bai"):
        print(bamfile + ".bai")
        subprocess.run(["samtools", "index", bamfile])
    proc = subprocess.Popen(
        ["samtools", "view", "-H", bamfile], stdout=subprocess.PIPE, text=True
    )
    existchromosomes = []
    for line in proc.stdout:
        if line.startswith("@SQ"):
            fields = line.strip().split("\t")
            chrom = fields[1].replace("SN:", "")
            existchromosomes.append(chrom)
    proc.stdout.close()
    proc.wait()

    validgenegtf = os.path.join(outputfolder, "validgene.gtf")
    with open(gtfin, "r") as gtf, open(validgenegtf, "w") as tm:
        for line in gtf:
            line = line.strip()
            # print(line)
            parts = line.split("\t")
            # print(parts)
            # sys.exit(1)
            genename = parts[8]
            if genename in genes2checkdict:
                chrom = parts[0]
                start = int(parts[3]) + 1 - hangout
                end = int(parts[4]) + hangout
                bedlen = end - start + 1

                bkstart = start - math.floor(bedlen / 2)
                bkend = end + math.ceil(bedlen / 2)
                bkstart0 = start - hangout
                bkend0 = end + hangout
                countbegin, countend, _, _, bkcountbegin0, bkcountend0 = (
                    count_reads_end_in_region(
                        bamfile, chrom, start, end, bkstart, bkend, bkstart0, bkend0
                    )
                )
                ratiobegin = countbegin / bedlen
                ratiobkbegin0 = bkcountbegin0 / (start - bkstart0)

                ratioend = countend / bedlen
                ratiobkend0 = bkcountend0 / (bkend0 - end)
                if ratiobegin > 2 * ratiobkbegin0:
                    res = test_poisson_2indep(
                        count1=countbegin + 1e-7,
                        exposure1=bedlen,
                        count2=bkcountbegin0 + 1e-7,
                        exposure2=start - bkstart0,
                        compare="diff",
                        method="score",
                        alternative="larger",
                    )
                    p_value = res.pvalue
                    if p_value < 0.05:
                        ratioend = countend / bedlen
                        ratiobkend0 = bkcountend0 / (bkend0 - end)
                        if ratioend > 2 * ratiobkend0:
                            res = test_poisson_2indep(
                                count1=countend + 1e-7,
                                exposure1=bedlen,
                                count2=bkcountend0 + 1e-7,
                                exposure2=bkend0 - end,
                                compare="diff",
                                method="score",
                                alternative="larger",
                            )
                            p_value = res.pvalue
                            if p_value < 0.05:
                                tm.write(line + "\n")
                                continue
            else:
                tm.write(line + "\n")
    return validgenegtf


def getvalidedgtf_worker(line, genes2checkdict, bamfile, existchromosomes, hangout):
    line = line.strip()
    parts = line.split("\t")
    genename = parts[8]
    chrom = parts[0]
    if genename in genes2checkdict:
        if chrom not in existchromosomes:
            return None
        start = int(parts[3]) + 1 - hangout
        end = int(parts[4]) + hangout
        bedlen = end - start + 1

        bkstart = start - math.floor(bedlen / 2)
        bkend = end + math.ceil(bedlen / 2)
        bkstart0 = start - hangout
        bkend0 = end + hangout
        countbegin, countend, _, _, bkcountbegin0, bkcountend0 = (
            count_reads_end_in_region(
                bamfile, chrom, start, end, bkstart, bkend, bkstart0, bkend0
            )
        )
        ratiobegin = countbegin / bedlen
        ratiobkbegin0 = bkcountbegin0 / (start - bkstart0)

        ratioend = countend / bedlen
        ratiobkend0 = bkcountend0 / (bkend0 - end)
        if ratiobegin > 2 * ratiobkbegin0:
            res = test_poisson_2indep(
                count1=countbegin + 1e-7,
                exposure1=bedlen,
                count2=bkcountbegin0 + 1e-7,
                exposure2=start - bkstart0,
                compare="diff",
                method="score",
                alternative="larger",
            )
            p_value = res.pvalue
            if p_value < 0.05:
                ratioend = countend / bedlen
                ratiobkend0 = bkcountend0 / (bkend0 - end)
                if ratioend > 2 * ratiobkend0:
                    res = test_poisson_2indep(
                        count1=countend + 1e-7,
                        exposure1=bedlen,
                        count2=bkcountend0 + 1e-7,
                        exposure2=bkend0 - end,
                        compare="diff",
                        method="score",
                        alternative="larger",
                    )
                    p_value = res.pvalue
                    if p_value < 0.05:
                        return line
    else:
        return line


def getvalidedgtf_parallel(gtfin, outputfolder, genes2check, hangout=5, threadsnum=4):
    """
    Validate suspicious genes in GTF annotation using statistical analysis of read alignments.

    This function performs advanced filtering of potentially problematic genes (e.g., miRNAs, piRNAs)
    by analyzing their read coverage patterns compared to flanking control regions. Uses statistical
    tests to determine if genes show anomalous enrichment that might indicate spurious mapping.

    Args:
        gtfin (str): Path to input GTF annotation file
        outputfolder (str): Output directory containing BAM files
        genes2check (str): Path to file listing gene names/IDs to validate
        hangout (int, optional): Size of flanking control regions (kb). Defaults to 5.
        threadsnum (int, optional): Number of parallel threads. Defaults to 4.

    Returns:
        None: Creates filtered GTF files in the output directory

    Creates:
        - Modified GTF files with suspicious genes removed
        - Statistical analysis reports for gene validation

    Note:
        Requires indexed BAM file (tempfiltered.bam) from previous alignment step.
        Gene names in genes2check file must match GTF annotation entries.
    """
    genes2checkdict = {}
    with open(genes2check, "r") as bf:
        for line in bf:
            line = line.strip()
            genes2checkdict[line] = 1

    sortedbam = os.path.join(outputfolder, "STAR/tempfiltered.sorted.bam")
    bamfile   = os.path.join(outputfolder, "STAR/tempfiltered.bam")

    if not os.path.isfile(sortedbam):
        subprocess.run(["samtools", "sort", "-@", str(threadsnum), "-o", sortedbam, bamfile], stdout=subprocess.DEVNULL, check=True)

    if not os.path.exists(sortedbam+".bai"):
        subprocess.run(["samtools", "index", sortedbam], stdout=subprocess.DEVNULL)
    proc = subprocess.Popen(
        ["samtools", "view", "-H", sortedbam],
        stdout=subprocess.PIPE,
        text=True
    )
    existchromosomes = []
    for line in proc.stdout:
        if line.startswith("@SQ"):
            fields = line.strip().split("\t")
            chrom = fields[1].replace("SN:", "")
            existchromosomes.append(chrom)
    proc.stdout.close()
    proc.wait()

    with open(gtfin, "r") as f:
        lines = f.readlines()

    pool = mp.Pool(int(threadsnum))
    results = pool.starmap(
        getvalidedgtf_worker,
        [(line, genes2checkdict, bamfile, existchromosomes, hangout) for line in lines],
    )
    pool.close()
    pool.join()

    validgenegtf = os.path.join(outputfolder, "validgene.gtf")
    with open(validgenegtf, "w") as outfile:
        for res in results:
            if res:
                outfile.write(res + "\n")
    return validgenegtf


if __name__ == "__main__":
    main()
