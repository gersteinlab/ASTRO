#!/usr/bin/env python3
import sys
import subprocess
import re
import numpy as np
import math
from statsmodels.stats.rates import test_poisson_2indep

#python3 /home/dz287/mark/SHARE/scripts/validgene.py hsamirhigh.bed ~/mark/SHARE/data/LM0623P100/tempfiltered.bam  LM0623a
def alignment_end(pos, cigar):
    ref_len = 0
    matches = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    for length_str, op in matches:
        length = int(length_str)
        if op in ['M', 'D', 'N', '=', 'X']:
            ref_len += length
    return pos + ref_len - 1


def count_reads_end_in_region(bam, chrom, start, end, bkstart, bkend, bkstart0, bkend0):
    region = f"{chrom}:{start}-{end}"
    proc = subprocess.Popen(
        ["samtools", "view", bam, region],
        stdout=subprocess.PIPE,
        text=True
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
        print("Usage: python count_bam_end_in_region.py <bedfile> <bamfile>", file=sys.stderr)
        sys.exit(1)

    bedfile = sys.argv[1]
    bamfile = sys.argv[2]
    outputprefix = sys.argv[3]
    if (len(sys.argv) == 5):
        hangout = int(sys.argv[4])
    else:
        hangout = 5
    outputtrue = outputprefix + '.map'
    outputfalse = outputprefix + '.unmap'
    with open(bedfile, "r") as bf, open(outputtrue, "w") as tm, open(outputfalse, "w") as fm:
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
            countbegin, countend, _, _, bkcountbegin0, bkcountend0 = count_reads_end_in_region(bamfile, chrom, start, end, bkstart, bkend, bkstart0, bkend0)
            ratiobegin = countbegin / bedlen
            ratiobkbegin0 = bkcountbegin0 / (start - bkstart0)
            #print(line)
            ratioend = countend / bedlen
            ratiobkend0 = bkcountend0 / (bkend0 - end)
            #print((ratiobegin, ratioend, ratiobkbegin0, ratiobkend0, countbegin, countend, bkcountbegin0, bkcountend0))
            if ratiobegin > 2 * ratiobkbegin0:
                #print(line)
                #print((countbegin, bedlen,bkcountbegin0, start - bkstart0))
                #print((countend, bedlen,bkcountbegin0, bkend0 - end))
                res = test_poisson_2indep(count1=countbegin+1e-7, exposure1=bedlen,count2=bkcountbegin0+1e-7, exposure2=start - bkstart0, compare='diff', method = 'score', alternative='larger')
                p_value = res.pvalue
                #print(p_value)
                if p_value < 0.05:
                    ratioend = countend / bedlen
                    ratiobkend0 = bkcountend0 / (bkend0 - end)
                    if ratioend > 2 * ratiobkend0:
                        res = test_poisson_2indep(count1=countend+1e-7, exposure1=bedlen,count2=bkcountend0+1e-7, exposure2=bkend0 - end, compare='diff', method = 'score', alternative='larger')
                        p_value = res.pvalue
                        #print(p_value)
                        if p_value < 0.05:
                            tm.write(line + '\n')
                            continue
            #ratiobkbegin = bkcountbegin / (start - bkstart  + bkend - end)
            #if ratiobegin > 1.5 * ratiobkbegin:
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
            fm.write(line + '\n')
            


if __name__ == "__main__":
    main()