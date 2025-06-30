import os
import sys
from collections import defaultdict
from statistics import variance
from math import log2
from .countfeature import genemat2tsv


def sum_by_group(values, groups):
    """
    Calculate sum of values grouped by category.

    Groups numeric values by their corresponding group identifiers and
    returns the sum for each group.

    Args:
        values (iterable): Numeric values to be summed
        groups (iterable): Group identifiers corresponding to each value

    Returns:
        dict_values: Collection of group sums
    """
    sums = defaultdict(float)
    for val, grp in zip(values, groups):
        sums[grp] += val
    return sums.values()


def removedims(inputtsv, outputtsv, filterlogratio):
    """
    Remove genes and barcodes with abnormal variance from expression matrix.

    Filters out genes (rows) and barcodes (columns) that show extreme variance
    compared to the overall distribution. Uses log2 ratio thresholding to
    identify and remove outliers.

    Args:
        inputtsv (str): Path to input expression matrix (TSV format)
        outputtsv (str): Path to output filtered expression matrix
        filterlogratio (str or float): Log2 ratio threshold for filtering
                                     (e.g., "2" means remove if variance > 4x median)

    Returns:
        None: Creates filtered output file and a .rt file with removal statistics

    Creates:
        - outputtsv: Filtered expression matrix
        - outputtsv.rt: Report of removed genes/barcodes
    """
    filterlogratio = float(filterlogratio)
    with open(inputtsv, "r") as f, open(outputtsv, "w") as fo, open(
        outputtsv + ".rt", "w"
    ) as fr:
        line1 = next(f)
        fo.write(line1)
        header = line1.rstrip("\n").split("\t")
        header = [item.split("x") for item in header[1:]]
        rowi, coli = zip(*[(int(item[0]), int(item[1])) for item in header])
        for line in f:
            line = line.rstrip("\n")
            cols = line.split("\t")
            values = list(map(int, cols[1:]))
            colvs = sum_by_group(values, coli)
            rowvs = sum_by_group(values, rowi)
            colvar = variance(colvs)
            rowvar = variance(rowvs)
            if colvar == 0 or rowvar == 0:
                continue
            else:
                logratio = abs(log2(rowvar / colvar))
                fr.write(cols[0] + "\t" + str(logratio) + "\n")
                if logratio < filterlogratio:
                    fo.write(line + "\n")


def mergetsv(inputtsv1, inputtsv2, outputtsv):
    dictA = {}
    header = []
    with open(inputtsv2, "r") as f:
        line1 = next(f)
        header = line1.rstrip("\n").split("\t")
        for line in f:
            line = line.rstrip("\n")
            cols = line.split("\t")
            key = cols[0]
            values = list(map(int, cols[1:]))
            dictA[key] = values
    with open(inputtsv1, "r") as f, open(outputtsv, "w") as fo:
        large_header = f.readline().rstrip("\n").split("\t")
        if large_header != header:
            sys.stderr.write(
                f"Error, {inputtsv1} and {inputtsv2} have different column names\n"
            )
        fo.write("\t".join(header) + "\n")
        for line in f:
            line = line.rstrip("\n")
            cols = line.split("\t")
            key = cols[0]
            values = list(map(int, cols[1:]))
            if key in dictA:
                old_values = dictA.pop(key)
                values = [ov + nv for ov, nv in zip(old_values, values)]
                values_str = list(map(str, values))
                line_out = key + "\t" + "\t".join(values_str) + "\n"
                fo.write(line_out)
            else:
                values_str = list(map(str, values))
                line_out = key + "\t" + "\t".join(values_str) + "\n"
                fo.write(line_out)
        for key, values in dictA.items():
            values_str = list(map(str, values))
            line_out = key + "\t" + "\t".join(values_str) + "\n"
            fo.write(line_out)


def filtMATbyRT(expmatgood, expmatbad, finalexpmat, filterlogratio=2):
    expmatbad2 = expmatbad + "-2"
    removedims(expmatbad, expmatbad2, filterlogratio)
    mergetsv(expmatgood, expmatbad2, finalexpmat)


def featurefilter(gtffile, options, barcodes_file, filterlogratio, outputfolder):
    """
    Apply feature filtering to gene expression matrix.

    Performs quality control filtering on the gene expression matrix by:
    1. Processing excluded features (if any)
    2. Merging standard and excluded expression data
    3. Removing genes/barcodes with abnormal variance patterns
    4. Generating final filtered expression matrices

    Args:
        gtffile (str): Path to GTF gene annotation file
        options (str): Processing options:
            - 'H': Hard mode - affects gene assignment strategy
        barcodes_file (str): Path to barcode coordinate mapping file
        filterlogratio (str or float): Log2 ratio threshold for variance filtering
        outputfolder (str): Output directory path

    Returns:
        None: Creates filtered expression matrices in the output directory

    Creates:
        - finalexpmat.tsv: Final filtered gene expression matrix
        - expmat.tsv.excl: Expression matrix for excluded features
        - Various intermediate processing files
    """
    do_easy_mode = 1
    if "H" in options:
        do_easy_mode = 0
    expmatbedexcl = os.path.join(outputfolder, "expmat.bed.excl")
    expmattsvexcel = os.path.join(outputfolder, "expmat.tsv.excl")
    expmattsv = os.path.join(outputfolder, "expmat.tsv")
    finalexpmattsv = os.path.join(outputfolder, "finalexpmat.tsv")
    genemat2tsv(
        input_file=expmatbedexcl,
        output_file=expmattsvexcel,
        barcodes_file=barcodes_file,
        gtf_file=gtffile,
        filter_str="0:0",
        easy_mode=do_easy_mode,
    )
    filtMATbyRT(expmattsv, expmattsvexcel, finalexpmattsv, filterlogratio)
    if "auto_barcode.tsv" in barcodes_file:

        auto_barcode_path = os.path.join(outputfolder, "temps", "auto_barcode.tsv")
        if os.path.exists(auto_barcode_path):
            i2bc = {}
            with open(auto_barcode_path, "r") as bf:
                for line in bf:
                    line = line.strip()
                    if not line:
                        continue
                    bc_str, i_str, _ = line.split("\t")
                    i2bc[i_str] = bc_str

            tmp_renamed = finalexpmattsv + ".renamed"
            with open(finalexpmattsv, "r") as fin, open(tmp_renamed, "w") as fout:
                lines = fin.readlines()
                if not lines:
                    return
                header = lines[0].rstrip("\n").split("\t")
                new_header = []
                for col in header:
                    if col == "gene":
                        new_header.append(col)
                    else:
                        part = col.split("x")
                        if len(part) == 2 and part[0] == part[1] and part[0] in i2bc:
                            new_header.append(i2bc[part[0]])
                        else:
                            new_header.append(col)
                fout.write("\t".join(new_header) + "\n")
                for line in lines[1:]:
                    fout.write(line)

            os.remove(finalexpmattsv)
            os.rename(tmp_renamed, finalexpmattsv)
