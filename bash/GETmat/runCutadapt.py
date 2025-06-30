import re
import argparse
import os
import gzip
import subprocess


def simple_fastq_iterator(handle):
    """
    Parse FASTQ format files and yield individual sequence records.

    Args:
        handle: File handle object for reading FASTQ data

    Yields:
        tuple: (title, sequence, quality_string) for each FASTQ record
            - title (str): Sequence identifier without '@' prefix
            - sequence (str): DNA/RNA sequence string
            - quality_string (str): Quality scores string
    """
    while True:
        title_line = handle.readline().rstrip()
        seq_string = handle.readline().rstrip()
        sep_line = handle.readline().rstrip()
        qual_string = handle.readline().rstrip()
        if not (title_line or seq_string or sep_line or qual_string):
            break
        yield title_line[1:], seq_string, qual_string


def cutfastq(filein, forward, length):
    """
    Trim FASTQ sequences to a specified length from either 5' or 3' end.

    Creates a temporary file during processing and replaces the original file.
    Handles both gzip-compressed and uncompressed FASTQ files.

    Args:
        filein (str): Path to input FASTQ file (can be .gz compressed)
        forward (bool): If True, trim from 5' end (keep first 'length' bases)
                       If False, trim from 3' end (keep last 'length' bases)
        length (int): Number of bases to keep after trimming

    Returns:
        None: Modifies the input file in place
    """
    with (
        gzip.open(filein, "rt") if filein.endswith(".gz") else open(filein, "rt")
    ) as in_handle:
        with open(filein + ".temp", "w") as out_handle:
            for title, seq, qual in simple_fastq_iterator(in_handle):
                if forward:
                    new_seq = seq[0:length]
                    new_qual = qual[0:length]
                    out_handle.write("@%s\n%s\n+\n%s\n" % (title, new_seq, new_qual))
                else:
                    new_seq = seq[-length : len(seq)]
                    new_qual = qual[-length : len(qual)]
                    out_handle.write("@%s\n%s\n+\n%s\n" % (title, new_seq, new_qual))
    os.remove(filein)
    os.rename(filein + ".temp", filein)


def combineFqs(outputfile, fqs):
    """
    Concatenate sequences from multiple FASTQ files into a single output file.

    Reads corresponding entries from multiple FASTQ files simultaneously and
    concatenates their sequences and quality scores. All input files must have
    the same number of reads in the same order.

    Args:
        outputfile (str): Path for the output combined FASTQ file (unused in current implementation)
        fqs (list): List of paths to input FASTQ files to be combined

    Returns:
        None: Creates 'index.fastq' file in the global outputfolder directory

    Note:
        The function creates a file named 'index.fastq' in the outputfolder directory,
        not using the outputfile parameter.
    """
    iterators0 = [
        gzip.open(file_path, "rt") if file_path.endswith(".gz") else open(file_path)
        for file_path in fqs
    ]
    iterators = [simple_fastq_iterator(file) for file in iterators0]
    with open(outputfolder + "index.fastq", "w") as out_handle:
        for entries in zip(*iterators):
            newseq = ""
            newqual = ""
            newtitle = ""
            for title, seq, qual in entries:
                newseq += seq
                newqual += qual
            out_handle.write("@%s\n%s\n+\n%s\n" % (title, newseq, newqual))


ap = argparse.ArgumentParser()
ap.add_argument("-b", "--barcodestr", required=True, help="Workpath with the BAM file")
ap.add_argument("-o", "--outputfolder", required=True, help="Path to input BAM file")
ap.add_argument("-i", "--inputfa", required=True, help="Path to input BAM file")
ap.add_argument(
    "-t", "--thread", required=False, default="8", type=str, help="Number of threads"
)
args = vars(ap.parse_args())


barcodestr = args["barcodestr"]
outputfolder = args["outputfolder"]
remainfa = args["inputfa"]
threadnum = args["thread"]

outputfolder = re.sub("/+$", "", outputfolder) + "/"
array4input = re.split("_", barcodestr)

matchstr = [bool(re.search(r"^[ACGTacgtnN]+$", s)) for s in array4input]
matchnum = [bool(re.search(r"^([0-9]*|b)$", s)) for s in array4input]
if (sum(matchnum) + sum(matchstr)) != len(matchstr):
    exit("Error: wrong barcode format: unexpected code")

matchstri = [i for i, x in enumerate(matchstr) if x]
matchnumi = [i for i, x in enumerate(matchnum) if x]

matchstrdiff = [matchstri[i + 1] - matchstri[i] for i in range(len(matchstri) - 1)]
matchnumdiff = [matchnumi[i + 1] - matchnumi[i] for i in range(len(matchnumi) - 1)]

if any(element != 2 for element in matchstrdiff):
    exit("Error: wrong barcode format: order error")
if any(element != 2 for element in matchnumdiff):
    exit("Error: wrong barcode format: order error")


customsetting = "-j " + threadnum

tmpi = 0
finalfqs = []
for stri in matchstri:
    thestr = [array4input[stri]]
    if (stri - 1) in matchnumi:
        tmpfaname = outputfolder + "temp_" + str(tmpi) + ".fastq"
        outputcommand = (
            ["cutadapt"]
            + [customsetting]
            + ["-a"]
            + thestr
            + ["-o"]
            + [tmpfaname]
            + [remainfa]
        )
        finalfqs.append(tmpfaname)
        subprocess.run(" ".join(outputcommand), shell=True)
        tmpi = tmpi + 1
    tmpfaname = outputfolder + "temp_" + str(tmpi) + ".fastq"
    outputcommand = (
        ["cutadapt"]
        + [customsetting]
        + ["-g"]
        + thestr
        + ["-o"]
        + [tmpfaname]
        + [remainfa]
    )
    remainfa = tmpfaname
    tmpi = tmpi + 1
    subprocess.run(" ".join(outputcommand), shell=True)


if matchnumi[-1]:
    finalfqs.append(tmpfaname)
print("Barcode Fastqs")
print(finalfqs)

for i, x in enumerate(matchnumi):
    if array4input[x] != "b":
        designedlength = int(array4input[x])
        forwardvalue = False
        if x == len(array4input) - 1:
            forwardvalue = True
            if designedlength != [0]:
                cutfastq(finalfqs[i], forwardvalue, designedlength)
        elif x == 0:
            if designedlength != [0]:
                cutfastq(finalfqs[i], forwardvalue, designedlength)

all_reads = set()


for fastq_file in finalfqs:
    current_reads = set()
    with (
        gzip.open(fastq_file, "rt") if fastq_file.endswith(".gz") else open(fastq_file)
    ) as f:
        for line in f:
            if line.startswith("@") or line.startswith(">"):
                read_name = line.split()[0][1:]
                if read_name in current_reads:
                    raise ValueError(
                        f"Duplicate read name found: {read_name} in file {fastq_file}"
                    )
                current_reads.add(read_name)
        if len(all_reads) == 0:
            all_reads = current_reads
        else:
            all_reads = all_reads & current_reads


for fastq_file in finalfqs:
    with (
        gzip.open(fastq_file, "rt") if fastq_file.endswith(".gz") else open(fastq_file)
    ) as f:
        with open(fastq_file + ".temp", "w") as out_handle:
            for title, seq, qual in simple_fastq_iterator(f):
                title = title.split()[0]
                if title in all_reads:
                    out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    os.rename(fastq_file + ".temp", fastq_file)


del all_reads


linkfqs = []
umifqs = []
for i, x in enumerate(matchnumi):
    if array4input[x] == "b":
        linkfqs.append(finalfqs[i])
    else:
        umifqs.append(finalfqs[i])

if len(linkfqs) >= 2:
    combineFqs(outputfolder + "index.fastq", linkfqs)
else:
    if len(linkfqs) == 1:
        os.rename(linkfqs[0], outputfolder + "index.fastq")

if len(umifqs) >= 2:
    combineFqs(outputfolder + "UMI.fastq", umifqs)
else:
    if len(umifqs) == 1:
        os.rename(umifqs[0], outputfolder + "UMI.fastq")
