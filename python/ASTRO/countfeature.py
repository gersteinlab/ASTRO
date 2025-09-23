import os
import subprocess
import re
import math
import multiprocessing
import logging
from .validgene import getvalidedgtf_parallel

def inter_bed2geneFile(
    inputbed, outputfile, gtffile, keep_temp=False, log_function=print
):
    """
    Intersect BED intervals with GTF annotations to assign reads to genes.

    Uses bedtools intersect to find overlaps between genomic intervals (from BED file)
    and gene annotations (from GTF file). For each interval, assigns the gene with
    the largest overlap or highest score based on the annotation attributes.

    Args:
        inputbed (str): Path to input BED file with genomic intervals
        outputfile (str): Path to output file for gene assignments
        gtffile (str): Path to GTF annotation file
        keep_temp (bool, optional): Whether to keep temporary files. Defaults to False.
        log_function (callable, optional): Function for logging messages. Defaults to print.

    Returns:
        None: Creates output file with gene assignments

    Output format:
        Tab-separated file with columns for interval coordinates, gene assignments,
        and overlap information.
    """

    cmd_intersect = ["bedtools", "intersect", "-a", inputbed, "-b", gtffile, "-wao"]
    proc = subprocess.Popen(cmd_intersect, stdout=subprocess.PIPE, text=True)

    with open(outputfile, "w") as f_out:
        infovector0 = ""
        biggestvalue = 0.0
        geneoutputall = ""

        for line in proc.stdout:
            cols = line.rstrip("\n").split("\t")
            if cols[-1] == "0":
                continue
            parts = [
                cols[0],
                cols[1],
                cols[2],
                cols[3],
                cols[5],
                cols[-7],
                cols[-6],
                cols[-2],
                cols[-1],
            ]
            # sys.exit(parts)
            if not line:
                continue
            if len(parts) < 9:
                continue

            chrname = parts[0]
            read_begin = int(parts[1])
            read_end = int(parts[2])
            readname = parts[3]
            strand = parts[4]
            gtf_start = int(parts[5])
            gtf_end = int(parts[6])
            gtf_name = parts[7]
            overlap_length = int(parts[8])

            read_length = read_end - read_begin
            gtf_length = gtf_end - gtf_start + 1
            valuetempit = (overlap_length - (read_length - overlap_length)) / float(
                gtf_length
            )
            # = (2 * overlap_length - read_length) / gtf_length

            match_obj = re.match(r"^(\d+)_(\d+)_(\d+):(.+)$", readname)
            if not match_obj:
                raise RuntimeError(f"wrong read name format: {readname}")

            readID = match_obj.group(1)
            Coord1 = match_obj.group(2)
            Coord2 = match_obj.group(3)
            the_rest = match_obj.group(4)

            splitted = the_rest.split(":")
            if len(splitted) > 0:
                splitted.pop(0)  # shift
            if len(splitted) >= 2:
                ASscore = splitted[0]
                fastq_read_length = splitted[1]
            elif len(splitted) == 1:
                ASscore = splitted[0]
                fastq_read_length = "."
            else:
                ASscore = "."
                fastq_read_length = "."

            if infovector0 == "":
                # first line
                biggestvalue = valuetempit
                geneoutputall = (
                    f"{Coord1}\t{Coord2}\t{gtf_name}\t{ASscore}\t{fastq_read_length}\t"
                    f"{overlap_length}\t{gtf_length}\t{read_length}\t{readID}"
                )
                infovector0 = readID
            elif readID != infovector0:
                # previous line
                f_out.write(geneoutputall + "\n")
                biggestvalue = valuetempit
                geneoutputall = (
                    f"{Coord1}\t{Coord2}\t{gtf_name}\t{ASscore}\t{fastq_read_length}\t"
                    f"{overlap_length}\t{gtf_length}\t{read_length}\t{readID}"
                )
                infovector0 = readID
            else:
                # same readID
                if valuetempit > biggestvalue:
                    biggestvalue = valuetempit
                    geneoutputall = (
                        f"{Coord1}\t{Coord2}\t{gtf_name}\t{ASscore}\t{fastq_read_length}\t"
                        f"{overlap_length}\t{gtf_length}\t{read_length}\t{readID}"
                    )
                elif abs(valuetempit - biggestvalue) < 1e-12:
                    geneout = (
                        f"{Coord1}\t{Coord2}\t{gtf_name}\t{ASscore}\t{fastq_read_length}\t"
                        f"{overlap_length}\t{gtf_length}\t{read_length}"
                    )
                    geneoutputall += f"\t{geneout}"

        if geneoutputall:
            f_out.write(geneoutputall + "\n")

    # if not keep_temp:
    #    try:
    #        os.remove(interbed)
    #    except OSError:
    #        pass


def bam_to_bed(input_bam, output_bed):
    with open(output_bed, "w") as outf:
        ret = subprocess.run(
            ["bedtools", "bamtobed", "-i", input_bam, "-split"], stdout=outf, check=True
        )
    if ret.returncode != 0:
        raise RuntimeError(
            f"[bedtools_bamtobed] bedtools bamtobed failed on {input_bam}"
        )


def collate_bam(input_bam, output_bam, threads=1):
    result = subprocess.run(
        ["samtools", "collate", "-@", str(threads), "-o", output_bam, input_bam],
        check=True,
    )
    if result.returncode != 0:
        raise RuntimeError("[samtools_collate] samtools collate failed.")


def load_barcodes(barcodes_file):
    locations = []
    with open(barcodes_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            x = parts[1]
            y = parts[2]
            locations.append(f"{x}_{y}")
    return locations


def load_genes_from_gtf(gtf_file):
    gene2num = {}
    num2gene = {}
    idx = 1
    with open(gtf_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            # GTF column 9 is gene
            gene = fields[8]
            if gene not in gene2num:
                gene2num[gene] = idx
                num2gene[idx] = gene
                idx += 1
    return gene2num, num2gene


def load_genes_from_input(input_file):
    gene_set = set()
    with open(input_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            fieldsin = line.split("\t")
            if len(fieldsin) >= 9:
                fieldsin.pop(8)

            i = 0
            while i + 7 < len(fieldsin):
                genename = fieldsin[i + 2]
                gene_set.add(genename)
                i += 8

    # assign gene2num, num2gene
    gene2num = {}
    num2gene = {}
    idx = 1
    for g in sorted(gene_set):
        gene2num[g] = idx
        num2gene[idx] = g
        idx += 1
    return gene2num, num2gene


def parse_input_build_matrix_exactly(
    input_file, gene2num, num2gene, as_threshold, ratio_threshold, do_filter, easy_mode
):
    from collections import defaultdict

    matrix = defaultdict(int)

    if do_filter:
        fo = open(input_file + ".excl", "w")
    with open(input_file, "r") as fh:

        for line in fh:
            line = line.strip()
            if not line:
                continue
            fieldsin = line.split("\t")
            # if(@fieldsin >=9){ splice(@fieldsin,8,1);}
            if len(fieldsin) >= 9:
                fieldsin.pop(8)

            genelist = []
            last_x = None
            last_y = None

            i = 0
            while i + 7 < len(fieldsin):
                x = fieldsin[i]
                y = fieldsin[i + 1]
                genename = fieldsin[i + 2]
                ASscore = fieldsin[i + 3]
                readlen = fieldsin[i + 4]
                overlaplen = fieldsin[i + 5]
                genelen = fieldsin[i + 6]
                read_length = fieldsin[i + 7]
                i += 8

                skip = False
                if do_filter and ASscore != ".":
                    try:
                        as_val = float(ASscore)
                        glen_val = float(genelen)
                        if (
                            as_val <= as_threshold
                            and as_val <= ratio_threshold * glen_val
                        ):
                            skip = True
                    except ValueError:
                        pass

                if skip:
                    fo.write(line + "\n")
                    continue

                # gene -> ID
                if genename in gene2num:
                    gid = gene2num[genename]
                    genelist.append(gid)

                last_x = x
                last_y = y

            if not genelist:
                continue

            # judge single / multi
            if len(genelist) == 1:
                g1 = genelist[0]
                key = f"{last_x}_{last_y}_{g1}"
                matrix[key] += 1
            else:
                genelist = sorted(set(genelist))
                if easy_mode:
                    # select first one
                    g1 = genelist[0]
                    key = f"{last_x}_{last_y}_{g1}"
                    matrix[key] += 1
                else:
                    # connect with '-'
                    joined = "-".join(str(g) for g in genelist)
                    key = f"{last_x}_{last_y}_{joined}"
                    matrix[key] += 1
    if do_filter:
        fo.close()
    return matrix
def summarize_matrix_with_per_pixel(matrix):
    genes_per_pixel = {}
    count_per_pixel = {}

    for key, cnt in matrix.items():
        pos, genes = key.rsplit('_', 1)
        count_per_pixel[pos] = count_per_pixel.get(pos, 0) + cnt
        if '-' not in genes:
            if pos not in genes_per_pixel:
                genes_per_pixel[pos] = set()
            genes_per_pixel[pos].add(genes)

    genes_per_stats = {pos: len(gset) for pos, gset in genes_per_pixel.items()}
    return genes_per_stats, count_per_pixel
def featurebed2mattsv(
    input_file, 
    output_file, 
    barcodes_file, 
    gtf_file=None, 
    filter_str="25:0.75", 
    easy_mode=False
):

    do_filter = False
    as_threshold = 25.0
    ratio_threshold = 0.75
    if filter_str:
        parts = filter_str.split(":")
        if len(parts) == 2:
            try:
                as_threshold = float(parts[0])
                ratio_threshold = float(parts[1])
                do_filter = True
            except ValueError:
                pass
        if as_threshold == 0 and ratio_threshold == 0:
            do_filter = False

    locations = load_barcodes(barcodes_file)

    if gtf_file:
        gene2num, num2gene = load_genes_from_gtf(gtf_file)
    else:
        gene2num, num2gene = load_genes_from_input(input_file)

    matrix = parse_input_build_matrix_exactly(
        input_file=input_file,
        gene2num=gene2num,
        num2gene=num2gene,
        as_threshold=as_threshold,
        ratio_threshold=ratio_threshold,
        do_filter=do_filter,
        easy_mode=easy_mode,
    )

    import statistics as stats
    genes_per_stats, count_per_pixel = summarize_matrix_with_per_pixel(matrix)

    gene_counts = list(genes_per_stats.values())
    read_counts = list(count_per_pixel.values())

    maxgene = max(gene_counts)
    mingene = min(gene_counts)
    meangene = stats.mean(gene_counts)
    stdgene = stats.stdev(gene_counts) if len(gene_counts) > 1 else 0

    maxcount = max(read_counts)
    mincount = min(read_counts)
    meancount = stats.mean(read_counts)
    stdcount = stats.stdev(read_counts) if len(read_counts) > 1 else 0

    logging.info(f"Max number of gene features over all positions: {maxgene}.")
    logging.info(f"Min number of gene features over all positions: {mingene}.")
    logging.info(f"Average number of gene features over all positions: {meangene}.")
    logging.info(f"STD of gene features over all positions: {stdgene}.\n")

    logging.info(f"Max number of UMI counts over all positions: {maxcount}.")
    logging.info(f"Min number of UMI counts over all positions: {mincount}.")
    logging.info(f"Average number of UMI counts over all positions: {meancount}.")
    logging.info(f"STD of UMI counts over all positions: {stdcount}.\n")
   

    used_gene_keys = set()
    for k in matrix.keys():
        arr = k.split("_")
        if len(arr) < 3:
            continue
        gene_part = arr[-1]
        if "-" in gene_part:
            used_gene_keys.add(gene_part)
        else:
            try:
                used_gene_keys.add(int(gene_part))
            except ValueError:
                used_gene_keys.add(gene_part)

    single_gid = []
    multi_gid = []

    for item in used_gene_keys:
        if isinstance(item, int):
            single_gid.append(item)
        elif isinstance(item, str):
            if "-" in item:
                multi_gid.append(item)
            else:
                try:
                    val = int(item)
                    single_gid.append(val)
                except ValueError:
                    multi_gid.append(item)

    single_gid.sort()

    def parse_multi_ids(mg_str):
        return [int(x) for x in mg_str.split("-")]

    multi_gid.sort(key=lambda x: parse_multi_ids(x))

    with open(output_file, 'w') as fo:
        header_elems = ["gene"] + [loc.replace('_', 'x') for loc in locations]

        

        fo.write("\t".join(header_elems) + "\n")

        for gid in single_gid:
            if gid in num2gene:
                gene_name = num2gene[gid]
            else:
                gene_name = str(gid)
            row_counts = []
            for loc in locations:
                key = f"{loc}_{gid}"
                cnt = matrix.get(key, 0)
                row_counts.append(str(cnt))
            fo.write(gene_name + "\t" + "\t".join(row_counts) + "\n")

        for mg in multi_gid:
            sub_ids = mg.split("-")
            sub_names = []
            for sid_str in sub_ids:
                try:
                    sid = int(sid_str)
                    if sid in num2gene:
                        sub_names.append(num2gene[sid])
                    else:
                        sub_names.append(sid_str)
                except ValueError:
                    sub_names.append(sid_str)
            mg_name = "+".join(sub_names)
            row_counts = []
            for loc in locations:
                key = f"{loc}_{mg}"
                cnt = matrix.get(key, 0)
                row_counts.append(str(cnt))
            fo.write(mg_name + "\t" + "\t".join(row_counts) + "\n")

    
    


def estimate_total_lines(input_bam):
    p_view = subprocess.Popen(
        ["samtools", "view", input_bam], stdout=subprocess.PIPE, text=True
    )
    total_lines = 0
    for _ in p_view.stdout:
        total_lines += 1
    p_view.wait()
    return total_lines


def collate_and_split(input_bam, out_prefix, threads):
    target_size = math.floor(estimate_total_lines(input_bam) / int(threads))

    cmd_collate = [
        "samtools",
        "collate",
        "-O",
        "--output-fmt",
        "SAM",
        "-@",
        str(threads),
        input_bam,
    ]
    proc_collate = subprocess.Popen(cmd_collate, stdout=subprocess.PIPE, text=True)

    chunk_files = []
    chunk_idx = 0
    proc_write = None
    lines_in_chunk = 0

    def start_new_chunk(idx):
        chunk_bam = f"{out_prefix}.chunk{idx}.bam"
        cmd_write = ["samtools", "view", "-b", "-o", chunk_bam, "-"]
        p = subprocess.Popen(cmd_write, stdin=subprocess.PIPE, text=True)
        chunk_files.append(chunk_bam)
        # print(f"start_new_chunk: {chunk_bam}")
        return p

    header_lines = []
    for line in proc_collate.stdout:
        line = line.rstrip("\n")
        if not line.startswith("@"):
            break
        header_lines.append(line)

    prev_qname = None
    lines_in_chunk = 1
    chunk_index = 0
    proc_write = start_new_chunk(chunk_index)
    for header_line in header_lines:
        proc_write.stdin.write(header_line + "\n")
    proc_write.stdin.write(line + "\n")

    for line in proc_collate.stdout:
        line = line.rstrip("\n")
        if lines_in_chunk >= target_size:
            qname = line.split("\t")[0]
            if prev_qname and prev_qname != qname:
                # print(f"chunk {prev_qname} {qname}")
                # print(f"{line}")
                proc_write.stdin.close()
                proc_write.wait()
                chunk_index += 1
                proc_write = start_new_chunk(chunk_index)
                for header_line in header_lines:
                    proc_write.stdin.write(header_line + "\n")
                lines_in_chunk = 1
                prev_qname = None
            elif prev_qname == None:
                prev_qname = qname

        lines_in_chunk += 1
        proc_write.stdin.write(line + "\n")

    if proc_write:
        proc_write.stdin.close()
        proc_write.wait()

    retc = proc_collate.wait()
    if retc != 0:
        raise RuntimeError("samtools collate failed")
    return chunk_files


def worker(chunk_bam, gtffile):
    chunk_bed = chunk_bam.replace(".bam", ".bed")
    chunk_expmatbed = chunk_bam.replace(".bam", ".expmat.bed")
    bam_to_bed(chunk_bam, chunk_bed)
    inter_bed2geneFile(chunk_bed, chunk_expmatbed, gtffile)
    os.remove(chunk_bed)
    return chunk_expmatbed

def parse_star_log(filepath):
    metrics = {}
    with open(filepath, "r") as f:
        for line in f:
            if "|" not in line:
                continue
            key, val = [x.strip() for x in line.split("|", 1)]
            clean_val = val.rstrip("%").strip()
            try:
                value = float(clean_val) if "." in clean_val else int(clean_val)
            except ValueError:
                value = clean_val
            metrics[key] = value
    return metrics


def countfeature(
    gtffile, threadnum, options, barcodes_file, outputfolder, qualityfilter, genes2check=False
):
    """
    Generate gene expression count matrix from aligned BAM files.

    This function processes aligned reads to create a gene expression matrix by:
    1. Converting BAM alignments to BED intervals
    2. Intersecting intervals with GTF gene annotations
    3. Counting reads per gene per barcode
    4. Applying quality filters and generating final expression matrix

    Args:
        gtffile (str): Path to GTF gene annotation file
        threadnum (int): Number of threads for parallel processing
        options (str): Processing options:
            - 'H': Hard mode - report multi-gene alignments as hyphenated entries
            - Other options for different counting strategies
        barcodes_file (str): Path to barcode coordinate mapping file
        outputfolder (str): Output directory path
        qualityfilter (str): Quality filter specification (e.g., "25:0.75")
                           Format: "min_alignment_score:min_coverage_fraction"

    Returns:
        None: Creates expression matrix files in the output directory

    Creates:
        - expmat.bed: BED file with gene assignments
        - expmat.tsv: Gene expression count matrix
        - Intermediate processing files in STAR/ directory
    """

    for h in logging.root.handlers[:]:
        logging.root.removeHandler(h)
        
    logfilename = os.path.join(outputfolder, ".logs/countfeature.log")
    os.makedirs(os.path.dirname(logfilename), exist_ok=True)
    logging.basicConfig(filename=logfilename, filemode="w", level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")

    logging.info("countfeature step starts.\n")

    if genes2check:
        logging.info("gene feature check step starts.\n")
        usedgtf = getvalidedgtf_parallel(gtfin=gtffile, outputfolder=outputfolder, genes2check=genes2check, hangout=5, threadsnum=threadnum)
        def count_lines(fname):
            with open(fname, "r") as f:
                return sum(1 for _ in f)

        n1 = count_lines(gtffile)
        n2 = count_lines(usedgtf)
        logging.info(f"Get valid GTF: {usedgtf}. check step ends.\n")
        logging.info("gene feature check step ends.\n")
        gtffile = usedgtf

    with open(os.path.join(outputfolder, "STAR/tempLog.final.out"), "r") as f:
        content = f.read()
        logging.info(f"\nload STAR genome mapping result...\nSTAR genome mapping stat:\n{content}\n")

    STAR_Log_final_out = parse_star_log(os.path.join(outputfolder, "STAR/tempLog.final.out"))
    input_read_number = STAR_Log_final_out['Number of input reads']
    mapped_read_number = STAR_Log_final_out['Uniquely mapped reads number'] + STAR_Log_final_out['Number of reads mapped to multiple loci']
    mapped_read_percent = round(100*mapped_read_number/input_read_number,2)
    logging.info(f"STAR Mapped Read Number: {mapped_read_number}.\n")
    logging.info(f"STAR alignment rate: {mapped_read_percent}%.\n")
    
    
    do_easy_mode = 1
    if "H" in options:
        do_easy_mode = 0
    tempfilteredbam = os.path.join(outputfolder, "STAR/tempfiltered.bam")
    tempfilteredsortedbam = os.path.join(outputfolder, "STAR/tempfiltered.sorted.bam")
    if os.path.isfile(tempfilteredsortedbam):
        tempfilteredbam = tempfilteredsortedbam
    collatedbams = os.path.join(outputfolder, "STAR/tempfiltered.collated")
    # collatedbed = os.path.join(outputfolder, "STAR/tempfiltered.bed")
    expmatbed = os.path.join(outputfolder, "expmat.bed")
    expmattsv = os.path.join(outputfolder, "expmat.tsv")

    #collate_bam(input_bam = tempfilteredbam, output_bam = collatedbam, threads = threadnum)
    #bam_to_bed(input_bam=collatedbam, output_bed=collatedbed)
    #inter_bed2geneFile(inputbed=collatedbed, outputfile=expmatbed, gtffile=gtffile)
    #featurebed2mattsv(input_file=expmatbed, output_file=expmattsv, barcodes_file=barcodes_file, gtf_file=gtffile, filter_str="25:0.75", easy_mode=do_easy_mode)
    #collate_and_split(tempfilteredbam)
    chunk_bams = collate_and_split(tempfilteredbam, collatedbams, threadnum)

    logging.info("overlapping feature with reads.\n")
    pool = multiprocessing.Pool(processes=int(threadnum))
    results = pool.starmap(worker, [(cbam, gtffile) for cbam in chunk_bams])
    pool.close()
    pool.join()
    with open(expmatbed, 'w') as fout:
        for cef in results:
            with open(cef, "r") as fin:
                for line in fin:
                    fout.write(line)
            os.remove(cef)
    
    logging.info("remove temp overlap files.\n")
    for cef in chunk_bams:
        os.remove(cef)

    
    with open(expmatbed, "r") as f:
        readwithfeat = sum(1 for _ in f)
    readwithoutfeat = mapped_read_number - readwithfeat
    unfeat_read_percent = round(100*(readwithoutfeat/mapped_read_number),2)
    feat_read_percent = round(100*(readwithfeat/mapped_read_number),2)
    logging.info(f"Total dedup reads overlapped with feautre: {readwithfeat}.\n")
    logging.info(f"Reads aligned to features (count & rate): {readwithfeat} [{feat_read_percent}%].\n")
    logging.info(f"Intergenic alignment (count & rate): {readwithoutfeat} [{unfeat_read_percent}%].\n")


    logging.info(f"feature bed to mat tsv with quality filter: {qualityfilter}.")
    logging.info(f"\n\n\nMatrix Stat of expmat.tsv:\n")
    if qualityfilter == 'NA' or qualityfilter == '0:0':
        featurebed2mattsv(
            input_file=expmatbed, 
            output_file=expmattsv, 
            barcodes_file=barcodes_file, 
            gtf_file=gtffile, 
            filter_str="0:0", 
            easy_mode=do_easy_mode
        )
    else:
        featurebed2mattsv(
            input_file=expmatbed, 
            output_file=expmattsv, 
            barcodes_file=barcodes_file, 
            gtf_file=gtffile, 
            filter_str=qualityfilter, 
            easy_mode=do_easy_mode
        )
    logging.info(f"countfeature step ends.\n")
    return gtffile
    
    
    
    
