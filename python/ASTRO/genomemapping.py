#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import multiprocessing as mp
import tempfile
import logging

def star_align(genome_dir, read_files_in, out_prefix, run_thread_n=16, gtf_file=None, extra_params=None):

    if extra_params is None:
        extra_params = []
    if isinstance(extra_params, str):
        extra_params = extra_params.split()

    gtf_part = []
    if gtf_file is not None and gtf_file != "NA":
        gtf_part = ["--sjdbGTFfile", gtf_file]

    cmd = [
        "STAR",
        "--genomeDir", genome_dir,
        "--readFilesIn", read_files_in,
        "--outFileNamePrefix", out_prefix,
        "--runThreadN", str(run_thread_n),
    ] + gtf_part + extra_params

    result4out = subprocess.run(cmd, text=True, capture_output=True, check=True)
    logging.info("STAR for genome mapping:\n%s", result4out.stdout)


def dedup_bam_samtools_markdup_2step(input_bam, output_bam, threads=1):
    markdup_bam = output_bam + ".markdup.bam"
    regex = r'.+\|:_:\|(.+)$'
    markdup_cmd = [
        "samtools", "markdup",
        "-r",
        "--barcode-rgx", regex,
        "-@", str(threads),
        "-O", "BAM",         
        input_bam,
        markdup_bam
    ]
    logging.info(f"[dedup_bam_2step] Step1 => {' '.join(markdup_cmd)}")
    ret = subprocess.run(markdup_cmd)
    if ret.returncode != 0:
        raise RuntimeError(f"[Error] samtools markdup failed on {input_bam}")

    temp_sam = output_bam + ".temp.sam"
    view_cmd_1 = ["samtools", "view", "-h", markdup_bam]
    logging.info(f"[dedup_bam_2step] Step2 => {' '.join(view_cmd_1)}")

    with open(temp_sam, "w") as out_sam:
        ret2 = subprocess.run(view_cmd_1, stdout=out_sam)
    if ret2.returncode != 0:
        raise RuntimeError("[Error] samtools view -h (reading markdup_bam) failed")

    final_sam = output_bam + ".final.sam"
    dedup_dict = {}
    rename_dict = {}
    global_index = 0

    with open(temp_sam, "r") as f_in, open(final_sam, "w") as f_out:
        for line in f_in:
            line = line.rstrip("\n")
            if line.startswith('@'):
                f_out.write(line + "\n")
                continue

            fields = line.split('\t')
            flag = int(fields[1])
            qname = fields[0]

            if (flag & 1024) != 0:
                splitted = qname.split("|:_:|", 1)
                old_readname = splitted[0] if splitted else qname
                dedup_dict[old_readname] = True
                continue

            if (flag & 4) != 0:
                continue

            splitted = qname.split("|:_:|", 1)
            if len(splitted) == 2:
                oldNamePart, posPart = splitted
            else:
                oldNamePart, posPart = ("","")

            if oldNamePart in dedup_dict:
                continue

            seq = fields[9]
            readlen = len(seq)
            ASv = -999999
            oq_idx = -1
            for i in range(11, len(fields)):
                if fields[i].startswith("AS:i:"):
                    ASv = int(fields[i][5:])
                elif fields[i].startswith("OQ:Z:"):
                    oq_idx = i

            if oq_idx != -1:
                fields.pop(oq_idx)

            if oldNamePart not in rename_dict:
                global_index += 1
                rename_dict[oldNamePart] = f"{global_index}_{posPart}"

            thename = rename_dict[oldNamePart]
            newQname = f"{thename}:{ASv}:{readlen}"
            fields[0] = newQname

            f_out.write("\t".join(fields) + "\n")

    view_cmd_2 = ["samtools", "view", "-b", "-o", output_bam, final_sam]
    logging.info(f"[dedup_bam_2step] Step3 => {' '.join(view_cmd_2)}")
    ret3 = subprocess.run(view_cmd_2)
    if ret3.returncode != 0:
        raise RuntimeError("[Error] samtools view final_sam => output_bam failed")

    for tmpf in [markdup_bam, temp_sam, final_sam]:
        if os.path.exists(tmpf):
            os.remove(tmpf)

    logging.info(f"[dedup_bam_2step] Done => {output_bam}")

def dedup_chunk(chunk_index, lines, temp_dir):    
    global_index = 0
    prev_qname   = ""
    records      = []
    best_AS      = -999999
    best_readlen = 0
    best_OQ      = ""

    def flush_group(output):
        nonlocal global_index, prev_qname, records, best_AS, best_readlen, best_OQ
        if not prev_qname:
            return

        global_index += 1
        new_prefix = f"{global_index}_{prev_qname}:{best_AS}:{best_readlen}"

        for rec_line in records:
            fds = rec_line.split('\t')
            oq_val = ""
            oq_tag_idx = -1
            for i in range(11, len(fds)):
                if fds[i].startswith("OQ:Z:"):
                    oq_val = fds[i][5:]
                    oq_tag_idx = i
                    break
            if oq_val == best_OQ:
                # rename QNAME
                fds[0] = new_prefix
                if oq_tag_idx != -1:
                    fds.pop(oq_tag_idx)
                output.write("\t".join(fds) + "\n")

        records.clear()
        best_AS      = -999999
        best_readlen = 0
        best_OQ      = ""
        prev_qname   = ""

    with tempfile.NamedTemporaryFile(mode='w', delete=False, dir=temp_dir, prefix=f"chunk_{chunk_index}_", suffix=".txt") as fw:
        for line in lines:
            line = line.rstrip("\n")
            if line.startswith('@'):
                fw.write(line + '\n')
                continue

            fds = line.split('\t')
            qname = fds[0]

            if qname != prev_qname:
                flush_group(fw)
                prev_qname = qname
                records = []

            records.append(line)

            # parse AS/readlen
            seq = fds[9]
            readlen = len(seq)
            ASv = -999999
            for tag in fds[11:]:
                if tag.startswith("AS:i:"):
                    ASv = int(tag[5:])

            if ASv > best_AS:
                best_AS      = ASv
                best_readlen = readlen
                OQv = ""
                for tag in fds[11:]:
                    if tag.startswith("OQ:Z:"):
                        OQv = tag[5:]
                        break
                best_OQ = OQv
        flush_group(fw)
    return fw.name


def parallel_dedup(cmd_collated_view, cmd_final_write, nproc=4, chunk_size=100000, temp_dir="/tmp"):

    pool = mp.Pool(processes=nproc)
    proc_in2 = subprocess.Popen(cmd_collated_view, stdout=subprocess.PIPE, text=True)

    async_results = []
    temp_files = []

    lines_buffer = []
    chunk_index = 0

    def flush_buffer_to_pool(buf, idx):
        return pool.apply_async(dedup_chunk, (idx, buf, temp_dir))

    prev_qname = ''
    for line in proc_in2.stdout:
        line = line.rstrip("\n")

        if len(lines_buffer) >= chunk_size:
            qname = line.split('\t')[0]
            if prev_qname != '' and prev_qname != qname:
                res = flush_buffer_to_pool(lines_buffer, chunk_index)
                async_results.append(res)
                lines_buffer = []
                chunk_index += 1
                prev_qname = ''
            elif prev_qname == '':
                prev_qname = qname
            lines_buffer.append(line)
        else:
            lines_buffer.append(line)

    if lines_buffer:
        res = flush_buffer_to_pool(lines_buffer, chunk_index)
        async_results.append(res)

    proc_in2.stdout.close()
    ret_view2 = proc_in2.wait()
    if ret_view2 != 0:
        raise RuntimeError("samtools view (collated_view) failed.")

    pool.close()
    for r in async_results:
        tf = r.get()
        temp_files.append(tf)
    pool.join()

    proc_out2 = subprocess.Popen(cmd_final_write, stdin=subprocess.PIPE, text=True)
    for tf in temp_files:
        with open(tf, 'r') as f:
            for line in f:
                proc_out2.stdin.write(line)
    proc_out2.stdin.close()
    ret_write2 = proc_out2.wait()
    if ret_write2 != 0:
        raise RuntimeError("samtools view (final_write) failed.")

    for tf in temp_files:
        try:
            os.remove(tf)
        except Exception as e:
            sys.stderr.write(f"Warning: Could not remove temp file {tf}: {e}\n")
    logging.info("[INFO] parallel_sam_view done.")


def dedup_bam_own(input_bam, output_bam, threads=1, temp_dir = '/tmp', chunk_size=100000):
    filtered_bam = f"{input_bam}.mappedOnly.bam"
    cmd_filter = [
        "samtools", "view",
        "-b", "-F", "4",
        input_bam,
        "-o", filtered_bam
    ]
    logging.info(f"[v2] Step0: Filtering unmapped => {' '.join(cmd_filter)}")
    ret = subprocess.run(cmd_filter)
    if ret.returncode != 0:
        raise RuntimeError(f"[Error] samtools view -F 4 failed on {input_bam}")

    temp_renamed_bam = f"{output_bam}.renamed.bam"
    cmd_view  = ["samtools", "view", "-h", filtered_bam]
    cmd_write = ["samtools", "view", "-b", '-@', str(int(threads)-1), "-o", temp_renamed_bam, "-"]

    logging.info("[v2] Step1: Renaming QNAME => bc+UMI, storing old in OQ:Z:")
    with subprocess.Popen(cmd_view, stdout=subprocess.PIPE, text=True) as proc_in, \
         subprocess.Popen(cmd_write, stdin=subprocess.PIPE, text=True) as proc_out:

        for line in proc_in.stdout:
            line = line.rstrip("\n")
            if line.startswith('@'):
                proc_out.stdin.write(line + "\n")
                continue

            fields = line.split('\t')
            old_qname = fields[0]

            if "|:_:|" in old_qname:
                splitted = old_qname.split("|:_:|", 1)
                oldName  = splitted[0]
                bcumi    = splitted[1]
                fields[0] = bcumi
                fields.append(f"OQ:Z:{oldName}")

            proc_out.stdin.write("\t".join(fields) + "\n")

        proc_out.stdin.close()
        ret_write = proc_out.wait()
    ret_view = proc_in.wait()
    if ret_view != 0 or ret_write != 0:
        raise RuntimeError("[v2] Step1 rename QNAME step failed")

    temp_collated_bam = f"{output_bam}.collated.bam"
    cmd_collate = [
        "samtools", "collate",
        "-@", str(threads),
        "-o", temp_collated_bam,
        temp_renamed_bam
    ]
    logging.info(f"[v2] Step2: Collate => {' '.join(cmd_collate)}")
    ret_collate = subprocess.run(cmd_collate)
    if ret_collate.returncode != 0:
        raise RuntimeError("[v2] samtools collate failed")

    logging.info("[v2] Step3: flush group => keep best AS => final rename => out.bam")

    cmd_collated_view = ["samtools", "view", "-h", temp_collated_bam]
    cmd_final_write   = ["samtools", "view", "-b", "-o", output_bam, "-"]
    
    parallel_dedup(cmd_collated_view, cmd_final_write, nproc=threads, chunk_size=chunk_size, temp_dir=temp_dir)

    os.remove(filtered_bam)
    os.remove(temp_renamed_bam)
    os.remove(temp_collated_bam)

    logging.info(f"[v2] Done => {output_bam}")

def collate_bam(input_bam, output_bam, threads=1):
    logging.info(f"[samtools_collate] Collating {input_bam} => {output_bam} (threads={threads})")
    result = subprocess.run([
        "samtools", "collate",
        "-@", str(threads),
        "-o", output_bam,
        input_bam
    ])
    if result.returncode != 0:
        raise RuntimeError("[samtools_collate] samtools collate failed.")
    logging.info("[samtools_collate] Done.")


def genomemapping(starref, gtffile, threadnum, options, outputfolder, STARparamfile='NA'):

    for h in logging.root.handlers[:]:
        logging.root.removeHandler(h)
        
    logfilename = os.path.join(outputfolder, ".logs/genomemapping.log")
    os.makedirs(os.path.dirname(logfilename), exist_ok=True)
    logging.basicConfig(filename=logfilename, filemode="w", level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


    logging.info("genomemapping step starts...\n")
    

    star_output_prefix = os.path.join(outputfolder, "STAR/temp")

    star_sorted_bam    = star_output_prefix + "Aligned.sortedByCoord.out.bam"
    star_unsorted_bam    = star_output_prefix + "Aligned.out.bam"
    

    if starref != "NA":
        if STARparamfile == "NA":
            if "M" in options:
                extra_star_params = [
                "--outSAMattributes", "NH", "HI", "AS", "nM", "NM",
                "--genomeLoad", "NoSharedMemory",
                "--limitOutSAMoneReadBytes", "200000000",
                "--outFilterMultimapNmax", "-1",
                "--outFilterMultimapScoreRange", "0",
                "--readMatesLengthsIn", "NotEqual",
                "--limitBAMsortRAM", "0",
                "--outMultimapperOrder", "Random",
                "--outSAMtype", "BAM", "SortedByCoordinate",
                "--outSAMunmapped", "Within",
                "--outSAMorder", "Paired",
                "--outSAMprimaryFlag", "AllBestScore",
                "--outSAMmultNmax", "-1",
                "--outFilterType", "Normal",
                "--outFilterScoreMinOverLread", "0",
                "--alignSJDBoverhangMin", "30",
                "--outFilterMatchNmin", "15",
                "--outFilterMatchNminOverLread", "0",
                "--outFilterMismatchNoverLmax", "0.1",
                "--outFilterMismatchNoverReadLmax", "0.15",
                "--alignIntronMin", "20",
                "--alignIntronMax", "1000000",
                "--alignEndsType", "Local",
                "--outBAMsortingBinsN", "200"
                ]
            else:
                extra_star_params = [
                "--outSAMattributes", "NH", "HI", "AS", "nM", "NM",
                "--genomeLoad", "NoSharedMemory",
                "--limitOutSAMoneReadBytes", "200000000",
                "--outFilterMultimapNmax", "-1",
                "--outFilterMultimapScoreRange", "0",
                "--readMatesLengthsIn", "NotEqual",
                "--outMultimapperOrder", "Random",
                "--outSAMtype", "BAM", "Unsorted",
                "--outSAMunmapped", "Within",
                "--outSAMorder", "Paired",
                "--outSAMprimaryFlag", "AllBestScore",
                "--outSAMmultNmax", "-1",
                "--outFilterType", "Normal",
                "--outFilterScoreMinOverLread", "0",
                "--alignSJDBoverhangMin", "30",
                "--outFilterMatchNmin", "15",
                "--outFilterMatchNminOverLread", "0",
                "--outFilterMismatchNoverLmax", "0.1",
                "--outFilterMismatchNoverReadLmax", "0.15",
                "--alignIntronMin", "20",
                "--alignIntronMax", "1000000",
                "--alignEndsType", "Local"]
        else:
            with open(STARparamfile, 'r') as sf:
                lines = [line.strip() for line in sf if line.strip()]
            extra_star_params = []
            for line in lines:
                extra_star_params.extend(line.split())

        logging.info("STAR genome mapping log:\n")
        star_align(
            genome_dir=starref,
            read_files_in=os.path.join(outputfolder, "combine.fq"), 
            out_prefix=star_output_prefix,
            run_thread_n=int(threadnum),
            gtf_file=gtffile,
            extra_params=extra_star_params
        )

    filtered_bam = os.path.join(outputfolder, "STAR/tempfiltered.bam")

    

    if "M" in options:
        logging.info("genomemapping Using two-step markdup => dedup_bam_samtools_markdup_2step()")
        dedup_bam_samtools_markdup_2step(
            input_bam = star_sorted_bam,
            output_bam= filtered_bam,
            threads   = int(threadnum),
        )
    else:
        logging.info("genomemapping Using alignment performance => dedup_bam_own()")
        temp_dir = os.path.join(outputfolder, "temps/")
        dedup_bam_own(
            input_bam = star_unsorted_bam,
            output_bam= filtered_bam,
            threads   = int(threadnum),
            temp_dir  = temp_dir,
        )

    logging.info(f"[genomemapping] Done. final => {filtered_bam}")
    logging.info(f"genomemapping step ends\n")
    
    return filtered_bam
