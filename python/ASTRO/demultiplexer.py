#!/usr/bin/env python3
import collections
import re
import os
import gzip
import subprocess
import multiprocessing as mp
import sys
import tempfile
from collections import defaultdict, Counter
import logging


def singleCutadapt(barcodestr,outputfile,remainfa,threadnum):


    def simple_fastq_iterator(handle):
        while True:
            title_line = handle.readline().rstrip()
            seq_string = handle.readline().rstrip()
            sep_line = handle.readline().rstrip()
            qual_string = handle.readline().rstrip()
            if not (title_line or seq_string or sep_line or qual_string):
                break
            yield title_line[1:], seq_string, qual_string
    
    def cutfastq(filein,forward,length):
        with (gzip.open(filein, "rt") if filein.endswith('.gz') else open(filein, "rt")) as in_handle: 
            with open(filein+'.temp', "w") as out_handle:
                for title, seq, qual in simple_fastq_iterator(in_handle):
                    if forward:
                        new_seq = seq[0:length]
                        new_qual = qual[0:length]
                        out_handle.write("@%s\n%s\n+\n%s\n" % (title, new_seq, new_qual))
                    else:
                        new_seq = seq[-length:len(seq)]
                        new_qual = qual[-length:len(qual)]
                        out_handle.write("@%s\n%s\n+\n%s\n" % (title, new_seq, new_qual))
        os.remove(filein)
        os.rename(filein+'.temp', filein)
    def cutfastq_abs(srcfile, outfile, startpos, length):
        with (gzip.open(srcfile, "rt") if srcfile.endswith('.gz') else open(srcfile, "rt")) as in_handle:
            with open(outfile, "w") as out_handle:
                for title, seq, qual in simple_fastq_iterator(in_handle):
                    s = int(startpos); k = int(length); e = s + k
                    if s >= len(seq):
                        new_seq  = ""
                        new_qual = ""
                    else:
                        new_seq  = seq[s:e]
                        new_qual = qual[s:e]
                    out_handle.write(f"@{title}\n{new_seq}\n+\n{new_qual}\n")
    def combineFqs(outputfile,fqs):
        iterators0 = [gzip.open(file_path, "rt") if file_path.endswith('.gz') else open(file_path) for file_path in fqs]
        iterators = [simple_fastq_iterator(file) for file in iterators0]
        with open(outputfile, 'w') as out_handle:
            for entries in zip(*iterators):
                newseq = ""
                newqual = ""
                newtitle = ''
                for title, seq, qual in entries:
                    newseq += seq
                    newqual += qual
                    if newtitle != '' and newtitle != title:
                        print(fqs)
                        print(newtitle)
                        print(title)
                        exit('Error: fastqs have different read name')
                out_handle.write("@%s\n%s\n+\n%s\n" % (title, newseq, newqual))




    customsetting = '-j ' + threadnum
    array4barcodes = re.split(':',barcodestr)

    ii = 0
    linkfqs = []
    for stri in array4barcodes:
        array4input = re.split('_', stri)
        tempout = outputfile + '.temp' + str(ii)
        if len(array4input) == 2:
            if re.fullmatch(r'\d+', array4input[0]) and re.fullmatch(r'\d+', array4input[1]):
                cutfastq_abs(remainfa, tempout, int(array4input[0])-1, int(array4input[1]))
                linkfqs.append(tempout)

            # original behaviors kept: N_ADAPTER or ADAPTER_N
            elif re.fullmatch(r'\d+', array4input[0]):
                outputcommand = ['cutadapt'] + [customsetting] + ['-a'] + [array4input[1]] + ['-e'] + ['0.25'] + ['-o'] + [tempout] + [remainfa]
                result4out = subprocess.run(' '.join(outputcommand), shell=True, text=True, capture_output=True, check=True)
                logging.info("\n%s", result4out.stdout)
                cutfastq(tempout, False, int(array4input[0]))
                linkfqs.append(tempout)

            elif re.fullmatch(r'\d+', array4input[1]):
                outputcommand = ['cutadapt'] + [customsetting] + ['-g'] + [array4input[0]] + ['-e'] + ['0.25'] + ['-o'] + [tempout] + [remainfa]
                result4out = subprocess.run(' '.join(outputcommand), shell=True, text=True, capture_output=True, check=True)
                logging.info("\n%s", result4out.stdout)
                cutfastq(tempout, True, int(array4input[1]))
                linkfqs.append(tempout)

            else:
                print(stri)
                exit('Error: wrong barcode format: order error')

        elif len(array4input) == 1:
            outputcommand = ['cutadapt'] + [customsetting] + ['-g'] + [array4input[0]] + ['-e'] + ['0.25'] + ['-o'] + [tempout] + [remainfa]
            result4out = subprocess.run(' '.join(outputcommand), shell=True, text=True, capture_output=True, check=True)
            logging.info("\n%s", result4out.stdout)
            linkfqs.append(tempout)
        else:
            print(stri)
            exit('Error: wrong barcode format: order error')
        ii = ii + 1

    if len(linkfqs) >= 2:
        combineFqs(outputfile,linkfqs)
        for fqi in linkfqs:
            os.remove(fqi)
    else:
        if len(linkfqs) == 1:
            os.rename(linkfqs[0], outputfile)




def merge_chunk_to_tempfile(chunk_index, linesA, linesB, linesC, temp_dir):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, dir=temp_dir, prefix=f"chunk_{chunk_index}_", suffix=".fastq") as tf:
            ii = 0
            for i in range(0, len(linesA), 4):
                a1 = linesA[i].rstrip('\n')   # @readA
                a2 = linesA[i+1].rstrip('\n') # SEQ
                a3 = linesA[i+2].rstrip('\n') # +
                a4 = linesA[i+3].rstrip('\n') # QUAL


                b1 = linesB[i].rstrip('\n')
                b2 = linesB[i+1].rstrip('\n')
                b3 = linesB[i+2].rstrip('\n')
                b4 = linesB[i+3].rstrip('\n')


                c1 = linesC[i].rstrip('\n')
                c2 = linesC[i+1].rstrip('\n')
                c3 = linesC[i+2].rstrip('\n')
                c4 = linesC[i+3].rstrip('\n')

                a1 = a1.split(' ')[0]
                b1 = b1.split(' ')[0]
                c1 = c1.split(' ')[0]
                if a1 != b1 or a1 != c1:
                    raise ValueError(f"Read names do not match: {a1} {b1} {c1}")

                # filter out empty reads
                if len(a2) == 0:
                    continue
                
                new_read_name = f"@{chunk_index}_{ii}---{b2}---{c2}---{b4}{c4}"
                ii += 1

                new_line1 = new_read_name
                new_line2 = a2
                new_line3 = a3
                new_line4 = a4

                tf.write(f"{new_line1}\n{new_line2}\n{new_line3}\n{new_line4}\n")

        return tf.name

def parse4barcodeposition(barcodeposition):
    startpos, endpos = 0, 0

    if m := re.fullmatch(r'([0-9_]+)b', barcodeposition):
        parts = m.group(1).split('_')
        startpos = int(parts[0])-1 if len(parts) > 1 else 0
        endpos = int(parts[-1])+startpos
    elif m := re.fullmatch(r'b([0-9_]+)', barcodeposition):
        parts = m.group(1).split('_')
        startpos = -int(parts[0])
        endpos = int(parts[-1])+startpos if len(parts) > 1 else None
        if endpos is not None and endpos >= 0:
            endpos = None
    else:
        raise ValueError(f"wrong format in barcodeposition: {barcodeposition}")
    return startpos, endpos
def merge_chunk_to_tempfile_to2files(chunk_index, linesA, linesB, linesC, temp_dir, barcodepositionvector, barcodelengthmin, barcodelengthmax, barcode_dict):
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, dir=temp_dir, prefix=f"chunk_{chunk_index}_", suffix=".fastq") as tf, tempfile.NamedTemporaryFile(mode='w', delete=False, dir=temp_dir, prefix=f"chunk_{chunk_index}_", suffix=".fastq") as tf2:
            ii = 0
            for i in range(0, len(linesA), 4):
                a1 = linesA[i].rstrip('\n')   # @readA
                a2 = linesA[i+1].rstrip('\n') # SEQ
                if not (int(barcodelengthmin) <= len(a2) <= int(barcodelengthmax)):
                    continue
                shouldbebarcode = a2[barcodepositionvector[0]:barcodepositionvector[1]]
                
                if shouldbebarcode in barcode_dict:
                    spatialcode = barcode_dict[shouldbebarcode]
                    b1, b2, b3, b4 = [line.rstrip('\n') for line in linesB[i:i+4]]
                    if len(b2) < 10:
                        continue
                    c1, c2, c3, c4 = [line.rstrip('\n') for line in linesC[i:i+4]]
                    a1 = a1.split(' ')[0]; b1 = b1.split(' ')[0]; c1 = c1.split(' ')[0]
                    if a1 != b1 or a1 != c1:
                        raise ValueError(f"Read names do not match: {a1} {b1} {c1}")
                    read_name = '@' + f"{chunk_index}_{ii}" + '|:_:|' + spatialcode + ':' + c2
                    
                    new_line1, new_line2, new_line3, new_line4 = read_name, b2, b3, b4
                    tf2.write(f"{new_line1}\n{new_line2}\n{new_line3}\n{new_line4}\n")
                    
                    ii += 1
                    continue

                a3 = linesA[i+2].rstrip('\n') # +
                a4 = linesA[i+3].rstrip('\n') # QUAL


                b1 = linesB[i].rstrip('\n')
                b2 = linesB[i+1].rstrip('\n')
                b3 = linesB[i+2].rstrip('\n')
                b4 = linesB[i+3].rstrip('\n')


                c1 = linesC[i].rstrip('\n')
                c2 = linesC[i+1].rstrip('\n')
                c3 = linesC[i+2].rstrip('\n')
                c4 = linesC[i+3].rstrip('\n')

                a1 = a1.split(' ')[0]
                b1 = b1.split(' ')[0]
                c1 = c1.split(' ')[0]
                if a1 != b1 or a1 != c1:
                    raise ValueError(f"Read names do not match: {a1} {b1} {c1}")

                # filter out empty reads
                if len(a2) == 0:
                    continue
                
                new_read_name = f"@{chunk_index}_{ii}---{b2}---{c2}---{b4}{c4}"
                ii += 1

                new_line1 = new_read_name
                new_line2 = a2
                new_line3 = a3
                new_line4 = a4

                tf.write(f"{new_line1}\n{new_line2}\n{new_line3}\n{new_line4}\n")

        return tf.name, tf2.name
def read_in_chunks(fA, fB, fC, chunk_num_reads):
    chunk_index = 0
    while True:
        linesA = []
        linesB = []
        linesC = []
        for _ in range(chunk_num_reads * 4):
            lineA = fA.readline()
            lineB = fB.readline()
            lineC = fC.readline()
            if not lineA or not lineB or not lineC:
                break
            linesA.append(lineA)
            linesB.append(lineB)
            linesC.append(lineC)
        if not linesA:
            break
        yield (chunk_index, linesA, linesB, linesC)
        chunk_index += 1
def Fqs2_1fq(fa, fb, fc, out, nproc, chunk_size, temp_dir, barcodeposition='NA',barcodelengthrange='NA',barcode_dict={}):

        pool = mp.Pool(processes=nproc)
        if barcodeposition == 'NA':
            temp_files = []   
            with open(fa, 'r') as fA, open(fb, 'r') as fB, open(fc, 'r') as fC:
                async_results = []
                for (chunk_index, linesA, linesB, linesC) in read_in_chunks(fA, fB, fC, chunk_num_reads=chunk_size):
                    res = pool.apply_async(merge_chunk_to_tempfile, (chunk_index, linesA, linesB, linesC, temp_dir))
                    async_results.append(res)
            pool.close()

            for r in async_results:
                tempfile_path = r.get()
                temp_files.append(tempfile_path)

            pool.join()
        else:
            if barcodelengthrange and barcodelengthrange != 'NA':
                barcodelengthmin, barcodelengthmax = barcodelengthrange.split('_')[0:2]
            else:
                barcodelengthmin = 1
                barcodelengthmax = 10000
            temp_files = []
            temp_files_final = []
            barcodeposition = parse4barcodeposition(barcodeposition)
            with open(fa, 'r') as fA, open(fb, 'r') as fB, open(fc, 'r') as fC:
                async_results = []
                for (chunk_index, linesA, linesB, linesC) in read_in_chunks(fA, fB, fC, chunk_num_reads=chunk_size):
                    res = pool.apply_async(merge_chunk_to_tempfile_to2files, (chunk_index, linesA, linesB, linesC, temp_dir, barcodeposition, barcodelengthmin, barcodelengthmax, barcode_dict))
                    async_results.append(res)
            pool.close()

            for r in async_results:
                tempfile_path1, tempfile_path2 = r.get()
                temp_files.append(tempfile_path1)
                temp_files_final.append(tempfile_path2)

            pool.join()

            with open(out+'0', 'w') as fw:
               for tf in temp_files_final:
                    with open(tf, 'r') as f:
                        for line in f:
                            fw.write(line)
        
            for tf in temp_files_final:
                try:
                    os.remove(tf)
                except Exception as e:
                    sys.stderr.write(f"Warning: Could not remove temp file {tf}: {e}\n")
        with open(out, 'w') as fw:
            for tf in temp_files:
                with open(tf, 'r') as f:
                    for line in f:
                        fw.write(line)

        for tf in temp_files:
            try:
                os.remove(tf)
            except Exception as e:
                sys.stderr.write(f"Warning: Could not remove temp file {tf}: {e}\n")
        
        




def filter_sam_NH(input_sam, output_fq):
    with open(input_sam, 'r') as infile, open(output_fq, 'w') as outfile:
        for line in infile:
            if not line.startswith('@'):
                line = line.split('\t')
                spatialcode = line[2]
                if spatialcode == '*':
                    continue
                if (int(line[1]) & 16) != 0: 
                    continue
                for tag in line[11:]:
                        if tag.startswith('NH:i:'):
                            nh_value = int(tag.split(':', 2)[2])
                        break
                if nh_value >= 2:
                    continue
                arrays = line[0].split('---', maxsplit=3)
                read_name = arrays[0]
                read1seq  = arrays[1]
                lenread1 = len(read1seq)
                if lenread1 <  10:
                    continue
                UMIseq = arrays[2]
                read1qual = arrays[3][:lenread1]
                read_name2 = '@' + read_name + '|:_:|' + spatialcode + ':' + UMIseq
                line = read_name2 + '\n' + read1seq + '\n' + '+' + '\n' + read1qual + '\n'
                outfile.write(line)

def filter_sam_nbhd(input_sam, output_fq):
    previousname = ''
    previousloca = ''
    wrongmultiple = 0
    rdyline = ''
    with open(input_sam, 'r') as infile, open(output_fq, 'w') as outfile:
        for line in infile:
            if not line.startswith('@'):
                line = line.split('\t')
                spatialcode = line[2]
                arrays = line[0].split('---', maxsplit=3)
                read_name = arrays[0]
                read2seq  = arrays[1]
                lenread1 = len(read2seq)
                if spatialcode == '*'  or (int(line[1]) & 16) != 0 or lenread1 <  10:
                    thisline = ''
                else:
                    UMIseq = arrays[2]
                    read2qual = arrays[3][:lenread1]
                    read_name2 = '@' + read_name + '|:_:|' + spatialcode + ':' + UMIseq
                    thisline = read_name2 + '\n' + read2seq + '\n' + '+' + '\n' + read2qual + '\n'
                if previousname != read_name:
                    if wrongmultiple == 0:
                        outfile.write(rdyline)
                    rdyline = thisline
                    wrongmultiple = 0
                    previousname = read_name
                    previousloca = spatialcode
                else:
                    #if previousloca != spatialcode:
                        wrongmultiple = 1
        if wrongmultiple == 0:
            outfile.write(rdyline)
def process_chunk(chunk_lines):
    read_count = set()
    read_delete = set()  
    line = chunk_lines[0] 
    read_name0 = line
    for line in chunk_lines[1:]:
        read_name = line
        if read_name == read_name0:
            read_delete.add(read_name)
        else:
            read_count.add(read_name0)
            read_name0 = read_name
    if read_name0 not in read_delete:
        read_count.add(read_name0)
    return read_count, read_delete
def total_delete(all_stats):
    total_read_count = set()
    total_read_delete = set()
    for read_count, read_delete in all_stats:
        for read_name in read_count:
            if read_name in total_read_count:
                total_read_delete.add(read_name)
            else:
                total_read_count.add(read_name)
        for read_name in read_delete:
            total_read_delete.add(read_name)
    return total_read_delete
def read_in_chunks_from_sam(filename, chunk_size):
    chunk = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('@'):
                continue
            lines = line.split('\t')
            if lines[2] == '*':
                continue
            #if (int(lines[1]) & 16) != 0: 
            #    continue
            read_name = line.split('\t')[0].split('---')[0]
            chunk.append(read_name)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        if chunk:
            yield chunk
def filter_sam0(input_sam, reads_to_remove, output_fq):
    with open(input_sam, 'r') as infile, open(output_fq, 'w') as outfile:
        for line in infile:
            if not line.startswith('@'):
                lines = line.split('\t')
                spatialcode = lines[2]
                if spatialcode == '*':
                    continue
                if (int(lines[1]) & 16) != 0: 
                    continue
                arrays = lines[0].split('---', maxsplit=3)
                read_name = arrays[0]
                if read_name not in reads_to_remove:
                    read1seq  = arrays[1]
                    lenread1 = len(read1seq)
                    if lenread1 <  10:
                        continue
                    UMIseq = arrays[2]
                    read1qual = arrays[3][:lenread1]
                    read_name2 = '@' + read_name + '|:_:|' + spatialcode + ':' + UMIseq
                    line = read_name2 + '\n' + read1seq + '\n' + '+' + '\n' + read1qual + '\n'
                    outfile.write(line)
def filter_sam(input_sam, output_fq, num_processes, chunk_size):
    pool = mp.Pool(processes=int(num_processes))
    chunks = read_in_chunks_from_sam(input_sam, chunk_size=chunk_size)
    results = pool.map(process_chunk, chunks)
    pool.close()
    pool.join()

    reads_to_remove = total_delete(results)
    filter_sam0(input_sam, reads_to_remove, output_fq)
    



    #def filter_sam_chunk(chunk_lines, chunk_index, reads_to_remove, temp_dir):
    #    with tempfile.NamedTemporaryFile(mode='w', delete=False, dir=temp_dir, prefix=f"chunk_{chunk_index}_", suffix=".fastq") as tf:
    #        for line in chunk_lines:
    #            if not line.startswith('@'):
    #                spatialcode = line.split('\t')[2]
    #                arrays = line.split('\t')[0].split('---', maxsplit=3)
    #                read_name = arrays[0]
    #                if read_name not in reads_to_remove:
    #                    read1seq  = arrays[1]
    #                    UMIseq = arrays[2]
    #                    read1qual = arrays[3][:len(read1seq)]
    #                    read_name2 = '@' + read_name + '|:_:|' + spatialcode + ':' + UMIseq
    #                    line = read_name2 + '\n' + read1seq + '\n' + '+' + '\n' + read1qual + '\n'
    #                    tf.write(line)
    #def worker_init(temp_dir):
    #    global TEMP_DIR
    #    TEMP_DIR = temp_dir 
    
    #    temp_dir = 'temp_dir/'
    #    if not temp_dir:
    #        temp_dir = tempfile.gettempdir()
    #    pool = mp.Pool(processes=num_processes, initializer=worker_init, initargs=(temp_dir,))
    #    async_results = []
    #    for (chunk_index, chunk_lines) in read_in_chunks_from_sam2(input_sam, chunk_size=chunk_size):
    #        res = pool.apply_async(filter_sam_chunk, (chunk_lines, chunk_index, reads_to_remove, temp_dir) )
    #    async_results.append(res)
    #    temp_files = []    
    #    pool.close()
    #    for r in async_results:
    #            tempfile_path = r.get()
    #            temp_files.append(tempfile_path)
    #    pool.join()
    #
    #    with open(output_fq, 'w') as fw:
    #        for tf in temp_files:
    #            with open(tf, 'r') as f:
    #                for line in f:
    #                    fw.write(line)
    #
    #    for tf in temp_files:
    #        try:
    #            os.remove(tf)
    #        except Exception as e:
    #            sys.stderr.write(f"Warning: Could not remove temp file {tf}: {e}\n")
     #def read_in_chunks_from_sam2(filename, chunk_size):
    #    chunk_index = 0
    #    chunk = []
    #    with open(filename, 'r') as file:
    #        for line in file:
    #            if line.startswith('@'):
    #                continue
    #            if (int(line.split('\t')[1]) & 16) != 0: 
    #                continue
    #            read_name = line.split('\t')[0].split('---')[0]
    #            chunk.append(read_name)
    #            if len(chunk) >= chunk_size:
    #                yield (chunk_index, chunk)
    #                chunk_index += 1
    #                chunk = []
    #        if chunk:
    #            yield (chunk_index, chunk)

def demultiplexing(read1,read2, barcode_file, PrimerStructure, StructureUMI, StructureBarcode, threadnum, outputfolder, limitOutSAMoneReadBytes4barcodeMapping, barcodeposition='NA', barcodelengthrange='NA'):
    
    os.makedirs(os.path.join(outputfolder, 'temps'), exist_ok=True)
    CleanFq1 = os.path.join(outputfolder, "temps/CleanFq1.fq")
    CleanFq2 = os.path.join(outputfolder, "temps/CleanFq2.fq")
    CombineFq = os.path.join(outputfolder, "combine.fq")
    barcode_db_fa = os.path.join(outputfolder, "temps/barcode_xy.fasta")
    barcode_db_path = os.path.join(outputfolder, "temps/barcode_db")
    index_fq = os.path.join(outputfolder, "temps/index.fastq")
    UMI_fq = os.path.join(outputfolder, "temps/UMI.fastq")
    temps_path = os.path.join(outputfolder, "temps/")

    for h in logging.root.handlers[:]:
        logging.root.removeHandler(h)

    logfilename = os.path.join(outputfolder, ".logs/demultiplexing.log")
    os.makedirs(os.path.dirname(logfilename), exist_ok=True)
    logging.basicConfig(filename=logfilename, filemode="w", level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    
    logging.info(f"demultiplexing step starts\n")

    threadnum = str(threadnum)
    if PrimerStructure != "NA":
        prefixread1 = PrimerStructure.split('_', 1)[0]
        suffixread1 = PrimerStructure.rsplit('_', 1)[-1]
        result4out = subprocess.run([
            "cutadapt", "-e", "0.25",
            "-a", suffixread1, "--times", "4",
            "-g", prefixread1,
            "-j", threadnum,
            "-o", CleanFq2, "-p", CleanFq1,
            read2, read1], text=True, capture_output=True, check=True)
        logging.info("cutadapt for mRNA read log:\n%s", result4out.stdout)
    else:
        def compress_one(src, dst):
            if not os.path.exists(str(src)):
                return
            tmp = dst + ".tmp"
            if which("pigz"):
                with open(tmp, "wb") as out:
                    subprocess.run(["pigz", "-p", str(threads), "-c", str(src)], check=True, stdout=out)
            elif which("bgzip"):
                with open(tmp, "wb") as out:
                    subprocess.run(["bgzip", "-@", str(threads), "-c", str(src)], check=True, stdout=out)
            else:
                with open(src, "rb") as f_in, gzip.open(tmp, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.replace(tmp, dst)   
            os.remove(src)
        def is_gzip(path):
            with open(path, "rb") as f:
                return f.read(2) == b"\x1f\x8b"
        if is_gzip(CleanFq1):
            compress_one(read1, CleanFq1)
        else:
            CleanFq1 = read1
        if is_gzip(CleanFq2):
            compress_one(read2, CleanFq2)
        else:
            CleanFq2 = read2
    
    logging.info("cutadapt for UMI log:\n")
    singleCutadapt(StructureUMI,UMI_fq,CleanFq1,threadnum)
    logging.info("cutadapt for barcode log:\n")
    singleCutadapt(StructureBarcode,index_fq,CleanFq1,threadnum)
    
    barcode_dict = {}
    with open(barcode_file, 'r') as barcodes_in, open(barcode_db_fa, 'w') as barcode_db_file:
      for line in barcodes_in:
          fields = line.strip().split('\t')
          header = f"{fields[1]}_{fields[2]}"
          sequence = fields[0]
          barcode_db_file.write(f">{header}\n{sequence}\n")
          barcode_dict[sequence] = header
    
    Fqs2_1fq(index_fq,CleanFq2,UMI_fq,CombineFq,16,1000000,temps_path,barcodeposition,barcodelengthrange,barcode_dict)
    result4out = subprocess.run([ "STAR", "--runMode", "genomeGenerate", "--runThreadN", threadnum, "--genomeDir", barcode_db_path, "--genomeFastaFiles", barcode_db_fa, "--genomeSAindexNbases", "7" ,"--limitGenomeGenerateRAM", "60000000000"], text=True, capture_output=True, check=True)
    logging.info("STAR genomeGenerate for barcode reference:\n%s", result4out.stdout)

    if limitOutSAMoneReadBytes4barcodeMapping != 'NA':
        limitOutSAMoneReadBytes4barcodeMapping = ['--limitOutSAMoneReadBytes', str(limitOutSAMoneReadBytes4barcodeMapping)]
    else:
        limitOutSAMoneReadBytes4barcodeMapping=[]

    result4out = subprocess.run([
    "STAR",
    "--runThreadN", threadnum,
    "--genomeDir", barcode_db_path,
    "--readFilesIn", CombineFq,
    "--readNameSeparator", "space",
    "--outFileNamePrefix", os.path.join(outputfolder, "temps/barcodeMapping/temp"),
    "--outSAMtype", "SAM",
    "--outFilterMismatchNmax", "3",
    "--outFilterMatchNmin", "13",
    "--alignEndsType", "Local",
    "--scoreGapNoncan", "-1000",
    "--scoreGapATAC", "-1000",
    "--alignIntronMax", "1",
    "--outFilterMultimapNmax", "-1",
    "--outSAMunmapped", "Within",
    "--outFilterMultimapScoreRange", "0",
    "--seedSearchStartLmax", "12",
    "--outFilterScoreMinOverLread", "0",
    "--outFilterMatchNminOverLread", "0",
    "--outFilterMismatchNoverLmax", "0.7",
    "--outFilterMismatchNoverReadLmax", "0.7"
    ]+limitOutSAMoneReadBytes4barcodeMapping, text=True, capture_output=True, check=True)
    logging.info("STAR genomeGenerate for barcode mapping:\n%s", result4out.stdout)
    filter_sam_nbhd(os.path.join(outputfolder, "temps/barcodeMapping/tempAligned.out.sam"), CombineFq)
    #filter_sam(os.path.join(outputfolder, "temps/barcodeMapping/tempAligned.out.sam"), CombineFq, int(threadnum), 1000000)

    if os.path.isfile(CombineFq + '0'):
        import shutil
        tmp_file = CombineFq + ".tmp"
        with open(tmp_file, 'wb') as out, open(CombineFq, 'rb') as f1, open(CombineFq + '0', 'rb') as f2:
            shutil.copyfileobj(f1, out)
            shutil.copyfileobj(f2, out)
        os.replace(tmp_file, CombineFq)
        os.remove(CombineFq + '0') 
    
    with open(os.path.join(outputfolder, "temps/barcodeMapping/tempLog.final.out"), "r") as f:
        content4log = f.read()
    logging.info("STAR Log.final.out for barcode mapping:\n%s", content4log)

    barcode_mapping_dir = os.path.join(outputfolder, "temps/barcodeMapping/")
    for file_name in os.listdir(barcode_mapping_dir):
        file_path = os.path.join(barcode_mapping_dir, file_name)
        if os.path.isfile(file_path):
            os.remove(file_path)  
    os.rmdir(barcode_mapping_dir)
    
    logging.info(f"demultiplexing step ends\n")


def get_barcode_for_single_cell(read1, read2, barcode_file, PrimerStructure, StructureUMI, StructureBarcode, threadnum, outputfolder,barcode_threshold=100,barcodelength=0):
    barcode_threshold = int(barcode_threshold)
    os.makedirs(os.path.join(outputfolder, 'temps'), exist_ok=True)
    CleanFq1 = os.path.join(outputfolder, "temps/cleanfq1.fq")
    CleanFq2 = os.path.join(outputfolder, "temps/cleanfq2.fq")
    CombineFq = os.path.join(outputfolder, "combine.fq")
    barcode_db_fa = os.path.join(outputfolder, "temps/barcode_xy.fasta")
    barcode_db_path = os.path.join(outputfolder, "temps/barcode_db")
    index_fq = os.path.join(outputfolder, "temps/index.fastq")
    UMI_fq = os.path.join(outputfolder, "temps/UMI.fastq")
    temps_path = os.path.join(outputfolder, "temps/") 
    auto_bc_path = os.path.join(outputfolder, "temps", "auto_barcode.tsv")
    
    prefixread1 = PrimerStructure.split('_', 1)[0]
    suffixread1 = PrimerStructure.rsplit('_', 1)[-1]
    
    threadnum = str(threadnum)
    subprocess.run([ "cutadapt", "-e", "0.25", "-a", suffixread1, "--times", "4", "-g", prefixread1, "-j", threadnum, "-o", CleanFq2, "-p", CleanFq1, read2, read1])
    singleCutadapt(StructureUMI,UMI_fq,CleanFq1,threadnum)
    singleCutadapt(StructureBarcode,index_fq,CleanFq1,threadnum)

    bar_counter = Counter()

    def fastq_iter(fq):
        while True:
            l1 = fq.readline()
            if not l1:
                break
            seq_line  = fq.readline()
            plus_line = fq.readline()
            qual_line = fq.readline()
            if not qual_line:
                break
            yield seq_line.strip()

    with open(index_fq, 'r') as fin:
        for seq in fastq_iter(fin):
            bar_counter[seq] += 1

    if barcode_file != "notavailable":
        whitelist = set()
        with open(barcode_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line: 
                    continue
                arr = line.split('\t')
                user_bc = arr[0]
                whitelist.add(user_bc)

        filtered = []
        for bc, ct in bar_counter.most_common():
            if bc not in whitelist:
                continue
            if ct < barcode_threshold:
                break
            if barcodelength>0 and len(bc)!=barcodelength:
                continue
            filtered.append(bc)
    else:
        filtered = []
        for bc, ct in bar_counter.most_common():
            if ct < barcode_threshold:
                break
            if barcodelength>0 and len(bc)!=barcodelength:
                continue
            filtered.append(bc)
    
    with open(auto_bc_path, 'w') as bf:
        i = 1
        for bc in filtered:
            bf.write(f"{bc}\t{i}\t{i}\n")
            i += 1

    return auto_bc_path
