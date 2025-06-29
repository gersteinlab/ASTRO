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
        array4input = re.split('_',stri)
        tempout = outputfile + '.temp' + str(ii)
        if len(array4input) == 2:
            if re.search(r'^[0-9]*$', array4input[0]):
                outputcommand = ['cutadapt'] + [customsetting] + ['-a'] + [array4input[1]] + ['-e'] + ['0.25'] + ['-o'] + [tempout] + [remainfa]
                subprocess.run(' '.join(outputcommand), shell=True)
                cutfastq(tempout,False,int(array4input[0]))
                linkfqs.append(tempout)
            elif re.search(r'^[0-9]*$', array4input[1]):
                outputcommand = ['cutadapt'] + [customsetting] + ['-g'] + [array4input[0]] + ['-e'] + ['0.25'] + ['-o'] + [tempout] + [remainfa]
                subprocess.run(' '.join(outputcommand), shell=True)
                cutfastq(tempout,True,int(array4input[1]))
                linkfqs.append(tempout)
            else:
                print(stri)
                exit('Error: wrong barcode format: order error')
        elif len(array4input) == 1:
            outputcommand = ['cutadapt'] + [customsetting] + ['-g'] + [array4input[0]] + ['-e'] + ['0.25'] + ['-o'] + [tempout] + [remainfa]
            subprocess.run(' '.join(outputcommand), shell=True)
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


                new_read_name = f"@{chunk_index}_{ii}---{b2}---{c2}---{b4}{c4}"
                ii += 1

                new_line1 = new_read_name
                new_line2 = a2
                new_line3 = a3
                new_line4 = a4

                tf.write(f"{new_line1}\n{new_line2}\n{new_line3}\n{new_line4}\n")

        return tf.name
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
def Fqs2_1fq(fa, fb, fc, out, nproc, chunk_size, temp_dir):

        pool = mp.Pool(processes=nproc)

        temp_files = []  
        with open(fa, 'r') as fA, open(fb, 'r') as fB, open(fc, 'r') as fC:
            async_results = []
            for (chunk_index, linesA, linesB, linesC) in read_in_chunks(fA, fB, fC, chunk_num_reads=chunk_size):
                res = pool.apply_async(
                    merge_chunk_to_tempfile,
                    (chunk_index, linesA, linesB, linesC, temp_dir)
                )
                async_results.append(res)

            pool.close()

            for r in async_results:
                tempfile_path = r.get()
                temp_files.append(tempfile_path)

            pool.join()

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
                read1seq  = arrays[1]
                lenread1 = len(read1seq)
                if spatialcode == '*'  or (int(line[1]) & 16) != 0 or lenread1 <  10:
                    thisline = ''
                else:
                    UMIseq = arrays[2]
                    read1qual = arrays[3][:lenread1]
                    read_name2 = '@' + read_name + '|:_:|' + spatialcode + ':' + UMIseq
                    thisline = read_name2 + '\n' + read1seq + '\n' + '+' + '\n' + read1qual + '\n'
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

def demultiplexing(R1, R2, barcode_file, PrimerStructure1, StructureUMI, StructureBarcode, threadnum, outputfolder, limitOutSAMoneReadBytes4barcodeMapping):
    os.makedirs(os.path.join(outputfolder, 'temps'), exist_ok=True)
    Cleanr1Fq1 = os.path.join(outputfolder, "temps/cleanr1fq1.fq")
    Cleanr1Fq2 = os.path.join(outputfolder, "temps/cleanr1fq2.fq")
    CombineFq = os.path.join(outputfolder, "combine.fq")
    barcode_db_fa = os.path.join(outputfolder, "temps/barcode_xy.fasta")
    barcode_db_path = os.path.join(outputfolder, "temps/barcode_db")
    index_fq = os.path.join(outputfolder, "temps/index.fastq")
    UMI_fq = os.path.join(outputfolder, "temps/UMI.fastq")
    temps_path = os.path.join(outputfolder, "temps/") 
    
    prefixread1 = PrimerStructure1.split('_', 1)[0]
    suffixread1 = PrimerStructure1.rsplit('_', 1)[-1]
    
    threadnum = str(threadnum)
    subprocess.run([ "cutadapt", "-e", "0.25", "-a", suffixread1, "--times", "4", "-g", prefixread1, "-j", threadnum, "-o", Cleanr1Fq1, "-p", Cleanr1Fq2, R1, R2])
    singleCutadapt(StructureUMI,UMI_fq,Cleanr1Fq2,threadnum)
    singleCutadapt(StructureBarcode,index_fq,Cleanr1Fq2,threadnum)
    with open(barcode_file, 'r') as barcodes_in, open(barcode_db_fa, 'w') as barcode_db_file:
      for line in barcodes_in:
          fields = line.strip().split('\t')
          header = f">{fields[1]}_{fields[2]}"
          sequence = fields[0]
          barcode_db_file.write(f"{header}\n{sequence}\n")
    subprocess.run([ "STAR", "--runMode", "genomeGenerate", "--runThreadN", threadnum, "--genomeDir", barcode_db_path, "--genomeFastaFiles", barcode_db_fa, "--genomeSAindexNbases", "7" ,"--limitGenomeGenerateRAM", "50000000000"])
    
    
    Fqs2_1fq(index_fq,Cleanr1Fq1,UMI_fq,CombineFq,16,1000000,temps_path)
    if limitOutSAMoneReadBytes4barcodeMapping != 'NA':
        limitOutSAMoneReadBytes4barcodeMapping = ['--limitOutSAMoneReadBytes', str(limitOutSAMoneReadBytes4barcodeMapping)]
    else:
        limitOutSAMoneReadBytes4barcodeMapping=[]
    subprocess.run([
    "STAR",
    "--runThreadN", threadnum,
    "--genomeDir", barcode_db_path,
    "--readFilesIn", CombineFq,
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
    ]+limitOutSAMoneReadBytes4barcodeMapping)
    
    filter_sam_nbhd(os.path.join(outputfolder, "temps/barcodeMapping/tempAligned.out.sam"), CombineFq)
    #filter_sam(os.path.join(outputfolder, "temps/barcodeMapping/tempAligned.out.sam"), CombineFq, int(threadnum), 1000000)
    
    os.replace(os.path.join(outputfolder, "temps/barcodeMapping/tempLog.final.out"),os.path.join(outputfolder, "barcodeMapping.out"))
    barcode_mapping_dir = os.path.join(outputfolder, "temps/barcodeMapping/")
    for file_name in os.listdir(barcode_mapping_dir):
        file_path = os.path.join(barcode_mapping_dir, file_name)
        if os.path.isfile(file_path):
            os.remove(file_path)  
    os.rmdir(barcode_mapping_dir) 


def get_barcode_for_single_cell(R1, R2, barcode_file, PrimerStructure1, StructureUMI, StructureBarcode, threadnum, outputfolder,barcode_threshold=100,barcodelength=0):
    barcode_threshold = int(barcode_threshold)
    os.makedirs(os.path.join(outputfolder, 'temps'), exist_ok=True)
    Cleanr1Fq1 = os.path.join(outputfolder, "temps/cleanr1fq1.fq")
    Cleanr1Fq2 = os.path.join(outputfolder, "temps/cleanr1fq2.fq")
    CombineFq = os.path.join(outputfolder, "combine.fq")
    barcode_db_fa = os.path.join(outputfolder, "temps/barcode_xy.fasta")
    barcode_db_path = os.path.join(outputfolder, "temps/barcode_db")
    index_fq = os.path.join(outputfolder, "temps/index.fastq")
    UMI_fq = os.path.join(outputfolder, "temps/UMI.fastq")
    temps_path = os.path.join(outputfolder, "temps/") 
    auto_bc_path = os.path.join(outputfolder, "temps", "auto_barcode.tsv")
    
    prefixread1 = PrimerStructure1.split('_', 1)[0]
    suffixread1 = PrimerStructure1.rsplit('_', 1)[-1]
    
    threadnum = str(threadnum)
    subprocess.run([ "cutadapt", "-e", "0.25", "-a", suffixread1, "--times", "4", "-g", prefixread1, "-j", threadnum, "-o", Cleanr1Fq1, "-p", Cleanr1Fq2, R1, R2])
    singleCutadapt(StructureUMI,UMI_fq,Cleanr1Fq2,threadnum)
    singleCutadapt(StructureBarcode,index_fq,Cleanr1Fq2,threadnum)

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
