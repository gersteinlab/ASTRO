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


def singleCutadapt(barcodestr,outputfile,remainfa,threadnum,keepi=1,keepfilename='NA'):


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


    if keepfilename != 'NA':
        os.rename(linkfqs[keepi], keepfilename)
        del linkfqs[keepi]


    if len(linkfqs) >= 2:
        combineFqs(outputfile,linkfqs)
        for fqi in linkfqs:
            os.remove(fqi)
    else:
        if len(linkfqs) == 1:
            os.rename(linkfqs[0], outputfile)
    


def process_chunk(chunk_lines):
    read_count = set()
    read_delete = set()  
    line = chunk_lines[0] 
    read_name0 = line
    read_count.add(read_name0)
    for line in chunk_lines[1:]:
        read_name = line
        if read_name == read_name0:
            read_delete.add(read_name)
    read_count.add(read_name)
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
            if (int(line.split('\t')[1]) & 16) != 0: 
                continue
            read_name = line.split('\t')[0].split('--')[0]
            chunk.append(read_name)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        if chunk:
            yield chunk
def filter_samPair0(input_sam, reads_to_remove, output_fq1, output_fq2):
    with open(input_sam, 'r') as infile, open(output_fq1, 'w') as outfile1, open(output_fq2, 'w') as outfile2:
        for line in infile:
            if not line.startswith('@'):
                spatialcode = line.split('\t')[2]
                if spatialcode == '*':
                    continue
                arrays = line.split('\t')[0].split('--', maxsplit=3)
                read_name = arrays[0]
                if read_name not in reads_to_remove:
                    read1seq, read2seq   = arrays[1].split('_', maxsplit=1)
                    lenread1 = len(read1seq)
                    if lenread1 <  10:
                        continue
                    UMIseq = arrays[2]
                    read1qual = arrays[3][:lenread1]
                    read2qual = arrays[3][lenread1:]
                    read_name2 = '@' + read_name + '|:_:|' + spatialcode + ':' + UMIseq
                    line = read_name2 + '\n' + read1seq + '\n' + '+' + '\n' + read1qual + '\n'
                    outfile1.write(line)
                    line2 = read_name2 + '\n' + read2seq + '\n' + '+' + '\n' + read2qual + '\n'
                    outfile2.write(line2)
def filter_samPair(input_sam, output_fq1, output_fq2, num_processes, chunk_size):
    pool = mp.Pool(processes=int(num_processes))
    chunks = read_in_chunks_from_sam(input_sam, chunk_size=chunk_size)
    results = pool.map(process_chunk, chunks)
    pool.close()
    pool.join()

    reads_to_remove = total_delete(results) # multiple mapping is too annoying
    filter_samPair0(input_sam, reads_to_remove, output_fq1, output_fq2)


def merge_chunk_to_tempfilePair(chunk_index, linesA, linesB, linesB2, linesC, lowlength, highlength, temp_dir):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, dir=temp_dir, prefix=f"chunk_{chunk_index}_", suffix=".fastq") as tf:
            ii = 0
            for i in range(0, len(linesA), 4):
                a1 = linesA[i].rstrip('\n')   # @readA
                a2 = linesA[i+1].rstrip('\n') # SEQ
                a3 = linesA[i+2].rstrip('\n') # +
                a4 = linesA[i+3].rstrip('\n') # QUAL

                if len(a2)  == 0:
                    continue


                b1 = linesB[i].rstrip('\n')
                b2 = linesB[i+1].rstrip('\n')
                b3 = linesB[i+2].rstrip('\n')
                b4 = linesB[i+3].rstrip('\n')

                b21 = linesB2[i].rstrip('\n')
                b22 = linesB2[i+1].rstrip('\n')
                b23 = linesB2[i+2].rstrip('\n')
                b24 = linesB2[i+3].rstrip('\n')


                c1 = linesC[i].rstrip('\n')
                c2 = linesC[i+1].rstrip('\n')
                c3 = linesC[i+2].rstrip('\n')
                c4 = linesC[i+3].rstrip('\n')

                a1 = a1.split(' ')[0]
                b1 = b1.split(' ')[0]
                c1 = c1.split(' ')[0]
                if a1 != b1 or a1 != c1:
                    raise ValueError(f"Read names do not match: {a1} {b1} {c1}")


                new_read_name = f"@{chunk_index}_{ii}--{b2}_{b22}--{c2}--{b4}{b24}"
                ii += 1

                new_line1 = new_read_name
                new_line2 = a2
                new_line3 = a3
                new_line4 = a4

                tf.write(f"{new_line1}\n{new_line2}\n{new_line3}\n{new_line4}\n")
                    

        return tf.name
def read_in_chunksPair(fA, fB, fB2, fC, chunk_num_reads):
    chunk_index = 0
    while True:
        linesA = []
        linesB = []
        linesB2 = []
        linesC = []
        for _ in range(chunk_num_reads * 4):
            lineA = fA.readline()
            lineB = fB.readline()
            lineB2 = fB2.readline()
            lineC = fC.readline()
            if not lineA or not lineB or not lineC:
                break
            linesA.append(lineA)
            linesB.append(lineB)
            linesB2.append(lineB2)
            linesC.append(lineC)
        if not linesA:
            break
        yield (chunk_index, linesA, linesB, linesB2, linesC)
        chunk_index += 1
def Fqs2_2fq(fa, fb, fb2, fc, out, nproc, chunk_size, temp_dir, lowlength = 20, highlength = 60):
        import gzip

        def smart_open(path):
            with open(path, 'rb') as f:
                magic = f.read(2)
                if magic == b'\x1f\x8b':
                    return gzip.open(path, 'rt', encoding='utf-8')
                else:
                    return open(path, 'r')
        nproc = int(nproc)
        pool = mp.Pool(processes=nproc)

        temp_files = []  
        with open(fa, 'r') as fA, smart_open(fb) as fB, open(fb2, 'r') as fB2, open(fc, 'r') as fC:
            async_results = []
            
            for (chunk_index, linesA, linesB, linesB2, linesC) in read_in_chunksPair(fA, fB, fB2, fC, chunk_num_reads=chunk_size):
                res = pool.apply_async(
                    merge_chunk_to_tempfilePair,
                    (chunk_index, linesA, linesB, linesB2, linesC, lowlength, highlength, temp_dir)
                )
                async_results.append(res)

            pool.close()

            for r in async_results:
                tempfile_path = r.get()
                print(tempfile_path)
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



    


def demultiplexingPair(R1, R2, barcode_file, PrimerStructure1, StructureUMI, StructureBarcode, threadnum, outputfolder):
    os.makedirs(os.path.join(outputfolder, 'temps'), exist_ok=True)
    Cleanr1Fq1 = os.path.join(outputfolder, "temps/cleanr1fq1.fq")
    Cleanr1Fq2 = os.path.join(outputfolder, "temps/cleanr1fq2.fq")
    RemainFq2 = os.path.join(outputfolder, "temps/remainfq2.fq")
    CombineFq1 = os.path.join(outputfolder, "combine1.fq")
    CombineFq2 = os.path.join(outputfolder, "combine2.fq")
    barcode_db_fa = os.path.join(outputfolder, "temps/barcode_xy.fasta")
    barcode_db_path = os.path.join(outputfolder, "temps/barcode_db")
    index_fq = os.path.join(outputfolder, "temps/index.fastq")
    UMI_fq = os.path.join(outputfolder, "temps/UMI.fastq")
    temps_path = os.path.join(outputfolder, "temps/") 
    
    
    threadnum = str(threadnum)

    if PrimerStructure1 != "NA":
        prefixread1 = PrimerStructure1.split('_', 1)[0]
        suffixread1 = PrimerStructure1.rsplit('_', 1)[-1]
        subprocess.run([ "cutadapt", "-e", "0.25", "-a", suffixread1, "--times", "4", "-g", prefixread1, "-j", threadnum, "-o", Cleanr1Fq1, "-p", Cleanr1Fq2, R1, R2])
    else:
        Cleanr1Fq1 = R1
        Cleanr1Fq2 = R2


    
    singleCutadapt(StructureUMI, UMI_fq, Cleanr1Fq2, threadnum)
    singleCutadapt(StructureBarcode, index_fq, Cleanr1Fq2, threadnum, 2, RemainFq2)
    Fqs2_2fq(index_fq, Cleanr1Fq1, RemainFq2, UMI_fq, CombineFq1, threadnum, 1000000, temps_path)

    with open(barcode_file, 'r') as barcodes_in, open(barcode_db_fa, 'w') as barcode_db_file:
      for line in barcodes_in:
          fields = line.strip().split('\t')
          header = f">{fields[1]}_{fields[2]}"
          sequence = fields[0]
          barcode_db_file.write(f"{header}\n{sequence}\n")
    subprocess.run([ "STAR", "--runMode", "genomeGenerate", "--runThreadN", threadnum, "--genomeDir", barcode_db_path, "--genomeFastaFiles", barcode_db_fa, "--genomeSAindexNbases", "7" ,"--limitGenomeGenerateRAM", "50000000000"])
    
    
    subprocess.run([
    "STAR",
    "--runThreadN", threadnum,
    "--genomeDir", barcode_db_path,
    "--readFilesIn", CombineFq1,
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
    ])
    
    filter_samPair(os.path.join(outputfolder, "temps/barcodeMapping/tempAligned.out.sam"), CombineFq1, CombineFq2, int(threadnum), 1000000)
    os.replace(os.path.join(outputfolder, "temps/barcodeMapping/tempLog.final.out"),os.path.join(outputfolder, "barcodeMapping.out"))
    barcode_mapping_dir = os.path.join(outputfolder, "temps/barcodeMapping/")
    for file_name in os.listdir(barcode_mapping_dir):
        file_path = os.path.join(barcode_mapping_dir, file_name)
        if os.path.isfile(file_path):
            os.remove(file_path)  
    os.rmdir(barcode_mapping_dir) 