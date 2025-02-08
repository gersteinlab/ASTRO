import re
import argparse
import gzip

def simple_fastq_iterator(handle):
    while True:
        title_line = handle.readline().rstrip()
        seq_string = handle.readline().rstrip()
        sep_line = handle.readline().rstrip()
        qual_string = handle.readline().rstrip()
        if not (title_line and seq_string and sep_line and qual_string):
            break
        yield title_line[1:], seq_string, qual_string



ap = argparse.ArgumentParser()
ap.add_argument("-o", "--outputfolder", required=True, help="Path to input BAM file")
ap.add_argument("-i", "--inputfa", required=True, help="Path to input BAM file")
ap.add_argument("-r", "--inputfa2", required=True, help="Path to input BAM file")
args = vars(ap.parse_args())

outputfolder = args["outputfolder"]
outputfolder = re.sub('/+$','',outputfolder)+'/'
inputfa = args["inputfa"]
inputfa2 = args["inputfa2"]


inputfas  = re.split(':',inputfa)
inputfas2 = re.split(':',inputfa2)

all_reads = set()

for fastq_file in inputfas+inputfas2:
    current_reads = set()
    with (gzip.open(fastq_file, "rt") if fastq_file.endswith('.gz') else open(fastq_file)) as f:
        for line in f:
            if line.startswith('@') or line.startswith('>'):
                read_name = line.split()[0][1:]
                if read_name in current_reads:
                    raise ValueError(f"Duplicate read name found: {read_name} in file {fastq_file}")
                current_reads.add(read_name)
        if len(all_reads) == 0:
            all_reads = current_reads
        else:
            all_reads = all_reads & current_reads
            
for fastq_file in inputfas:
    with (gzip.open(fastq_file, "rt") if fastq_file.endswith('.gz') else open(fastq_file)) as f:
        with open(outputfolder+'Flt_'+fastq_file, "w") as out_handle:
            for title, seq, qual in simple_fastq_iterator(in_handle):
                if title in all_reads:
                    out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))


if len(inputfas2) > 0:
    iterators0 = [gzip.open(file_path, "rt") if file_path.endswith('.gz') else open(file_path) for file_path in inputfas2]
    #iterators0 = [open(file_path, 'r') for file_path in inputfas2]
    iterators = [simple_fastq_iterator(file) for file in iterators0]

    with open(outputfolder+'barcodes.fastq', 'w') as out_handle:
        for entries in zip(*iterators):
            newseq = ""
            newqual = ""
            newtitle = ''
            for title, seq, qual in entries:
                if not newseq:
                    newtitle = title
                newseq += seq
                newqual += qual
            out_handle.write("@%s\n%s\n+\n%s\n" % (title, newseq, newqual))


for file in iterators0:
    file.close()