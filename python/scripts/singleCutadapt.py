import re
import argparse
import os
import gzip
import subprocess

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



ap = argparse.ArgumentParser()
ap.add_argument("-b", "--barcodestr", required=True, help="Workpath with the BAM file")
ap.add_argument("-o", "--outputfile", required=True, help="Path to input BAM file")
ap.add_argument("-i", "--inputfa", required=True, help="Path to input BAM file")
ap.add_argument("-t", "--thread", required=False, default="8", type=str, help="Number of threads")
args = vars(ap.parse_args())


barcodestr = args["barcodestr"]
outputfile = args["outputfile"]
remainfa = args["inputfa"]
threadnum = args["thread"]


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

