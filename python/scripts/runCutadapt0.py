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



ap = argparse.ArgumentParser()
ap.add_argument("-b", "--barcodestr", required=True, help="Workpath with the BAM file")
ap.add_argument("-o", "--outputfolder", required=True, help="Path to input BAM file")
ap.add_argument("-i", "--inputfa", required=True, help="Path to input BAM file")
ap.add_argument("-t", "--thread", required=False, default="8", type=str, help="Number of threads")
args = vars(ap.parse_args())


barcodestr = args["barcodestr"]
outputfolder = args["outputfolder"]
remainfa = args["inputfa"]
threadnum = args["thread"]

outputfolder = re.sub('/+$','',outputfolder)+'/'
array4input = re.split('_',barcodestr)

matchstr = [bool(re.search(r'^[ACGTacgtnN]+$', s)) for s in array4input]
matchnum = [bool(re.search(r'^([0-9]*|b)$', s)) for s in array4input]
if (sum(matchnum)+sum(matchstr)) != len(matchstr):
    exit('Error: wrong barcode format: unexpected code')

matchstri = [i for i, x in enumerate(matchstr) if x]
matchnumi = [i for i, x in enumerate(matchnum) if x]

matchstrdiff = [matchstri[i+1] - matchstri[i] for i in range(len(matchstri)-1) ]
matchnumdiff = [matchnumi[i+1] - matchnumi[i] for i in range(len(matchnumi)-1) ]

if any(element != 2 for element in matchstrdiff):
    exit('Error: wrong barcode format: order error')
if any(element != 2 for element in matchnumdiff):
    exit('Error: wrong barcode format: order error')


customsetting = '-j ' + threadnum

tmpi = 0
finalfqs = []
for stri in matchstri:
    thestr = [array4input[stri]]
    if (stri - 1) in matchnumi:
        tmpfaname = outputfolder + 'temp_' + str(tmpi) + '.fastq'
        outputcommand = ['cutadapt'] + [customsetting] + ['-a'] + thestr + ['-o'] + [tmpfaname] + [remainfa]
        finalfqs.append(tmpfaname)
        subprocess.run(' '.join(outputcommand), shell=True)
        tmpi = tmpi+1
    tmpfaname = outputfolder + 'temp_' + str(tmpi) + '.fastq'
    outputcommand = ['cutadapt'] + [customsetting] + ['-g'] + thestr + ['-o'] + [tmpfaname] + [remainfa]
    remainfa = tmpfaname
    tmpi = tmpi+1
    subprocess.run(' '.join(outputcommand), shell=True)


if matchnumi[-1]:
    finalfqs.append(tmpfaname)
print('Barcode Fastqs')    
print(finalfqs)

for i,x in enumerate(matchnumi):
    if array4input[x] != 'b':
        designedlength = int(array4input[x])
        forwardvalue = False
        if x == len(array4input)-1:
            forwardvalue = True
            if designedlength != [0]:
                cutfastq(finalfqs[i],forwardvalue,designedlength)
        elif x == 0:
            if designedlength != [0]:
                cutfastq(finalfqs[i],forwardvalue,designedlength)
                
                
linkfqs = []
for i,x in enumerate(matchnumi):
    if array4input[x] == 'b':
        linkfqs.append(finalfqs[i])


iterators0 = [gzip.open(file_path, "rt") if file_path.endswith('.gz') else open(file_path) for file_path in linkfqs]
iterators = [simple_fastq_iterator(file) for file in iterators0]

with open(outputfolder + 'index.fastq', 'w') as out_handle:
	for entries in zip(*iterators):
	        newseq = ""
	        newqual = ""
	        newtitle = ''
	        for title, seq, qual in entries:
	            if not newtitle:
	                newtitle = title
	            newseq += seq
	            newqual += qual
	        out_handle.write("@%s\n%s\n+\n%s\n" % (title, newseq, newqual))
