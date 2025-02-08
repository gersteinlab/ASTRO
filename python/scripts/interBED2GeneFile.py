import pandas as pd
import re
import io
import subprocess
import argparse
import os


ap = argparse.ArgumentParser()
ap.add_argument("-b", "--inputbam", required=True, help="Workpath with the BAM file")
ap.add_argument("-o", "--outputfile", required=True, help="Path to input BAM file")
ap.add_argument("-C", "--chunk_size", required=False, default=100000, help="read line number 1 time")
ap.add_argument("-g", "--gtffile", required=True,  help="Gene type name")
args = vars(ap.parse_args())


inputbam = args["inputbam"]
outputfile = args["outputfile"]
chunk_size = args["chunk_size"]
gtffile  = args["gtffile"]


inputbampath = os.path.dirname(inputbam)
tempbed = inputbampath + '/bamtemp.bed'
interbed = inputbampath + '/tempinter.bed'


command = "bedtools bamtobed -i " + inputbam + " -split > " + tempbed
subprocess.check_output(command, shell=True, universal_newlines=True)
command = 'bedtools intersect -a ' + tempbed + " -b " + gtffile + ' -wao | awk \'BEGIN{OFS="\t"} { if ($NF != 0) { print $1, $3-$2, $4, $6, $(NF-5)-$(NF-6)+1, $(NF-1), $NF }}\' > ' + interbed 
subprocess.check_output(command, shell=True, universal_newlines=True)




def process_chunk(chunk):
    chunk['valuetempit'] = (chunk.iloc[:,6] - (chunk.iloc[:,1] - chunk.iloc[:,6]))  / chunk.iloc[:, 4]
    #chunk = chunk[chunk.iloc[:,7] > 0].copy()
    #ranks = chunk.groupby(chunk.columns[2])['valuetempit'].rank(method='dense', ascending=False)
    #return chunk[ranks == 1]
    inforvector = [item.split(":")[0] for item in chunk.iloc[:,2]]
    chunk.iloc[:,2] = inforvector
    return chunk.loc[chunk.groupby(chunk.columns[2])['valuetempit'].idxmax()]



chunks = []
for chunk in pd.read_csv(interbed, sep='\t', header=None, chunksize=chunk_size):
    processed_chunk = process_chunk(chunk)
    processed_chunk = processed_chunk.iloc[:,[2,5,7]]
    chunks.append(processed_chunk)




final_result = pd.concat(chunks)
final_result = final_result.loc[final_result.groupby(final_result.columns[0])['valuetempit'].idxmax()]


split_data = [item.split("_") for item in final_result.iloc[:,0]]
vector1 = [int(pair[1]) for pair in split_data]
vector2 = [int(pair[2]) for pair in split_data]


genemat = {
    'X': vector1,
    'Y': vector2,
    'Gene': final_result.iloc[:,1]
}
genemat = pd.DataFrame(genemat)
genemat.to_csv(outputfile, sep='\t', index=False, header=False)

