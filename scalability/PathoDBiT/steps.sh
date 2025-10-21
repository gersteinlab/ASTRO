# prepare gtf files, they are also in this Github repository
zcat ../../Built_GTFs/mmu.mod.gtf.gz > mmu.mod.gtf 
zcat ../../Built_GTFs/mmu.mod.gtf.gz | grep '__miRNA' | cut -f 9 | sort | uniq > mmu.mir2check.txt
zcat ../../Built_GTFs/hsa.mod.gtf.gz > hsa.mod.gtf 
zcat ../../Built_GTFs/hsa.mod.gtf.gz | grep '__miRNA' | cut -f 9 | sort | uniq > hsa.mir2check.txt
# please build STAR reference folder for mm39 and GRCh38
# just like GRCh38_StarIndex/   or mm39_STARref/


# MALT sample 
# get fastq files from GSM8454082
# organized files to get MALTinputR1.fq.gz and MALTinputR2.fq.gz
ASTRO MALT.json

# healthy donor lymph node
# get fastq files from GSM8454083
# organized files to get healthyinputR1.fq.gz and healthyinputR2.fq.gz

ASTRO healthy_donor_lymph_node.json


# mouse embryo 1
# get fastq files from GSM8454077
# organized files to get mmembryo1inputR1.fq.gz and mmembryo1inputR2.fq.gz
ASTRO mmEmb1.json


# mouse embryo 2
# get fastq files from GSM8454077
# organized files to get mmembryo2inputR1.fq.gz and mmembryo2inputR2.fq.gz
ASTRO mmEmb2.json


