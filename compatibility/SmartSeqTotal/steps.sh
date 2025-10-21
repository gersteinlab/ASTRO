# prepare gtf files, they are also in this Github repository
zcat ../../Built_GTFs/hsa.mod.gtf.gz > hsa.mod.gtf 
# please build STAR reference folder for GRCh38: GRCh38_StarIndex


#########################
#########################
# download data from GSE151334
# data/GSE151334set.humanmat.txt is downloaded from GEO database, however, it seems that previous link does not work.
# simple processing geo data
    cut -f10 data/GSE151334set.humanmat.txt | tail -n +2 > cache/srx.list
    prefetch --output-directory ./sra --option-file cache/srx.list
    
    python3 scripts/addsrr.py
    fasterq-dump --split-files --threads 8 --outdir ./fastq ./sra/*.sra
    
    for dir in ./sra/SRR*; do
        fasterq-dump --split-files --threads 1 --outdir ./fastq $dir/*.sra
    done
    python3 scripts/mergeFqs21.py
    cut -f13 cache/GSE151334set.humanmat_with_srr.txt | tail -n +2 | awk 'BEGIN{OFS="\t"}{ print $1, NR, NR }' > cache/SMARTposition.txt
    
    rm fastq/*
    rmdir fastq/
    rm sra/*
    rmdir sra/




#########################
#########################
# or you can directly use my download or I give another pipeline based on SRA Run Selector  https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA635600
# click the button Metadata to download	SraRunTable.csv
# put it in data/
# simple processing geo data

    cat data/SraRunTable.csv | sed 's/,/\t/g' | egrep '(cell_line|MCF7 cell|primary cell|HEK293T cell)' | egrep '(LibraryLayout|NovaSeq)' > cache/SraRunTable_human_smartseq.tsv
    cut -f1 cache/SraRunTable_human_smartseq.tsv | tail -n +2 > cache/sra.list
    prefetch --output-directory ./sra --option-file cache/sra.list
    
    fasterq-dump --split-files --threads 8 --outdir ./fastq ./sra/*.sra
    
    for dir in ./sra/SRR*; do
        fasterq-dump --split-files --threads 1 --outdir ./fastq $dir/*.sra
    done
    python3 scripts/mergeFqs21_SraRunTable.py
    cut -f1 cache/SraRunTable_human_smartseq.tsv | tail -n +2 | awk 'BEGIN{OFS="\t"}{ print $1, NR, NR }' > cache/SMARTposition.txt 

    rm fastq/*
    rmdir fastq/
    rm sra/*
    rmdir sra/



# Then run ASTRO on SMART-seq-all fastq file
ASTRO parameter.json


# download expression matrix from original paper
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151334/suppl/GSE151334%5Fcounts.human.tsv.gz
mv GSE151334_counts.human.tsv.gz data/


