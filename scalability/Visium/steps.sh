# prepare cell barcode to position files.
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4565nnn/GSM4565823/suppl/GSM4565823%5FP4%5Frep1%5Ftissue%5Fpositions%5Flist.csv.gz
zcat GSM4565823_P4_rep1_tissue_positions_list.csv.gz   | awk -F',' 'BEGIN { OFS="\t" } { print $1, $4, $3 }'   |  sed 's/-1//' > data/GSM4565823_P4_rep1_tissue_positions_list.txt
rm GSM4565823_P4_rep1_tissue_positions_list.csv.gz


# prepare gtf files, they are also in this Github repository
zcat ../../Built_GTFs/hsa.mod.gtf.gz > hsa.mod.gtf 
# please build STAR reference folder for GRCh38: GRCh38_StarIndex


# download, merge and trim the input fastqs

MERGED="cache/P4_B1_S1_R2.merged.fastq.gz"
MERGED_R1="data/P4_B1_S1_R1.merged.fastq.gz"
TRIMMED="data/P4_B1_S1_R2.trimmed.fastq.gz"
CUTADAPT_TXT="cache/cutadapt_cSCC.txt"
CUTADAPT_JSON="cache/cutadapt_cSCC.json"
THREAD=4

for SRR in SRR11832354 SRR11832355 SRR11832356 SRR11832357; do
    fasterq-dump $SRR --split-files -e $THREAD --progress -O cache/
    pigz -p $THREAD cache/${SRR}*.fastq
done

mkdir -p data
cat cache/SRR11832354_1.fastq.gz cache/SRR11832355_1.fastq.gz cache/SRR11832356_1.fastq.gz cache/SRR11832357_1.fastq.gz > "${MERGED}"
cat cache/SRR11832354_2.fastq.gz cache/SRR11832355_2.fastq.gz cache/SRR11832356_2.fastq.gz cache/SRR11832357_2.fastq.gz > "${MERGED_R1}"


########################
# trimed by cutadapet, this is legecy code, it has two problems:
# 1.ASTRO could do the trim by itself.
# 2. cutadapter may lead R1 and R2 have different results.
# But it works so we do not want to change all the codes
    TSO_SEQ="AAGCAGTGGTATCAACGCAGAGTACATGGG"
    MIN_OVL=5
    ERR=0.1
    POLY_A_ADAPT="A{10}"
    cutadapt -j $THREAD \
    -e "${ERR}" -O "${MIN_OVL}" \
    -g "^${TSO_SEQ}" \
    -a "${POLY_A_ADAPT}" \
    --json "${CUTADAPT_JSON}" \
    -o - "${MERGED}" 2> "${CUTADAPT_TXT}" | pigz -c -p $THREAD > "${TRIMMED}"

    # run ASTRO on fastq files
    ASTRO parameter.json
    ########################
    #
    #
    #
    #
    #
    rm cache/*
    rmdir cache/

########################
# run ASTRO which do the trim step
    mv ${MERGED} data/
    ASTRO parameter2.json
    ########################
    #
    #
    #
    #
    #
    rm cache/*
    rmdir cache/
    
