#processing RNAcentral gtf:
########################################
########################################
# make sure bedtools is available for your envirnoment, it is needed for some scripts, if not, there will be error.
# If there is no error without bedtools, it means the scripts you uesd do not need bedtools


### this download maybe change with RNAcentral version, you can see the version at RNAcentral/release_notes.txt
wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/gff3/homo_sapiens.GRCh38.gff3.gz
###
mkdir -p RNAcentral
mv homo_sapiens.GRCh38.gff3.gz RNAcentral/homo_sapiens.GRCh38.gff3.gz
gunzip RNAcentral/homo_sapiens.GRCh38.gff3.gz
cat RNAcentral/homo_sapiens.GRCh38.gff3 | grep $'\t'transcript$'\t' | egrep -v '^[^0-9XY]+[\t]' | sed 's/^/chr/' > RNAcentral/homo_sapiens.GRCh38-addchr_findtrans.gff3
python3 scripts/div_clpsGTF.py -f RNAcentral/ -i RNAcentral/homo_sapiens.GRCh38-addchr_findtrans.gff3 -l vault_RNA:Y_RNA
########################################
########################################







#processing gencode gtf:
########################################
########################################
mkdir -p gencode/
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz
gunzip gencode.v45.annotation.gtf.gz
mv gencode.v45.annotation.gtf gencode/
cat gencode/gencode.v45.annotation.gtf | egrep -v '^#' | awk '/gene_name "Y_RNA"/ {gsub(/gene_type/,  "old_type"); gsub(/gene_name/,  "gene_type"); gsub(/gene_id/,  "gene_name"); print; next} 1' | grep -v 'gene_type "miRNA";' > gencode/gencode.v45.formod.gtf
########################################
########################################




#processing mirna gtf:
########################################
########################################
mkdir -p mirbase/
wget https://www.mirbase.org/download/hsa.gff3
mv hsa.gff3 mirbase/
cat mirbase/hsa.gff3 | grep -v '#' | sed s/$'\t'miRNA_primary_transcript$'\t'/$'\t'exon$'\t'/ | sed s/$'\t'miRNA$'\t'/$'\t'exon$'\t'/ | sed -e 's/ID=[^;]*;/ID=miRNA;/' > mirbase/hsa.formod.gtf 
########################################
########################################




#processing GtRNAdb gtf:
# download GtRNAdb/hg38-tRNAs.fa as GtRNAdb/Screenshot 2024-05-01 at 9.54.45 PM.png ,that is,
# but sometimes the wget does not work so just download it form browser manually
########################################
########################################
#wget --no-check-certificate https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa
#mkdir -p GtRNAdb/
#mv hg38-tRNAs.fa GtRNAdb/hg38-tRNAs.fa.txt
cat GtRNAdb/hg38-tRNAs.fa.txt | grep '>' | awk 'BEGIN{OFS="\t"} {split($(NF-1), parts, "[:-]"); gsub(/[()]/, "", $NF); gsub(/>/, "", $1); print parts[1], "ENSEMBL", "exon", parts[2], parts[3], ".", $NF, ".", $1"__tRNA"  }' > GtRNAdb/hg38-tRNAs.mod.gtf
########################################
########################################






#processing piRBase gtf:
########################################
########################################
mkdir -p piRBase/
wget http://bigdata.ibp.ac.cn/piRBase/download/v3.0/fasta/hsa.gold.fa.gz
wget http://bigdata.ibp.ac.cn/piRBase/download/v3.0/bed/hsa.align.bed.gz
gunzip hsa.align.bed.gz 
gunzip hsa.gold.fa.gz 
mv hsa.align.bed piRBase/
mv hsa.gold.fa piRBase/
python3 scripts/filter_pi2gtf.py -o piRBase/hsa.piRNA.gtf -r piRBase/hsa.gold.fa -i piRBase/hsa.align.bed
########################################
########################################








#modifying gtfs:
########################################
########################################
python3 scripts/modGTF.py -i gencode/gencode.v45.formod.gtf   -T gene_type -G gene_name -f gencode/gencode.v45.mod.gtf
python3 scripts/modGTF.py -i RNAcentral/vault_RNA.gtf   -T gene_type -G gene_name -f RNAcentral/vault_RNA.mod.gtf
python3 scripts/modGTF.py -i RNAcentral/Y_RNA.gtf   -T gene_type -G gene_name -f RNAcentral/Y_RNA.mod.gtf
python3 scripts/modGTF.py -i mirbase/hsa.formod.gtf   -T ID -G Name -f mirbase/hsa.mod.gtf
python3 scripts/modGTF.py -i piRBase/hsa.piRNA.gtf -T gene_type -G gene_id -f piRBase/hsa.piRNA.mod.gtf
########################################
########################################







#final gtfs
########################################
########################################
python3 scripts/clpsGTF.py -i mirbase/hsa.mod.gtf:GtRNAdb/hg38-tRNAs.mod.gtf:RNAcentral/vault_RNA.mod.gtf:RNAcentral/Y_RNA.mod.gtf:gencode/gencode.v45.mod.gtf:piRBase/hsa.piRNA.mod.gtf -o hsa.with_piRNA.gtf
python3 scripts/clpsGTF.py -i mirbase/hsa.mod.gtf:GtRNAdb/hg38-tRNAs.mod.gtf:RNAcentral/vault_RNA.mod.gtf:RNAcentral/Y_RNA.mod.gtf:gencode/gencode.v45.mod.gtf -o hsa.mod.gtf
########################################
########################################

#only reserve right chrs.
awk '$1 !~ /[_\.]/' GRCh38/hsa.with_piRNA.gtf | sed 's/^MT/M/' > GRCh38/hsa.with_piRNA.gtf2
awk '$1 !~ /[_\.]/' GRCh38/hsa.mod.gtf | sed 's/^MT/M/' > GRCh38/hsa.mod.gtf2
mv GRCh38/hsa.with_piRNA.gtf2 GRCh38/hsa.with_piRNA.gtf
mv GRCh38/hsa.mod.gtf2 GRCh38/hsa.mod.gtf

