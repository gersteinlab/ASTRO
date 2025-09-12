#processing RNAcentral gtf:
########################################
########################################
### this download maybe change with RNAcentral version, you can see the version at RNAcentral/release_notes.txt
wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/gff3/mus_musculus.GRCm39.gff3.gz
####
mkdir -p RNAcentral
mv mus_musculus.GRCm39.gff3.gz RNAcentral/mus_musculus.GRCm39.gff3.gz
gunzip RNAcentral/mus_musculus.GRCm39.gff3.gz
cat RNAcentral/mus_musculus.GRCm39.gff3 | grep $'\t'transcript$'\t' | egrep -v '^[^0-9XY]+[\t]' | sed 's/^/chr/' > RNAcentral/mus_musculus.GRCm39-addchr_findtrans.gff3
python3 scripts/div_clpsGTF.py -f RNAcentral/ -i RNAcentral/mus_musculus.GRCm39-addchr_findtrans.gff3 -l vault_RNA:Y_RNA
########################################
########################################







#processing gencode gtf:
########################################
########################################
mkdir -p gencode
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M34/gencode.vM34.annotation.gtf.gz
gunzip gencode.vM34.annotation.gtf.gz
mv gencode.vM34.annotation.gtf.gz gencode
cat gencode/gencode.vM34.annotation.gtf | egrep -v '^#' | awk '/Gm55767|Gm56322|Gm56480|Gm56181|Gm55795|Gm56246|Gm56393|Gm55481|Gm54376|Gm54851/ {gsub(/misc_RNA/,  "Y_RNA"); print; next} 1' | grep -v 'gene_type "miRNA";' > gencode/gencode.vM34.formod.gtf

########################################
########################################





#processing mirna gtf:
########################################
########################################
mkdir -p mirbase/
wget https://www.mirbase.org/download/mmu.gff3
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToMm39.over.chain.gz
mv mm10ToMm39.over.chain.gz mirbase/
mv mmu.gff3 mirbase/
liftOver -gff mirbase/mmu.gff3 mirbase/mm10ToMm39.over.chain.gz mirbase/mmu.mm39.gff3 mirbase/mmu.mm39unmapped.gff3 # warning does not influence this task
cat mirbase/mmu.mm39.gff3 | grep -v '#' | sed s/$'\t'miRNA_primary_transcript$'\t'/$'\t'exon$'\t'/ | sed s/$'\t'miRNA$'\t'/$'\t'exon$'\t'/ | sed -e 's/ID=[^;]*;/ID=miRNA;/' > mirbase/mmu.formod.gtf 
########################################
########################################






#processing GtRNAdb gtf:
# download GtRNAdb/mm39-tRNAs.fa as GtRNAdb/Screenshot 2024-08-11 at 7.01.37 PM.png ,that is,
#wget https://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc39/mm39-tRNAs.fa
# mkdir -p GtRNAdb/
# mv mm39-tRNAs.fa GtRNAdb/mm39-tRNAs.fa.txt
# but sometimes the wget does not work so just download it form browser manually

########################################
########################################
cat GtRNAdb/mm39-tRNAs.fa.txt | grep '>' | awk 'BEGIN{OFS="\t"} {split($(NF-1), parts, "[:-]"); gsub(/[()]/, "", $NF); gsub(/>/, "", $1); print parts[1], "ENSEMBL", "exon", parts[2], parts[3], ".", $NF, ".", $1"__tRNA"  }' > GtRNAdb/mm39-tRNAs.mod.gtf
########################################
########################################







#processing piRBase gtf:
########################################
########################################
mkdir -p piRBase/
wget http://bigdata.ibp.ac.cn/piRBase/download/v3.0/fasta/mmu.gold.fa.gz
wget http://bigdata.ibp.ac.cn/piRBase/download/v3.0/bed/mmu.align.bed.gz
gunzip mmu.align.bed.gz 
gunzip mmu.gold.fa.gz 
mv mmu.align.bed piRBase/
mv mmu.gold.fa piRBase/
python3 scripts/filter_pi2gtf.py -o piRBase/hsmmua.piRNA.gtf -r piRBase/mmu.gold.fa -i piRBase/mmu.align.bed
bedtools getfasta -fi mm39.fa  -fo piRBase/mmu.piRNA_mm39.fa -bed piRBase/mmu.piRNA_mm39.bed -name -s
# use UCSC Genome browser (liftover  set "Minimum ratio of bases that must remap" as 1) change the genome version form GRCm38 to GRCm39/mm39  and the file mmu.piRNA_mm39.bed is gotten
python3 scripts/examineFAs.py -i piRBase/mmu.piRNA_mm39.fa -r piRBase/mmu.gold.fa -o piRBase/mmu.piRNA_mm39.gtf
########################################
########################################







#modifying gtfs:
########################################
########################################
python3 scripts/modGTF.py -i gencode/gencode.vM34.formod.gtf  -T gene_type -G gene_name -f gencode/gencode.vM34.mod.gtf
python3 scripts/modGTF.py -i RNAcentral/vault_RNA.gtf   -T gene_type -G gene_name -f RNAcentral/vault_RNA.mod.gtf
python3 scripts/modGTF.py -i RNAcentral/Y_RNA.gtf   -T gene_type -G gene_name -f RNAcentral/Y_RNA.mod.gtf
python3 scripts/modGTF.py -i mirbase/mmu.formod.gtf   -T ID -G Name -f mirbase/mmu.mod.gtf
python3 scripts/modGTF.py -i piRBase/mmu.piRNA_mm39.gtf        -T gene_type -G gene_id -f piRBase/mmu.piRNA.mod.gtf
########################################
########################################







#final gtfs
########################################
########################################
python3 scripts/clpsGTF.py -i mirbase/mmu.mod.gtf:GtRNAdb/mm39-tRNAs.mod.gtf:RNAcentral/vault_RNA.mod.gtf:RNAcentral/Y_RNA.mod.gtf:gencode/gencode.vM34.mod.gtf:piRBase/mmu.piRNA.mod.gtf -o mmu.with_piRNA.gtf -C
python3 scripts/clpsGTF.py -i mirbase/mmu.mod.gtf:GtRNAdb/mm39-tRNAs.mod.gtf:RNAcentral/vault_RNA.mod.gtf:RNAcentral/Y_RNA.mod.gtf:gencode/gencode.vM34.mod.gtf -o mmu.mod.gtf -C
########################################
########################################




