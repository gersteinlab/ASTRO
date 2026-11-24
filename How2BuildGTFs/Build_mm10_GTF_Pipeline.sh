#processing gencode gtf:
########################################
########################################
mkdir -p gencode
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz
gunzip gencode.vM23.annotation.gtf.gz
mv gencode.vM23.annotation.gtf gencode
cat gencode/gencode.vM23.annotation.gtf | egrep -v '^#' | awk '/Gm55767|Gm56322|Gm56480|Gm56181|Gm55795|Gm56246|Gm56393|Gm55481|Gm54376|Gm54851/ {gsub(/misc_RNA/,  "Y_RNA"); print; next} 1' | grep -v 'gene_type "miRNA";' > gencode/gencode.vM23.formod.gtf

########################################
########################################





#processing mirna gtf:
########################################
########################################
mkdir -p mirbase/
wget https://www.mirbase.org/download/mmu.gff3
mv mmu.gff3 mirbase/
cat mirbase/mmu.gff3 | grep -v '#' | sed s/$'\t'miRNA_primary_transcript$'\t'/$'\t'exon$'\t'/ | sed s/$'\t'miRNA$'\t'/$'\t'exon$'\t'/ | sed -e 's/ID=[^;]*;/ID=miRNA;/' > mirbase/mmu.formod.gtf 
########################################
########################################







#modifying gtfs:
########################################
########################################
python3 scripts/modGTF.py -i gencode/gencode.vM34.formod.gtf  -T gene_type -G gene_name -f gencode/gencode.vM34.mod.gtf
python3 scripts/modGTF.py -i mirbase/mmu.formod.gtf   -T ID -G Name -f mirbase/mmu.mod.gtf
########################################
########################################







#final gtfs
########################################
########################################
python3 scripts/clpsGTF.py -i mirbase/mmu.mod.gtf:GtRNAdb/mm39-tRNAs.mod.gtf:RNAcentral/vault_RNA.mod.gtf:RNAcentral/Y_RNA.mod.gtf:gencode/gencode.vM34.mod.gtf:piRBase/mmu.piRNA.mod.gtf -o mmu.with_piRNA.gtf -C
python3 scripts/clpsGTF.py -i mirbase/mmu.mod.gtf:GtRNAdb/mm39-tRNAs.mod.gtf:RNAcentral/vault_RNA.mod.gtf:RNAcentral/Y_RNA.mod.gtf:gencode/gencode.vM34.mod.gtf -o mmu.mod.gtf -C
########################################
########################################



