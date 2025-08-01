OUTPUT=$4
OUTPUT="${OUTPUT%/}"
mkdir -p $OUTPUT
mkdir -p $OUTPUT/temps
logf=$OUTPUT/log.txt

R1=$1 && echo $R1 > $logf
R2=$2 && echo $R2 >> $logf
barcodes=$3 && echo $barcodes >> $logf
echo $OUTPUT >> $logf # this is also $4
starref=${5:-data/reference/human}  && echo $starref >> $logf
gtffile=${6:-data/reference/RNAcentral_hsa}  && echo $gtffile >> $logf
PrimerStructure1=${7:-AAGCAGTGGTATCAACGCAGAGTGAATGGG_b_A{10\}N{150\}}  && echo $PrimerStructure1 >> $logf
StructureUMI=${8:-b_ATCCACGTGCTTGAGAGGCCAGAGCATTCG_b_GTGGCCGATGTTTCGCATCGGCGTACGACT_10}  && echo $StructureUMI >> $logf
StructureBarcode=${9:-b_ATCCACGTGCTTGAGAGGCCAGAGCATTCG_b_GTGGCCGATGTTTCGCATCGGCGTACGACT_10}  && echo $StructureBarcode >> $logf
barcodeposition=${10:-NA}  && echo $barcodeposition >> $logf #for access files in scripts/GETmat/
barcodelengthrange=${11:-NA}  && echo $barcodelengthrange >> $logf #for access files in scripts/GETmat/
scriptFolder=${12:-scripts/}  && echo $scriptFolder >> $logf #for access files in scripts/GETmat/
STARparamfile=${13:-NA} && echo $STARparamfile >> $logf
options=${14:-NA}


easymap=' -e '
if [[ $options =~ 'H' ]]
then
  easymap=' '
fi

scriptFolder="${scriptFolder%/}"
scriptFolder=$scriptFolder'/'
###intermediate
Cleanr1Fq1=$OUTPUT/temps/cleanr1fq1.fq
Cleanr1Fq2=$OUTPUT/temps/cleanr1fq2.fq
CombineFq=$OUTPUT/combine.fq
barcode_db_fa=$OUTPUT/temps/barcode_xy.fasta
barcode_db_path=$OUTPUT/temps/barcode_db
index_fq=$OUTPUT/temps/index.fastq
UMI_fq=$OUTPUT/temps/UMI.fastq
index0_out=$OUTPUT/temps/index0.out
index_out=$OUTPUT/temps/index.out
unmapfq=$OUTPUT/temps/unmap.fastq
###
STARoutput=$OUTPUT/STAR/temp
STARbam=${STARoutput}Aligned.sortedByCoord.out.bam
pureSTARbam=${STARoutput}filtered.bam
###
expmatbed=$OUTPUT/expmat.bed
expmattsv=$OUTPUT/expmat.tsv



if [ "$R1" != "NA" ]; then
  prefixread1=${PrimerStructure1%%_*}
  suffixread1=${PrimerStructure1##*_}
  cutadapt -e 0.25 -a $suffixread1 --times 4 -g $prefixread1 -j 16 -o $Cleanr1Fq1 -p $Cleanr1Fq2 $R1 $R2
  python3 ${scriptFolder}GETmat/singleCutadapt.py -i $Cleanr1Fq2 -o $UMI_fq -b $StructureUMI  -t ${SLURM_JOB_CPUS_PER_NODE}
  python3 ${scriptFolder}GETmat/singleCutadapt.py -i $Cleanr1Fq2 -o $index_fq -b $StructureBarcode  -t ${SLURM_JOB_CPUS_PER_NODE}
  paste <(paste -d '_' <(cut $barcodes -f 2) <(cut $barcodes -f 3)) <(cut $barcodes -f 1) | sed 's/^/>/g' | sed 's/\t/\n/g' > $barcode_db_fa
  bowtie2-build $barcode_db_fa $barcode_db_path
  if [ "$barcodeposition" == "NA" ]; then
    python3 ${scriptFolder}GETmat/runBowtie.py -i $index_fq -o $index_out -r $barcode_db_path -t 16
    perl ${scriptFolder}GETmat/barcodedFq.pl -r $OUTPUT/temps/index.out -u $UMI_fq -o $CombineFq -i $Cleanr1Fq1
  else
    perl ${scriptFolder}GETmat/ez_barcode.pl  -u $unmapfq -o $index0_out -r $barcodes -i $index_fq -b $barcodeposition -l $barcodelengthrange
    python3 ${scriptFolder}GETmat/runBowtie.py -i $unmapfq -o $index_out -r $barcode_db_path -t 16
    cat $index0_out $index_out > ${index_out}.temp && mv ${index_out}.temp $index_out
    perl ${scriptFolder}GETmat/barcodedFq.pl -r $index_out -u $UMI_fq -o $CombineFq -i $Cleanr1Fq1
  fi 
fi



if [ "$STARparamfile" == "NA" ]; then
  STARparamfile=${scriptFolder}GETmat/starpara.txt
fi
STARparams=$(cat $STARparamfile)

if [ "$starref" != "NA" ]; then
  if [ "$starref" != "NANA" ]; then
    STAR --genomeDir $starref --readFilesIn $CombineFq --outFileNamePrefix $STARoutput --runThreadN ${SLURM_JOB_CPUS_PER_NODE} --sjdbGTFfile $gtffile $STARparams
  fi
   perl ${scriptFolder}GETmat/redupBAM.pl -i $STARbam -o $pureSTARbam -p ${SLURM_JOB_CPUS_PER_NODE}
fi




singleconvert () {
    chrname=$1
    inputbam=$2
    prefix=$3
    gtfin=$4
    scriptFolder2=$5
    tempgtf=${prefix}_${chrname}.gtf
    tempbed=${prefix}_${chrname}.bed
    tempoutbed=${prefix}_${chrname}._tempout.bed
    echo "cat $gtfin | awk -v var="$chrname" '$1 == var' > $tempgtf"
    cat $gtfin | awk -v var="$chrname" '$1 == var' > $tempgtf
    bedtools bamtobed -i $inputbam -split | awk -v var="$chrname" '$1 == var' > $tempbed
    perl ${scriptFolder2}GETmat/interBED2GeneFile.pl -i $tempbed -o $tempoutbed -g $tempgtf 
}
export -f singleconvert
samtools view -H $pureSTARbam | grep '@SQ' | sed 's/.*SN://' | sed 's/\s.*//' | parallel -j $SLURM_CPUS_PER_TASK singleconvert {} $pureSTARbam $OUTPUT/temps/count $gtffile $scriptFolder
cat $OUTPUT/temps/count*_tempout.bed > ${expmatbed}.temp
perl ${scriptFolder}GETmat/mergeReadname.pl -i ${expmatbed}.temp -o $expmatbed
#rm ${expmatbed}.temp
#rm $OUTPUT/temps/count*
perl ${scriptFolder}GETmat/genemat2tsv.pl -o $expmattsv -d $barcodes -i $expmatbed -g $gtffile $easymap
