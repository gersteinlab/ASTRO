# The script that is called when the workflow is set to “old.” Below is the explanation.

**1. CombFqs_AlignReads.py**

   Merges multiple FASTQ files and aligns them, typically consolidating reads before alignment.

**2. barcodedFq.pl**

   Inserts parsed barcode and UMI information into read headers to produce new FASTQ files containing full barcode/UMI details.

**3. ez_barcode.pl**

   Aggregates barcodes by coordinate positions and outputs an expression matrix indexed by (x, y), mainly used in spatial workflows.

**4. genemat2tsv.pl**

   Merge annotation and barcode data into the final gene expression matrix (TSV).

**5. interBED2GeneFile.py**

   Uses bedtools intersect to annotate BED data with GTF features, labeling each read with its corresponding gene and producing an annotated output.

**6. runBowtie.py**

   Wrapping Bowtie2 alignment in Python, similarly outputting read alignment information used for barcode or other matching tasks.

**7. redupBAM.pl**
   
   Deduplicates the STAR-generated BAM file, keeping only the best alignment records.

**8. singleCutadapt.py**

   A core Python script for barcode/UMI extraction, capable of multiple sequential calls to Cutadapt to trim specified segments and then merge them into one or more FASTQ files; it has been updated to improve parallel processing and merging logic.
