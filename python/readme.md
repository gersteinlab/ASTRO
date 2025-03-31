# Descriptions of the functions (or purposes) of each script:
**1. ASTRO_core.py**
   Defines the main ASTRO function, which reads parameters and sequentially calls core workflows such as demultiplexing, genome mapping, and feature counting.
**2. ASTRO_run.py**
   Provides a command-line entry point (main(), etc.), parsing arguments and then calling ASTRO or legacy flows to unify the CLI interface.
**3. init.py**
   Initializes the ASTRO package and provides common utility functions (like inter_bed2geneFile, bam_to_bed, etc.) for gene annotation and file format conversion.
**4. countfeature.py**
   Converts genome alignment results (BAM files) to BED format, annotates them to genes, and aggregates the data into an expression matrix (TSV).
**5. demultiplexer.py**
   Auto-generates or filters barcodes (single-cell or spatial) via Cutadapt + multi-process merging, aligns them with STAR, and outputs final FASTQs containing barcodes and UMIs.
**7. featurefilter.py**
   Filters and merges the expression matrix (e.g., via filtMATbyRT), and runs additional feature filtering for spatial mode to remove unreliable genes or barcodes.
**8. genomemapping.py**
   Runs STAR alignment to map reads to the genome, then performs deduplication (with markdup or a custom approach), producing a final filtered BAM file.
**9. olddriver.py**
   A driver script for the older pipeline (run_old_pipeline), containing legacy-style logic compatible with earlier workflows.
