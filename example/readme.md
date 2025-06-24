# ASTRO Demo Guide

This directory contains example data and instructions for running the ASTRO pipeline demonstration.

## Prerequisites

Before running the demo, ensure you have:
- ASTRO installed and configured
- A STAR reference genome for mouse (mm39)
- Required bioinformatics tools (STAR, bedtools, samtools, cutadapt)

## Step 1: Uncompress the GTF file

Extract the compressed GTF annotation file:

```bash
tar -xJvf mmu.all.gtf.tar.xz
```

## Step 2: Configure STAR reference

The STAR reference is too large to be included in this repository. You need to:

1. Download or build a STAR reference for mouse genome (version: mm39)
2. Update the `starref` variable in the `test.json` file to point to your mm39 STAR index path

## Step 3: Run the pipeline

Execute the ASTRO pipeline with the test configuration:

```bash
ASTRO test.json
```

## Step 4: View example output

To examine the expected output format, extract the provided output example:

```bash
tar -xJvf output.tar.xz
```

**Note**: Some large files have been removed from the example output to reduce file size. The following files were removed:

```bash
rm output/STAR/tempLog.out
rm -rf output/STAR/temp_STARgenome/
rm -rf output/temps/barcode_db/*
```
