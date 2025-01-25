#!/usr/bin/env python
from ASTRO.ASTRO_core import ASTRO
from ASTRO.featurefilter import filtMATbyRT
import argparse

def main():
    parser = argparse.ArgumentParser(description="get information")
    
    parser.add_argument("json_file_path1", nargs="?", default=None, help="json file for the input")
    parser.add_argument("--json_file_path", required = False, help="json file for the input")
    parser.add_argument("--R1", help="fastq files containing input RNA")
    parser.add_argument("--R2", help="fastq files including barcode information")
    parser.add_argument("--barcode_file", help="files including spatial barcodes")
    parser.add_argument("--outputfolder", help="output folder")
    parser.add_argument("--starref", help="STAR referebce folder")
    parser.add_argument("--gtffile", help="gtf file")
    parser.add_argument("--PrimerStructure1", help="structure for R1, like AAGCAGTGGTATCAACGCAGAGTGAATGGG_b_A{10}N{150}") 
    parser.add_argument("--StructureUMI", help="structure for UMI, like CAAGCGTTGGCTTCTCGCATCT_10") 
    parser.add_argument("--StructureBarcode", help="structure for UMI, like 20_ATCCACGTGCTTGAGAGGCCAGAGCATTCG:ATCCACGTGCTTGAGAGGCCAGAGCATTCG...GTGGCCGATGTTTCGCATCGGCGTACGACT")
    parser.add_argument("--barcodeposition", default="NA")
    parser.add_argument("--barcodelengthrange", default="NA")
    parser.add_argument("--threadnum", required = False)
    parser.add_argument("--options", default="", help="H:hardmode for gene2tsv; M: samtools markdup for redup")
    parser.add_argument("--steps", help="1 => demultiplexing; 2 => genomemapping; 4 => feature counting")
    parser.add_argument("--STARparamfile4genome", default="NA", help="whether change the input of STAR for genome mapping")
    parser.add_argument("--qualityfilter", help="quality filter for reads")
    parser.add_argument("--removeByDim", help="remove wrong features by dimension")
    parser.add_argument("--filterlogratio", help="exclude extreme genes")
    args = parser.parse_args()

    args = vars(args)

    ASTRO(**args)

def filtmatbyrt():
    parser = argparse.ArgumentParser(description="get information")

    parser.add_argument("pos_expmatgood", nargs='?', default=None, help="Positional: expmatgood")
    parser.add_argument("pos_expmatbad", nargs='?', default=None, help="Positional: expmatbad")
    parser.add_argument("pos_finalexpmat", nargs='?', default=None, help="Positional: finalexpmat")
    parser.add_argument("pos_filterlogratio", nargs='?', default=None, help="Positional: filterlogratio")

    parser.add_argument("--expmatgood", required = False, help="json file for the input")
    parser.add_argument("--expmatbad", help="fastq files containing input RNA")
    parser.add_argument("--finalexpmat", help="fastq files including barcode information")
    parser.add_argument("--filterlogratio", default=2, help="files including spatial barcodes")
    args = parser.parse_args()

    expmatgood = args.expmatgood if args.expmatgood else args.pos_expmatgood
    expmatbad = args.expmatbad if args.expmatbad else args.pos_expmatbad
    finalexpmat = args.finalexpmat if args.finalexpmat else args.pos_finalexpmat
    filterlogratio = args.filterlogratio if args.filterlogratio else args.pos_filterlogratio

    filtMATbyRT(expmatgood, expmatbad, finalexpmat, filterlogratio)
