#!/usr/bin/env python
from ASTRO.ASTRO_core import ASTRO
from ASTRO.featurefilter import filtMATbyRT
import argparse
import json

def main():
    parser = argparse.ArgumentParser(description="get information")
    
    parser.add_argument("json_file_path1", nargs="?", default=None, help="json file for the input")
    parser.add_argument("--json_file_path", required = False, help="json file for the input")
    parser.add_argument("--R1", '--barcode_read', dest='R1', help="barcode read, fastq files containing input RNA")
    parser.add_argument("--R2", '--transcript_read', dest='R2', help="transcript read, fastq files including barcode information")
    parser.add_argument("--barcode_file", help="files including spatial barcodes")
    parser.add_argument("--outputfolder", help="output folder")
    parser.add_argument("--starref", help="STAR referebce folder")
    parser.add_argument("--gtffile", help="gtf file")
    parser.add_argument("--PrimerStructure1", help="structure for R1, like AAGCAGTGGTATCAACGCAGAGTGAATGGG_b_A{10}N{150}") 
    parser.add_argument("--StructureUMI", help="structure for UMI, like CAAGCGTTGGCTTCTCGCATCT_10") 
    parser.add_argument("--StructureBarcode", help="structure for UMI, like 20_ATCCACGTGCTTGAGAGGCCAGAGCATTCG:ATCCACGTGCTTGAGAGGCCAGAGCATTCG...GTGGCCGATGTTTCGCATCGGCGTACGACT")
    parser.add_argument("--barcodemode", choices=["singlecell", "spatial"], help="Barcode processing mode: 'singlecell' or 'spatial'")
    parser.add_argument("--genes2check", help="path to a file listing targets to validate (each line must equal GTF col9)") 
    parser.add_argument("--barcodeposition", default="NA")
    parser.add_argument("--barcodelengthrange", default="NA")
    parser.add_argument("--threadnum", required = False)
    parser.add_argument("--options", default="", help="H:hardmode for gene2tsv; M: samtools markdup for redup")
    parser.add_argument("--steps", help="1 => demultiplexing; 2 => genomemapping; 4 => feature counting")
    parser.add_argument("--STARparamfile4genome", default="NA", help="whether change the input of STAR for genome mapping")
    parser.add_argument("--qualityfilter", help="quality filter for reads")
    parser.add_argument("--removeByDim", help="remove wrong features by dimension")
    parser.add_argument("--filterlogratio", help="exclude extreme genes")
    parser.add_argument("--workflow", default="new", help="which workflow to run, old or new")
    parser.add_argument("--ReadLayout", default="singleend", help="which Read Layout, singleend or pairedend")
    parser.add_argument("--limitOutSAMoneReadBytes4barcodeMapping", default="NA", help="limitOutSAMoneReadBytes for barcode mapping")
    parser.add_argument("--not_organize_result", action="store_true", help="not try to organize outputfolder by removing tmp,  compresing files and moving important intermediate files to interim")

    args = parser.parse_args()
    workflow = args.workflow
    if args.json_file_path1:
        with open(args.json_file_path1, "r") as f:
            config = json.load(f)
        for k, v in config.items():
            setattr(args, k, v)
    if workflow == "old":
        from .olddriver import run_old_pipeline
        all_args = dict(vars(args))
        allowed_keys = {
            "R1", "R2", "barcode_file", "outputfolder", "starref", "gtffile",
            "PrimerStructure1", "StructureUMI", "StructureBarcode",
            "scriptFolder", "barcodeposition", "barcodelengthrange",
            "threadnum",
        }
        
        call_args = {}
        for k,v in all_args.items():
            if k in allowed_keys:
                call_args[k] = v
        
        run_old_pipeline(**call_args)
    else:
        from .ASTRO_core import ASTRO
        ASTRO(**vars(args))


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
