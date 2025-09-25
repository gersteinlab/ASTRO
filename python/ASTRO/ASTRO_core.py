#!/usr/bin/env python3

import json
import os
import sys
from .genomemapping import genomemapping
from .countfeature import countfeature
from .featurefilter import featurefilter
from .demultiplexer import demultiplexing, get_barcode_for_single_cell
import subprocess



def auto_set_barcodes(structure_barcode_in: str):
    parts = structure_barcode_in.split(':')

    out_structs = []
    out_lens = []
    start_position = 1
    n = len(parts)

    for i, part in enumerate(parts, start=1):
        elems = part.split('_')

        if len(elems) == 2:
            a, b = elems
            if a.isdigit() and b.isdigit():
                out_lens.append(int(b))
                outi = f"{int(a)}_{int(b)}"
                out_structs.append(outi)
            elif a.isdigit():
                out_lens.append(int(a))
                if i == 1:
                    outi = f"{int(a) + 4}_{b}"
                    out_structs.append(outi)
                    start_position = 5
                else:
                    out_structs.append(part)
            elif b.isdigit():
                out_lens.append(int(b))
                if i == n:
                    outi = f"{a}_{int(b) + 4}"
                    out_structs.append(outi)
                else:
                    out_structs.append(part)
            else:
                raise ValueError(f"1. wrong barcode format at '{part}': expected digits around '_'")
        
        elif len(elems) == 3:
            left, mid, right = elems
            if not mid.isdigit():
                raise ValueError(f"2. wrong barcode length '{mid}' in '{part}'")
            out_structs.append(f"{left}...{right}")
            out_lens.append(int(mid))
        
        else:
            raise ValueError(f"3. wrong barcode format at '{part}'")

    barcode_length = sum(out_lens)
    structure_barcode = ':'.join(out_structs)
    barcode_position = f"{start_position}_{barcode_length}b"
    barcode_length_range = f"{barcode_length - 4}_{barcode_length + 12}"

    return structure_barcode, barcode_position, barcode_length_range


def ASTRO(**kwargs):
    """
    Main ASTRO pipeline function for spatial transcriptomics data processing.

    This function orchestrates the complete ASTRO workflow including:
    1. Demultiplexing: Adapter trimming, UMI, and barcode splitting
    2. Genome Mapping: STAR alignment and duplicate removal
    3. Feature Counting: Gene expression quantification from GTF annotations
    4. Feature Filtering: Quality control and filtering of genes/barcodes

    Args:
        **kwargs: Configuration parameters that can include:
            json_file_path (str): Path to JSON configuration file
            json_file_path1 (str): Alternative JSON file path parameter
            options (str): Pipeline options (e.g., 'H' for hard mode, 'M' for markdup)
            threadnum (int): Number of threads for parallel processing (default: 16)
            steps (int): Bitwise integer specifying which steps to run:
                - 1: Demultiplexing
                - 2: Genome Mapping
                - 4: Feature Counting and Filtering
                - 7: All steps (1+2+4)
            outputfolder (str): Output directory path
            R1 (str): Path to R1 FASTQ file (main RNA sequences)
            R2 (str): Path to R2 FASTQ file (barcodes and UMIs)
            barcode_file (str): Path to barcode coordinate file
            starref (str): Path to STAR genome index directory
            gtffile (str): Path to GTF annotation file
            PrimerStructure1 (str): R1 primer structure definition
            StructureUMI (str): UMI structure definition
            StructureBarcode (str): Spatial barcode structure definition
            workflow (str): Pipeline version ('new' or 'old', default: 'new')
            barcodemode (str): Mode for barcode processing ('spatial' or 'singlecell')
            And many other optional parameters...

    Returns:
        None: Results are written to the specified output directory

    Raises:
        SystemExit: If required parameters are missing or invalid
    """

    args = kwargs

    json_file_path = args["json_file_path"] or args["json_file_path1"]
    if json_file_path:
        with open(json_file_path, "r") as file:
            data = json.load(file)
    else:
        data = {}

    
    
    for k, v in data.items():
        if args.get(k) is None:
            args[k] = v

    WHITELIST = {'options', 'threadnum', 'barcodemode', 'R1', 'R2', 'PrimerStructure', 'StructureUMI', 'StructureBarcode', 'filterlogratio', 'starref', 'barcode_read', 'transcript_read', 'qualityfilter', 'workflow', 'removeByDim', 'barcode_file', 'json_file_path', 'json_file_path1',
    'gtffile', 'steps','outputfolder', 'STARparamfile4genome','genes2check', 'barcode_threshold', 'barcodelength', 'barcodeposition', 'barcodelengthrange', 'ReadLayout', 'limitOutSAMoneReadBytes4barcodeMapping',  'not_organize_result', 'manually_set_barcode_details'}
    invalid_keys = set(args.keys()) - WHITELIST
    if invalid_keys:
        sys.exit(f"Error: invalid argument(s) detected: {', '.join(invalid_keys)}")

    args['options'] = args.get('options') or data.get('options') or ""
    args['threadnum'] = args.get('threadnum') or data.get('threadnum') or 16
    args['steps'] = args.get('steps') or data.get('steps') or 7
    args['steps'] = int(args['steps'])
    args["outputfolder"] = (
        args.get("outputfolder")
        or data.get("outputfolder")
        or sys.exit("outputfolder is not specified in both parser and json")
    )
    args["STARparamfile4genome"] = (
        args.get("STARparamfile4genome") or data.get("STARparamfile4genome") or "NA"
    )
    args['genes2check'] = args.get("genes2check", data.get("genes2check", False))
    args["barcode_threshold"] = int(
        args.get("barcode_threshold", data.get("barcode_threshold", 100))
    )
    args["barcodelength"] = int(
        args.get("barcodelength", data.get("barcodelength", 0))
    )
    
    args['barcodeposition']  = args.get('barcodeposition') or data.get('barcodeposition') or "NA"
    args['barcodelengthrange']  = args.get('barcodelengthrange') or data.get('barcodelengthrange') or "NA"
    args['ReadLayout']  = args.get('ReadLayout') or data.get('ReadLayout') or "singleend"
    args['limitOutSAMoneReadBytes4barcodeMapping']  = args.get('limitOutSAMoneReadBytes4barcodeMapping') or data.get('limitOutSAMoneReadBytes4barcodeMapping') or "NA"
    args['not_organize_result'] = args.get("not_organize_result", data.get("not_organize_result", False))
    args['manually_set_barcode_details'] = args.get("manually_set_barcode_details", data.get("manually_set_barcode_details", False))
    
    os.makedirs(args['outputfolder'], exist_ok=True)

    args["barcodemode"] = (
        "singlecell" if args.get("singlecell") or data.get("singlecell") else "spatial"
    )

    if args['steps'] & 1:
        if args.get('barcode_read') and args.get('R1'):
            sys.exit("barcode_read and R1 are same things, duplicated settings")
        if args.get('barcode_read'):
            args['R1'] = args.pop('barcode_read')

        if args.get('transcript_read') and args.get('R2'):
            sys.exit("transcript_read and R2 are same things, duplicated settings")
        if args.get('transcript_read'):
            args['R2'] = args.pop('transcript_read')
        
        if not args["R1"]:
            sys.exit("R1 is not specified in both parser and json")
        if not args["R2"]:
            sys.exit("R2 is not specified in both parser and json")
        
        args["PrimerStructure"] = (
            args.get("PrimerStructure") or data.get("PrimerStructure") or "NA"
        )
        
        args["StructureUMI"] = (
            args.get("StructureUMI")
            or data.get("StructureUMI")
            or sys.exit("StructureUMI is not specified in both parser and json")
        )
        args["StructureBarcode"] = (
            args.get("StructureBarcode")
            or data.get("StructureBarcode")
            or sys.exit("StructureBarcode is not specified in both parser and json")
        )
        
        if not args['manually_set_barcode_details']:
            args['StructureBarcode'], args['barcodeposition'], args['barcodelengthrange']  =  auto_set_barcodes(args['StructureBarcode'])

        if args['barcodemode'] == "singlecell":
            user_bc_file = (
                args.get("barcode_file") or data.get("barcode_file") or "notavailable"
            )

            final_bc_file = get_barcode_for_single_cell(
                read1=args['R1'],
                read2=args['R2'],
                barcode_file=user_bc_file,
                PrimerStructure=args['PrimerStructure'],
                StructureUMI=args['StructureUMI'],
                StructureBarcode=args['StructureBarcode'],
                threadnum=args['threadnum'],
                outputfolder=args['outputfolder'],
                barcode_threshold=args['barcode_threshold'],
                barcodelength=args['barcodelength']
            )

            if args["ReadLayout"] == "pairedend":
                from .demultiplexer2 import demultiplexingPair

                demultiplexingPair(
                    read1=args['R1'],
                    read2=args['R2'],
                    barcode_file=final_bc_file, 
                    PrimerStructure=args['PrimerStructure'],
                    StructureUMI=args['StructureUMI'],
                    StructureBarcode=args['StructureBarcode'],
                    threadnum=args['threadnum'], 
                    outputfolder=args['outputfolder'],
                    limitOutSAMoneReadBytes4barcodeMapping=args['limitOutSAMoneReadBytes4barcodeMapping'],
                )
            else:
                demultiplexing(
                    read1=args['R1'],
                    read2=args['R2'],
                    barcode_file=final_bc_file, 
                    PrimerStructure=args['PrimerStructure'],
                    StructureUMI=args['StructureUMI'],
                    StructureBarcode=args['StructureBarcode'],
                    threadnum=args['threadnum'], 
                    outputfolder=args['outputfolder'],
                    limitOutSAMoneReadBytes4barcodeMapping=args['limitOutSAMoneReadBytes4barcodeMapping'],
                    barcodeposition=args['barcodeposition'],
                    barcodelengthrange=args['barcodelengthrange']
                )

            args["barcode_file"] = final_bc_file

        else:
            bcfile = (
                args.get("barcode_file")
                or data.get("barcode_file")
                or sys.exit("barcode_file missing for spatial")
            )
            if args["ReadLayout"] == "pairedend":
                from .demultiplexer2 import demultiplexingPair

                demultiplexingPair(
                    read1=args['R1'],
                    read2=args['R2'],
                    barcode_file=bcfile, 
                    PrimerStructure=args['PrimerStructure'],
                    StructureUMI=args['StructureUMI'],
                    StructureBarcode=args['StructureBarcode'],
                    threadnum=args['threadnum'], 
                    outputfolder=args['outputfolder'],
                    limitOutSAMoneReadBytes4barcodeMapping=args['limitOutSAMoneReadBytes4barcodeMapping']
                )
            else:
                demultiplexing(
                    read1=args['R1'],
                    read2=args['R2'],
                    barcode_file=bcfile, 
                    PrimerStructure=args['PrimerStructure'],
                    StructureUMI=args['StructureUMI'],
                    StructureBarcode=args['StructureBarcode'],
                    threadnum=args['threadnum'], 
                    outputfolder=args['outputfolder'],
                    barcodeposition=args['barcodeposition'],
                    barcodelengthrange=args['barcodelengthrange'],
                    limitOutSAMoneReadBytes4barcodeMapping=args['limitOutSAMoneReadBytes4barcodeMapping']
                )
            args["barcode_file"] = bcfile

    if args["steps"] & 2:
        args["starref"] = (
            args.get("starref")
            or data.get("starref")
            or sys.exit("starref is not specified in both parser and json")
        )
        args["gtffile"] = (
            args.get("gtffile")
            or data.get("gtffile")
            or sys.exit("gtffile is not specified in both parser and json")
        )
        genomemapping(
            args["starref"],
            args["gtffile"],
            args["threadnum"],
            args["options"],
            args["outputfolder"],
            args["STARparamfile4genome"],
        )
        target_dir = os.path.join(args['outputfolder'], "STAR/")
        for root, dirs, files in os.walk(target_dir):
            for d in dirs:
                os.chmod(os.path.join(root, d), 0o755)
            for f in files:
                os.chmod(os.path.join(root, f), 0o644)
    if args["steps"] & 4:
        from .countfeature import countfeature

        bcfile = (
            args.get("barcode_file")
            or data.get("barcode_file")
            or sys.exit("barcode_file missing in Step4")
        )
        gtffile = (
            args.get("gtffile")
            or data.get("gtffile")
            or sys.exit("gtffile missing in Step4")
        )
        qualityfilter = (
            args.get("qualityfilter") or data.get("qualityfilter") or "25:0.75"
        )

        usedgtf = countfeature(args['gtffile'], args['threadnum'], args['options'], bcfile, args['outputfolder'], qualityfilter, args['genes2check'])
        
        if qualityfilter not in ['0:0','NA']:
            os.makedirs(os.path.join(args['outputfolder'], "filteredout/"), exist_ok=True)
            expmatbedexcl = os.path.join(args['outputfolder'], "filteredout/expmat.bed.excl")
            os.rename(os.path.join(args['outputfolder'], "expmat.bed.excl"), expmatbedexcl)

        if args["barcodemode"] == "singlecell":
            final_bc_file = bcfile
            i2bc = {}
            with open(final_bc_file, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    arr = line.split("\t")
                    if len(arr) < 3:
                        continue
                    bc_seq, i_str, i_str2 = arr[0], arr[1], arr[2]
                    i2bc[i_str] = bc_seq

            expmat_tsv = os.path.join(args["outputfolder"], "expmat.tsv")
            if os.path.exists(expmat_tsv):
                tmp_renamed = expmat_tsv + ".renamed"
                with open(expmat_tsv, "r") as fin, open(tmp_renamed, "w") as fout:
                    lines = fin.readlines()
                    if lines:
                        header = lines[0].rstrip("\n").split("\t")
                        new_header = []
                        for col in header:
                            if col == "gene":
                                new_header.append(col)
                            else:
                                part = col.split("x")
                                if len(part) == 2 and part[0] == part[1]:
                                    i_val = part[0]
                                    if i_val in i2bc:
                                        new_header.append(i2bc[i_val])
                                    else:
                                        new_header.append(col)
                                else:
                                    new_header.append(col)
                        fout.write("\t".join(new_header) + "\n")
                        for line in lines[1:]:
                            fout.write(line)
                os.remove(expmat_tsv)
                os.rename(tmp_renamed, expmat_tsv)

        else:
            if qualityfilter not in ["0:0", "NA"]:
                removeByDim = args.get("removeByDim", data.get("removeByDim", False))
                if removeByDim:
                    from .featurefilter import featurefilter
                    args['filterlogratio'] = (
                        args.get('filterlogratio') or data.get('filterlogratio') or 2
                    )
                    #gtffile = args.get('gtffile') or data.get('gtffile')
                    #featurefilter(args['gtffile'], args['options'], args['barcode_file'], args['filterlogratio'], args['outputfolder'])
                    featurefilter(
                        usedgtf,
                        args['options'],
                        args['barcode_file'],
                        args['filterlogratio'],
                        args['outputfolder']
                    )
                    
    if not args['not_organize_result']:
        interimfolder = os.path.join(args['outputfolder'], "interim/")
        os.makedirs(interimfolder, exist_ok=True)
        import gzip
        import shutil
        from shutil import which

        expmatbed = os.path.join(args['outputfolder'], "expmat.bed")
        combinefq = os.path.join(args['outputfolder'], "combine.fq")
        expmatbed2 = os.path.join(args['outputfolder'], "interim/expmat.bed.gz")
        combinefq2 = os.path.join(args['outputfolder'], "interim/combine.fq.gz")
        

        threads = args['threadnum']
        def compress_one(src, dst):
            if not os.path.exists(str(src)):
                return
            tmp = dst + ".tmp"
            if which("pigz"):
                
                with open(tmp, "wb") as out:
                    subprocess.run(["pigz", "-p", str(threads), "-c", str(src)], check=True, stdout=out)
            elif which("bgzip"):
                with open(tmp, "wb") as out:
                    subprocess.run(["bgzip", "-@", str(threads), "-c", str(src)], check=True, stdout=out)
            else:
                with open(src, "rb") as f_in, gzip.open(tmp, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)

            os.replace(tmp, dst)   
            os.remove(src)
        compress_one(expmatbed,  expmatbed2)
        compress_one(combinefq,  combinefq2)


        barcode_db = os.path.join(args['outputfolder'], "temps/barcode_db/")
        temps = os.path.join(args['outputfolder'], "temps/")
        if os.path.isdir(barcode_db):
            for file_name in os.listdir(barcode_db):
                file_path = os.path.join(barcode_db, file_name)
                if os.path.isfile(file_path):
                    os.remove(file_path)
            os.rmdir(barcode_db)
        if os.path.isdir(temps):
            for file_name in os.listdir(temps):
                file_path = os.path.join(temps, file_name)
                if os.path.isfile(file_path):
                    os.remove(file_path)
            os.rmdir(temps)
    
    

    with open(os.path.join(args['outputfolder'],"ASTRO.log.out"),"w") as logout:
        logout.write("input json information:\n")
        json.dump(args, logout, indent=4)

        demultiplexerlog = os.path.join(args['outputfolder'], ".logs/demultiplexing.log")
        genomemappiglog = os.path.join(args['outputfolder'], ".logs/genomemapping.log")
        countfeaturelog = os.path.join(args['outputfolder'], ".logs/countfeature.log")
        featurefilterlog = os.path.join(args['outputfolder'], ".logs/featurefilter.log")


        if os.path.isfile(demultiplexerlog):
            logout.write(f"\n###############################\n###############################\n###############################\n[STEP] demultiplexer step start\n###############################\n###############################\n###############################\n")
            with open(demultiplexerlog, "r") as f:
                content = f.read()
                logout.write(f"\n{content}\n")
        if os.path.isfile(genomemappiglog):
            logout.write(f"\n###############################\n###############################\n###############################\n[STEP] genomemappig step start\n###############################\n###############################\n###############################\n")
            with open(genomemappiglog, "r") as f:
                content = f.read()
                logout.write(f"\n{content}\n")
        if os.path.isfile(countfeaturelog):
            logout.write(f"\n###############################\n###############################\n###############################\n[STEP] countfeature step start\n###############################\n###############################\n###############################\n")
            with open(countfeaturelog, "r") as f:
                content = f.read()
                logout.write(f"\n{content}\n")
        if os.path.isfile(featurefilterlog):
            logout.write(f"\n###############################\n###############################\n###############################\n[STEP] featurefilter step start\n###############################\n###############################\n###############################\n")
            with open(featurefilterlog, "r") as f:
                content = f.read()
                logout.write(f"\n{content}\n")
        
