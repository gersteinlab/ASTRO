#!/usr/bin/env python3

import json
import os
import sys
from .genomemapping import genomemapping
from .countfeature import countfeature
from .featurefilter import featurefilter
from .demultiplexer import demultiplexing, get_barcode_for_single_cell
from .validgene import getvalidedgtf_parallel


def ASTRO (**kwargs):

    args = kwargs
    
    json_file_path = args['json_file_path'] or args['json_file_path1']
    if json_file_path:
        with open(json_file_path, "r") as file:
            data = json.load(file)
    else:
        data = {}

    if data.get('barcode_read') and data.get('R1'):
        sys.exit("barcode_read and R1 are same things, duplicated settings")
    if data.get('barcode_read'):
        data['R1'] = data.pop('barcode_read')

    if data.get('transcript_read') and data.get('R2'):
        sys.exit("transcript_read and R2 are same things, duplicated settings")
    if data.get('transcript_read'):
        data['R2'] = data.pop('transcript_read')

    
    args['options'] = args.get('options') or data.get('options') or ""
    args['threadnum'] = args.get('threadnum') or data.get('threadnum') or 16
    args['steps'] = args.get('steps') or data.get('steps') or 7
    args['steps'] = int(args['steps'])
    args['outputfolder'] = args.get('outputfolder') or data.get('outputfolder') or sys.exit("outputfolder is not specified in both parser and json")
    args['STARparamfile4genome'] = args.get('STARparamfile4genome') or data.get('STARparamfile4genome') or 'NA'
    args['genes2check'] = args.get('genes2check') or data.get('genes2check') or False
    args['barcode_threshold'] = int(args.get('barcode_threshold') or data.get('barcode_threshold') or 100)
    args['barcodelength'] = int(args.get('barcodelength') or data.get('barcodelength') or 0) 
    args['ReadLayout']  = args.get('ReadLayout') or data.get('ReadLayout') or "singleend"
    args['limitOutSAMoneReadBytes4barcodeMapping']  = args.get('limitOutSAMoneReadBytes4barcodeMapping') or data.get('limitOutSAMoneReadBytes4barcodeMapping') or "NA"

    os.makedirs(args['outputfolder'], exist_ok=True)

    args['barcodemode'] = "singlecell" if args.get('singlecell') or data.get('singlecell') else "spatial"

    
    if args['steps'] & 1:
        args['R1'] = args.get('R1') or data.get('R1') or sys.exit("R1 is not specified in both parser and json")
        args['R2'] = args.get('R2') or data.get('R2') or sys.exit("R2 is not specified in both parser and json")
        args['PrimerStructure1'] = args.get('PrimerStructure1') or data.get('PrimerStructure1') or 'NA'
        
        args['StructureUMI'] = args.get('StructureUMI') or data.get('StructureUMI') or sys.exit("StructureUMI is not specified in both parser and json")
        args['StructureBarcode'] = args.get('StructureBarcode') or data.get('StructureBarcode') or sys.exit("StructureBarcode is not specified in both parser and json")

        if args['barcodemode'] == "singlecell":
            user_bc_file = args.get('barcode_file') or data.get('barcode_file') or "notavailable"

            final_bc_file = get_barcode_for_single_cell(read1=args['R1'],read2=args['R2'],barcode_file=user_bc_file,PrimerStructure1=args['PrimerStructure1'],StructureUMI=args['StructureUMI'],StructureBarcode=args['StructureBarcode'],
                threadnum=args['threadnum'],
                outputfolder=args['outputfolder'],
                barcode_threshold=args['barcode_threshold'],
                barcodelength=args['barcodelength']
            )
            if args['ReadLayout'] == "pairedend":
                from .demultiplexer2 import demultiplexingPair
                demultiplexingPair(
                    read1=args['R1'],read2=args['R2'], barcode_file=final_bc_file, 
                    PrimerStructure1=args['PrimerStructure1'],
                    StructureUMI=args['StructureUMI'],
                    StructureBarcode=args['StructureBarcode'],
                    threadnum=args['threadnum'], 
                    outputfolder=args['outputfolder'],
                    limitOutSAMoneReadBytes4barcodeMapping=args['limitOutSAMoneReadBytes4barcodeMapping']
                )
            else:
                demultiplexing(
                    read1=args['R1'],read2=args['R2'], barcode_file=final_bc_file, 
                    PrimerStructure1=args['PrimerStructure1'],
                    StructureUMI=args['StructureUMI'],
                    StructureBarcode=args['StructureBarcode'],
                    threadnum=args['threadnum'], 
                    outputfolder=args['outputfolder'],
                    limitOutSAMoneReadBytes4barcodeMapping=args['limitOutSAMoneReadBytes4barcodeMapping']
                )
            
           
            args['barcode_file'] = final_bc_file

        else:
            bcfile = args.get('barcode_file') or data.get('barcode_file') or sys.exit("barcode_file missing for spatial")
            if args['ReadLayout'] == "pairedend":
                from .demultiplexer2 import demultiplexingPair
                demultiplexingPair(
                    read1=args['R1'],read2=args['R2'], barcode_file=bcfile, 
                    PrimerStructure1=args['PrimerStructure1'],
                    StructureUMI=args['StructureUMI'],
                    StructureBarcode=args['StructureBarcode'],
                    threadnum=args['threadnum'], 
                    outputfolder=args['outputfolder'],
                    limitOutSAMoneReadBytes4barcodeMapping=args['limitOutSAMoneReadBytes4barcodeMapping']
                )
            else:
                demultiplexing(
                    read1=args['R1'],read2=args['R2'], barcode_file=bcfile, 
                    PrimerStructure1=args['PrimerStructure1'],
                    StructureUMI=args['StructureUMI'],
                    StructureBarcode=args['StructureBarcode'],
                    threadnum=args['threadnum'], 
                    outputfolder=args['outputfolder'],
                    limitOutSAMoneReadBytes4barcodeMapping=args['limitOutSAMoneReadBytes4barcodeMapping']
                )
            args['barcode_file'] = bcfile

    if args['steps'] & 2:
        args['starref'] = args.get('starref') or data.get('starref') or sys.exit("starref is not specified in both parser and json")
        args['gtffile'] = args.get('gtffile') or data.get('gtffile') or sys.exit("gtffile is not specified in both parser and json")
        genomemapping(args['starref'], args['gtffile'], args['threadnum'], args['options'], args['outputfolder'],args['STARparamfile4genome'])
    if args['steps'] & 4:
        from .countfeature import countfeature
        bcfile = args.get('barcode_file') or data.get('barcode_file') or sys.exit("barcode_file missing in Step4")
        gtffile = args.get('gtffile') or data.get('gtffile') or sys.exit("gtffile missing in Step4")
        qualityfilter = args.get('qualityfilter') or data.get('qualityfilter') or "25:0.75"
        usedgtf = args['gtffile']
        if args['genes2check']:
            usedgtf = getvalidedgtf_parallel(gtfin=gtffile,outputfolder=args['outputfolder'],genes2check=args['genes2check'],hangout=5,threadsnum=args['threadnum'])
        countfeature(usedgtf, args['threadnum'], args['options'], bcfile, args['outputfolder'], qualityfilter)
        
        if args['barcodemode'] == "singlecell":
            final_bc_file = bcfile  
            i2bc = {}
            with open(final_bc_file, 'r') as f:
                for line in f:
                    line=line.strip()
                    if not line:
                        continue
                    arr=line.split('\t')
                    if len(arr)<3:
                        continue
                    bc_seq, i_str, i_str2 = arr[0], arr[1], arr[2]
                    i2bc[i_str] = bc_seq
            
            expmat_tsv = os.path.join(args['outputfolder'], "expmat.tsv")
            if os.path.exists(expmat_tsv):
                tmp_renamed = expmat_tsv + ".renamed"
                with open(expmat_tsv, 'r') as fin, open(tmp_renamed, 'w') as fout:
                    lines = fin.readlines()
                    if lines:
                        header = lines[0].rstrip('\n').split('\t')
                        new_header = []
                        for col in header:
                            if col == "gene":
                                new_header.append(col)
                            else:
                                part = col.split('x')
                                if len(part)==2 and part[0]==part[1]:
                                    i_val = part[0]
                                    if i_val in i2bc:
                                        new_header.append(i2bc[i_val])
                                    else:
                                        new_header.append(col)
                                else:
                                    new_header.append(col)
                        fout.write('\t'.join(new_header) + '\n')
                        for line in lines[1:]:
                            fout.write(line)
                os.remove(expmat_tsv)
                os.rename(tmp_renamed, expmat_tsv)

        else:
            if qualityfilter not in ['0:0','NA']:
                removeByDim = args.get('removeByDim') or data.get('removeByDim') or True
                if removeByDim:
                    from .featurefilter import featurefilter
                    args['filterlogratio'] = args.get('filterlogratio') or data.get('filterlogratio') or 2
                    gtffile = args.get('gtffile') or data.get('gtffile')
                    featurefilter(args['gtffile'], args['options'], args['barcode_file'], args['filterlogratio'], args['outputfolder'])

    with open(os.path.join(args['outputfolder'],"input.json"),"w") as jf:
        json.dump(args, jf, indent=4)
