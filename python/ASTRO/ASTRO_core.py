#!/usr/bin/env python3

import json
import os
import sys
from .demultiplexer import demultiplexing
from .genomemapping import genomemapping
from .countfeature import countfeature
from .featurefilter import featurefilter




def ASTRO (**kwargs):

    args = kwargs
    
    json_file_path = args['json_file_path'] or args['json_file_path1']
    if json_file_path:
        with open(json_file_path, "r") as file:
            data = json.load(file)
    else:
        data = {}

            
    
    
    args['options'] = args.get('options') or data.get('options') or ""
    args['threadnum'] = args.get('threadnum') or data.get('threadnum') or 16
    args['steps'] = args.get('steps') or data.get('steps') or 7
    args['steps'] = int(args['steps'])
    args['outputfolder'] = args.get('outputfolder') or data.get('outputfolder') or sys.exit("outputfolder is not specified in both parser and json")
    args['outputfolder'] = args.get('outputfolder') or data.get('outputfolder') or sys.exit("outputfolder is not specified in both parser and json")
    args['STARparamfile4genome'] = args.get('STARparamfile4genome') or data.get('STARparamfile4genome') or 'NA'
    os.makedirs(args['outputfolder'], exist_ok=True)

    
    
    if args['steps'] & 1:
        args['R1'] = args.get('R1') or data.get('R1') or sys.exit("R1 is not specified in both parser and json")
        args['R2'] = args.get('R2') or data.get('R2') or sys.exit("R2 is not specified in both parser and json")
        args['barcode_file'] = args.get('barcode_file') or data.get('barcode_file') or sys.exit("barcode_file is not specified in both parser and json")
        args['PrimerStructure1'] = args.get('PrimerStructure1') or data.get('PrimerStructure1') or sys.exit("PrimerStructure1 is not specified in both parser and json")
        args['StructureUMI'] = args.get('StructureUMI') or data.get('StructureUMI') or sys.exit("StructureUMI is not specified in both parser and json")
        args['StructureBarcode'] = args.get('StructureBarcode') or data.get('StructureBarcode') or sys.exit("StructureBarcode is not specified in both parser and json")
        demultiplexing(args['R1'], args['R2'], args['barcode_file'], args['PrimerStructure1'], args['StructureUMI'], args['StructureBarcode'], args['threadnum'], args['outputfolder'])
    if args['steps'] & 2:
        args['starref'] = args.get('starref') or data.get('starref') or sys.exit("starref is not specified in both parser and json")
        args['gtffile'] = args.get('gtffile') or data.get('gtffile') or sys.exit("gtffile is not specified in both parser and json")
        genomemapping(args['starref'], args['gtffile'], args['threadnum'], args['options'], args['outputfolder'],args['STARparamfile4genome'])
    if args['steps'] & 4:
        args['gtffile'] = args.get('gtffile') or data.get('gtffile') or sys.exit("gtffile is not specified in both parser and json")
        args['barcode_file'] = args.get('barcode_file') or data.get('barcode_file') or sys.exit("barcode_file is not specified in both parser and json")
        args['qualityfilter'] = args.get('qualityfilter') or data.get('qualityfilter') or '25:0.75'
        args['genes2check'] = args.get('genes2check') or data.get('genes2check') or False
        usedgtf = args['gtffile']
        if args['genes2check']:
            from .validgene import getvalidedgtf_parallel
            usedgtf = getvalidedgtf_parallel(usedgtf, args['outputfolder'], args['genes2check'], args['threadnum'])
        countfeature(usedgtf, args['threadnum'], args['options'], args['barcode_file'], args['outputfolder'], args['qualityfilter'])
        if args['qualityfilter'] != '0:0' or args['qualityfilter'] != 'NA':
            args['removeByDim'] = args.get('removeByDim') or data.get('removeByDim') or True
            args['filterlogratio'] = args.get('filterlogratio') or data.get('filterlogratio') or 2
            if args['removeByDim']:
                featurefilter(usedgtf, args['options'], args['barcode_file'], args['filterlogratio'], args['outputfolder'])

    with open(args['outputfolder'] + "/input.json", "w") as json_file:
        json.dump(args, json_file, indent=4)
    
