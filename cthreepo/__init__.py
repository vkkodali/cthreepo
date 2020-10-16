#!/usr/bin/env python3

import os
import csv
import sys
import argparse
from collections import namedtuple
from .fetch_assm_report import *
from .converters import *

csv.field_size_limit(5000000)

# See http://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

def processargs(args):
    Args = namedtuple('Args', ['fi', 'fo', 'mapfile', 'accn', 'id_from',
            'id_to', 'ku', 'p', 'conv_func', 'col'])

    mapfile_dict = {
                'h38': 'mapfiles/h38.map',
                'h37': 'mapfiles/h37.map',
                'm38': 'mapfiles/m38.map',
                'm37': 'mapfiles/m37.map'
                }
    id_dict = {
                'ens': 0,
                'gb': 4,
                'rs': 6,
                'uc': 9,
                'ensembl': 0,
                'genbank': 4,
                'refseq': 6,
                'ucsc': 9,
                }
    format_dict = {
                'gff3'     : convgxf,
                'gtf'      : convgxf,
                'vcf'      : convgxf,
                'bed'      : convbed,
                'bedgraph' : convbed,
                'tsv'      : convtsv,
                'sam'      : convsam,
                'wig'      : convwig
                }
    ## check in and out files
    if args.infile:
        fi = open(args.infile, 'rt')
    else:
        fi = sys.stdin

    if args.outfile:
        fo = open(args.outfile, 'wt')
    else:
        fo = sys.stdout

    ## check mapfile
    if args.mapfile:
        if args.mapfile in mapfile_dict:
            # Determine the directory of the current module
            this_dir, _ = os.path.split(__file__)
            # Derive relative directory
            mapfile = os.path.join(this_dir, mapfile_dict[args.mapfile])
        else:
            mapfile = args.mapfile
    else:
        mapfile = None

    ## check assm_acc
    accn = args.accn if args.accn else None    

    ## check id_from and id_to
    id_from = id_dict[args.id_from.lower().strip()]
    id_to = id_dict[args.id_to.lower().strip()]

    ku = 'T' if args.keep_unmapped else 'F'

    p = 'T' if args.primary else 'F'

    ## check format
    format = args.format.lower().strip()
    conv_func = format_dict[format]

    ## check column
    col = 'NA' # default value in case column not provided
    if args.format == 'tsv':
        if not args.column:
            print("ERROR: `col` required for `tsv` format", file = sys.stderr)
            sys.exit()
        else:
            col = args.column - 1 # user provides col in 1-based number

    ## return args
    return Args(fi, fo, mapfile, accn, id_from, id_to, ku, p, conv_func, col)

def create_maptbl(mapfile):
    with open(mapfile, 'r') as f:
        maptbl = f.readlines()
    return maptbl

def chrnamedict(maptbl, id_from, id_to, p):
    """
    create a mapping dict that will be used swap seq-ids in the input file
    ## params
    mapfile : path to NCBI assembly_report.txt format file
    id_from : column number of the seq-id format in input file
    id_to   : column number of the seq-id format to convert to
    ## returns
    chrmap  : a dict object with id_from:id_to
    """
    chrmap = {}
    tbl = csv.reader(maptbl, delimiter = '\t')
    if p == 'F' and id_from == 0:
        for line in tbl:
            if not line[0].startswith('#') and line[id_to] != 'na':
                chrmap[line[id_from]] = line[id_to]
                # to deal with ens using gb seq-ids in their GTF
                chrmap[line[4]] = line[id_to]
                # to deal with ens prepending CHR to their patches, etc
                chrmap['CHR_'+line[id_from]] = line[id_to]
    if p == 'F' and id_from != 0:
        for line in tbl:
            if not line[0].startswith('#') and line[id_to] != 'na':
                    chrmap[line[id_from]] = line[id_to]
    elif p == 'T':
        for line in tbl:
            if not line[0].startswith('#') and line[0] == '1':
                au = line[7]
                print(au, file=sys.stderr)
                chrmap[line[id_from]] = line[id_to]
                break
        for line in tbl:
            if not line[0].startswith('#') and line[7] == au:
                chrmap[line[id_from]] = line[id_to]
    return chrmap

def main():
    parser = argparse.ArgumentParser(description ="""This script parses input
                file and converts the seq-id name from one kind to the other""")
    parser.add_argument('-i', '--infile', help="input file")
    parser.add_argument('-o', '--outfile', help="output file")
    mapgroup = parser.add_mutually_exclusive_group(required=True)
    mapgroup.add_argument('-m', '--mapfile',
                        help = "NCBI style assembly_report file for mapping")
    mapgroup.add_argument('-a', '--accn',
                        help = "NCBI Assembly Accession with version")
    parser.add_argument('-if', '--id_from', default = 'uc',
                        choices = ['uc', 'rs', 'gb', 'ens'],
                        type = lambda s: s.lower().strip(),
                        help = "seq-id format in the input file; can be \
                            `ens`, `uc`, `gb`, or `rs`; default is `uc`")
    parser.add_argument('-it', '--id_to', default = 'rs', 
                        choices = ['uc', 'rs', 'gb', 'ens'],
                        type = lambda s: s.lower().strip(),
                        help = "seq-id format in the output file; can be \
                            `ens`, `uc`, `gb`, or `rs`; default is `rs`")
    parser.add_argument('-ku', '--keep_unmapped', action='store_true',
                        help = "keep lines that don't have seq-id matches")
    parser.add_argument('-p', '--primary', action='store_true',
                        help = "restrict to primary assembly only")
    parser.add_argument('-f', '--format', default = 'gff3', 
                        choices = ['gff3', 'gtf', 'bedgraph', 'bed',
                            'sam', 'vcf', 'wig', 'tsv'],
                        type = lambda s: s.lower().strip(),
                        help = "input file format; can be `gff3`, `gtf`, \
                            `bedgraph`, `bed`, `sam`, `vcf`, `wig` or `tsv`; \
                            default is `gff3`")
    parser.add_argument('-c', '--column', type = int,
                        help = "column where the seq-id is located; \
                            required for `tsv` format")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    A = processargs(args)
    
    if A.mapfile:
        maptbl = create_maptbl(A.mapfile)
    elif A.accn:
        esearch_result = perform_esearch(A.accn)
        esummary_result = perform_esummary(esearch_result)
        maptbl = fetch_asm_rpt(esummary_result)
    
    chrmap = chrnamedict(maptbl, A.id_from, A.id_to, A.p)

    if A.conv_func == convtsv:
        A.conv_func(A.fi, A.fo, chrmap, A.ku, A.col)
    else:
        A.conv_func(A.fi, A.fo, chrmap, A.ku)
