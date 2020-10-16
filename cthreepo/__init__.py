#!/usr/bin/env python3

import os
import re
import csv
import sys
import argparse
import collections
from .fetch_assm_report import *

csv.field_size_limit(5000000)

# See http://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

def processargs(args):
    Args = collections.namedtuple('Args', ['fi', 'fo', 'mapfile', 'accn', 'id_from',
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

def convgxf(fi, fo, chrmap, ku):
    """
    converts seq-ids in the `infile` to the desired format and writes output
    to `outfile`
    ## params
    infile  : input file, gff3 or gtf format
    outfile : output file with swapped seq-ids; same format as input
    chrmap  : dict object with id_from:id_to
    ## returns
    unmapped : list of seq-ids that were unmapped
    """
    tblin = csv.reader( fi,
                        delimiter = '\t'
                        )
    tblout = csv.writer(
                        fo,
                        delimiter = '\t',
                        # doublequote = False,
                        # quoting = csv.QUOTE_NONE,
                        quotechar = "\xb6",
                        # escapechar = '\\',
                        lineterminator=os.linesep
                        )
    all_lines = 0
    um_lines = 0
    um_acc = set()
    for line in tblin:
        if not line[0].startswith("#"):
            all_lines = all_lines + 1
            if line[0] in chrmap:
                chrom = chrmap[line[0]]
                newline = [chrom] + line[1:]
                tblout.writerow(newline)
            elif ku == 'T':
                um_lines = um_lines + 1
                um_acc.add(line[0])
                tblout.writerow(line)
            else:
                um_lines = um_lines + 1
                um_acc.add(line[0])
        elif line[0].startswith('##sequence-region'):
            line = line[0].split(' ')
            if line[1] in chrmap:
                chrom = chrmap[line[1]]
                newline = ' '.join([line[0]] + [chrom] + line[2:])
                tblout.writerow([newline])
            elif ku == 'T':
                um_lines = um_lines + 1
                um_acc.add(line[1])
                tblout.writerow(line)
            else:
                um_lines = um_lines + 1
                um_acc.add(line[1])
        else:
            tblout.writerow(line)
    if len(um_acc) > 0 and ku == 'F':
        print(
        "WARNING: {} accessions were not present in the mapfile; they are "
        "dropped in the output file. {} of {} lines were dropped. "
        "Use `-ku` option to keep them instead."
        .format(len(um_acc), um_lines, all_lines),
        file = sys.stderr)
    fi.close()
    fo.close()

def convbed(fi, fo, chrmap, ku):
    """
    converts seq-ids in the `infile` to the desired format and writes output
    to `outfile`
    ## params
    infile  : input file, bed, bigbed or bedgraph formats
    outfile : output file with swapped seq-ids; same format as input
    chrmap  : dict object with id_from:id_to
    ## returns
    unmapped : list of seq-ids that were unmapped
    """
    tblin = csv.reader(fi, delimiter = '\t')
    tblout = csv.writer(
                        fo,
                        delimiter = '\t',
                        quotechar = "\xb6",
                        lineterminator=os.linesep
                        )
    all_lines = 0
    um_lines = 0
    um_acc = set()
    for line in tblin:
        if not (
            line[0].startswith('#') or
            line[0].startswith('track') or
            line[0].startswith('browser')
            ):
            all_lines = all_lines + 1
            if len(line) == 1:
                line = re.sub(' +',' ', line[0])
                line = [item for item in line.split()]
            if line[0] in chrmap:
                chrom = chrmap[line[0]]
                newline = [chrom] + line[1:]
                tblout.writerow(newline)
            elif ku == 'T':
                um_lines = um_lines + 1
                um_acc.add(line[0])
                tblout.writerow(line)
            else:
                um_lines = um_lines + 1
                um_acc.add(line[0])
        elif 'browser position' in line[0]:
            line = re.sub(' +',' ', line[0])
            line = [item for item in line.split()]
            f_seqid = line[2].split(':')[0]
            t_seqid = chrmap[line[2].split(':')[0]]
            line[2] = re.sub(f_seqid, t_seqid, line[2])
            line = [' '.join(line)]
            tblout.writerow(line)
        else:
            tblout.writerow(line)
    if len(um_acc) > 0 and ku == 'F':
        print(
        "WARNING: {} accessions were not present in the mapfile; they are "
        "dropped in the output file. {} of {} lines were dropped. "
        "Use `-ku` option to keep them instead."
        .format(len(um_acc), um_lines, all_lines),
        file = sys.stderr)
    fi.close()
    fo.close()

def convwig(fi, fo, chrmap, ku):
    """
    Dropping lines with unmapped seq-ids does not make sense for wig files
    because it essentially breaks the entire wig file. I am keeping it for
    now, with the intention to remove it at a later point.
    """
    all_lines = 0
    um_lines = 0
    um_acc = set()
    for line in fi:
        if (
            line.startswith('variableStep') or
            line.startswith('fixedStep')
            ):
            all_lines = all_lines + 1
            line = line.split()
            f_seqid = line[1].split('=')[1]
            if f_seqid in chrmap:
                chrom = chrmap[f_seqid]
                line[1] = re.sub(f_seqid, chrom, line[1])
                fo.write(' '.join(line)+'\n')
            elif ku == 'T':
                um_lines = um_lines + 1
                um_acc.add(f_seqid)
                fo.write(' '.join(line)+'\n')
            else:
                um_lines = um_lines + 1
                um_acc.add(f_seqid)
        else:
            fo.write(line)
    if len(um_acc) > 0 and ku == 'F':
        print(
        "WARNING: {} accessions were not present in the mapfile; they are "
        "dropped in the output file. {} of {} lines were dropped. "
        "Use `-ku` option to keep them instead."
        .format(len(um_acc), um_lines, all_lines),
        file = sys.stderr)
    fi.close()
    fo.close()

def convsam(fi, fo, chrmap, ku):
    tblin = csv.reader(fi, delimiter = '\t')
    tblout = csv.writer(
                        fo,
                        delimiter = '\t',
                        quotechar = "\xb6",
                        lineterminator=os.linesep
                        )
    all_lines = 0
    um_lines = 0
    um_acc = set()
    for line in tblin:
        if not line[0].startswith('@'):
            all_lines = all_lines + 1
            if line[2] in chrmap:
                chrom = chrmap[line[2]]
                newline = line[:2] + [chrom] + line[3:]
                tblout.writerow(newline)
            elif ku == 'T':
                um_lines = um_lines + 1
                um_acc.add(line[2])
                tblout.writerow(line)
            else:
                um_lines = um_lines + 1
                um_acc.add(line[2])
        elif line[0] == '@SQ':
            f_seqid = line[1].split(':')[1]
            if f_seqid in chrmap:
                chrom = chrmap[f_seqid]
                newline = [line[0]] + ['SN:' + chrom] + line[2:]
                tblout.writerow(newline)
            elif ku == 'T':
                um_acc.add(f_seqid)
                tblout.writerow(line)
            else:
                um_acc.add(f_seqid)
        else:
            tblout.writerow(line)
    if len(um_acc) > 0 and ku == 'F':
        print(
        "WARNING: {} accessions were not present in the mapfile; they are "
        "dropped in the output file. {} of {} lines were dropped. "
        "Use `-ku` option to keep them instead."
        .format(len(um_acc), um_lines, all_lines),
        file = sys.stderr)
    fi.close()
    fo.close()

def convtsv(fi, fo, chrmap, ku, col):
    """
    converts seq-ids in the `infile` to the desired format and writes output
    to `outfile`
    ## params
    infile  : input file, some kind of tab-delimited format
    outfile : output file with swapped seq-ids; same format as input
    col     : column number where the seq-id is located (1-based)
    chrmap  : dict object with id_from:id_to
    ## returns
    unmapped : list of seq-ids that were unmapped
    """
    tblin = csv.reader(fi, delimiter = '\t')
    tblout = csv.writer(
                        fo,
                        delimiter = '\t',
                        quotechar = "\xb6",
                        lineterminator=os.linesep
                        )
    all_lines = 0
    um_lines = 0
    um_acc = set()
    for line in tblin:
        if not (
            line[0].startswith('#') or
            line[0].startswith('track') or
            line[0].startswith('browser')
            ):
            all_lines = all_lines + 1
            if len(line) == 1:
                line = re.sub(' +',' ', line[0])
                line = [item for item in line.split()]
            if line[col] in chrmap:
                chrom = chrmap[line[col]]
                newline = line[:col] + [chrom] + line[col+1:]
                tblout.writerow(newline)
            elif ku == 'T':
                um_lines = um_lines + 1
                um_acc.add(line[col])
                tblout.writerow(line)
            else:
                um_lines = um_lines + 1
                um_acc.add(line[col])
        elif 'browser position' in line[0]:
            line = re.sub(' +',' ', line[0])
            line = [item for item in line.split()]
            f_seqid = line[2].split(':')[0]
            t_seqid = chrmap[line[2].split(':')[0]]
            line[2] = re.sub(f_seqid, t_seqid, line[2])
            line = [' '.join(line)]
            tblout.writerow(line)
        else:
            tblout.writerow(line)
    if len(um_acc) > 0 and ku == 'F':
        print(
        "WARNING: {} accessions were not present in the mapfile; they are "
        "dropped in the output file. {} of {} lines were dropped. "
        "Use `-ku` option to keep them instead."
        .format(len(um_acc), um_lines, all_lines),
        file = sys.stderr)
    fi.close()
    fo.close()


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
