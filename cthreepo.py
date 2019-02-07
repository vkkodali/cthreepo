#!/usr/bin/env python

import os
import re
import csv
import sys
import argparse
import collections
# import multiprocessing as mp

csv.field_size_limit(5000000)

# See http://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

def processargs(args):
    Args = collections.namedtuple('Args', ['fi', 'fo', 'mapfile', 'id_from',
            'id_to', 'ku', 'p', 'conv_func'])

    mapfile_dict = {
                'h38' : '/mapfiles/h38.map',
                'h37' : '/mapfiles/h37.map',
                'm38' : '/mapfiles/m38.map',
                'm37' : '/mapfiles/m37.map'
                }
    id_dict = {
                'ens': 0,
                'gb': 4,
                'rs': 6,
                'uc': 9
                }
    format_dict = {
                'gff3'     : convgxf,
                'gtf'      : convgxf,
                'vcf'      : convgxf,
                'bed'      : convbed,
                'bedgraph' : convbed,
                # 'psl'      : convpsl,
                'sam'      : convsam,
                'wig'      : convwig
                }
    ## check in and out files
    if args.infile:
        fi = open(args.infile, 'r')
    else:
        fi = sys.stdin

    if args.outfile:
        fo = open(args.outfile, 'w')
    else:
        fo = sys.stdout

    ## check mapfile
    if not args.mapfile:
        print(
        "ERROR: mapfile required. Can be an NCBI assembly_report file or "
        "one of `h37`, `h38`, `m37` and `m38` for preloaded lists"
        )
        sys.exit()
    elif args.mapfile in mapfile_dict:
        mapfile = os.path.abspath(os.path.dirname(sys.argv[0])) + mapfile_dict[args.mapfile]
    else:
        mapfile = args.mapfile

    ## check id_from and id_to
    if args.id_from not in id_dict or args.id_to not in id_dict:
        print(
        "ERROR: id_from and id_to can only be one of the following:"
        "`ens`, `gb`, `rs` or `uc`")
        sys.exit()
    else:
        id_from = id_dict[args.id_from]
        id_to = id_dict[args.id_to]

    if args.keep_unmapped:
        ku = 'T'
    else:
        ku = 'F'

    if args.primary:
        p = 'T'
    else:
        p = 'F'

    ## check format
    if args.format not in format_dict:
        print(
            "ERROR: invalid format. Choose from: {}"
            .format(list(format_dict.keys()))
            )
        sys.exit()
    else:
        conv_func = format_dict[args.format ]

    ## return args
    return Args(fi, fo, mapfile, id_from, id_to, ku, p, conv_func)

def chrnamedict(mapfile, id_from, id_to, p):
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
    with open(mapfile, 'r')  as f:
        tbl = csv.reader(f, delimiter = '\t')
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
    tblin = csv.reader(fi, delimiter = '\t')
    tblout = csv.writer(
                        fo,
                        delimiter = '\t',
                        quotechar = "`" ,
                        escapechar = '\\',
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
        .format(len(um_acc), um_lines, all_lines)
        )
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
        .format(len(um_acc), um_lines, all_lines)
        )
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
        .format(len(um_acc), um_lines, all_lines)
        )
    fi.close()
    fo.close()

def convsam(fi, fo, chrmap, ku):
    tblin = csv.reader(fi, delimiter = '\t')
    tblout = csv.writer(
                        fo,
                        delimiter = '\t',
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
        .format(len(um_acc), um_lines, all_lines)
        )
    fi.close()
    fo.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ="""This script parses input
                file and converts the seq-id name from one kind to the other""")
    parser.add_argument('-i', '--infile', help="input file")
    parser.add_argument('-o', '--outfile', help="output file")
    parser.add_argument('-m', '--mapfile',
                        help = "NCBI style assembly_report file for mapping")
    parser.add_argument('-if', '--id_from', default = 'uc', help = "seq-id \
                        format in the input file; can be `ens`, `uc`, \
                        `gb`, or `rs`; default is `uc`")
    parser.add_argument('-it', '--id_to', default = 'rs', help = "seq-id \
                        format in the output file; can be `ens`, `uc`, \
                        `gb`, or `rs`; default is `rs`")
    parser.add_argument('-ku', '--keep_unmapped', action='store_true',
                        help = "keep lines that don't have seq-id matches")
    parser.add_argument('-p', '--primary', action='store_true',
                        help = "restrict to primary assembly only")
    parser.add_argument('-f', '--format', default = 'gff3', help = "input \
                        file format; can be `gff3`, `gtf`, `bedgraph` \
                        `bed`, `sam`, `vcf` or `wig`; default is `gff3`")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    A = processargs(args)
    chrmap = chrnamedict(A.mapfile, A.id_from, A.id_to, A.p)
    A.conv_func(A.fi, A.fo, chrmap, A.ku)
