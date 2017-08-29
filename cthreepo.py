#!/home/kodalivk/envs/vkenv3/bin/python

import csv
import sys
import argparse
import collections
# import multiprocessing as mp

# See http://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

def processargs(args):
    Args = collections.namedtuple('Args', ['fi', 'fo', 'mapfile', 'id_from',
            'id_to', 'ku', 'p'])

    mapfile_dict = {
                'h38' : 'mapfiles/h38.map',
                'h37' : 'mapfiles/h37.map',
                'm38' : 'mapfiles/m38.map',
                'm37' : 'mapfiles/m37.map'
                }
    id_dict = {
                'ens': 0,
                'gb': 4,
                'rs': 6,
                'uc': 9
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
        mapfile = mapfile_dict[args.mapfile]
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
    ## return args
    return Args(fi, fo, mapfile, id_from, id_to, ku, p)

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
        if p == 'F':
            for line in tbl:
                if not line[0].startswith('#'):
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

def convgff3(fi, fo, chrmap, ku):
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
    tblout = csv.writer(fo, delimiter = '\t', quotechar = "'" , escapechar = '\\')
    all_lines = 0
    um_lines = 0
    um_acc = set()
    for line in tblin:
        if not line[0].startswith("#"):
            all_lines = all_lines + 1
            if line[0] in chrmap:
                chrom = chrmap[line[0]]
                newline = [chrom] + line[1:]
                x = tblout.writerow(newline)
            elif ku == 'T':
                um_lines = um_lines + 1
                um_acc.add(line[0])
                x = tblout.writerow(line)
            else:
                um_lines = um_lines + 1
                um_acc.add(line[0])
        else:
            x = tblout.writerow(line)
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
    parser = argparse.ArgumentParser(description ="""This script parses gtf/gff3
                file and converts the seq-id name from one kind to the other""")
    parser.add_argument('-i', '--infile', help="input gff3 file")
    parser.add_argument('-o', '--outfile', help="output gff3 file")
    parser.add_argument('-m', '--mapfile',
                        help = "NCBI style assembly_report file for mapping")
    parser.add_argument('-if', '--id_from', default = 'uc', help = "seq-id \
                        format in the input gff3 file; default is `uc`")
    parser.add_argument('-it', '--id_to', default = 'rs', help = "seq-id \
                        format in the output gff3 file; default is `rs`")
    parser.add_argument('-ku', '--keep_unmapped', action='store_true',
                        help = "keep lines that don't have seq-id matches")
    parser.add_argument('-p', '--primary', action='store_true',
                        help = "restrict to primary assembly only")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    A = processargs(args)
    chrmap = chrnamedict(A.mapfile, A.id_from, A.id_to, A.p)
    convgff3(A.fi, A.fo, chrmap, A.ku)
