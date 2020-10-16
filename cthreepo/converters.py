import csv
import re 
import sys
import os 

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