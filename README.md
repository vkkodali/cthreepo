# cthreepo
A python script to interconvert seq-ids in gff3, gtf, bed and other files.
---
## Quick start for the impatient
1. Clone the repository
2. Run the following to install: 
```
python3 setup.py install
```
3. Execute as follows:
```
## convert seq-ids in <input.gff3> from refseq format (NC_000001.11)
## to UCSC format (chr1) using the Human GRCh38 mapping dictionary
cthreepo -i <input.gff3> -if rs -it uc -f gff3 -m h38 -o <output.gff3>
```
---
## Introduction
NCBI RefSeq, UCSC and Ensembl use different identifiers for chromosomes in annotation and other files such as GFF3, GTF, etc. Users interested in using a mix of files downloaded from different sources and use them in a single pipeline may end up with seq-id mismatch related errors. This script converts seq-ids from one style to the other in order to make the files compatible with each other.

## Installation and Usage
Python3 is required for this script to work. With that requirement satisfied, download/clone the repository, install and run the script `cthreepo.py` as shown below.
```
## installation
python3 setup.py install

## help
cthreepo --help 

## usage
## convert seq-ids in <input.gff3> from refseq format (NC_000001.11)
## to UCSC format (chr1) using the Human GRCh38 mapping dictionary
cthreepo \
    --infile <input.gff3> \
    --id_from rs \
    --id_to uc \
    --format gff3 \
    --mapfile h38 \
    --outfile <output.gff3>
```

## File formats supported
1. GFF3 (default)
2. GTF
3. BedGraph
4. BED
5. SAM
6. VCF
7. WIG
8. TSV

## Mapping files
`cthreepo` expects a `mapfile` that it uses to figure out how seq-ids map from one style to the other. For human and mouse assemblies, one can use the built-in shortcuts but for all other organisms, an NCBI assembly report file needs to be provided. For a given assembly it can be downloaded from the [NCBI Assembly](https://www.ncbi.nlm.nih.gov/assembly) website. If the 'Download' button is used, this file is called 'Assembly structure report'. On the NCBI Genomes FTP site, these files have the suffix `assembly_report.txt`. 
