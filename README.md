# cthreepo
[![PyPI version](https://badge.fury.io/py/cthreepo.svg)](https://badge.fury.io/py/cthreepo)
[![Conda](https://img.shields.io/conda/dn/bioconda/cthreepo?label=bioconda-install&style=flat)](https://anaconda.org/bioconda/cthreepo)

A python script to interconvert seq-ids in gff3, gtf, bed and other files.

---
## Quick start for the impatient
1. Install using conda 
```
conda install -c bioconda cthreepo 
```
2. Execute as follows:
```
## convert seq-ids in <input.gff3> from refseq format (NC_000001.11)
## to UCSC format (chr1) using the Human GRCh38 mapping dictionary
cthreepo -i <input.gff3> -if rs -it uc -f gff3 -m h38 -o <output.gff3>
```
---

## Introduction
NCBI RefSeq, UCSC and Ensembl use different identifiers for chromosomes in annotation and other files such as GFF3, GTF, etc. Users interested in using a mix of files downloaded from different sources and use them in a single pipeline may end up with seq-id mismatch related errors. This script converts seq-ids from one style to the other in order to make the files compatible with each other.

## Installation and Usage
Python3 is required for this script to work. With that requirement satisfied, you can install as shown below:
### Install using conda 
```
conda install -c bioconda cthreepo 
```
### Install using pip
```
pip install cthreepo
```
### Install from this repository
First, download/clone the repository. Then run: 
```
python3 setup.py install
```
### Usage
```
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
`cthreepo` needs a `mapfile` that it uses to figure out how seq-ids map from one style to the other. 
* Use the built-in shortcuts -- `h38`, `h37`, `m38` and `m37` for GRCh38/hg38, GRCh37/hg19, MGSCv37/mm9 and GRCm38/mm10 respectively. I try to keep these files up-to-date but if they don't work as expected, I suggest using the latest file by following one of the two options described below.
* Provide NCBI assembly accession using the `-a` parameter. A complete, legal accession.version such as GCF_000001405.39 should be provided. 
* Provide an NCBI assembly report file. For a given assembly it can be downloaded from the [NCBI Assembly](https://www.ncbi.nlm.nih.gov/assembly) website. If the 'Download' button is used, this file is called 'Assembly structure report'. On the [NCBI Genomes FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/) site, these files have the suffix `assembly_report.txt`. 

