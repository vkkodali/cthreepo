# cthreepo
A python script to interconvert seq-ids in gff3, gtf, bed and other files.
---
## Quick start for the impatient
1. Clone the repository
2. Execute as follows:
```
## convert seq-ids in <input.gff3> from refseq format (NC_000001.11)
## to UCSC format (chr1) using the Human GRCh38 mapping dictionary
cthreepo.py \
    --infile <input.gff3> \
    --id_from rs \
    --id_to uc \
    --format gff3 \
    --mapfile h38 \
    --outfile <output.gff3>
```
---
## Introduction
NCBI RefSeq, UCSC and Ensembl use different identifiers for chromosomes in annotation and other files such as GFF3, GTF, etc. Users interested in using a mix of files downloaded from different sources and use them in a single pipeline may end up with seq-id mismatch related errors. This script converts seq-ids from one style to the other in order to make the files compatible with each other.

## File formats supported
1. GFF3 (default)
2. GTF
3. BedGraph
4. BED
5. SAM
6. VCF
7. WIG

## Installation and Usage
Python3 is required for this script to work. With that requirement satisfied, download/clone the repository and run the script `cthreepo.py` as shown in the example above.
