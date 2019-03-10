# cthreepo
A python script to interconvert seq-ids in gff3, gtf, bed and other files. 

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
