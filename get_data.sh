#!/bin/bash

set -e
set -u
set -o
set -x

mkdir -p data/

# annotation
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gff3.gz > data/annotations.gff3.gz

# variant
curl https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > data/variants.vcf.gz

# reference
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz > data/references.fasta.gz
