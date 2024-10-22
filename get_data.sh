#!/bin/bash

set -e
set -u
set -o
set -x

mkdir -p data/

# annotation
curl https://ftp.ensembl.org/pub/release-92/gff3/homo_sapiens/Homo_sapiens.GRCh38.92.gff3.gz > data/annotations.gff3.gz

# variant
curl https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz > data/variants.vcf.gz

# reference
curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz > data/references.fasta.gz
