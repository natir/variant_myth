#!/usr/bin/env python

# std import
import argparse
import copy
import csv
import gzip
import logging
import os
import pathlib
import random
import re
import string
import sys
import typing
import urllib.request

from collections import defaultdict

# 3rd party import

# project import


GFF_ID_REGEX = re.compile(r"ID=(?P<id>[^;]+)")
GFF_PARENT_REGEX = re.compile(r"Parent=(?P<parent>[^;]+)")


def is_gz(filepath: pathlib.Path):
    """Return true if first two bit match gzip magic number."""
    with open(filepath, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def open_file(filepath: pathlib.Path) -> typing.IO:
    """Open file if it's gzip or not."""
    if is_gz(filepath):
        return gzip.open(filepath, mode="rt")
    else:
        return open(filepath, mode="rt")


def random_seq(length: int) -> str:
    """Generate a random sequence."""
    return "".join(random.choices("actg", k=length))


def edit_record(record: list[str], chrom: str, new_start: int) -> list[str]:
    """Edit gff record."""
    record = copy.deepcopy(record)
    prev_begin = int(record[3])
    prev_end = int(record[4])
    new_end = new_start + prev_end - prev_begin

    record[0] = chrom
    record[3] = str(new_start)
    record[4] = str(new_end)

    return record


def download() -> int:
    """Download  example input file."""
    os.makedirs("data", exist_ok=True)

    logging.info("Start download annotations")
    urllib.request.urlretrieve(
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gff3.gz",
        "data/annotations.gff3.gz",
    )
    logging.info("End download annotations")

    logging.info("Start download variants")
    urllib.request.urlretrieve(
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        "data/variants.vcf.gz",
    )
    logging.info("End download variants")

    logging.info("Start download sequence")
    urllib.request.urlretrieve(
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz",
        "data/references.fasta.gz",
    )
    logging.info("End download sequence")

    return 0


def subsample(
    annotations_input: pathlib.Path,
    references_input: pathlib.Path,
    annotations_output: pathlib.Path,
    references_output: pathlib.Path,
    variants_output: pathlib.Path,
    seed: int,
) -> int:
    """Extract some information in input to generate a fake genome."""
    random.seed(seed)

    ###############
    # Set parameter
    ###############
    logging.info("Start generate target value")
    nb_gene = 20
    nb_transcript_range = (1, 5)
    gene2nb_transcript = [random.randint(*nb_transcript_range) for _ in range(nb_gene)]

    nb_chromosome = 3
    chr_name = [f"chr{a}" for a in string.ascii_uppercase[:nb_chromosome]]
    gene_dist_range = (50, 500)

    nb_snv = 20
    nb_mnv = 20
    mnv_range = (2, 5)
    nb_struct = 10
    struct_type = ["INS", "DEL", "DUP", "INV", "CNV"]
    struct_len = (50, 1000)
    logging.info("End generate target value")

    #################
    # Read annotation
    #################
    logging.info("Start read annotation input")
    id2record: dict[str, list[str]] = {}
    feature2id: dict[str, list[str]] = defaultdict(list)
    parent2child: dict[str: list[list[str]]] = defaultdict(list)
    with open_file(annotations_input) as fh:
        reader = csv.reader(filter(lambda row: row[0] != "#", fh), delimiter="\t")
        for record in reader:
            if id_match := GFF_ID_REGEX.search(record[8]):
                id2record[id_match.groupdict()["id"]] = record

                feature2id[record[2]].append(id_match.groupdict()["id"])
                if parent_match := GFF_PARENT_REGEX.search(record[8]):
                    parent2child[parent_match.groupdict()["parent"]].append(
                        id_match.groupdict()["id"]
                    )

    logging.info("End read annotation input")

    ###############
    # Read sequence
    ###############
    logging.info("Start read reference input")
    chr2seq_ = defaultdict(list)
    with open_file(references_input) as fh:
        for index, line in enumerate(fh):
            if line.startswith(">"):
                current_chr = line.strip().split()[0][1:]
            else:
                chr2seq_[current_chr].append(line.strip())

    chr2seq_in = {k: "".join(v) for k, v in chr2seq_.items()}
    logging.info("End read reference input")

    logging.info("Start extract and edit annotation")
    ####################
    # Extract annotation
    ####################
    random.shuffle(feature2id["gene"])

    nb_transcript2gene_id = defaultdict(list)
    for gene_id in feature2id["gene"]:
        nb_transcript2gene_id[len(parent2child[gene_id])].append(gene_id)

    select_gene_id = [
        random.choices(nb_transcript2gene_id[nb_transcript])[0]
        for nb_transcript in gene2nb_transcript
    ]

    ######################################
    # Edit annotation and extract sequence
    ######################################
    chr2position = {name: 1 for name in chr_name}
    chr2seq_ = {name: [] for name in chr_name}

    select_record = []
    for gene_id in select_gene_id:
        # Get random value
        choose_chr = random.choice(chr_name)
        random_len = random.randint(*gene_dist_range)
        start_position = chr2position[choose_chr] + random_len
        chr2seq_[choose_chr].append(random_seq(random_len))

        # Get record
        record = id2record[gene_id]
        chr_origin = record[0]
        prev_begin = int(record[3])
        prev_end = int(record[4])

        # Edit record and child
        select_record.append(edit_record(record, choose_chr, start_position))
        for transcript_id in parent2child[gene_id]:
            transcript_record = id2record[transcript_id]
            new_start = start_position + int(transcript_record[3]) - prev_begin
            select_record.append(edit_record(transcript_record, choose_chr, new_start))

            for child_id in parent2child[transcript_id]:
                child_record = id2record[child_id]
                if child_record[2] == "CDS":
                    print(child_record)
                new_start = start_position + int(child_record[3]) - prev_begin
                select_record.append(edit_record(child_record, choose_chr, new_start))

        chr2position[choose_chr] = int(record[4])
        chr2seq_[choose_chr].append(chr2seq_in[chr_origin][prev_begin - 1 : prev_end])

    chr2seq_out = {k: "".join(v) for k, v in chr2seq_.items()}
    logging.info("End extract and edit annotation")

    ##################
    # Generate variant
    ##################
    logging.info("Start generate variant")
    variant_out = []
    variant_out.append("##fileformat=VCFv4.3")
    for name, seq in chr2seq_out.items():
        variant_out.append(f"##contig=<ID={name},length={len(seq)}>")
    variant_out.append(
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">'
    )
    variant_out.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    ####################################
    # Generate Single Nucleotive Variant
    ####################################
    for _ in range(nb_snv):
        chrom = random.choice(chr_name)
        pos = random.randint(1, len(chr2seq_out[chrom]))
        ref = chr2seq_out[chrom][pos]
        alt = random.choice([nuc for nuc in "actg" if nuc.upper() != ref.upper()])

        variant_out.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t99\tPASS\t.")

    ####################################
    # Generate Multiple Nucleotive Variant
    ####################################
    for _ in range(nb_mnv):
        chrom = random.choice(chr_name)
        pos = random.randint(1, len(chr2seq_out[chrom]))
        length = random.randint(mnv_range[0], mnv_range[1])

        if random.choice(["ins", "del"]) == "ins":
            ref = chr2seq_out[chrom][pos]
            alt = ref + random_seq(length)
        else:
            max_range = min(pos + length, len(chr2seq_out[chrom]))
            ref = chr2seq_out[chrom][pos:max_range]
            alt = chr2seq_out[chrom][pos]

        variant_out.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t99\tPASS\t.")

    ########################################
    # Generate Structural Nucleotive Variant
    ########################################
    for _ in range(nb_struct):
        sv_type = random.choice(struct_type)
        pos = random.randint(1, len(chr2seq_out[chrom]))
        ref = chr2seq_out[chrom][pos]
        length = random.randint(struct_len[0], struct_len[1])

        variant_out.append(
            f"{chrom}\t{pos}\t.\t{ref}\t<{sv_type}>\t99\tPASS\tSVLEN={length}"
        )
    logging.info("End generate variant")

    ##################
    # Write annotation
    ##################
    logging.info("Start write annotation")
    with open(annotations_output, "wt") as fh:
        for record in select_record:
            print("\t".join(record), file=fh)
    logging.info("End write annotation")

    ################
    # Write sequence
    ################
    logging.info("Start write reference")
    with open(references_output, "wt") as fh:
        for chrom, seq in chr2seq_out.items():
            print(f">{chrom}\n{seq}", file=fh)
    logging.info("End write reference")

    ###############
    # Write variant
    ###############
    logging.info("Start write variant")
    with open(variants_output, "wt") as fh:
        for record in variant_out:
            print(record, file=fh)
    logging.info("End write variant")

    return 0


if __name__ == "__main__":
    logging.basicConfig(
        format="%(asctime)s - %(levelname)s: %(message)s",
        level=logging.DEBUG,
        encoding="utf-8",
        stream=sys.stderr,
    )

    parser = argparse.ArgumentParser(
        description="Script to get data to test variant_myth"
    )
    subparser = parser.add_subparsers()

    dl_parser = subparser.add_parser("download", help="Download data.")
    dl_parser.set_defaults(func=download)

    subsample_parser = subparser.add_parser("subsample", help="Subsample data.")
    subsample_parser.set_defaults(func=subsample)

    subsample_parser.add_argument(
        "-a", "--annotations-input", type=pathlib.Path, help="annotations input path"
    )
    subsample_parser.add_argument(
        "-r",
        "--references-input",
        type=pathlib.Path,
        help="reference sequence input path",
    )
    subsample_parser.add_argument(
        "-A",
        "--annotations-output",
        type=pathlib.Path,
        help="annotations output path",
    )
    subsample_parser.add_argument(
        "-R",
        "--references-output",
        type=pathlib.Path,
        help="reference output path",
    )
    subsample_parser.add_argument(
        "-V",
        "--variants-output",
        type=pathlib.Path,
        help="variant output path",
    )
    subsample_parser.add_argument(
        "-s", "--seed", type=int, help="random seed", default=42
    )

    argument = vars(parser.parse_args())

    sys.exit(argument.pop("func")(**argument))
