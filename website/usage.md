# Usage

To perform annotation run :
```bash
variant_myth -i variants.vcf -r reference.fasta -a annotations.gff \
-t translate_table.txt parquet -o output.parquet
```

All input could be compress in gzip, bzip2, bgzip or xz format.

If you install variant_myth with `parallel` feature you can add `--threads` option.

## Get test data

You could run `get_data.py download` to create `data` directory and download reference sequence, annotation and some variant from human public data (runtime ~5 min)
