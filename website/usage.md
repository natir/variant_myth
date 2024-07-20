# Usage

To perform annotation run :
```bash
variant_myth -i variants.vcf -r reference.fasta -a annotations.gff \
-t translate_table.txt -o output.json
```

All input could be compress in gzip, bzip2, bgzip or xz format.

If you install variant_myth with `parallel` feature you can add `--threads` option.
