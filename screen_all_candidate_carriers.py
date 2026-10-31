#!/usr/bin/env python3

import csv
import gzip
import os
import sys

VARIANT_LIST = "CADD15_variant_list.tsv"
VCF_DIR = "/home/snondzasa/pd_qc_project/topmed_results_geno01"
OUTPUT = "CADD15_all_carriers.tsv"

# Read candidate variants
variants = []

with open(VARIANT_LIST, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        variants.append(row)

# Group variants by chromosome
variants_by_chr = {}

for v in variants:
    variants_by_chr.setdefault(v["Chr"], []).append(v)

with open(OUTPUT, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "Sample_ID",
        "Gene",
        "Chromosome",
        "Position",
        "REF",
        "ALT",
        "Genotype",
        "Zygosity",
        "Dosage"
    ])

    total_carriers = 0

    for chrom in sorted(variants_by_chr):
        vcf = f"{VCF_DIR}/{chrom}.dose.vcf.gz"

        if not os.path.exists(vcf):
            print(f"WARNING: VCF not found: {vcf}")
            continue

        # Match variants using chromosome + position + REF + ALT
        targets = {
            (v["Start"], v["Ref"], v["Alt"]): v
            for v in variants_by_chr[chrom]
        }

        print(f"Scanning {chrom}: {len(targets)} candidate variant(s)")

        with gzip.open(vcf, "rt") as f:

            sample_ids = []

            for line in f:
                if line.startswith("#CHROM"):
                    fields = line.rstrip("\n").split("\t")
                    sample_ids = fields[9:]
                    break

            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.rstrip("\n").split("\t")

                pos = fields[1]
                ref = fields[3]
                alt = fields[4]

                key = (pos, ref, alt)

                if key not in targets:
                    continue

                variant = targets[key]
                genotypes = fields[9:]

                for sample, genotype_field in zip(sample_ids, genotypes):

                    parts = genotype_field.split(":")
                    gt = parts[0]

                    if len(parts) > 3:
                        dosage = parts[3]
                    else:
                        dosage = "."

                    if gt in ("0|1", "1|0", "0/1", "1/0"):
                        zygosity = "Heterozygous"
                    elif gt in ("1|1", "1/1"):
                        zygosity = "Homozygous_alt"
                    else:
                        continue

                    writer.writerow([
                        sample,
                        variant["Gene"],
                        chrom,
                        pos,
                        ref,
                        alt,
                        gt,
                        zygosity,
                        dosage
                    ])

                    total_carriers += 1

print(f"\nCarrier screen complete.")
print(f"Total carrier records: {total_carriers}")
print(f"Output: {OUTPUT}")
