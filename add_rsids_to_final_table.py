import csv
import glob
import os
import re

base = "annovar_output_geno01/final_candidates_cadd15"
input_file = "final_case_table_geno01/CADD15_final_cases_R2gt0.8.tsv"
output_file = "final_case_table_geno01/CADD15_final_cases_R2gt0.8_with_rsIDs.tsv"

def norm_chr(value):
    value = value.strip().lower()
    return value[3:] if value.startswith("chr") else value

def key(chrom, pos, ref, alt):
    return (norm_chr(chrom), str(pos).strip(), ref.strip().upper(), alt.strip().upper())

rsids = {}
for path in glob.glob(os.path.join(base, "*_final_cadd15.tsv")):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            match = re.search(r"\brs[0-9]+\b", "\t".join(row.values()))
            if not match:
                continue
            variant_key = key(row["Chr"], row["Start"], row["Ref"], row["Alt"])
            rsid = match.group(0)
            rsids.setdefault(variant_key, set()).add(rsid)

with open(input_file, newline="") as src:
    reader = csv.DictReader(src, delimiter="\t")
    columns = reader.fieldnames
    rows = list(reader)

updated = 0
ambiguous = []
with open(output_file, "w", newline="") as dst:
    writer = csv.DictWriter(dst, fieldnames=columns, delimiter="\t")
    writer.writeheader()

    for row in rows:
        variant_key = key(
            row["Chromosome"], row["Position"], row["REF"], row["ALT"]
        )
        matches = rsids.get(variant_key, set())

        if len(matches) == 1:
            row["rsID"] = next(iter(matches))
            updated += 1
        elif len(matches) > 1:
            ambiguous.append((variant_key, sorted(matches)))
            row["rsID"] = "NA"
        else:
            row["rsID"] = "NA"

        writer.writerow(row)

print(f"Case-carrier records: {len(rows)}")
print(f"Records assigned an rsID: {updated}")
print(f"Records without a unique rsID: {len(rows) - updated}")
print(f"Output: {output_file}")

if ambiguous:
    print("Ambiguous mappings (left as NA):")
    for variant_key, matches in ambiguous:
        print(variant_key, ",".join(matches))
