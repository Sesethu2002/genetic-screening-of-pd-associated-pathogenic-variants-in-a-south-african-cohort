import csv
from pathlib import Path

BASE = Path("/home/snondzasa/pd_qc_project")

CLINICAL = BASE / "CADD15_clinical_carrier_detailed.tsv"
IMPUTATION = (
    BASE / "annovar_output_geno01/final_candidates_cadd15/"
    "CADD15_imputation_metrics_with_status.tsv"
)

OUTPUT_DIR = BASE / "final_case_table_geno01"
OUTPUT_DIR.mkdir(exist_ok=True)

OUTPUT = OUTPUT_DIR / "CADD15_final_cases_R2gt0.8.tsv"
AUDIT = OUTPUT_DIR / "CADD15_case_variants_R2_filter_audit.tsv"

def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def variant_key(row, chr_col, pos_col, ref_col, alt_col):
    chrom = row[chr_col].strip()
    if not chrom.lower().startswith("chr"):
        chrom = "chr" + chrom
    return (
        chrom.lower(),
        row[pos_col].strip(),
        row[ref_col].strip().upper(),
        row[alt_col].strip().upper()
    )

clinical_rows = read_tsv(CLINICAL)
imputation_rows = read_tsv(IMPUTATION)

# Match variants by chromosome, position, REF and ALT.
r2_by_variant = {}
for row in imputation_rows:
    key = variant_key(row, "Chr", "Start", "Ref", "Alt")
    r2_by_variant[key] = row.get("R2", "").strip()

columns = [
    "Sample_ID",
    "Gene",
    "Variant",
    "Consequence",
    "Amino_acid_change",
    "Chromosome",
    "Position",
    "REF",
    "ALT",
    "rsID",
    "Genotype",
    "Zygosity",
    "Sex",
    "Age_at_sample",
    "AAO",
    "Family_History",
    "gnomAD_FAF95_max",
    "CADD_Phd",
    "REVEL",
    "ClinVar",
    "R2",
]

# The detailed clinical table has no rsID column.
# Keep rsID explicitly marked as NA until a verified rsID source is joined.
kept = []
audit = []

for row in clinical_rows:
    # Case-only output. Do not infer status from genotype files.
    if row.get("Case_Control", "").strip() != "Case":
        continue

    key = variant_key(row, "Chromosome", "Position", "REF", "ALT")
    r2_text = r2_by_variant.get(key, "")

    try:
        r2 = float(r2_text)
    except (TypeError, ValueError):
        r2 = None

    record = {
        "Sample_ID": row.get("Sample_ID", ""),
        "Gene": row.get("Gene", ""),
        "Variant": row.get("Variant", ""),
        "Consequence": row.get("Consequence", ""),
        "Amino_acid_change": row.get("Amino_acid_change", ""),
        "Chromosome": row.get("Chromosome", ""),
        "Position": row.get("Position", ""),
        "REF": row.get("REF", ""),
        "ALT": row.get("ALT", ""),
        "rsID": "NA",
        "Genotype": row.get("Genotype", ""),
        "Zygosity": row.get("Zygosity", ""),
        "Sex": row.get("Sex", ""),
        "Age_at_sample": row.get("Age_at_sample", ""),
        "AAO": row.get("AAO", ""),
        "Family_History": row.get("Family_History", ""),
        "gnomAD_FAF95_max": row.get("gnomAD_FAF95_max", ""),
        "CADD_Phd": row.get("CADD_Phd", ""),
        "REVEL": row.get("REVEL", ""),
        "ClinVar": row.get("ClinVar", ""),
        "R2": r2_text if r2_text else "NA",
    }

    if r2 is not None and r2 > 0.8:
        kept.append(record)
        audit_status = "PASS_R2_gt_0.8"
    else:
        audit_status = (
            "FAIL_R2_not_gt_0.8" if r2 is not None
            else "MISSING_R2"
        )

    audit_record = dict(record)
    audit_record["R2_filter_status"] = audit_status
    audit.append(audit_record)

def write_tsv(path, fieldnames, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)

write_tsv(OUTPUT, columns, kept)
write_tsv(AUDIT, columns + ["R2_filter_status"], audit)

print(f"Case carrier records in detailed input: "
      f"{sum(r.get('Case_Control', '').strip() == 'Case' for r in clinical_rows)}")
print(f"Case records passing R2 > 0.8: {len(kept)}")
print(f"Case records failing or missing R2: {len(audit) - len(kept)}")
print(f"Final table: {OUTPUT}")
print(f"Case-only R2 audit: {AUDIT}")
print("Note: rsID is NA pending a verified rsID annotation source.")
