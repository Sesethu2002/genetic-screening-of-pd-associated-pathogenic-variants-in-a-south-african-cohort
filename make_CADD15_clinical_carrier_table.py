#!/usr/bin/env python3

import csv

CARRIER_FILE = "annovar_output_geno01/final_candidates_cadd15/CADD15_all_carriers.tsv"
MASTERKEY = "updated_R12_master_key_subset_STELLENBOS.txt"
OUTPUT = "CADD15_clinical_carrier_table.tsv"

# ---------------------------------------------------------
# Read Masterkey
# ---------------------------------------------------------

masterkey = {}

with open(MASTERKEY, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        sample_id = row["FID"]

        diagnosis = row["diagnosis"]

        if diagnosis == "PD":
            case_control = "Case"
        elif diagnosis == "Control":
            case_control = "Control"
        elif diagnosis == "Unaffected family member":
            case_control = "Unaffected family member"
        else:
            case_control = diagnosis

        masterkey[sample_id] = {
            "Case_Control": case_control,
            "Sex": row["sex"],
            "Age_at_sample": row["age_at_sample_collection"],
            "AAO": row["age_of_onset"],
            "Family_History": row["family_history_for_qc"],
            "Race": row["race_for_qc"]
        }

# ---------------------------------------------------------
# Join carrier information to Masterkey
# ---------------------------------------------------------

columns = [
    "Sample_ID",
    "Gene",
    "Chromosome",
    "Position",
    "REF",
    "ALT",
    "Genotype",
    "Zygosity",
    "Dosage",
    "Case_Control",
    "Sex",
    "Age_at_sample",
    "AAO",
    "Family_History",
    "Race"
]

missing = set()
records = []

with open(CARRIER_FILE, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        sample_id = row["Sample_ID"]

        if sample_id not in masterkey:
            missing.add(sample_id)
            continue

        clinical = masterkey[sample_id]

        records.append([
            sample_id,
            row["Gene"],
            row["Chromosome"],
            row["Position"],
            row["REF"],
            row["ALT"],
            row["Genotype"],
            row["Zygosity"],
            row["Dosage"],
            clinical["Case_Control"],
            clinical["Sex"],
            clinical["Age_at_sample"],
            clinical["AAO"],
            clinical["Family_History"],
            clinical["Race"]
        ])

# ---------------------------------------------------------
# Write output
# ---------------------------------------------------------

with open(OUTPUT, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(columns)
    writer.writerows(records)

# ---------------------------------------------------------
# Quality-control checks
# ---------------------------------------------------------

unique_ids = set(r[0] for r in records)

print(f"Carrier records written: {len(records)}")
print(f"Unique carrier individuals: {len(unique_ids)}")
print(f"Masterkey IDs missing: {len(missing)}")

if missing:
    print("Missing IDs:")
    for sample_id in sorted(missing):
        print(sample_id)

# Check diagnosis categories represented in the carrier table
diagnosis_counts = {}

for r in records:
    diagnosis = r[9]
    diagnosis_counts[diagnosis] = diagnosis_counts.get(diagnosis, 0) + 1

print("\nCarrier records by clinical status:")

for diagnosis, count in sorted(diagnosis_counts.items()):
    print(f"{diagnosis}: {count}")

