import pandas as pd
import re

# ---------------------------------------------------------
# INPUT FILES
# ---------------------------------------------------------

CARRIERS = "/home/snondzasa/pd_qc_project/annovar_output_geno01/final_candidates_cadd15/CADD15_all_carriers.tsv"

ANNOTATION = "/home/snondzasa/pd_qc_project/annovar_output_geno01/final_candidates_cadd15/CADD15_master_candidates.tsv"

MASTERKEY = "/home/snondzasa/pd_qc_project/updated_R12_master_key_subset_STELLENBOS.txt"

OUTPUT = "/home/snondzasa/pd_qc_project/CADD15_clinical_carrier_detailed.tsv"


# ---------------------------------------------------------
# FUNCTION: READ FILE WITH CR LINE ENDINGS
# ---------------------------------------------------------

def read_cr_file(path):
    with open(path, "r", newline="") as f:
        text = f.read()

    # Convert carriage returns to normal line breaks
    text = text.replace("\r\n", "\n")
    text = text.replace("\r", "\n")

    from io import StringIO

    return pd.read_csv(
        StringIO(text),
        sep="\t",
        dtype=str
    )


# ---------------------------------------------------------
# 1. READ ORIGINAL CARRIER FILE
# ---------------------------------------------------------

carriers = read_cr_file(CARRIERS)

carriers.columns = carriers.columns.str.strip()

print(f"Carrier records read: {len(carriers)}")


# ---------------------------------------------------------
# 2. READ MASTERKEY
# ---------------------------------------------------------

masterkey = read_cr_file(MASTERKEY)

masterkey.columns = masterkey.columns.str.strip()

print(f"Masterkey records read: {len(masterkey)}")


# ---------------------------------------------------------
# 3. READ VARIANT ANNOTATION
# ---------------------------------------------------------

annotation = pd.read_csv(
    ANNOTATION,
    sep="\t",
    dtype=str
)

annotation.columns = annotation.columns.str.strip()


# ---------------------------------------------------------
# 4. CLEAN MATCHING COLUMNS
# ---------------------------------------------------------

for col in ["Gene", "Chromosome", "Position"]:
    carriers[col] = carriers[col].astype(str).str.strip()

carriers["Chromosome"] = carriers["Chromosome"].str.replace(
    r"^chr",
    "",
    regex=True
)

annotation["Chr"] = annotation["Chr"].astype(str).str.strip()
annotation["Chr"] = annotation["Chr"].str.replace(
    r"^chr",
    "",
    regex=True
)

annotation["Start"] = annotation["Start"].astype(str).str.strip()

annotation["Gene.refGene"] = (
    annotation["Gene.refGene"]
    .astype(str)
    .str.strip()
)


# ---------------------------------------------------------
# 5. CALCULATE GNOMAD FAF95 MAX
# ---------------------------------------------------------

def max_faf95(row):

    values = []

    for col in [
        "gnomad41_genome_fafmax_faf95_max",
        "gnomad41_exome_fafmax_faf95_max"
    ]:

        value = row.get(col, ".")

        if pd.notna(value) and str(value).strip() not in ["", "."]:

            try:
                values.append(float(value))
            except ValueError:
                pass

    if values:
        return max(values)

    return "."


annotation["gnomAD_FAF95_max"] = annotation.apply(
    max_faf95,
    axis=1
)


# ---------------------------------------------------------
# 6. EXTRACT AMINO ACID CHANGE
# ---------------------------------------------------------

def extract_protein_change(value):

    if pd.isna(value):
        return "."

    matches = re.findall(
        r"p\.[A-Za-z*?0-9]+",
        str(value)
    )

    if matches:
        return matches[0]

    return "."


annotation["Amino_acid_change"] = annotation[
    "AAChange.refGene"
].apply(extract_protein_change)


# ---------------------------------------------------------
# 7. PREPARE ANNOTATION TABLE
# ---------------------------------------------------------

annotation_small = annotation[
    [
        "Gene.refGene",
        "Chr",
        "Start",
        "ExonicFunc.refGene",
        "Amino_acid_change",
        "gnomAD_FAF95_max",
        "CADD_phred",
        "REVEL_score",
        "CLNSIG"
    ]
].copy()


annotation_small = annotation_small.rename(
    columns={
        "Gene.refGene": "Gene",
        "Chr": "Chromosome",
        "Start": "Position",
        "ExonicFunc.refGene": "Consequence",
        "CADD_phred": "CADD_Phd",
        "REVEL_score": "REVEL",
        "CLNSIG": "ClinVar"
    }
)

annotation_small["Variant"] = annotation_small["Amino_acid_change"]

# ---------------------------------------------------------
# 8. MATCH CARRIERS TO VARIANT ANNOTATION
# ---------------------------------------------------------

merged = carriers.merge(
    annotation_small,
    on=["Gene", "Chromosome", "Position"],
    how="left",
    validate="many_to_one"
)


# ---------------------------------------------------------
# 9. CHECK ANNOTATION MATCHING
# ---------------------------------------------------------

missing_annotation = merged["Amino_acid_change"].isna().sum()

if missing_annotation > 0:

    print(
        f"WARNING: {missing_annotation} carrier records "
        "did not match annotation."
    )

else:

    print(
        "All carrier records matched to variant annotation."
    )


# ---------------------------------------------------------
# 10. PREPARE MASTERKEY
# ---------------------------------------------------------

masterkey["GP2ID"] = masterkey["GP2ID"].astype(str).str.strip()

# Keep only clinical information needed for the table
masterkey_small = masterkey[
    [
        "GP2ID",
        "diagnosis",
        "sex",
        "age_at_sample_collection",
        "age_of_onset",
        "family_history_for_qc",
        "race_for_qc"
    ]
].copy()


# Rename Masterkey fields
masterkey_small = masterkey_small.rename(
    columns={
        "GP2ID": "Sample_ID",
        "age_at_sample_collection": "Age_at_sample",
        "age_of_onset": "AAO",
        "family_history_for_qc": "Family_History",
        "race_for_qc": "Race"
    }
)


# ---------------------------------------------------------
# 11. CONVERT DIAGNOSIS TO CASE/CONTROL LABEL
# ---------------------------------------------------------

def diagnosis_to_status(value):

    value = str(value).strip()

    if value == "PD":
        return "Case"

    if value == "Control":
        return "Control"

    if value == "Unaffected family member":
        return "Unaffected family member"

    return value


masterkey_small["Case_Control"] = masterkey_small[
    "diagnosis"
].apply(diagnosis_to_status)


# ---------------------------------------------------------
# 12. MERGE MASTERKEY CLINICAL INFORMATION
# ---------------------------------------------------------

merged = merged.merge(
    masterkey_small[
        [
            "Sample_ID",
            "Case_Control",
            "sex",
            "Age_at_sample",
            "AAO",
            "Family_History",
            "Race"
        ]
    ],
    on="Sample_ID",
    how="left",
    validate="many_to_one"
)


# ---------------------------------------------------------
# 13. CHECK MASTERKEY MATCHING
# ---------------------------------------------------------

missing_masterkey = merged["Case_Control"].isna().sum()

if missing_masterkey > 0:

    print(
        f"WARNING: {missing_masterkey} carrier records "
        "did not match the Masterkey."
    )

else:

    print(
        "All carrier records matched to Masterkey clinical information."
    )


# ---------------------------------------------------------
# 14. FINAL COLUMN ORDER
# ---------------------------------------------------------

final_columns = [
    "Sample_ID",
    "Gene",
    "Variant",
    "Consequence",
    "Amino_acid_change",
    "gnomAD_FAF95_max",
    "CADD_Phd",
    "REVEL",
    "ClinVar",
    "Chromosome",
    "Position",
    "REF",
    "ALT",
    "Genotype",
    "Zygosity",
    "Dosage",
    "Case_Control",
    "sex",
    "Age_at_sample",
    "AAO",
    "Family_History",
    "Race"
]

merged = merged[final_columns]


# ---------------------------------------------------------
# 15. RENAME SEX COLUMN FOR THESIS TABLE
# ---------------------------------------------------------

merged = merged.rename(
    columns={
        "sex": "Sex"
    }
)


# ---------------------------------------------------------
# 16. WRITE CLEAN TSV
# ---------------------------------------------------------

merged.to_csv(
    OUTPUT,
    sep="\t",
    index=False,
    lineterminator="\n"
)


# ---------------------------------------------------------
# 17. FINAL VALIDATION
# ---------------------------------------------------------

print()
print("======================================")
print("DETAILED CADD15 TABLE CREATED")
print("======================================")
print(f"Output: {OUTPUT}")
print(f"Carrier records: {len(merged)}")
print(f"Columns: {len(merged.columns)}")

if len(merged) == 77:
    print("PASS: 77 carrier records retained.")
else:
    print("WARNING: Expected 77 carrier records.")

print("======================================")
