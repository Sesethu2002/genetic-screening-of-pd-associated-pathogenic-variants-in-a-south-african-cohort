import pandas as pd

BASE = "/home/snondzasa/pd_qc_project"

CARRIERS = f"{BASE}/CADD15_clinical_carrier_detailed.tsv"
SAMPLES = f"{BASE}/final_QC_1355_sampleIDs.txt"
MASTERKEY = f"{BASE}/updated_R12_master_key_subset_STELLENBOS.txt"

OUTPUT = f"{BASE}/CADD15_ancestry_carrier_frequencies.tsv"

# Read final QC sample IDs
with open(SAMPLES, "r") as f:
    final_ids = {line.strip() for line in f if line.strip()}

# Read Masterkey
masterkey = pd.read_csv(
    MASTERKEY,
    sep="\t",
    dtype=str,
    engine="python"
)

# Remove carriage returns from column names and values
masterkey.columns = masterkey.columns.str.replace("\r", "", regex=False).str.strip()

for col in masterkey.columns:
    masterkey[col] = (
        masterkey[col]
        .astype(str)
        .str.replace("\r", "", regex=False)
        .str.strip()
    )

# Keep only the final 1,355 QC individuals
masterkey_final = masterkey[
    masterkey["GP2ID"].isin(final_ids)
].copy()

print(f"Final QC individuals matched to Masterkey: {len(masterkey_final)}")

# Verify ancestry denominators
denominators = (
    masterkey_final["race_for_qc"]
    .value_counts()
    .to_dict()
)

print("\nAncestry denominators:")
for ancestry, count in sorted(denominators.items()):
    print(f"{ancestry}: {count}")

# Read carrier table
carriers = pd.read_csv(
    CARRIERS,
    sep="\t",
    dtype=str
)

carriers.columns = carriers.columns.str.strip()

# Clean IDs
carriers["Sample_ID"] = carriers["Sample_ID"].str.strip()

# Match carrier records to ancestry
carrier_ancestry = carriers.merge(
    masterkey_final[["GP2ID", "race_for_qc"]],
    left_on="Sample_ID",
    right_on="GP2ID",
    how="left"
)

# Check that all carriers matched
missing = carrier_ancestry["race_for_qc"].isna().sum()

if missing:
    print(f"\nWARNING: {missing} carrier records did not match ancestry metadata.")
else:
    print("\nAll carrier records matched to ancestry metadata.")

# Count carriers by variant and ancestry
counts = (
    carrier_ancestry
    .groupby(
        ["Gene", "Variant", "race_for_qc"],
        dropna=False
    )
    .size()
    .reset_index(name="Carriers")
)

# Ancestry categories in the desired order
ancestries = [
    "Asian",
    "Black or African American",
    "Multi-racial",
    "Other",
    "White"
]

rows = []

for _, row in counts.iterrows():
    pass

# Create one row per variant
variants = (
    carrier_ancestry[["Gene", "Variant"]]
    .drop_duplicates()
    .sort_values(["Gene", "Variant"])
)

for _, variant in variants.iterrows():

    gene = variant["Gene"]
    protein = variant["Variant"]

    subset = carrier_ancestry[
        (carrier_ancestry["Gene"] == gene) &
        (carrier_ancestry["Variant"] == protein)
    ]

    result = {
        "Gene": gene,
        "Variant": protein
    }

    total = 0

    for ancestry in ancestries:
        n = (subset["race_for_qc"] == ancestry).sum()
        denominator = denominators.get(ancestry, 0)

        percentage = (
            (n / denominator) * 100
            if denominator > 0
            else 0
        )

        result[f"{ancestry}_n"] = n
        result[f"{ancestry}_percent"] = round(percentage, 4)

        total += n

    result["Total_carriers"] = total

    rows.append(result)

output = pd.DataFrame(rows)

output.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)

print(f"\nOutput written to:")
print(OUTPUT)

print(f"\nVariants analysed: {len(output)}")
print(f"Carrier records analysed: {len(carrier_ancestry)}")
