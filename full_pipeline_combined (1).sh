#!/bin/bash
########################################################################################
# PD PATHOGENIC VARIANT SCREENING - FULL PIPELINE (combined)
# South African cohort (GP2 Stellenbosch), 21 monogenic PD genes
# (Blauwendraat et al. 2020, Lancet Neurol, Table 1)
#
# This file concatenates all pipeline steps in order for reference/GitHub.
#
# STEP TYPES
#   Bash steps   : 1, 2, 3, 4, 5, 7, 9, 12
#   Python steps : 8, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21
#                  (all embedded below via `python3 << 'PYEOF'` heredocs)
#   PBS step     : 6 - PBS array job; submit it SEPARATELY with
#                  `qsub 06_annovar_annotation.pbs`, do not run it inline with the rest.
#
# DIRECTORY LAYOUT (all set once in the CONFIGURATION block below)
#   PROJECT_DIR   working directory: QC outputs, clinical and ancestry tables, final_case_table/
#   ANNOVAR_DIR   annotated chromosomes, gene_extracts/, final_candidates/, final_candidates_cadd15/
#   TOPMED_DIR    TOPMed imputation results (chr*.dose.vcf.gz)
#
# INPUTS NOT BUILT IN THIS SCRIPT
#   ${ANNOVAR_DIR}/final_candidates_cadd15/CADD15_variant_list.tsv                  (Step 14)
#   ${ANNOVAR_DIR}/final_candidates_cadd15/CADD15_master_candidates.tsv             (Step 15)
#   ${ANNOVAR_DIR}/final_candidates_cadd15/CADD15_imputation_metrics_with_status.tsv (Step 17)
#   ${PROJECT_DIR}/cohort_ancestry.tsv                                              (Steps 20-21)
#
# Python steps 8-21 need pandas; the virtual environment is activated in Step 8.
#
# Must run inside an interactive/batch HPC session, not the login node.
########################################################################################
set -euo pipefail

######################################[ CONFIGURATION ]###################################
# Edit the paths in this block only. Python steps read them from the environment.
# Run this script from the project working directory (it becomes PROJECT_DIR).
##########################################################################################

export PROJECT_DIR="$(pwd)"
export ANNOVAR_DIR="/path/to/annovar_output"
export TOPMED_DIR="/path/to/topmed_results"
export MASTERKEY="/path/to/updated_R12_master_key_subset_STELLENBOS.txt"

REF_FASTA="/path/to/Homo_sapiens_assembly38.fasta"
GWASQC_DIR="/path/to/GWASQC"
NATORA_DIR="/path/to/NAToRA_Public"

#####################################[ STEP 1 ]#########################################
# Type: Bash
# Strand/reference flip, BEFORE QC. Flipping first lets bcftools +fixref resolve
# ambiguous-strand variants before the QC "ATCG" step would otherwise discard them.
##########################################################################################

RAW_INPUT="GP2_NBA_raw_STELLENBOS"

module load plink/2.00a bcftools samtools

samtools faidx "${REF_FASTA}"
plink2 --pfile "${RAW_INPUT}" --recode vcf-iid bgz --out "${RAW_INPUT}_preflip" --output-chr chrMT
bcftools index "${RAW_INPUT}_preflip.vcf.gz"
bcftools +fixref "${RAW_INPUT}_preflip.vcf.gz" -Oz -o "${RAW_INPUT}_flipped.vcf.gz" -- \
    -f "${REF_FASTA}" -m flip -d
bcftools index "${RAW_INPUT}_flipped.vcf.gz" -f
plink2 --vcf "${RAW_INPUT}_flipped.vcf.gz" --vcf-half-call missing \
    --make-bed --out "${RAW_INPUT}_flipped_bed"

#####################################[ STEP 2 ]#########################################
# Type: Bash
# Build sex/status covariate file. Join key: IID/GP2ID (unique) - NOT FID (shared by
# relatives). Status source: Masterkey "GP2_phenotype_for_qc" column.
##########################################################################################

RAW_INPUT_BED="${RAW_INPUT}_flipped_bed"

awk -F'\t' '
    BEGIN { print "ID,SEX,STATUS" }
    NR==FNR {
        if (FNR==1) next
        gp2id = $3
        sex_raw = $19
        status_raw = $16
        sex = (sex_raw=="Male" ? "M" : (sex_raw=="Female" ? "F" : "Unknown"))
        status = (status_raw=="PD" ? "Case" : (status_raw=="Control" ? "Control" : "Unknown"))
        mk_sex[gp2id] = sex
        mk_status[gp2id] = status
        next
    }
    { iid=$2; if (iid in mk_sex) print iid","mk_sex[iid]","mk_status[iid]; else print iid",Unknown,Unknown" }
' "${MASTERKEY}" "${RAW_INPUT_BED}.fam" > info_file.csv

awk -F',' 'NR==1 || $3!="Unknown"' info_file.csv > info_file_ccOnly.csv
sed 's/,/\t/g' info_file_ccOnly.csv > info_file_ccOnly.tsv
awk -F',' '$3=="Unknown"{print $1,$1}' info_file.csv > remove_unknown.txt

plink2 --bfile "${RAW_INPUT_BED}" --remove remove_unknown.txt \
    --make-bed --out "${RAW_INPUT_BED}_ccOnly"

echo "Status breakdown:"
tail -n +2 info_file_ccOnly.csv | cut -d',' -f3 | sort | uniq -c

#####################################[ STEP 3 ]#########################################
# Type: Bash
# QC Run 1 - related individuals EXCLUDED (population structure only, not used downstream)
# HWE deliberately excluded from both runs. qc_steps.txt MUST be tab-delimited.
##########################################################################################

INPUT="${RAW_INPUT_BED}_ccOnly"

module load plink/1.90 python/3.11.9

printf "sex-check\nATCG\nmind\t0.05\ngeno\t0.01\nduplicate\nheterozygosity\t3\nrelationship\t0.0884\n" > qc_steps.txt

python "${GWASQC_DIR}/main.py" \
    -i "${INPUT}" -o nba_qced -O qced_data \
    -I info_file_ccOnly.tsv -s qc_steps.txt \
    --NAToRA "${NATORA_DIR}/NAToRA_Public.py" \
    2>&1 | tee run1_full_log.txt

#####################################[ STEP 4 ]#########################################
# Type: Bash
# QC Run 2 - related individuals INCLUDED (MAIN ANALYSIS DATASET)
# Identical to Run 1 except the "relationship" step is omitted.
##########################################################################################

printf "sex-check\nATCG\nmind\t0.05\ngeno\t0.01\nduplicate\nheterozygosity\t3\n" > qc_steps_pathogenic.txt

python "${GWASQC_DIR}/main.py" \
    -i "${INPUT}" -o nba_qced_pathogenic -O qced_data_pathogenic \
    -I info_file_ccOnly.tsv -s qc_steps_pathogenic.txt \
    --NAToRA "${NATORA_DIR}/NAToRA_Public.py" \
    2>&1 | tee run2_full_log.txt

# Final QC Run 2 output: qced_data_pathogenic/nba_qced_pathogenic_Heterozygosity.pgen/.pvar/.psam

#####################################[ STEP 5 ]#########################################
# Type: Bash
# Second flip pass + split by chromosome, BEFORE imputation submission (safety check -
# this pass found 0% mismatch, confirming Step 1's flip already worked)
##########################################################################################

QC2_FINAL="qced_data_pathogenic/nba_qced_pathogenic_Heterozygosity"
OUT_PREFIX="nba_pathogenic_preimpute"

plink2 --pfile "${QC2_FINAL}" --recode vcf-iid bgz --out "${OUT_PREFIX}" --output-chr chrMT
bcftools index "${OUT_PREFIX}.vcf.gz"
bcftools +fixref "${OUT_PREFIX}.vcf.gz" -Oz -o "${OUT_PREFIX}_flipped.vcf.gz" -- \
    -f "${REF_FASTA}" -m flip -d
bcftools index "${OUT_PREFIX}_flipped.vcf.gz" -f

for chr in $(seq 1 22); do
    bcftools view --regions chr${chr} -Oz -o "${OUT_PREFIX}_chr${chr}.vcf.gz" "${OUT_PREFIX}_flipped.vcf.gz"
done

# Upload the 22 *_chr*.vcf.gz files to TOPMed Imputation Server:
#   Array build: GRCh38 | Reference panel: TOPMed-r3 (hg38) | Rsq filter: 0.3
#   Phasing Engine: Eagle v2.4 | Mode: Quality Control & Imputation

#####################################[ STEP 6 ]#########################################
# Type: PBS array job (submit separately)
# ANNOVAR annotation of all 22 imputed chromosomes.
# *** SUBMIT SEPARATELY, DO NOT RUN INLINE: ***
#     qsub 06_annovar_annotation.pbs
# (PBS array job #PBS -J 1-22; see 06_annovar_annotation.pbs for the full job script -
#  protocol: refGene,clinvar_20250721,dbnsfp47a,gnomad41_genome,gnomad41_exome)
##########################################################################################

#####################################[ STEP 7 ]#########################################
# Type: Bash
# Fast per-gene extraction from the annotated chromosome files (avoids loading full
# multi-GB chromosome files into pandas directly - that took >8h and was killed twice)
##########################################################################################

cd "${ANNOVAR_DIR}"
mkdir -p gene_extracts

extract_chrom() {
    local chrom=$1; shift; local genes=("$@")
    local wanted=$(printf ",%s" "${genes[@]}"); wanted=${wanted:1}
    echo "Extracting from chr${chrom}: ${genes[*]}"
    awk -F'\t' -v OFS='\t' -v genelist="$wanted" -v outdir="gene_extracts" '
        BEGIN { n = split(genelist, g, ","); for (i=1;i<=n;i++) wanted[g[i]]=1 }
        NR==1 { header=$0; for (i=1;i<=n;i++) print header > (outdir"/"g[i]"_variants.tsv"); next }
        ($7 in wanted) { print > (outdir"/"$7"_variants.tsv") }
    ' "annotated_chr${chrom}.hg38_multianno.txt"
}

extract_chrom 1 PARK7 PINK1 ATP13A2 GBA DNAJC6
extract_chrom 2 HTRA2 GIGYF2
extract_chrom 3 EIF4G1 DNAJC13
extract_chrom 4 SNCA UCHL1
extract_chrom 6 PRKN
extract_chrom 12 LRRK2
extract_chrom 14 LRP10
extract_chrom 15 POLG VPS13C
extract_chrom 16 VPS35
extract_chrom 20 TMEM230
extract_chrom 21 SYNJ1
extract_chrom 22 FBXO7 PLA2G6

wc -l gene_extracts/*.tsv

#####################################[ STEP 8 ]#########################################
# Type: Python (embedded heredoc)
# Apply variant filters per gene (MAF<=1%, exonic, CADD>=20, REVEL>=0.7 for missense)
# and build the per-gene summary table.
##########################################################################################

module load python/3.11.9
source ~/py311env/bin/activate   # venv with pandas installed

python3 << 'PYEOF'
import os
import pandas as pd

ANNOVAR_DIR = os.environ["ANNOVAR_DIR"]
WORK_DIR = f"{ANNOVAR_DIR}/gene_extracts"
OUT_DIR = ANNOVAR_DIR

genes = [
    "SNCA", "PRKN", "PARK7", "LRRK2", "PINK1", "ATP13A2", "FBXO7", "GBA",
    "PLA2G6", "VPS35", "POLG", "DNAJC6", "SYNJ1", "VPS13C", "UCHL1",
    "HTRA2", "GIGYF2", "EIF4G1", "DNAJC13", "TMEM230", "LRP10",
]
MAF_COLS = ["gnomad41_genome_fafmax_faf95_max", "gnomad41_exome_fafmax_faf95_max"]
CADD_COL = "CADD_phred"
REVEL_COL = "REVEL_score"
summary = {}

for gene in genes:
    print(f"WORKING ON: {gene}")
    df = pd.read_csv(f"{WORK_DIR}/{gene}_variants.tsv", sep="\t", low_memory=False)
    num_cols = MAF_COLS + [CADD_COL, REVEL_COL]
    df[num_cols] = df[num_cols].apply(pd.to_numeric, errors="coerce")
    total = len(df)

    func_counts = {}
    for func in ["intronic", "upstream", "downstream", "UTR3", "UTR5", "splicing", "exonic"]:
        func_counts[func] = len(df[df["Func.refGene"] == func])

    exonic = df[df["Func.refGene"] == "exonic"]
    exonic_func_counts = {}
    for ef in ["stopgain", "stoploss", "startloss", "frameshift deletion", "frameshift insertion",
               "nonframeshift deletion", "nonframeshift insertion", "synonymous SNV", "nonsynonymous SNV"]:
        exonic_func_counts[ef] = len(exonic[exonic["ExonicFunc.refGene"] == ef])

    maf_mask = (
        (df[MAF_COLS[0]].isna() | (df[MAF_COLS[0]] <= 0.01)) &
        (df[MAF_COLS[1]].isna() | (df[MAF_COLS[1]] <= 0.01))
    )
    maf_01pct = len(df[
        (df[MAF_COLS[0]].isna() | (df[MAF_COLS[0]] <= 0.001)) &
        (df[MAF_COLS[1]].isna() | (df[MAF_COLS[1]] <= 0.001))
    ])
    cadd_pass = len(df[df[CADD_COL] >= 20])
    revel_pass = len(df[df[REVEL_COL] >= 0.7])

    is_missense = df["ExonicFunc.refGene"] == "nonsynonymous SNV"
    revel_ok = (~is_missense) | (df[REVEL_COL] >= 0.7)
    final_mask = maf_mask & (df["Func.refGene"] == "exonic") & (df[CADD_COL] >= 20) & revel_ok
    final = df[final_mask]
    final.to_csv(f"{OUT_DIR}/{gene}_filtered_variants.tsv", sep="\t", index=False)

    counts = {"Gene": gene, "Total variants": total}
    counts.update({k.capitalize(): v for k, v in func_counts.items()})
    counts.update({k.capitalize(): v for k, v in exonic_func_counts.items()})
    counts["MAF <=1%"] = int(maf_mask.sum())
    counts["MAF <=0.1%"] = maf_01pct
    counts["CADD >=20"] = cadd_pass
    counts["REVEL >=0.7"] = revel_pass
    counts["Passing all filters"] = len(final)
    summary[gene] = counts

pd.DataFrame(summary).T.to_csv(f"{OUT_DIR}/PD_genes_variant_summary.tsv", sep="\t", index=False)
print("Saved summary table to PD_genes_variant_summary.tsv")
PYEOF

#####################################[ STEP 9 ]#########################################
# Type: Bash
# Exclude synonymous variants - only protein-altering variants are real candidates
##########################################################################################

mkdir -p final_candidates
for f in *_filtered_variants.tsv; do
    gene=$(basename "$f" _filtered_variants.tsv)
    awk -F'\t' 'NR==1 || $9 != "synonymous SNV"' "$f" > "final_candidates/${gene}_final.tsv"
    count=$(tail -n +2 "final_candidates/${gene}_final.tsv" | wc -l)
    echo "$gene: $count final candidate(s)"
done

#####################################[ STEP 10 ]########################################
# Type: Python (embedded heredoc)
# Identify which individuals carry each final candidate variant
##########################################################################################

python3 << 'PYEOF'
import csv
import gzip
import os

ANNOVAR_DIR = os.environ["ANNOVAR_DIR"]
FINAL_DIR = f"{ANNOVAR_DIR}/final_candidates"
VCF_DIR = os.environ["TOPMED_DIR"]
MASTERKEY = os.environ["MASTERKEY"]
OUT_FILE = f"{ANNOVAR_DIR}/candidate_variant_carriers.tsv"

genes = {
    "SNCA": 4, "PRKN": 6, "PARK7": 1, "LRRK2": 12, "PINK1": 1,
    "ATP13A2": 1, "FBXO7": 22, "GBA": 1, "PLA2G6": 22, "VPS35": 16,
    "POLG": 15, "DNAJC6": 1, "SYNJ1": 21, "VPS13C": 15,
    "UCHL1": 4, "HTRA2": 2, "GIGYF2": 2, "EIF4G1": 3, "DNAJC13": 3,
    "TMEM230": 20, "LRP10": 14,
}
FIRST_SAMPLE_COL = 222

print("Loading Masterkey...")
mk = {}
with open(MASTERKEY, encoding="utf-8", errors="replace") as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader)
    for row in reader:
        gp2id = row[2]
        status_raw = row[15]
        status = "Case" if status_raw == "PD" else ("Control" if status_raw == "Control" else "Unknown")
        mk[gp2id] = {"sex": row[18], "status": status, "aao": row[20], "famhist": row[25]}

sample_order_cache = {}
def get_sample_order(chrom):
    if chrom in sample_order_cache:
        return sample_order_cache[chrom]
    with gzip.open(f"{VCF_DIR}/chr{chrom}.dose.vcf.gz", "rt") as f:
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                sample_order_cache[chrom] = samples
                return samples
    raise ValueError(f"No #CHROM line in chr{chrom}")

results = []
for gene, chrom in genes.items():
    try:
        with open(f"{FINAL_DIR}/{gene}_final.tsv", encoding="utf-8", errors="replace") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)
            rows = list(reader)
    except FileNotFoundError:
        continue
    if not rows:
        continue
    print(f"{gene}: {len(rows)} final candidate variant(s)")
    samples = get_sample_order(chrom)
    for row in rows:
        chr_col, pos, ref, alt = row[0], row[1], row[3], row[4]
        exonicfunc = row[8] if len(row) > 8 else ""
        aachange = row[9] if len(row) > 9 else ""
        clnsig = row[14] if len(row) > 14 else ""
        genotype_cells = row[FIRST_SAMPLE_COL - 1:]
        if len(genotype_cells) != len(samples):
            print(f"  WARNING: genotype column mismatch for {gene} {chr_col}:{pos} -- skipping")
            continue
        for sample_id, cell in zip(samples, genotype_cells):
            gt = (cell.split(":")[0] if cell else "").replace("|", "/")
            if gt in ("0/0", "", "./.", ".|."):
                continue
            zygosity = "Het" if gt in ("0/1", "1/0") else ("Hom" if gt == "1/1" else f"Other({gt})")
            c = mk.get(sample_id, {"sex": "NA", "status": "NA", "aao": "NA", "famhist": "NA"})
            results.append({
                "Gene": gene, "Chr": chr_col, "Pos": pos, "Ref": ref, "Alt": alt,
                "ExonicFunc": exonicfunc, "AAChange": aachange, "CLNSIG": clnsig,
                "Sample": sample_id, "Zygosity": zygosity,
                "Sex": c["sex"], "Status": c["status"], "AAO": c["aao"], "FamilyHistory": c["famhist"],
            })

print(f"\nTotal carrier records: {len(results)}")
with open(OUT_FILE, "w", newline="", encoding="utf-8") as f:
    fieldnames = ["Gene","Chr","Pos","Ref","Alt","ExonicFunc","AAChange","CLNSIG",
                  "Sample","Zygosity","Sex","Status","AAO","FamilyHistory"]
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(results)
print(f"Saved: {OUT_FILE}")
PYEOF


#####################################[ STEP 11 ]########################################
# Type: Python (embedded heredoc)
# Sensitivity check: repeat the candidate filter at CADD_phred >= 15
# (result: 29 candidates vs 28; added LRRK2 p.E899D).
##########################################################################################

python3 << 'PYEOF'
"""Sensitivity check: repeat the candidate filter at CADD_phred >= 15 (result: 29 candidates vs 28; added LRRK2 p.E899D)."""
import os
import pandas as pd

ANNOVAR_DIR = os.environ["ANNOVAR_DIR"]
WORK_DIR = f"{ANNOVAR_DIR}/gene_extracts"
OUT_DIR = f"{ANNOVAR_DIR}/final_candidates_cadd15"
os.makedirs(OUT_DIR, exist_ok=True)

genes = ["SNCA", "PRKN", "PARK7", "LRRK2", "PINK1", "ATP13A2", "FBXO7", "GBA", "PLA2G6", "VPS35", "POLG",
         "DNAJC6", "SYNJ1", "VPS13C", "UCHL1", "HTRA2", "GIGYF2", "EIF4G1", "DNAJC13", "TMEM230", "LRP10"]
MAF_COLS = ["gnomad41_genome_fafmax_faf95_max", "gnomad41_exome_fafmax_faf95_max"]
CADD_COL, REVEL_COL = "CADD_phred", "REVEL_score"

total = 0
for gene in genes:
    df = pd.read_csv(f"{WORK_DIR}/{gene}_variants.tsv", sep="\t", low_memory=False)
    num = MAF_COLS + [CADD_COL, REVEL_COL]
    df[num] = df[num].apply(pd.to_numeric, errors="coerce")
    maf_ok = (df[MAF_COLS[0]].isna() | (df[MAF_COLS[0]] <= 0.01)) & (df[MAF_COLS[1]].isna() | (df[MAF_COLS[1]] <= 0.01))
    missense = df["ExonicFunc.refGene"] == "nonsynonymous SNV"
    revel_ok = (~missense) | (df[REVEL_COL] >= 0.7)
    keep = maf_ok & (df["Func.refGene"] == "exonic") & (df[CADD_COL] >= 15) & revel_ok & (df["ExonicFunc.refGene"] != "synonymous SNV")
    out = df[keep]
    out.to_csv(f"{OUT_DIR}/{gene}_final_cadd15.tsv", sep="\t", index=False)
    print(f"{gene}: {len(out)} candidate(s) at CADD>=15")
    total += len(out)
print("Total:", total)
PYEOF

#####################################[ STEP 12 ]########################################
# Type: Bash
# Demographics for the QC Run 2 sample (n = 1,355).
# Masterkey columns (1-indexed): 3 GP2ID, 16 GP2_phenotype_for_qc,
# 19 sex, 21 age_of_onset, 22 age_at_diagnosis, 25 race_for_qc, 26 family_history_for_qc.
# Note: the final QC .psam has IID in column 1.
##########################################################################################

cd "${PROJECT_DIR}"   # QC outputs (.psam) are relative to the project working directory
PSAM="qced_data_pathogenic/nba_qced_pathogenic_Heterozygosity.psam"

awk -F'\t' 'NR==FNR { if (FNR==1) next; sex[$3]=$19; aao[$3]=$21; fh[$3]=$26; st[$3]=$16; race[$3]=$25; next }
    FNR==1 { next }
    ($1 in sex) { print $1"\t"sex[$1]"\t"aao[$1]"\t"fh[$1]"\t"st[$1]"\t"race[$1] }' \
    "${MASTERKEY}" "${PSAM}" > demographics_run2.tsv

cut -f1 demographics_run2.tsv > final_QC_1355_sampleIDs.txt   # sample list used in Step 19
wc -l demographics_run2.tsv                                  # expect 1355
echo "Sex by status:";  cut -f5,2 demographics_run2.tsv | sort | uniq -c
echo "Ancestry (Masterkey race_for_qc):"; cut -f6 demographics_run2.tsv | sort | uniq -c
echo "Family history among cases:"; awk -F'\t' '$5=="PD"{print $4}' demographics_run2.tsv | sort | uniq -c
echo "Age at onset, cases (n, mean, SD):"
awk -F'\t' '$5=="PD" && $3!="" && $3!="NA"{s+=$3; q+=$3*$3; n++} END{m=s/n; printf "%d, %.1f, %.1f\n", n, m, sqrt(q/n-m*m)}' demographics_run2.tsv

#####################################[ STEP 13 ]########################################
# Type: Python (embedded heredoc)
# Per-gene CADD>=15 candidate filter (MAF<=1%, exonic, REVEL>=0.7 for missense, no synonymous).
# Standalone-script version of Step 11 (same filters, same outputs).
# Writes {gene}_final_cadd15.tsv.
##########################################################################################

python3 << 'PYEOF'
import os

import pandas as pd

ANNOVAR_DIR = os.environ["ANNOVAR_DIR"]
WORK_DIR = f"{ANNOVAR_DIR}/gene_extracts"
OUT_DIR = f"{ANNOVAR_DIR}/final_candidates_cadd15"
os.makedirs(OUT_DIR, exist_ok=True)

genes = [
    "SNCA", "PRKN", "PARK7", "LRRK2", "PINK1", "ATP13A2", "FBXO7", "GBA",
    "PLA2G6", "VPS35", "POLG", "DNAJC6", "SYNJ1", "VPS13C", "UCHL1",
    "HTRA2", "GIGYF2", "EIF4G1", "DNAJC13", "TMEM230", "LRP10",
]

MAF_COLS = ["gnomad41_genome_fafmax_faf95_max", "gnomad41_exome_fafmax_faf95_max"]
CADD_COL = "CADD_phred"
REVEL_COL = "REVEL_score"

summary = {}

for gene in genes:
    fname = f"{WORK_DIR}/{gene}_variants.tsv"
    df = pd.read_csv(fname, sep="\t", low_memory=False)
    num_cols = MAF_COLS + [CADD_COL, REVEL_COL]
    df[num_cols] = df[num_cols].apply(pd.to_numeric, errors="coerce")

    maf_mask = (
        (df[MAF_COLS[0]].isna() | (df[MAF_COLS[0]] <= 0.01)) &
        (df[MAF_COLS[1]].isna() | (df[MAF_COLS[1]] <= 0.01))
    )
    is_missense = df["ExonicFunc.refGene"] == "nonsynonymous SNV"
    revel_ok = (~is_missense) | (df[REVEL_COL] >= 0.7)
    not_synonymous = df["ExonicFunc.refGene"] != "synonymous SNV"

    final_mask = maf_mask & (df["Func.refGene"] == "exonic") & (df[CADD_COL] >= 15) & revel_ok & not_synonymous
    final = df[final_mask]

    print(f"{gene}: {len(final)} candidate(s) at CADD>=15")
    final.to_csv(f"{OUT_DIR}/{gene}_final_cadd15.tsv", sep="\t", index=False)
    summary[gene] = len(final)

print("\nTotal:", sum(summary.values()))
PYEOF

#####################################[ STEP 14 ]########################################
# Type: Python (embedded heredoc)
# Screen all 22 TOPMed-imputed chromosome VCFs for carriers of the CADD>=15 candidates.
# Requires: final_candidates_cadd15/CADD15_variant_list.tsv (columns Chr, Start, Ref, Alt, Gene).
# Output: final_candidates_cadd15/CADD15_all_carriers.tsv.
##########################################################################################

python3 << 'PYEOF'
#!/usr/bin/env python3

import csv
import gzip
import os
import sys

ANNOVAR_DIR = os.environ["ANNOVAR_DIR"]
CADD15_DIR = f"{ANNOVAR_DIR}/final_candidates_cadd15"

VARIANT_LIST = f"{CADD15_DIR}/CADD15_variant_list.tsv"
VCF_DIR = os.environ["TOPMED_DIR"]
OUTPUT = f"{CADD15_DIR}/CADD15_all_carriers.tsv"

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
PYEOF

#####################################[ STEP 15 ]########################################
# Type: Python (embedded heredoc)
# Expand the carrier table with variant annotation (aa change, gnomAD FAF95, CADD, REVEL, ClinVar)
# and Masterkey clinical data (case/control, sex, age, AAO, family history, ancestry).
# Requires: final_candidates_cadd15/CADD15_master_candidates.tsv.
# Output: ${PROJECT_DIR}/CADD15_clinical_carrier_detailed.tsv.
##########################################################################################

python3 << 'PYEOF'
import os
import re

import pandas as pd

# ---------------------------------------------------------
# INPUT FILES
# ---------------------------------------------------------

ANNOVAR_DIR = os.environ["ANNOVAR_DIR"]
PROJECT_DIR = os.environ["PROJECT_DIR"]
CADD15_DIR = f"{ANNOVAR_DIR}/final_candidates_cadd15"

CARRIERS = f"{CADD15_DIR}/CADD15_all_carriers.tsv"

ANNOTATION = f"{CADD15_DIR}/CADD15_master_candidates.tsv"

MASTERKEY = os.environ["MASTERKEY"]

OUTPUT = f"{PROJECT_DIR}/CADD15_clinical_carrier_detailed.tsv"


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
PYEOF

#####################################[ STEP 16 ]########################################
# Type: Python (embedded heredoc)
# Simple clinical carrier table: carrier records joined to Masterkey clinical fields.
# Output: ${PROJECT_DIR}/CADD15_clinical_carrier_table.tsv.
##########################################################################################

python3 << 'PYEOF'
#!/usr/bin/env python3

import csv
import os

ANNOVAR_DIR = os.environ["ANNOVAR_DIR"]
PROJECT_DIR = os.environ["PROJECT_DIR"]

CARRIER_FILE = f"{ANNOVAR_DIR}/final_candidates_cadd15/CADD15_all_carriers.tsv"
MASTERKEY = os.environ["MASTERKEY"]
OUTPUT = f"{PROJECT_DIR}/CADD15_clinical_carrier_table.tsv"

# ---------------------------------------------------------
# Read Masterkey
# ---------------------------------------------------------

masterkey = {}

with open(MASTERKEY, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        sample_id = row["GP2ID"]   # unique per individual; FID is shared by relatives

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

PYEOF

#####################################[ STEP 17 ]########################################
# Type: Python (embedded heredoc)
# Final case-only table, keeping variants with imputation R2 > 0.8 (plus an R2 audit table).
# Requires: final_candidates_cadd15/CADD15_imputation_metrics_with_status.tsv.
# Output: ${PROJECT_DIR}/final_case_table/.
##########################################################################################

python3 << 'PYEOF'
import csv
import os
from pathlib import Path

PROJECT_DIR = Path(os.environ["PROJECT_DIR"])
ANNOVAR_DIR = Path(os.environ["ANNOVAR_DIR"])

CLINICAL = PROJECT_DIR / "CADD15_clinical_carrier_detailed.tsv"
IMPUTATION = (
    ANNOVAR_DIR / "final_candidates_cadd15" /
    "CADD15_imputation_metrics_with_status.tsv"
)

OUTPUT_DIR = PROJECT_DIR / "final_case_table"
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
PYEOF

#####################################[ STEP 18 ]########################################
# Type: Python (embedded heredoc)
# Add rsIDs to the final case table (rsID assigned only when uniquely matched; otherwise NA).
# Output: ${PROJECT_DIR}/final_case_table/CADD15_final_cases_R2gt0.8_with_rsIDs.tsv.
##########################################################################################

python3 << 'PYEOF'
import csv
import glob
import os
import re

base = os.path.join(os.environ["ANNOVAR_DIR"], "final_candidates_cadd15")
case_table_dir = os.path.join(os.environ["PROJECT_DIR"], "final_case_table")
input_file = os.path.join(case_table_dir, "CADD15_final_cases_R2gt0.8.tsv")
output_file = os.path.join(case_table_dir, "CADD15_final_cases_R2gt0.8_with_rsIDs.tsv")

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
PYEOF

#####################################[ STEP 19 ]########################################
# Type: Python (embedded heredoc)
# Carrier counts and frequencies per variant by ancestry (race_for_qc), using the 1,355 QC samples.
# Requires: ${PROJECT_DIR}/final_QC_1355_sampleIDs.txt (created in Step 12).
# Output: ${PROJECT_DIR}/CADD15_ancestry_carrier_frequencies.tsv.
##########################################################################################

python3 << 'PYEOF'
import os

import pandas as pd

BASE = os.environ["PROJECT_DIR"]

CARRIERS = f"{BASE}/CADD15_clinical_carrier_detailed.tsv"
SAMPLES = f"{BASE}/final_QC_1355_sampleIDs.txt"
MASTERKEY = os.environ["MASTERKEY"]

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
PYEOF

#####################################[ STEP 20 ]########################################
# Type: Python (embedded heredoc)
# Global ancestry summary (mean, median, SD, min, max) - Overall / PD / Control.
# Requires: ${PROJECT_DIR}/cohort_ancestry.tsv. Output: ${PROJECT_DIR}/global_ancestry_summary.tsv.
##########################################################################################

python3 << 'PYEOF'
import csv
import os
import statistics

PROJECT_DIR = os.environ["PROJECT_DIR"]

components = ["NAMA", "EUR", "AFR", "EAS", "MAL"]
data = []

with open(f"{PROJECT_DIR}/cohort_ancestry.tsv", newline="") as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        if len(row) != 7:
            print("Skipping unexpected row:", row)
            continue

        sample_id, diagnosis = row[0], row[1]

        try:
            values = [float(x) * 100 for x in row[2:7]]
        except ValueError:
            print("Skipping non-numeric row:", row)
            continue

        data.append((diagnosis, values))

groups = {
    "Overall": data,
    "PD": [r for r in data if r[0] == "PD"],
    "Control": [r for r in data if r[0] == "Control"]
}

output = f"{PROJECT_DIR}/global_ancestry_summary.tsv"

with open(output, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow([
        "Group", "Ancestry", "N", "Mean_percent",
        "Median_percent", "SD_percent",
        "Minimum_percent", "Maximum_percent"
    ])

    for group_name, records in groups.items():
        for i, component in enumerate(components):
            vals = [r[1][i] for r in records]
            writer.writerow([
                group_name,
                component,
                len(vals),
                round(statistics.mean(vals), 4),
                round(statistics.median(vals), 4),
                round(statistics.stdev(vals), 4) if len(vals) > 1 else 0,
                round(min(vals), 4),
                round(max(vals), 4)
            ])

print("Summary saved to", output)
print("Individuals included:")
for name, records in groups.items():
    print(name + ":", len(records))
PYEOF

#####################################[ STEP 21 ]########################################
# Type: Python (embedded heredoc)
# Global ancestry median and IQR - Overall / PD / Control.
# Requires: ${PROJECT_DIR}/cohort_ancestry.tsv. Output: ${PROJECT_DIR}/global_ancestry_median_iqr.tsv.
##########################################################################################

python3 << 'PYEOF'
import csv
import math
import os
import statistics

PROJECT_DIR = os.environ["PROJECT_DIR"]

components = ["NAMA", "EUR", "AFR", "EAS", "MAL"]
data = []

def percentile(values, p):
    values = sorted(values)
    position = (len(values) - 1) * p
    lower = math.floor(position)
    upper = math.ceil(position)

    if lower == upper:
        return values[lower]

    return values[lower] + (values[upper] - values[lower]) * (position - lower)

with open(f"{PROJECT_DIR}/cohort_ancestry.tsv", newline="") as f:
    reader = csv.reader(f, delimiter="\t")
    for row in reader:
        if len(row) != 7:
            continue

        diagnosis = row[1]

        try:
            values = [float(x) * 100 for x in row[2:7]]
        except ValueError:
            continue

        data.append((diagnosis, values))

groups = {
    "Overall": data,
    "PD": [r for r in data if r[0] == "PD"],
    "Control": [r for r in data if r[0] == "Control"]
}

with open(f"{PROJECT_DIR}/global_ancestry_median_iqr.tsv", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow([
        "Group", "Ancestry", "N", "Median_percent",
        "Q1_percent", "Q3_percent", "IQR_percent"
    ])

    for group_name, records in groups.items():
        for i, component in enumerate(components):
            vals = [r[1][i] for r in records]

            q1 = percentile(vals, 0.25)
            median = statistics.median(vals)
            q3 = percentile(vals, 0.75)

            writer.writerow([
                group_name,
                component,
                len(vals),
                round(median, 4),
                round(q1, 4),
                round(q3, 4),
                round(q3 - q1, 4)
            ])

print("Summary saved to global_ancestry_median_iqr.tsv")
for name, records in groups.items():
    print(name + ":", len(records))
PYEOF

echo "PIPELINE COMPLETE"
