#!/bin/bash
########################################################################################
# PD PATHOGENIC VARIANT SCREENING - FULL PIPELINE (combined)
# South African cohort (GP2 Stellenbosch), 21 monogenic PD genes
# (Blauwendraat et al. 2020, Lancet Neurol, Table 1)
#
# This file concatenates all pipeline steps in order for reference/GitHub.
# Steps 1-5, 7, 9, 12 are bash; Steps 8, 10 and 11 are Python (embedded below via heredoc);
# Step 6 is a PBS array job - submit it SEPARATELY with `qsub 06_annovar_annotation.pbs`,
# do not run it inline with the rest.
#
# Must run inside an interactive/batch HPC session, not the login node.
########################################################################################
set -euo pipefail
START_DIR="$(pwd)"   # Step 12 returns here (Step 7 cd-s into annovar_output)

#####################################[ STEP 1 ]#########################################
# Strand/reference flip, BEFORE QC. Flipping first lets bcftools +fixref resolve
# ambiguous-strand variants before the QC "ATCG" step would otherwise discard them.
##########################################################################################

RAW_INPUT="GP2_NBA_raw_STELLENBOS"
REF_FASTA="/path/to/Homo_sapiens_assembly38.fasta"

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
# Build sex/status covariate file. Join key: IID/GP2ID (unique) - NOT FID (shared by
# relatives). Status source: Masterkey "GP2_phenotype_for_qc" column.
##########################################################################################

MASTERKEY="updated_R12_master_key_subset_STELLENBOS.txt"
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
# QC Run 1 - related individuals EXCLUDED (population structure only, not used downstream)
# HWE deliberately excluded from both runs. qc_steps.txt MUST be tab-delimited.
##########################################################################################

INPUT="${RAW_INPUT_BED}_ccOnly"
GWASQC_DIR="/path/to/GWASQC"
NATORA_DIR="/path/to/NAToRA_Public"

module load plink/1.90 python/3.11.9

printf "sex-check\nATCG\nmind\t0.05\ngeno\t0.01\nduplicate\nheterozygosity\t3\nrelationship\t0.0884\n" > qc_steps.txt

python "${GWASQC_DIR}/main.py" \
    -i "${INPUT}" -o nba_qced -O qced_data \
    -I info_file_ccOnly.tsv -s qc_steps.txt \
    --NAToRA "${NATORA_DIR}/NAToRA_Public.py" \
    2>&1 | tee run1_full_log.txt

#####################################[ STEP 4 ]#########################################
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
# ANNOVAR annotation of all 22 imputed chromosomes.
# *** SUBMIT SEPARATELY, DO NOT RUN INLINE: ***
#     qsub 06_annovar_annotation.pbs
# (PBS array job #PBS -J 1-22; see 06_annovar_annotation.pbs for the full job script -
#  protocol: refGene,clinvar_20250721,dbnsfp47a,gnomad41_genome,gnomad41_exome)
##########################################################################################

#####################################[ STEP 7 ]#########################################
# Fast per-gene extraction from the annotated chromosome files (avoids loading full
# multi-GB chromosome files into pandas directly - that took >8h and was killed twice)
##########################################################################################

cd /path/to/annovar_output
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
# Apply variant filters per gene (MAF<=1%, exonic, CADD>=20, REVEL>=0.7 for missense)
# and build the per-gene summary table. (Python - embedded via heredoc)
##########################################################################################

module load python/3.11.9
source ~/py311env/bin/activate   # venv with pandas installed

python3 << 'PYEOF'
import pandas as pd

WORK_DIR = "/path/to/annovar_output/gene_extracts"
OUT_DIR = "/path/to/annovar_output"

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
# Identify which individuals carry each final candidate variant (Python - embedded)
##########################################################################################

python3 << 'PYEOF'
import gzip, csv

ANNOVAR_DIR = "/path/to/annovar_output"
FINAL_DIR = f"{ANNOVAR_DIR}/final_candidates"
VCF_DIR = "/path/to/topmed_results"
MASTERKEY = "/path/to/updated_R12_master_key_subset_STELLENBOS.txt"
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
# Sensitivity check: repeat the candidate filter at CADD_phred >= 15
# (result: 29 candidates vs 28; added LRRK2 p.E899D). Python - embedded via heredoc.
##########################################################################################

python3 << 'PYEOF'
"""Sensitivity check: repeat the candidate filter at CADD_phred >= 15 (result: 29 candidates vs 28; added LRRK2 p.E899D)."""
import os
import pandas as pd

WORK_DIR = "/path/to/annovar_output/gene_extracts"
OUT_DIR = "/path/to/annovar_output/final_candidates_cadd15"
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
# Demographics for the QC Run 2 sample (n = 1,355).
# Masterkey columns (1-indexed): 3 GP2ID, 16 GP2_phenotype_for_qc,
# 19 sex, 21 age_of_onset, 22 age_at_diagnosis, 25 race_for_qc, 26 family_history_for_qc.
# Note: the final QC .psam has IID in column 1.
##########################################################################################

cd "${START_DIR}"   # QC outputs (.psam) are relative to the original working directory
MASTERKEY="/path/to/updated_R12_master_key_subset_STELLENBOS.txt"
PSAM="qced_data_pathogenic/nba_qced_pathogenic_Heterozygosity.psam"

awk -F'\t' 'NR==FNR { if (FNR==1) next; sex[$3]=$19; aao[$3]=$21; fh[$3]=$26; st[$3]=$16; race[$3]=$25; next }
    FNR==1 { next }
    ($1 in sex) { print $1"\t"sex[$1]"\t"aao[$1]"\t"fh[$1]"\t"st[$1]"\t"race[$1] }' \
    "${MASTERKEY}" "${PSAM}" > demographics_run2.tsv

wc -l demographics_run2.tsv                                  # expect 1355
echo "Sex by status:";  cut -f5,2 demographics_run2.tsv | sort | uniq -c
echo "Ancestry (Masterkey race_for_qc):"; cut -f6 demographics_run2.tsv | sort | uniq -c
echo "Family history among cases:"; awk -F'\t' '$5=="PD"{print $4}' demographics_run2.tsv | sort | uniq -c
echo "Age at onset, cases (n, mean, SD):"
awk -F'\t' '$5=="PD" && $3!="" && $3!="NA"{s+=$3; q+=$3*$3; n++} END{m=s/n; printf "%d, %.1f, %.1f\n", n, m, sqrt(q/n-m*m)}' demographics_run2.tsv

echo "PIPELINE COMPLETE"
