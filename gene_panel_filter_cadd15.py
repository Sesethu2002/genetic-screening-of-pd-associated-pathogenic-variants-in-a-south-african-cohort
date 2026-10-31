import pandas as pd

WORK_DIR = "/home/snondzasa/pd_qc_project/annovar_output_geno01/gene_extracts"
OUT_DIR = "/home/snondzasa/pd_qc_project/annovar_output_geno01/final_candidates_cadd15"
import os
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
