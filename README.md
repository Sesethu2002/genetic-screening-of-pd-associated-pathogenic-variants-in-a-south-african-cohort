**Genetic screening of Parkinson’s disease-associated pathogenic variants in a South African cohort.**

This research aims to identify and characterise candidate variants using the NeuroBooster array data, thereby addressing the underrepresentation of African populations in genetics studies.

This study utilised genotyping data from the South African PD Study Collection, comprising 691 individuals diagnosed with PD, 30 family members and 826 neurologically healthy controls. 

Overview 
QC -> Strand/reference correction -> TopMed imputation -> ANNOVAR annotation -> prioritisation (gnomAD FAF95_max ≤ 0.01, CADD ≥ 15, REVEL ≥ 0.7) -> Carrier identification.

Software and versions 
PLINK	 v1.9/ v2.00a3 (OC, file conversion)
bcftools (+fixref)	1.9 (htslib 1.9)	(Strand and REF correction)
samtools	version not recorded	(Reference indexing)
NAToRA	No version number: (files dated 22 July 2026; Leal et al., 2022)	(Relatedness)
TopMed Imputation Server	v2.0	(Imputation)
Eagle 	v2.4	(Phasing) 
ANNOVAR	2 March 2025	(Annotation) 
ClinVar	21 July 2025 release 	(Clinical classification)
gnomAD	v4.1	(Allele frequency (FAF95_max))
dbNSFP	v4.7a	(REVEL score)
CADD	v1.7 (via dbNSFP v4.7a; Schubach et al., 2024)	(Deleteriousness scores)
Python 	Python/3.11.9	(HPC) (Analysis and figures)


Data availability
Data can not be shared publicly due to it consisting of real participant information, individual-level genotype and clinical data. 

Running the pipeline 
Run the scripts in ‘scripts/‘ in numerical order on an CHPC/HPC session.
Edit the paths at the top of each script. 

Contact
Sesethu Nondzasa, 31083641@sun.ac.za

Acknowledgments
Prof. Petersen, Prof. Möller, Dr. Step, MJFF, ASAP, GP2, NIH, Stellenbosch University Department of Biomedical Sciences, MAGIC Lab and Parkinson’s Disease Research Group.
