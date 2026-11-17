import os
import pandas as pd

from pybiomart import Dataset
from pyfaidx import Fasta

INDIR = "scenicplus_input"
os.makedirs(INDIR, exist_ok=True)

# ---------------------------------------------------------
# Locate the Gallus_gallus.GRCg6a.dna.toplevel.fa FASTA
# ---------------------------------------------------------
fa_path = os.path.join(INDIR, "Gallus_gallus.GRCg6a.dna.toplevel.fa")

genome_fa = None

if os.path.isfile(fa_path):
    genome_fa = fa_path
elif os.path.isdir(fa_path):
    candidate = os.path.join(fa_path, "Gallus_gallus.GRCg6a.dna.toplevel.fa")
    if os.path.isfile(candidate):
        genome_fa = candidate

if genome_fa is None:
    raise FileNotFoundError(
        "Could not find Gallus_gallus.GRCg6a.dna.toplevel.fa.\n"
        f"Tried:\n  {fa_path}\n  {os.path.join(fa_path, 'Gallus_gallus.GRCg6a.dna.toplevel.fa')}"
    )

print(f"[make_genome_annots] Using genome FASTA: {genome_fa}")


# ---------------------------------------------------------
# 1) genome_annotation.tsv via Ensembl (pybiomart)
# ---------------------------------------------------------
print("[make_genome_annots] Fetching gene annotation from Ensembl (ggallus_gene_ensembl)...")

dataset = Dataset(
    name="ggallus_gene_ensembl",
    host="http://www.ensembl.org",
)

# NOTE: no duplicate attributes — the old script had "transcription_start_site" twice
# which caused pybiomart to emit a suffixed column name and break the rename dict.
df = dataset.query(
    attributes=[
        "ensembl_gene_id",
        "chromosome_name",
        "start_position",
        "end_position",
        "strand",
        "external_gene_name",
        "transcription_start_site",
        "gene_biotype",
    ]
)

print("[make_genome_annots] Raw columns returned by BioMart:", df.columns.tolist())

# ---------------------------------------------------------
# Robust column rename: detect actual column names returned
# by pybiomart (they vary slightly between Ensembl versions)
# ---------------------------------------------------------
col_map = {}
for col in df.columns:
    cl = col.lower().replace(" ", "_").replace("(", "").replace(")", "").replace("/", "_")
    if "ensembl_gene_id" in cl or "gene_stable_id" in cl or "stable_id" in cl:
        col_map[col] = "ensembl_gene_id"
    elif "chromosome" in cl or "scaffold" in cl:
        col_map[col] = "Chromosome"
    elif "start_position" in cl or "gene_start" in cl:
        col_map[col] = "Start"
    elif "end_position" in cl or "gene_end" in cl:
        col_map[col] = "End"
    elif cl == "strand":
        col_map[col] = "Strand"
    elif "external_gene_name" in cl or "gene_name" in cl:
        col_map[col] = "Gene"
    elif "transcription_start_site" in cl or "tss" in cl:
        # Only map the FIRST TSS column; ignore any _1 suffix duplicates
        if "Transcription_Start_Site" not in col_map.values():
            col_map[col] = "Transcription_Start_Site"
    elif "biotype" in cl or "gene_type" in cl:
        col_map[col] = "Transcript_type"

df = df.rename(columns=col_map)

# Drop any leftover duplicate TSS columns (suffixed by pybiomart)
tss_dupes = [c for c in df.columns if "transcription_start_site" in c.lower() and c != "Transcription_Start_Site"]
if tss_dupes:
    df = df.drop(columns=tss_dupes)
    print(f"[make_genome_annots] Dropped duplicate TSS columns: {tss_dupes}")

print("[make_genome_annots] Columns after rename:", df.columns.tolist())

needed_cols = [
    "Chromosome",
    "Start",
    "End",
    "Strand",
    "Gene",
    "ensembl_gene_id",
    "Transcription_Start_Site",
    "Transcript_type",
]

missing = [c for c in needed_cols if c not in df.columns]
if missing:
    raise RuntimeError(
        f"After renaming, these required columns are missing: {missing}\n"
        f"Available columns: {df.columns.tolist()}"
    )

genome_annot = df[needed_cols].copy()

# Fill blank Gene symbols with the Ensembl ID so SCENIC+ always has a non-empty Gene column
blank_gene = genome_annot["Gene"].isna() | (genome_annot["Gene"].astype(str).str.strip() == "")
genome_annot.loc[blank_gene, "Gene"] = genome_annot.loc[blank_gene, "ensembl_gene_id"]
print(f"[make_genome_annots] Filled {blank_gene.sum()} blank Gene symbols with Ensembl IDs")

out_genome = os.path.join(INDIR, "genome_annotation.tsv")
genome_annot.to_csv(out_genome, sep="\t", index=False)
print(f"[make_genome_annots] Wrote genome_annotation.tsv with {len(genome_annot)} rows to {out_genome}")

# Quick sanity: how many unique gene symbols vs Ensembl IDs?
n_symbols = genome_annot["Gene"].nunique()
n_ensembl = genome_annot["ensembl_gene_id"].nunique()
print(f"[make_genome_annots] Unique gene symbols: {n_symbols}  |  Unique Ensembl IDs: {n_ensembl}")


# ---------------------------------------------------------
# 2) chromsizes.tsv from FASTA
# ---------------------------------------------------------
print("[make_genome_annots] Computing chromsizes from FASTA...")
fa = Fasta(genome_fa)

rows = []
for chrom in fa.keys():
    length = len(fa[chrom])
    rows.append(
        {
            "Chromosome": chrom,
            "Start": 0,
            "End": int(length),
        }
    )

chromsizes = pd.DataFrame(rows)
out_chrom = os.path.join(INDIR, "chromsizes.tsv")
chromsizes.to_csv(out_chrom, sep="\t", index=False)
print(f"[make_genome_annots] Wrote chromsizes.tsv with {len(chromsizes)} chromosomes to {out_chrom}")

print("[make_genome_annots] Done.")