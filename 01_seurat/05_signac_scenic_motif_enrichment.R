#!/usr/bin/env Rscript
# signac_scenic_motif_enrichment.R
# =============================================================================
# Standalone Option 4: motif enrichment using JASPAR PWMs filtered to only
# TFs present in the SCENIC+ chicken cisTarget motif annotation table.
#
# This reconciles the Signac and SCENIC+ TF vocabularies so results are
# directly comparable.
#
# Inputs (edit paths below):
#   atac_with_motifs_and_stats.rds  — ATAC Seurat object from coaccessible_motif.Rmd
#   da_peaks.csv                    — DA peaks from combined_differential.Rmd
#   scenicplus_input/v10nr_clust_public/snapshots/motifs-v10-nr.chicken-m0.00001-o0.0.tbl
#
# Outputs:
#   enriched_motifs_scenic_reconciled.csv   — motif enrichment, TFs matched to SCENIC+
#   scenic_tf_coverage.csv                  — which SCENIC+ TFs got JASPAR coverage
#   motif_comparison_summary.csv            — side-by-side: old human-only vs new reconciled
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(JASPAR2020)
  library(TFBSTools)
  library(BSgenome.Ggallus.UCSC.galGal6)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(data.table)
  library(dplyr)
})

# ---------- USER CONFIG ----------
atac_rds        <- "atac_with_motifs_and_stats.rds"
da_peaks_csv    <- "da_peaks.csv"
motif_tbl_path  <- "scenicplus_input/v10nr_clust_public/snapshots/motifs-v10-nr.chicken-m0.00001-o0.0.tbl"
old_enriched_csv <- "enriched_motifs_sorted.csv"   # existing human-JASPAR result to compare against
outdir          <- "."

# DA peak filter thresholds (must match coaccessible_motif.Rmd)
DA_PADJ_THRESH  <- 0.05
DA_FC_THRESH    <- 0.5

# Vertebrate JASPAR species to search for PWMs
# (human, mouse, zebrafish — broadens coverage for conserved TF families)
JASPAR_SPECIES  <- c(9606L, 10090L, 7955L)
# ---------- END CONFIG ----------

cat("=================================================================\n")
cat("  SCENIC+-reconciled motif enrichment\n")
cat("=================================================================\n\n")

# ── 1. Load SCENIC+ motif annotation table ────────────────────────
cat("[1] Loading SCENIC+ motif annotation table...\n")
motif_tbl <- fread(motif_tbl_path, sep = "\t", data.table = FALSE)
colnames(motif_tbl)[1] <- "motif_id"   # rename #motif_id if leading # present

# Identify gene_name column
gene_col <- grep("^gene_name$|^TF$|^name$", colnames(motif_tbl),
                 ignore.case = TRUE, value = TRUE)[1]
if (is.na(gene_col)) {
  # Fallback: second non-id column
  gene_col <- colnames(motif_tbl)[2]
  warning("Could not identify gene_name column; using: ", gene_col)
}
cat("  Using TF column: '", gene_col, "'\n", sep = "")

scenic_tfs <- unique(na.omit(motif_tbl[[gene_col]]))
scenic_tfs <- scenic_tfs[nchar(trimws(scenic_tfs)) > 0]
cat("  Unique TF names in SCENIC+ motif table:", length(scenic_tfs), "\n")
cat("  Sample:", paste(head(scenic_tfs, 8), collapse = ", "), "\n\n")

# ── 2. Fetch JASPAR PWMs for those TFs ───────────────────────────
# getMatrixSet only accepts a single integer for species — loop and merge
cat("[2] Fetching JASPAR2020 vertebrate PWMs matching SCENIC+ TFs...\n")
pfm_all <- do.call(c, lapply(JASPAR_SPECIES, function(sp) {
  cat("  Fetching species", sp, "...\n")
  getMatrixSet(JASPAR2020, opts = list(species = sp, all_versions = FALSE))
}))
# Remove duplicates by matrix ID (same motif may appear under multiple species)
pfm_ids <- sapply(pfm_all, function(m) ID(m))
pfm_all <- pfm_all[!duplicated(pfm_ids)]
cat("  Total vertebrate PWMs after deduplication:", length(pfm_all), "\n")

# Match by TF name (case-insensitive)
scenic_tfs_upper <- toupper(scenic_tfs)
keep_mask <- sapply(pfm_all, function(m) toupper(name(m)) %in% scenic_tfs_upper)
pfm_scenic <- pfm_all[keep_mask]
cat("  PWMs matching SCENIC+ TFs:", length(pfm_scenic), "\n")

# Coverage report: which SCENIC+ TFs got a JASPAR match?
matched_tfs   <- toupper(sapply(pfm_scenic, name))
unmatched_tfs <- scenic_tfs[!scenic_tfs_upper %in% matched_tfs]
cat("  SCENIC+ TFs with JASPAR coverage:", sum(scenic_tfs_upper %in% matched_tfs),
    "of", length(scenic_tfs), "\n")
cat("  SCENIC+ TFs without JASPAR PWM:", length(unmatched_tfs), "\n")
if (length(unmatched_tfs) > 0)
  cat("  Sample unmatched:", paste(head(unmatched_tfs, 10), collapse = ", "), "\n")
cat("\n")

# Save coverage summary
coverage_df <- data.frame(
  TF            = scenic_tfs,
  in_JASPAR     = scenic_tfs_upper %in% matched_tfs,
  n_PWMs        = sapply(scenic_tfs_upper, function(tf)
    sum(sapply(pfm_scenic, function(m) toupper(name(m)) == tf)))
)
write.csv(coverage_df,
          file.path(outdir, "scenic_tf_coverage.csv"),
          row.names = FALSE)
cat("  Saved: scenic_tf_coverage.csv\n\n")

# ── 3. Load ATAC object ───────────────────────────────────────────
cat("[3] Loading ATAC object...\n")
atac <- readRDS(atac_rds)
DefaultAssay(atac) <- "peaks"
cat("  Cells:", ncol(atac), " | Peaks:", nrow(atac), "\n\n")

# ── 4. Drop old GC/length columns and recompute RegionStats ───────
cat("[4] Recomputing RegionStats for GC-matched background...\n")
mf <- atac@assays$peaks@meta.features
mf <- mf[, !colnames(mf) %in% c("GC.percent", "sequence.length",
                                "gc", "bias", "region.width"), drop = FALSE]
atac@assays$peaks@meta.features <- mf

atac <- RegionStats(
  object = atac,
  assay  = "peaks",
  genome = BSgenome.Ggallus.UCSC.galGal6
)

mf <- atac@assays$peaks@meta.features
mf$gc           <- mf$GC.percent
mf$bias         <- mf$sequence.length
mf$region.width <- mf$sequence.length
atac@assays$peaks@meta.features <- mf
cat("  RegionStats complete.\n\n")

# ── 5. Add SCENIC+-reconciled motifs ─────────────────────────────
cat("[5] Running AddMotifs with SCENIC+-reconciled PWMs...\n")
atac <- AddMotifs(
  object = atac,
  genome = BSgenome.Ggallus.UCSC.galGal6,
  pfm    = pfm_scenic
)
cat("  Motifs added:", nrow(atac@assays$peaks@motifs@data), "peak-motif entries\n\n")

# ── 6. Load DA peaks and build foreground set ─────────────────────
cat("[6] Building DA peak foreground set...\n")
da_peaks <- read.csv(da_peaks_csv, row.names = 1)
filtered_da <- da_peaks[
  !is.na(da_peaks$p_val_adj) &
    da_peaks$p_val_adj < DA_PADJ_THRESH &
    da_peaks$avg_log2FC > DA_FC_THRESH, ]
cat("  DA peaks (padj<", DA_PADJ_THRESH, ", FC>", DA_FC_THRESH, "):",
    nrow(filtered_da), "\n")

peak_names <- rownames(filtered_da)

# Handle chr prefix: try both with and without, keep whichever intersects
fg_with_chr    <- intersect(paste0("chr", peak_names), rownames(atac))
fg_without_chr <- intersect(peak_names, rownames(atac))

if (length(fg_with_chr) >= length(fg_without_chr)) {
  foreground_peaks <- fg_with_chr
  cat("  Using 'chr' prefix format\n")
} else {
  foreground_peaks <- fg_without_chr
  cat("  Using bare chromosome format (no 'chr' prefix)\n")
}
cat("  Foreground peaks matched to ATAC object:", length(foreground_peaks), "\n")

if (length(foreground_peaks) == 0)
  stop("Zero foreground peaks matched — check chr prefix format in da_peaks.csv vs rownames(atac)")
cat("\n")

# ── 7. Run FindMotifs ─────────────────────────────────────────────
cat("[7] Running FindMotifs (GC-matched background)...\n")
enriched_scenic <- FindMotifs(
  object   = atac,
  features = foreground_peaks
)
cat("  Motifs tested:", nrow(enriched_scenic), "\n")

# ── 8. Annotate with SCENIC+ TF names ───────────────────────────
cat("[8] Annotating with SCENIC+ TF names...\n")

# Build motif-ID → SCENIC+ TF map from the annotation table
motif_to_scenic_tf <- motif_tbl[, c("motif_id", gene_col), drop = FALSE]
motif_to_scenic_tf <- motif_to_scenic_tf[!is.na(motif_to_scenic_tf[[gene_col]]) &
                                           nchar(trimws(motif_to_scenic_tf[[gene_col]])) > 0, ]
colnames(motif_to_scenic_tf) <- c("motif_id", "scenic_TF")
motif_to_scenic_tf <- unique(motif_to_scenic_tf)

# JASPAR motif names map: motif.name in enriched_scenic = JASPAR matrix name (e.g. "ELK3")
# Map through: JASPAR name → scenic_TF (case-insensitive join on TF name)
enriched_scenic$jaspar_TF <- enriched_scenic$motif.name

# Direct TF name join (case-insensitive)
jaspar_upper    <- toupper(enriched_scenic$jaspar_TF)
scenic_tf_upper <- toupper(motif_to_scenic_tf$scenic_TF)

enriched_scenic$scenic_TF <- motif_to_scenic_tf$scenic_TF[
  match(jaspar_upper, scenic_tf_upper)
]

# Flag whether this TF is in the SCENIC+ motif universe
enriched_scenic$in_scenic_db <- !is.na(enriched_scenic$scenic_TF)

cat("  Enriched motifs with SCENIC+ TF match:",
    sum(enriched_scenic$in_scenic_db), "of", nrow(enriched_scenic), "\n\n")

# Sort by p-value
enriched_scenic <- enriched_scenic[order(enriched_scenic$pvalue), ]

# ── 9. Save reconciled enrichment ────────────────────────────────
out_reconciled <- file.path(outdir, "enriched_motifs_scenic_reconciled.csv")
write.csv(enriched_scenic, out_reconciled, row.names = FALSE)
cat("[9] Saved:", out_reconciled, "\n")

# ── 10. Side-by-side comparison with old human-JASPAR result ─────
cat("[10] Building comparison with original human-only enrichment...\n")
if (file.exists(old_enriched_csv)) {
  old_enriched <- read.csv(old_enriched_csv)
  
  cat("  Old enrichment columns:", paste(colnames(old_enriched), collapse = ", "), "\n")
  cat("  New enrichment columns:", paste(colnames(enriched_scenic), collapse = ", "), "\n")
  
  # motif.name must be present in both
  if (!"motif.name" %in% colnames(old_enriched))
    stop("old enrichment CSV has no 'motif.name' column; columns are: ",
         paste(colnames(old_enriched), collapse = ", "))
  
  # Defensively select only columns that actually exist in each dataframe
  old_keep <- intersect(c("motif.name", "pvalue", "motif.fc", "fold.enrichment"),
                        colnames(old_enriched))
  new_keep <- intersect(c("motif.name", "pvalue", "motif.fc", "fold.enrichment",
                          "in_scenic_db", "scenic_TF", "jaspar_TF"),
                        colnames(enriched_scenic))
  
  comp <- merge(
    old_enriched[, old_keep, drop = FALSE],
    enriched_scenic[, new_keep, drop = FALSE],
    by       = "motif.name",
    all      = TRUE,
    suffixes = c(".human_jaspar", ".scenic_reconciled")
  )
  comp <- comp[order(comp$pvalue.scenic_reconciled, na.last = TRUE), ]
  
  out_comp <- file.path(outdir, "motif_comparison_summary.csv")
  write.csv(comp, out_comp, row.names = FALSE)
  cat("  Saved:", out_comp, "\n")
  
  # Quick summary
  n_shared   <- sum(!is.na(comp$pvalue.human_jaspar) &
                      !is.na(comp$pvalue.scenic_reconciled))
  n_new_only <- sum(is.na(comp$pvalue.human_jaspar) &
                      !is.na(comp$pvalue.scenic_reconciled))
  n_old_only <- sum(!is.na(comp$pvalue.human_jaspar) &
                      is.na(comp$pvalue.scenic_reconciled))
  
  cat("\n  Comparison summary:\n")
  cat("    Motifs in both analyses:", n_shared, "\n")
  cat("    New (SCENIC+-reconciled only):", n_new_only, "\n")
  cat("    Dropped (human-JASPAR only, not in SCENIC+):", n_old_only, "\n")
  
  # Top motifs significant in both
  sig_both <- comp[
    !is.na(comp$pvalue.human_jaspar) & comp$pvalue.human_jaspar < 0.05 &
      !is.na(comp$pvalue.scenic_reconciled) & comp$pvalue.scenic_reconciled < 0.05, ]
  if (nrow(sig_both) > 0) {
    cat("\n  Top 10 motifs significant in BOTH analyses (by SCENIC+-reconciled p-value):\n")
    print(head(sig_both[order(sig_both$pvalue.scenic_reconciled), ], 10))
  }
  
  # Motifs significant only in reconciled (new findings)
  # isTRUE(comp$in_scenic_db) handles NA safely — avoids subsetting errors
  # when in_scenic_db is missing or NA for rows that came from the left join
  in_db_flag <- if ("in_scenic_db" %in% colnames(comp)) {
    isTRUE(comp$in_scenic_db) | (is.logical(comp$in_scenic_db) & !is.na(comp$in_scenic_db) & comp$in_scenic_db)
  } else {
    rep(TRUE, nrow(comp))   # column absent: don't filter on it
  }
  sig_new <- comp[
    (is.na(comp$pvalue.human_jaspar) | comp$pvalue.human_jaspar >= 0.05) &
      !is.na(comp$pvalue.scenic_reconciled) & comp$pvalue.scenic_reconciled < 0.05 &
      in_db_flag, ]
  if (nrow(sig_new) > 0) {
    cat("\n  Motifs significant ONLY in SCENIC+-reconciled (new, in SCENIC+ DB):\n")
    print(head(sig_new[order(sig_new$pvalue.scenic_reconciled), ], 10))
  }
} else {
  cat("  Old enrichment file not found at '", old_enriched_csv,
      "' — skipping comparison.\n", sep = "")
}

cat("\n=================================================================\n")
cat("  Done. Output files:\n")
cat("    enriched_motifs_scenic_reconciled.csv\n")
cat("    scenic_tf_coverage.csv\n")
if (file.exists(old_enriched_csv))
  cat("    motif_comparison_summary.csv\n")
cat("=================================================================\n")


rna <- readRDS("rna_prePlot.rds")
rna <- JoinLayers(rna)
avg <- AverageExpression(rna, group.by = "celltype", assay = "RNA")$RNA
write.csv(avg, "celltype_mean_expression.csv")


library(AnnotationHub)
ah <- AnnotationHub()
gg_ensdb <- ah[["AH113850"]]  # or whichever Ensembl 113 GRCg7b entry you used
gene_info <- genes(gg_ensdb, filter = GeneIdFilter("ENSGALG00010015567"))
# or use biomaRt:
library(biomaRt)
# Try human ortholog lookup specifically
getBM(
  attributes = c("ensembl_gene_id", "external_gene_name",
                 "hsapiens_homolog_ensembl_gene",
                 "hsapiens_homolog_associated_gene_name",
                 "hsapiens_homolog_orthology_type",
                 "hsapiens_homolog_perc_id"),
  filters    = "ensembl_gene_id",
  values     = c("ENSGALG00010015567","ENSGALG00010002949",
                 "ENSGALG00010006080","ENSGALG00010020018"),
  mart       = ensembl
)
