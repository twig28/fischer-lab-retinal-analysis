#!/usr/bin/env Rscript
# Exports SCENIC+ / pycisTopic-ready files from a combined Seurat v5 object.
# Also exports a combined UMAP (coords + selected metadata) for use in Python.

# ---------- USER CONFIG ----------
combined_rds <- "combined.rds"
outdir <- "scenicplus_input"
rna_assay_name  <- "RNA"
atac_assay_pref <- c("peaks","ATAC")
join_layers_on_export <- TRUE
export_per_treatment_lists <- TRUE                     # writes barcode lists per treatment
treatment_meta_keys <- c("treatment","condition","group") # first present will be used
umap_reduction_try <- c("wnn.umap","umap","ref.umap") # try in order
create_h5ad_via_reticulate <- FALSE                   # optional small RNA h5ad
species_ensembl_dataset <- "ggallus_gene_ensembl"     # Gallus gallus

# ---------- LIBS ----------
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(data.table)
  library(GenomicRanges)
  library(jsonlite)
  library(biomaRt)
})

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(combined_rds)) stop("combined.rds not found: ", combined_rds)
obj <- readRDS(combined_rds)

# Coerce factor meta->character
obj@meta.data[] <- lapply(obj@meta.data, function(x) if (is.factor(x)) as.character(x) else x)

# ---------- helpers ----------
get_best_layer_matrix <- function(seurat_obj, assay_name) {
  if (!assay_name %in% Assays(seurat_obj)) stop("Assay not present: ", assay_name)
  DefaultAssay(seurat_obj) <- assay_name
  layers <- Layers(seurat_obj[[assay_name]])
  if (length(layers) == 0) {
    # Seurat v5 Assay may still expose counts/data without named layers
    mat <- GetAssayData(seurat_obj, assay = assay_name, layer = NULL)
    if (!inherits(mat, "dgCMatrix")) mat <- as(mat, "dgCMatrix")
    return(list(mat = mat, layer = NA_character_))
  }
  # score layers
  score <- function(ln) {
    m <- tryCatch(GetAssayData(seurat_obj, assay = assay_name, layer = ln), error=function(e) NULL)
    if (is.null(m)) return(c(nnz=NA, total=NA))
    if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
    c(nnz=length(m@x), total=sum(m@x))
  }
  sums <- do.call(rbind, lapply(layers, score))
  rownames(sums) <- layers
  chosen <- if ("counts" %in% layers) {
    "counts"
  } else {
    cand <- grep("^counts", layers, value=TRUE)
    if (length(cand)) cand[which.max(sums[cand,"total"])] else if ("data" %in% layers) "data" else layers[which.max(sums[,"total"])]
  }
  m <- GetAssayData(seurat_obj, assay = assay_name, layer = chosen)
  if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
  list(mat=m, layer=chosen)
}

# ---------- optional layer join ----------
if (join_layers_on_export) {
  # RNA
  if (rna_assay_name %in% Assays(obj)) {
    suppressWarnings({
      tmp <- try(JoinLayers(obj, assay = rna_assay_name), silent=TRUE)
      if (!inherits(tmp, "try-error")) obj <- tmp
    })
  }
  # ATAC
  atac_assay_name <- atac_assay_pref[which(atac_assay_pref %in% Assays(obj))[1]]
  if (is.na(atac_assay_name)) atac_assay_name <- DefaultAssay(obj)
  suppressWarnings({
    tmp <- try(JoinLayers(obj, assay = atac_assay_name), silent=TRUE)
    if (!inherits(tmp, "try-error")) obj <- tmp
  })
}

# ---------- export RNA ----------
stopifnot(rna_assay_name %in% Assays(obj))
rna_res <- get_best_layer_matrix(obj, rna_assay_name)
rna_counts <- rna_res$mat
# features/genes
rna_features <- data.table(gene_id = rownames(rna_counts), gene_name = rownames(rna_counts))
# barcodes + metadata
rna_barcodes <- data.table(barcode = colnames(rna_counts))
rna_meta <- as.data.table(obj@meta.data, keep.rownames = "cell")

Matrix::writeMM(rna_counts, file.path(outdir, "rna_matrix.mtx"))
fwrite(rna_features, file.path(outdir, "rna_features.tsv"), sep="\t", col.names=FALSE)
fwrite(rna_barcodes, file.path(outdir, "rna_barcodes.tsv"), sep="\t", col.names=FALSE)
fwrite(rna_meta, file.path(outdir, "rna_cell_metadata.csv"))

# ---------- export ATAC ----------
# pick assay
atac_assay_name <- atac_assay_pref[which(atac_assay_pref %in% Assays(obj))[1]]
if (is.na(atac_assay_name)) atac_assay_name <- DefaultAssay(obj)
stopifnot(atac_assay_name %in% Assays(obj))

atac_res <- get_best_layer_matrix(obj, atac_assay_name)
atac_counts <- atac_res$mat

# peaks (prefer ranges from assay)
assay_obj <- obj@assays[[atac_assay_name]]
ranges_obj <- NULL
if (!is.null(assay_obj@ranges))    ranges_obj <- assay_obj@ranges
if (is.null(ranges_obj) && !is.null(assay_obj@rowRanges)) ranges_obj <- assay_obj@rowRanges

peak_names <- rownames(atac_counts)
if (is.null(peak_names)) stop("ATAC assay has no rownames (peaks)")

if (!is.null(ranges_obj) && length(ranges_obj) == length(peak_names)) {
  gr <- ranges_obj
  seq_chr <- as.character(seqnames(gr))
  seq_chr <- ifelse(grepl("^chr", seq_chr), seq_chr, paste0("chr", seq_chr))
  peaks_dt <- data.table(chr = seq_chr,
                         start = as.integer(start(gr) - 1L),
                         end   = as.integer(end(gr)),
                         name  = if (!is.null(names(gr))) names(gr) else peak_names)
} else {
  # fallback parser for rownames like 'chr1:100-200' / '1:100-200' / '1-100-200'
  parse_peak <- function(p) {
    p <- as.character(p)
    m <- regexpr("(chr?[A-Za-z0-9_.]+:\\d+-\\d+|[0-9XYM]+-\\d+-\\d+)", p, perl=TRUE)
    if (m[1] == -1) return(c(NA,NA,NA))
    tok <- regmatches(p, m)
    if (grepl(":", tok)) {
      chr <- sub(":.*","",tok); if (!grepl("^chr", chr)) chr <- paste0("chr", chr)
      coords <- sub(".*:","",tok); st <- sub("-.*","",coords); en <- sub(".*-","",coords)
      return(c(chr, st, en))
    } else {
      # 1-100-200
      parts <- strsplit(tok,"-",fixed=TRUE)[[1]]
      return(c(paste0("chr",parts[1]), parts[2], parts[3]))
    }
  }
  peaks_mat <- t(vapply(peak_names, parse_peak, FUN.VALUE=character(3)))
  peaks_dt <- data.table(chr=peaks_mat[,1], start=as.integer(peaks_mat[,2]), end=as.integer(peaks_mat[,3]), name=peak_names)
}

valid <- with(peaks_dt, !is.na(chr) & !is.na(start) & !is.na(end))
if (!all(valid)) {
  fwrite(peaks_dt, file.path(outdir, "peaks_with_names.tsv"), sep="\t")
  warning(sum(!valid), " peaks could not be parsed; saved to peaks_with_names.tsv for inspection.")
}
fwrite(peaks_dt[valid, .(chr,start,end)], file.path(outdir, "peaks.bed"), sep="\t", col.names=FALSE)

# matrices + barcodes + metadata for ATAC
Matrix::writeMM(atac_counts, file.path(outdir, "atac_matrix.mtx"))
fwrite(data.table(barcode = colnames(atac_counts)), file.path(outdir, "atac_barcodes.tsv"), sep="\t", col.names=FALSE)
fwrite(as.data.table(obj@meta.data, keep.rownames="cell"), file.path(outdir, "atac_cell_metadata.csv"))

# ---------- fragments discovery ----------
frag_paths <- character(0)
# try Fragments() on assay
try({
  fr <- NULL
  try({ fr <- Fragments(obj[[atac_assay_name]]) }, silent=TRUE)
  if (is.null(fr)) try({ fr <- Fragments(obj) }, silent=TRUE)
  add_path <- function(x) {
    p <- tryCatch({ if (!is.null(x@path)) x@path else NULL }, error=function(e) NULL)
    if (!is.null(p)) assign("frag_paths", unique(c(get("frag_paths", inherits=TRUE), p)), inherits=TRUE)
  }
  if (!is.null(fr)) {
    if (is.list(fr)) for (i in seq_along(fr)) add_path(fr[[i]]) else add_path(fr)
  }
}, silent=TRUE)
# fallback to assay slot
try({
  assay_frag <- obj@assays[[atac_assay_name]]@fragments
  add_path <- function(x) {
    p <- tryCatch({ if (!is.null(x@path)) x@path else NULL }, error=function(e) NULL)
    if (!is.null(p)) assign("frag_paths", unique(c(get("frag_paths", inherits=TRUE), p)), inherits=TRUE)
  }
  if (!is.null(assay_frag)) {
    if (is.list(assay_frag)) for (i in seq_along(assay_frag)) add_path(assay_frag[[i]]) else add_path(assay_frag)
  }
}, silent=TRUE)

write(toJSON(list(fragments = frag_paths), pretty=TRUE, auto_unbox=TRUE),
      file = file.path(outdir, "fragments_paths.json"))

# ---------- genes (convenience) ----------
genes_bed_path <- file.path(outdir, "genes.bed")
genes_table_path <- file.path(outdir, "genes_ensembl_table.tsv")
try({
  mart <- tryCatch(useEnsembl(biomart="genes", dataset=species_ensembl_dataset),
                   error=function(e) useMart("ENSEMBL_MART_ENSEMBL", dataset=species_ensembl_dataset))
  attrs <- c("external_gene_name","chromosome_name","start_position","end_position","strand","ensembl_gene_id")
  genes_res <- getBM(attributes=attrs, mart=mart)
  if (nrow(genes_res) > 0) {
    genes_bed <- data.table(
      chr = paste0("chr", genes_res$chromosome_name),
      start = pmax(0, genes_res$start_position - 1L),
      end = genes_res$end_position,
      gene_name = genes_res$external_gene_name,
      ensembl_id = genes_res$ensembl_gene_id
    )
    fwrite(genes_bed, genes_bed_path, sep="\t", col.names=FALSE)
    fwrite(genes_res, genes_table_path, sep="\t")
  }
}, silent=TRUE)

# ---------- common barcodes ----------
common_bcs <- intersect(colnames(rna_counts), colnames(atac_counts))
if (length(common_bcs) > 0) fwrite(data.table(barcode=common_bcs), file.path(outdir, "common_barcodes.tsv"), sep="\t", col.names=FALSE)

# ---------- per-treatment barcode lists ----------
if (export_per_treatment_lists) {
  key <- treatment_meta_keys[treatment_meta_keys %in% colnames(obj@meta.data)][1]
  if (!is.na(key)) {
    dt <- data.table(barcode = colnames(obj), treatment = obj@meta.data[[key]][colnames(obj)])
    dt <- unique(na.omit(dt))
    for (t in sort(unique(dt$treatment))) {
      fwrite(dt[treatment==t, .(barcode)], file.path(outdir, sprintf("barcodes_treatment_%s.tsv", gsub("[^A-Za-z0-9_.-]","_",t))), sep="\t", col.names=FALSE)
    }
  }
}

# ---------- export combined UMAP ----------
umap_name <- umap_reduction_try[umap_reduction_try %in% Reductions(obj)][1]
if (!is.na(umap_name)) {
  emb <- Embeddings(obj, reduction = umap_name)
  umap_dt <- data.table(cell = rownames(emb), UMAP_1 = emb[,1], UMAP_2 = emb[,2])
  # add key metadata columns commonly useful downstream
  add_cols <- intersect(c("treatment","condition","group","seurat_clusters","celltype","cell_type","cluster","sample"), colnames(obj@meta.data))
  for (ac in add_cols) umap_dt[[ac]] <- obj@meta.data[umap_dt$cell, ac, drop=TRUE]
  fwrite(umap_dt, file.path(outdir, "umap_combined.csv"))
} else {
  warning("No UMAP reduction found among: ", paste(umap_reduction_try, collapse=", "))
}

# ---------- optional: quick RNA h5ad ----------
if (create_h5ad_via_reticulate) {
  if (requireNamespace("reticulate", quietly=TRUE)) {
    library(reticulate)
    need <- c("anndata","scipy","numpy","pandas")
    avail <- sapply(need, py_module_available)
    if (!all(avail)) try(py_install(need[!avail], pip=TRUE), silent=TRUE)
    if (all(sapply(need, py_module_available))) {
      anndata <- import("anndata", convert=FALSE)
      scipy   <- import("scipy.sparse", convert=FALSE)
      np      <- import("numpy", convert=FALSE)
      r_mat <- rna_counts
      data_vals <- as.numeric(r_mat@x)
      indices <- as.integer(r_mat@i)
      indptr  <- as.integer(r_mat@p)
      shape   <- as.integer(c(nrow(r_mat), ncol(r_mat)))
      csr <- scipy$csr_matrix(tuple(list(np$array(data_vals),
                                         np$array(indices, dtype='int64'),
                                         np$array(indptr,  dtype='int64'))),
                              shape = tuple(shape))
      obs_df <- r_to_py(as.data.frame(obj@meta.data[ colnames(r_mat), , drop=FALSE ]), convert=TRUE)
      var_df <- r_to_py(data.frame(gene_id=rownames(r_mat), gene_name=rownames(r_mat)), convert=TRUE)
      ad <- anndata$AnnData(X=csr, obs=obs_df, var=var_df)
      ad$write_h5ad(file.path(outdir, "rna_for_scenicplus.h5ad"))
    }
  }
}

cat("Done. Files in:", outdir, "\n")
