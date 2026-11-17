# Load required libraries
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(pheatmap)
library(ggrepel)
library(Signac)
library(Seurat)
library(patchwork)
library(GenomicRanges)
library(future)
library(AnnotationHub)
library(BSgenome.Ggallus.UCSC.galGal6)
library(dplyr)
library(GenomeInfoDb)
library(JASPAR2020)
library(TFBSTools)
library(SeuratWrappers)
library(monocle3)
library(cicero)

#-----------------------------------------
# Cicero Co-accessibility Analysis on ATAC Object
#-----------------------------------------
#atac <- readRDS("atac_prePlot.rds")
#DefaultAssay(atac) <- "ATAC"
rna <-readRDS("rna_prePlot.rds")
DefaultAssay(rna) <- "RNA"
rna <- JoinLayers(rna)

#Function imported to get GRanges for a given gene
#————————————————————————————————————————————————————————
# getGeneGRanges_biomaRt.R
#————————————————————————————————————————————————————————

#’ Fetch a GRanges for a given chicken gene symbol (galGal6) via biomaRt
#’
#’ This uses Ensembl’s “ggallus_gene_ensembl” dataset to look up
#’ the chrom, start, end, strand for your gene, collapses across
#’ any transcripts, then optionally extends by upstream/downstream.
#’
#’ @param goi         Gene symbol (character), e.g. "EDA2R"
#’ @param upstream    Bases upstream to include (numeric, default 0)
#’ @param downstream  Bases downstream to include (numeric, default 0)
#’ @return A single GRanges object spanning the gene body ± flanks
#’ @importFrom GenomicRanges GRanges IRanges
#’ @importFrom biomaRt useEnsembl getBM
#’ @examples
#’ gr <- getGeneGRanges("EDA2R", upstream=1000, downstream=1000)
#’ CoveragePlot(..., region=gr, ...)
#————————————————————————————————————————————————————————
getGeneGRanges <- function(goi, upstream=0, downstream=0) {
  if (!requireNamespace("biomaRt", quietly=TRUE)) {
    stop("Please install the biomaRt package: BiocManager::install('biomaRt')")
  }
  library(biomaRt)
  library(GenomicRanges)
  
  # connect to Ensembl Genes for chicken
  ensembl <- useEnsembl(
    biomart="genes", 
    dataset="ggallus_gene_ensembl"
  )  # :contentReference[oaicite:0]{index=0}
  
  # fetch chrom, start, end, strand for all transcripts of goi
  attrs <- c("chromosome_name","start_position","end_position","strand")
  res <- getBM(
    attributes = attrs,
    filters    = "external_gene_name",
    values     = goi,
    mart       = ensembl
  )
  if (nrow(res)==0) {
    stop(sprintf("Gene ‘%s’ not found in Ensembl ggallus_gene_ensembl", goi))
  }
  
  # collapse across transcripts to full gene span
  seqs <- unique(res$chromosome_name)
  if (length(seqs)!=1) {
    stop("Gene spans multiple chromosomes; please check annotation.")
  }
  gr <- GRanges(
    seqnames = seqs,
    ranges   = IRanges(
      start = min(res$start_position),
      end   = max(res$end_position)
    ),
    strand   = if (all(res$strand==1)) "+" else if (all(res$strand==2)) "-" else "*"
  )
  
  # extend by flanking bases
  if (upstream>0 || downstream>0) {
    start(gr) <- pmax(1, start(gr) - upstream)
    end(gr)   <- end(gr)   + downstream
  }
  return(gr)
}


## THIS IS FOR RUNNING CICER ~3hrs FOR THIS OBJECT, BEWARE LOADING NAMESPACE IS A BIT BUGGY WITH MONOCLE
###########################################################################
# 2. Convert Seurat object to Monocle3 CellDataSet for Cicero
#atac.cds <- SeuratWrappers::as.cell_data_set(atac)

#atac.cds <- cluster_cells(atac.cds, 
#                          reduction_method = "UMAP",
#                          # you can tweak resolution to get more/fewer clusters:
#                          cluster_resolution = 1e-2)

#reduced_coords <- reducedDims(atac.cds)$UMAP

# 4. Create Cicero CDS
#atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = reduced_coords)

# 5. Prepare genomic coordinates (chromosome sizes)
#    Retrieve chromosome lengths from the Seurat object
#genome_sizes <- seqlengths(atac)

#if not want run whole genome
#genome <- genome[1]

#    Convert to data.frame
#genome_df <- data.frame(
# chr    = names(genome_sizes),
#  length = as.numeric(genome_sizes)
#)

# 6. Run Cicero to compute co-accessibility scores
#    Adjust sample_num and window parameters as needed
#conns <- run_cicero(
#  atac.cicero,
#  genomic_coords = genome_df,
#  sample_num      = 100
#)

#head(conns)

# 7. Generate cis-co-accessible networks (CCANs)
#ccans <- generate_ccans(conns)

#head(ccans)

# 8. Convert Cicero connections to Signac links
#links <- ConnectionsToLinks(conns = conns, ccans = ccans)

# 9. Add links to the ChromatinAssay in Seurat
#Links(atac) <- links

# 10. Save the updated object
#saveRDS(atac, file = "atac_with_cicero_links.rds")
##############################################################

#Load cicero object
atac <- readRDS("atac_with_cicero_links.rds")


# your gene lists
top_genes    <- c("EDA2R","FABP9","TUBB6","SYNJ2","ADCYAP1","CTSK",
                  "DBNDD2","FXYD6","HSPB1","KLF2","LGALS1A","MSX1")
bottom_genes <- c("CA2","ENSGALG00010001303","GLUL","FN1","GFAP",
                  "HSPB9","PADI3","PRLH","SPARCL1","WIF1")

# define how many bases to flank each gene
flank.up   <- 1000
flank.down <- 1000

top_tracks       <- list()
top_cov_plots    <- list()

bottom_tracks    <- list()
bottom_cov_plots <- list()

# Top genes
for (goi in top_genes) {
  cat("Generating plots for:", goi, "\n")
  try({
    region_gr <- getGeneGRanges(goi, upstream = flank.up, downstream = flank.down)
    
    expr_plot <- ExpressionPlot(
      object   = rna,
      features = goi,
      assay    = "RNA",
      group.by = "predicted_celltype"
    )
    
    cov_plot <- CoveragePlot(
      object     = atac,
      region     = region_gr,
      annotation = TRUE,
      peaks      = TRUE,
      links      = TRUE,
      assay      = "ATAC",
      group.by   = "celltype",
      extend.upstream   = flank.up,
      extend.downstream = flank.down
    )
    
    combined <- CombineTracks(
      plotlist        = list(cov_plot),
      expression.plot = expr_plot,
      heights         = c(3, 10)
    )
    
    top_tracks[[goi]]    <- combined
    top_cov_plots[[goi]] <- cov_plot
  }, silent = TRUE)
}

# Bottom genes
for (goi in bottom_genes) {
  cat("Generating plots for:", goi, "\n")
  try({
    region_gr <- getGeneGRanges(goi, upstream = flank.up, downstream = flank.down)
    
    expr_plot <- ExpressionPlot(
      object   = rna,
      features = goi,
      assay    = "RNA",
      group.by = "predicted_celltype"
    )
    
    cov_plot <- CoveragePlot(
      object     = atac,
      region     = region_gr,
      annotation = TRUE,
      peaks      = TRUE,
      links      = TRUE,
      assay      = "ATAC",
      group.by   = "celltype",
      extend.upstream   = flank.up,
      extend.downstream = flank.down
    )
    
    combined <- CombineTracks(
      plotlist        = list(cov_plot),
      expression.plot = expr_plot,
      heights         = c(3, 10),    
      widths          = c(3, 12)     
    )
    print(combined)
    
    
    bottom_tracks[[goi]]    <- combined
    bottom_cov_plots[[goi]] <- cov_plot
  }, silent = TRUE)
}

for (plot in top_tracks){
  print(plot)
}

for (plot in bottom_tracks){
  print(plot)
}

