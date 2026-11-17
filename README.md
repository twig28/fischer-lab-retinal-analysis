# Chicken Retina Multiome: Seurat + SCENIC+

Single-cell RNA + ATAC analysis of chicken (*Gallus gallus*, galGal6 / GRCg6a) retina from two samples (SS, NM). The main question is what separates **activated** from **resting Müller glia (MG)**, first through differential expression, accessibility and motif enrichment in Seurat/Signac, then through gene regulatory networks (eRegulons) inferred with SCENIC+.

The project runs in three stages, plus a separate side analysis:

1. **Seurat / Signac** (R): build the RNA and ATAC objects, label cell types, run differential analysis and motif enrichment.
2. **SCENIC+ prep**: export the data for Python and build a chicken cisTarget motif database.
3. **SCENIC+** (Python): run the SCENIC+ Snakemake pipeline and analyse the eRegulons for MG.
4. **SAMap** (side analysis): map chicken retina cells onto a mouse dataset.

---

## Repository layout

```
.
├── README.md
├── 01_seurat/
│   ├── 01_RNA_ATAC_Setup.Rmd
│   ├── 02_combined_differential.Rmd
│   ├── 03_cicero_coaccessibility.R             # atac_analysis.R
│   ├── 04_coaccessible_motif.Rmd
│   ├── 04b_coaccessible_motif_act_vs_rest.Rmd
│   └── 05_signac_scenic_motif_reconcile.R      # signac_scenic_extra_ports.R
├── 02_scenicplus_prep/
│   ├── 01_convert_combined_to_scenic.R
│   ├── 02_make_genome_annots.py
│   └── cistarget_db/
│       ├── steps.txt
│       └── run_cistarget_db_osc.sh
├── 03_scenicplus/
│   ├── 01_scenicplus_allcells.ipynb            # Scenicplus.ipynb
│   └── 02_scenicplus_MG_v2db.ipynb             # Scenicplusnewfeather.ipynb
├── 04_samap_chicken_mouse/                     # side analysis, separate from the pipeline
│   ├── 01_Seurat_To_Samap.R
│   ├── 02_Chicken_Mouse_Samap_Run.ipynb
│   ├── 02b_Chicken_Mouse_Samap_Run_Sub48Chicken.ipynb
│   ├── 03_Chicken_Mouse_Downstream.ipynb
│   ├── 03b_Chicken_Mouse_Downstream_Sub48Chicken.ipynb
│   └── 04_Samap_To_Seurat.R
├── data/                                       # not tracked: 10x outputs, genome FASTA, reference objects
├── objects/                                    # not tracked: .rds objects made by 01_seurat
├── scenicplus_input/                           # not tracked: exports, cisTopic objects, feather DBs
└── results/                                    # CSV tables, figures and the Cytoscape network
```

The numbers give the run order, and the comment shows a file's current name where it is renamed. Code lives in the numbered folders; data and generated objects live in the four folders at the bottom, which are best kept out of git with a `.gitignore` (except `results/` if you want the final tables tracked).

---

## Reading the diagrams

Each stage has its own flowchart. The shapes show what each box is:

```mermaid
flowchart LR
    i[/"Input you supply<br/>(raw data, reference, download)"/]:::input
    s["R or Python script"]:::script
    n(["Jupyter notebook"]):::notebook
    o[("Object or table<br/>made by the pipeline")]:::obj
    i --> s --> o --> n

    classDef input fill:#fef3c7,stroke:#b45309,color:#000
    classDef script fill:#dbeafe,stroke:#1d4ed8,color:#000
    classDef notebook fill:#ede9fe,stroke:#6d28d9,color:#000
    classDef obj fill:#f3f4f6,stroke:#6b7280,color:#000
```

Arrows point from what a step reads to what it writes. A small grey dot is a junction: everything flowing into it is used by every step flowing out of it. To keep the charts readable, each one draws only the main files a step builds on; the tables under each chart list every file a step reads.

---

## Stage 1: Seurat / Signac (`01_seurat/`)

```mermaid
flowchart TD
    raw[/"10x outputs per sample (SS, NM)<br/>RNA + ATAC .h5, fragments.tsv.gz"/]:::input
    ref[/"Labelled reference<br/>rna_Jack.rds"/]:::input
    ann[/"Ensembl 113 annotation<br/>(AnnotationHub)"/]:::input

    s1["01_RNA_ATAC_Setup.Rmd"]:::script
    raw --> s1
    ref --> s1
    ann --> s1
    s1 --> pre[("rna_prePlot.rds<br/>atac_prePlot.rds")]:::obj

    pre --> s2["02_combined_differential.Rmd"]:::script
    pre --> s3["03_cicero_coaccessibility.R"]:::script

    s2 --> comb[("combined.rds<br/>→ Stage 2")]:::obj
    s2 --> da[("da_peaks.csv, da_activity.csv,<br/>da_rna_full.csv, da_combined_full.csv")]:::obj
    s3 --> cic[("atac_with_cicero_links.rds")]:::obj

    da --- j(( )):::junction
    cic --- j
    j --> s4["04_coaccessible_motif.Rmd"]:::script
    j --> s4b["04b_coaccessible_motif_act_vs_rest.Rmd"]:::script

    s4 --> mot[("enriched_motifs_sorted.csv<br/>atac_with_motifs_and_stats.rds")]:::obj
    s4b --> motb[("enriched_motifs_resting_vs_activated_sorted.csv")]:::obj

    tbl[/"SCENIC+ chicken motif table<br/>(motif collection download)"/]:::input
    mot --> s5["05_signac_scenic_motif_reconcile.R"]:::script
    tbl --> s5
    s5 --> rec[("enriched_motifs_scenic_reconciled.csv<br/>motif_comparison_summary.csv<br/>celltype_mean_expression.csv<br/>→ Stage 3")]:::obj

    classDef input fill:#fef3c7,stroke:#b45309,color:#000
    classDef script fill:#dbeafe,stroke:#1d4ed8,color:#000
    classDef notebook fill:#ede9fe,stroke:#6d28d9,color:#000
    classDef obj fill:#f3f4f6,stroke:#6b7280,color:#000
    classDef junction fill:#6b7280,stroke:#6b7280
```

### Scripts

| File | What it does |
|---|---|
| `01_RNA_ATAC_Setup.Rmd` | Loads both samples, merges RNA, and transfers cell-type labels from the reference (`rna_Jack.rds`, which splits cell types into finer subsets than the earlier `NEW_clod_NMDA.rds` reference). Builds a shared ATAC peak set, runs QC (TSS enrichment, nucleosome signal, drops a low-quality cluster), and labels ATAC cells through gene activity. |
| `02_combined_differential.Rmd` | Imputes RNA into ATAC, merges both into one object with a joint UMAP, and runs activated vs resting MG tests: differential peaks, gene activity and RNA DE. |
| `03_cicero_coaccessibility.R` | Runs Cicero co-accessibility on the ATAC object (about 3 h) and makes early coverage plots. |
| `04_coaccessible_motif.Rmd` | Coverage + expression plots for the top and bottom DE genes with Cicero links, then JASPAR motif annotation, chromVAR, and GC-matched motif enrichment in activated-MG peaks. |
| `04b_coaccessible_motif_act_vs_rest.Rmd` | Variant that tests activated peaks against resting peaks directly, and the reverse. |
| `05_signac_scenic_motif_reconcile.R` | Reruns the Signac motif enrichment using only TFs that also appear in the SCENIC+ chicken motif table, so the two approaches can be compared. Also exports mean expression per cell type. |

### Objects

| Object | Made by | Used by |
|---|---|---|
| `rna_prePlot.rds` | 01 | 02, 04, 04b, 05 |
| `atac_prePlot.rds` | 01 | 02, 03, 04, 04b |
| `combined.rds` | 02 | Stage 2 export |
| `da_peaks.csv` | 02 | 04, 04b, 05 |
| `da_rna_full.csv`, `da_combined_full.csv` | 02 | Stage 3 notebooks |
| `atac_with_cicero_links.rds` | 03 | 04, 04b |
| `atac_with_motifs_and_stats.rds` | 04 | 05 |
| `enriched_motifs_sorted.csv` | 04 | 05 |
| `enriched_motifs_scenic_reconciled.csv`, `motif_comparison_summary.csv`, `celltype_mean_expression.csv` | 05 | Stage 3 notebooks |

---

## Stage 2: SCENIC+ prep (`02_scenicplus_prep/`)

The motif database is built from the regions in the cisTopic model, so it depends on the first part of the all-cells notebook (cell 0: RNA preprocessing and cisTopic). The rest of that notebook runs in Stage 3, once the database exists.

```mermaid
flowchart TD
    comb[("combined.rds<br/>from Stage 1")]:::obj
    p1["01_convert_combined_to_scenic.R"]:::script
    comb --> p1
    p1 --> exp[("scenicplus_input/<br/>rna_matrix.mtx, atac_matrix.mtx,<br/>barcodes, metadata, peaks.bed,<br/>fragments_paths.json, umap_combined.csv")]:::obj

    fasta[/"Gallus_gallus.GRCg6a.dna.toplevel.fa"/]:::input
    p2["02_make_genome_annots.py"]:::script
    fasta --> p2
    p2 --> gann[("genome_annotation.tsv<br/>chromsizes.tsv")]:::obj

    nbA(["01_scenicplus_allcells.ipynb<br/>cell 0: RNA + cisTopic"]):::notebook
    exp --> nbA
    nbA --> cisto[("rna_scanpy_preprocessed.h5ad<br/>atac_cistopic_selected_model_chrstr.pkl<br/>region_sets_topics_otsu.pkl")]:::obj

    gf["bedtools getfasta<br/>(steps.txt)"]:::script
    cisto --> gf
    fasta --> gf
    gf --> cfa[("consensus_regions.fa")]:::obj

    mcoll[/"v10nr_clust motif collection<br/>+ chicken snapshot .tbl"/]:::input
    mf["motif filtering<br/>(steps.txt + filter_motifs_by_expression.py)"]:::script
    mcoll --> mf
    mf --> mlist[("motifs_chicken_expr_for_db.lst")]:::obj

    db["run_cistarget_db_osc.sh<br/>(OSC HPC)"]:::script
    cfa --> db
    mlist --> db
    mcoll --> db
    db --> feather[("gallus_consensus_chicken_expr<br/>.rankings.feather / .scores.feather")]:::obj

    classDef input fill:#fef3c7,stroke:#b45309,color:#000
    classDef script fill:#dbeafe,stroke:#1d4ed8,color:#000
    classDef notebook fill:#ede9fe,stroke:#6d28d9,color:#000
    classDef obj fill:#f3f4f6,stroke:#6b7280,color:#000
```

| File | What it does |
|---|---|
| `01_convert_combined_to_scenic.R` | Exports `combined.rds` into `scenicplus_input/` as matrices, barcodes, metadata, `peaks.bed`, fragment paths and the UMAP. |
| `02_make_genome_annots.py` | Pulls the gene annotation from Ensembl BioMart and computes chromosome sizes from the GRCg6a FASTA. |
| `cistarget_db/steps.txt` | How the cisTarget database was built: filtering motifs to chicken, filtering by expression, building the consensus FASTA, the database command, and downloading the genome annotation for Snakemake. |
| `cistarget_db/run_cistarget_db_osc.sh` | SLURM job that builds the cisTarget motif database on OSC. |

---

## Stage 3: SCENIC+ (`03_scenicplus/`)

```mermaid
flowchart TD
    inp[("From Stages 1–2<br/>cisTopic model + RNA h5ad<br/>genome_annotation.tsv<br/>da_*.csv + reconciled motif tables")]:::obj
    v1[("cisTarget DB v1<br/>from Stage 2")]:::obj

    inp --> nb1(["01_scenicplus_allcells.ipynb<br/>remaining cells + Snakemake"]):::notebook
    v1 --> nb1
    nb1 --> out1[("scplusmdata.h5mu<br/>eRegulons_extended.tsv<br/>MG_TF_crosslayer_evidence.csv<br/>comparison tables, volcano plots")]:::obj

    v2[/"cisTarget DB v2<br/>motifs_chicken_for_db_v2_final.lst"/]:::input
    out1 --> nb2(["02_scenicplus_MG_v2db.ipynb<br/>MG-only cisTopic + Snakemake"]):::notebook
    inp --> nb2
    v2 --> nb2
    nb2 --> final[("results/<br/>Cytoscape network<br/>TF knockout summary<br/>top eRegulons + targets")]:::obj

    classDef input fill:#fef3c7,stroke:#b45309,color:#000
    classDef script fill:#dbeafe,stroke:#1d4ed8,color:#000
    classDef notebook fill:#ede9fe,stroke:#6d28d9,color:#000
    classDef obj fill:#f3f4f6,stroke:#6b7280,color:#000
```

Both notebooks run the SCENIC+ Snakemake pipeline from inside the notebook.

| File | What it does |
|---|---|
| `01_scenicplus_allcells.ipynb` | SCENIC+ run on all cell types. Cell 0 does the scanpy RNA preprocessing and cisTopic topic modelling (Stage 2). The remaining cells check and repair the feather database and region names, run the Snakemake pipeline, and compare activated vs resting MG across RNA, ATAC and eRegulon layers. |
| `02_scenicplus_MG_v2db.ipynb` | **Final analysis.** Reruns cisTopic on MG cells only, builds per-cell-type region sets, and runs SCENIC+ against the rebuilt v2 database. Produces the Cytoscape network, simulated TF knockouts (MYC, RUNX1, THRB, ELF3) and the top eRegulon/target tables. |

---

## Side analysis: chicken–mouse SAMap (`04_samap_chicken_mouse/`)

This maps chicken retina cells onto a mouse retina dataset with [SAMap](https://github.com/atarashansky/SAMap), which aligns species using BLAST gene homology instead of shared gene names. It is related to the main project but doesn't feed into, or depend on, the stages above. There are two parallel runs: one with the full chicken dataset and one with a smaller chicken subset (`sub48`). The `b` files are the subset versions of the same steps.

```mermaid
flowchart TD
    chrds[/"Chicken Seurat object (.rds)"/]:::input
    msrds[/"mouse.rds"/]:::input
    blast[/"BLAST gene maps<br/>chicken_and_mouse/chms/maps/"/]:::input

    a1["01_Seurat_To_Samap.R"]:::script
    chrds --> a1
    msrds --> a1
    a1 --> h5in[("samap_adata_sparse.h5ad (chicken)<br/>samap_adata_chicken_sub48.h5ad (subset)<br/>samap_adata_sparse_mouse.h5ad (mouse)")]:::obj

    a2(["02_Chicken_Mouse_Samap_Run.ipynb<br/>02b for subset"]):::notebook
    h5in --> a2
    blast --> a2
    a2 --> res[("samap_object_*.pkl<br/>samap_adata_*.pkl<br/>samap_results_h5ad_*.h5ad")]:::obj

    a3(["03_Chicken_Mouse_Downstream.ipynb<br/>03b for subset"]):::notebook
    res --> a3
    a3 --> dout[("Sankey + chord plots<br/>treatmentgenepairs.csv<br/>celltypegenepairs.csv")]:::obj

    a4["04_Samap_To_Seurat.R"]:::script
    res --> a4
    a4 --> rds[("SAMap results as Seurat .rds")]:::obj

    classDef input fill:#fef3c7,stroke:#b45309,color:#000
    classDef script fill:#dbeafe,stroke:#1d4ed8,color:#000
    classDef notebook fill:#ede9fe,stroke:#6d28d9,color:#000
    classDef obj fill:#f3f4f6,stroke:#6b7280,color:#000
```

| File | What it does |
|---|---|
| `01_Seurat_To_Samap.R` | Converts a Seurat object to `.h5ad` with sceasy so SAMap can read it, storing the cluster identities in `active_ident`. Set up for the mouse object, and rerun with the chicken object (or subset) to make the chicken input. |
| `02_Chicken_Mouse_Samap_Run.ipynb` | Loads the chicken and mouse `.h5ad` files and the BLAST homology maps, runs SAMap, shows the joint UMAP, and saves the SAMap object, the combined AnnData and a results `.h5ad`. |
| `02b_Chicken_Mouse_Samap_Run_Sub48Chicken.ipynb` | Same run using the chicken subset. |
| `03_Chicken_Mouse_Downstream.ipynb` | Loads the saved SAMap object and compares the species twice: chicken treatment against mouse cell types, then cell type against cell type. Each comparison gives mapping scores, Sankey and chord plots, and the gene pairs driving the alignment. |
| `03b_Chicken_Mouse_Downstream_Sub48Chicken.ipynb` | Same analysis on the subset run (groups chicken cells by `orig.ident` for the first comparison). |
| `04_Samap_To_Seurat.R` | Loads the SAMap pickles through reticulate and converts the results `.h5ad` back into a Seurat object for further analysis in R. |

The BLAST maps come from SAMap's `map_genes.sh` run on the chicken and mouse proteomes; that step isn't included here.

---

## Software

R: Seurat v5, Signac, SeuratWrappers, BSgenome.Ggallus.UCSC.galGal6, JASPAR2020, TFBSTools, monocle3, cicero, biomaRt, AnnotationHub.
Python: scanpy, pycisTopic, pycistarget, scenicplus (Snakemake pipeline), mudata, pyranges.
cisTarget database: `create_cisTarget_databases` with Cluster-Buster and the v10nr_clust motif collection.
SAMap: samap, sam-algorithm (Python 3.10), sceasy and reticulate for converting between Seurat and AnnData.
