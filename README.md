\# Chicken Retina Multiome: Seurat + SCENIC+



Single-cell RNA + ATAC analysis of chicken (\*Gallus gallus\*, galGal6 / GRCg6a) retina from two samples (SS, NM). The main question is what separates \*\*activated\*\* from \*\*resting Müller glia (MG)\*\*, first through differential expression, accessibility and motif enrichment in Seurat/Signac, then through gene regulatory networks (eRegulons) inferred with SCENIC+.



The project runs in three stages, plus a separate side analysis:



1\. \*\*Seurat / Signac\*\* (R): build the RNA and ATAC objects, label cell types, run differential analysis and motif enrichment.

2\. \*\*SCENIC+ prep\*\*: export the data for Python and build a chicken cisTarget motif database.

3\. \*\*SCENIC+\*\* (Python): run the SCENIC+ Snakemake pipeline and analyse the eRegulons for MG.

4\. \*\*SAMap\*\* (side analysis): map chicken retina cells onto a mouse dataset.



\---



\## Repository layout



```

.

├── README.md

├── 01\_seurat/

│   ├── 01\_RNA\_ATAC\_Setup.Rmd

│   ├── 02\_combined\_differential.Rmd

│   ├── 03\_cicero\_coaccessibility.R             # atac\_analysis.R

│   ├── 04\_coaccessible\_motif.Rmd

│   ├── 04b\_coaccessible\_motif\_act\_vs\_rest.Rmd

│   └── 05\_signac\_scenic\_motif\_reconcile.R      # signac\_scenic\_extra\_ports.R

├── 02\_scenicplus\_prep/

│   ├── 01\_convert\_combined\_to\_scenic.R

│   ├── 02\_make\_genome\_annots.py

│   └── cistarget\_db/

│       ├── steps.txt

│       └── run\_cistarget\_db\_osc.sh

├── 03\_scenicplus/

│   ├── 01\_scenicplus\_allcells.ipynb            # Scenicplus.ipynb

│   └── 02\_scenicplus\_MG\_v2db.ipynb             # Scenicplusnewfeather.ipynb

├── 04\_samap\_chicken\_mouse/                     # side analysis, separate from the pipeline

│   ├── 01\_Seurat\_To\_Samap.R

│   ├── 02\_Chicken\_Mouse\_Samap\_Run.ipynb

│   ├── 02b\_Chicken\_Mouse\_Samap\_Run\_Sub48Chicken.ipynb

│   ├── 03\_Chicken\_Mouse\_Downstream.ipynb

│   ├── 03b\_Chicken\_Mouse\_Downstream\_Sub48Chicken.ipynb

│   └── 04\_Samap\_To\_Seurat.R

├── data/                                       # not tracked: 10x outputs, genome FASTA, reference objects

├── objects/                                    # not tracked: .rds objects made by 01\_seurat

├── scenicplus\_input/                           # not tracked: exports, cisTopic objects, feather DBs

└── results/                                    # CSV tables, figures and the Cytoscape network

```



The numbers give the run order, and the comment shows a file's current name where it is renamed. Code lives in the numbered folders; data and generated objects live in the four folders at the bottom, which are best kept out of git with a `.gitignore` (except `results/` if you want the final tables tracked).



\---



\## Reading the diagrams



Each stage has its own flowchart. The shapes show what each box is:



```mermaid

flowchart LR

&#x20;   i\[/"Input you supply<br/>(raw data, reference, download)"/]:::input

&#x20;   s\["R or Python script"]:::script

&#x20;   n(\["Jupyter notebook"]):::notebook

&#x20;   o\[("Object or table<br/>made by the pipeline")]:::obj

&#x20;   i --> s --> o --> n



&#x20;   classDef input fill:#fef3c7,stroke:#b45309,color:#000

&#x20;   classDef script fill:#dbeafe,stroke:#1d4ed8,color:#000

&#x20;   classDef notebook fill:#ede9fe,stroke:#6d28d9,color:#000

&#x20;   classDef obj fill:#f3f4f6,stroke:#6b7280,color:#000

```



Arrows point from what a step reads to what it writes. A small grey dot is a junction: everything flowing into it is used by every step flowing out of it. To keep the charts readable, each one draws only the main files a step builds on; the tables under each chart list every file a step reads.



\---



\## Stage 1: Seurat / Signac (`01\_seurat/`)



```mermaid

flowchart TD

&#x20;   raw\[/"10x outputs per sample (SS, NM)<br/>RNA + ATAC .h5, fragments.tsv.gz"/]:::input

&#x20;   ref\[/"Labelled reference<br/>rna\_Jack.rds"/]:::input

&#x20;   ann\[/"Ensembl 113 annotation<br/>(AnnotationHub)"/]:::input



&#x20;   s1\["01\_RNA\_ATAC\_Setup.Rmd"]:::script

&#x20;   raw --> s1

&#x20;   ref --> s1

&#x20;   ann --> s1

&#x20;   s1 --> pre\[("rna\_prePlot.rds<br/>atac\_prePlot.rds")]:::obj



&#x20;   pre --> s2\["02\_combined\_differential.Rmd"]:::script

&#x20;   pre --> s3\["03\_cicero\_coaccessibility.R"]:::script



&#x20;   s2 --> comb\[("combined.rds<br/>→ Stage 2")]:::obj

&#x20;   s2 --> da\[("da\_peaks.csv, da\_activity.csv,<br/>da\_rna\_full.csv, da\_combined\_full.csv")]:::obj

&#x20;   s3 --> cic\[("atac\_with\_cicero\_links.rds")]:::obj



&#x20;   da --- j(( )):::junction

&#x20;   cic --- j

&#x20;   j --> s4\["04\_coaccessible\_motif.Rmd"]:::script

&#x20;   j --> s4b\["04b\_coaccessible\_motif\_act\_vs\_rest.Rmd"]:::script



&#x20;   s4 --> mot\[("enriched\_motifs\_sorted.csv<br/>atac\_with\_motifs\_and\_stats.rds")]:::obj

&#x20;   s4b --> motb\[("enriched\_motifs\_resting\_vs\_activated\_sorted.csv")]:::obj



&#x20;   tbl\[/"SCENIC+ chicken motif table<br/>(motif collection download)"/]:::input

&#x20;   mot --> s5\["05\_signac\_scenic\_motif\_reconcile.R"]:::script

&#x20;   tbl --> s5

&#x20;   s5 --> rec\[("enriched\_motifs\_scenic\_reconciled.csv<br/>motif\_comparison\_summary.csv<br/>celltype\_mean\_expression.csv<br/>→ Stage 3")]:::obj



&#x20;   classDef input fill:#fef3c7,stroke:#b45309,color:#000

&#x20;   classDef script fill:#dbeafe,stroke:#1d4ed8,color:#000

&#x20;   classDef notebook fill:#ede9fe,stroke:#6d28d9,color:#000

&#x20;   classDef obj fill:#f3f4f6,stroke:#6b7280,color:#000

&#x20;   classDef junction fill:#6b7280,stroke:#6b7280

```



\### Scripts



| File | What it does |

|---|---|

| `01\_RNA\_ATAC\_Setup.Rmd` | Loads both samples, merges RNA, and transfers cell-type labels from the reference (`rna\_Jack.rds`, which splits cell types into finer subsets than the earlier `NEW\_clod\_NMDA.rds` reference). Builds a shared ATAC peak set, runs QC (TSS enrichment, nucleosome signal, drops a low-quality cluster), and labels ATAC cells through gene activity. |

| `02\_combined\_differential.Rmd` | Imputes RNA into ATAC, merges both into one object with a joint UMAP, and runs activated vs resting MG tests: differential peaks, gene activity and RNA DE. |

| `03\_cicero\_coaccessibility.R` | Runs Cicero co-accessibility on the ATAC object (about 3 h) and makes early coverage plots. |

| `04\_coaccessible\_motif.Rmd` | Coverage + expression plots for the top and bottom DE genes with Cicero links, then JASPAR motif annotation, chromVAR, and GC-matched motif enrichment in activated-MG peaks. |

| `04b\_coaccessible\_motif\_act\_vs\_rest.Rmd` | Variant that tests activated peaks against resting peaks directly, and the reverse. |

| `05\_signac\_scenic\_motif\_reconcile.R` | Reruns the Signac motif enrichment using only TFs that also appear in the SCENIC+ chicken motif table, so the two approaches can be compared. Also exports mean expression per cell type. |



\### Objects



| Object | Made by | Used by |

|---|---|---|

| `rna\_prePlot.rds` | 01 | 02, 04, 04b, 05 |

| `atac\_prePlot.rds` | 01 | 02, 03, 04, 04b |

| `combined.rds` | 02 | Stage 2 export |

| `da\_peaks.csv` | 02 | 04, 04b, 05 |

| `da\_rna\_full.csv`, `da\_combined\_full.csv` | 02 | Stage 3 notebooks |

| `atac\_with\_cicero\_links.rds` | 03 | 04, 04b |

| `atac\_with\_motifs\_and\_stats.rds` | 04 | 05 |

| `enriched\_motifs\_sorted.csv` | 04 | 05 |

| `enriched\_motifs\_scenic\_reconciled.csv`, `motif\_comparison\_summary.csv`, `celltype\_mean\_expression.csv` | 05 | Stage 3 notebooks |



\---



\## Stage 2: SCENIC+ prep (`02\_scenicplus\_prep/`)



The motif database is built from the regions in the cisTopic model, so it depends on the first part of the all-cells notebook (cell 0: RNA preprocessing and cisTopic). The rest of that notebook runs in Stage 3, once the database exists.



```mermaid

flowchart TD

&#x20;   comb\[("combined.rds<br/>from Stage 1")]:::obj

&#x20;   p1\["01\_convert\_combined\_to\_scenic.R"]:::script

&#x20;   comb --> p1

&#x20;   p1 --> exp\[("scenicplus\_input/<br/>rna\_matrix.mtx, atac\_matrix.mtx,<br/>barcodes, metadata, peaks.bed,<br/>fragments\_paths.json, umap\_combined.csv")]:::obj



&#x20;   fasta\[/"Gallus\_gallus.GRCg6a.dna.toplevel.fa"/]:::input

&#x20;   p2\["02\_make\_genome\_annots.py"]:::script

&#x20;   fasta --> p2

&#x20;   p2 --> gann\[("genome\_annotation.tsv<br/>chromsizes.tsv")]:::obj



&#x20;   nbA(\["01\_scenicplus\_allcells.ipynb<br/>cell 0: RNA + cisTopic"]):::notebook

&#x20;   exp --> nbA

&#x20;   nbA --> cisto\[("rna\_scanpy\_preprocessed.h5ad<br/>atac\_cistopic\_selected\_model\_chrstr.pkl<br/>region\_sets\_topics\_otsu.pkl")]:::obj



&#x20;   gf\["bedtools getfasta<br/>(steps.txt)"]:::script

&#x20;   cisto --> gf

&#x20;   fasta --> gf

&#x20;   gf --> cfa\[("consensus\_regions.fa")]:::obj



&#x20;   mcoll\[/"v10nr\_clust motif collection<br/>+ chicken snapshot .tbl"/]:::input

&#x20;   mf\["motif filtering<br/>(steps.txt + filter\_motifs\_by\_expression.py)"]:::script

&#x20;   mcoll --> mf

&#x20;   mf --> mlist\[("motifs\_chicken\_expr\_for\_db.lst")]:::obj



&#x20;   db\["run\_cistarget\_db\_osc.sh<br/>(OSC HPC)"]:::script

&#x20;   cfa --> db

&#x20;   mlist --> db

&#x20;   mcoll --> db

&#x20;   db --> feather\[("gallus\_consensus\_chicken\_expr<br/>.rankings.feather / .scores.feather")]:::obj



&#x20;   classDef input fill:#fef3c7,stroke:#b45309,color:#000

&#x20;   classDef script fill:#dbeafe,stroke:#1d4ed8,color:#000

&#x20;   classDef notebook fill:#ede9fe,stroke:#6d28d9,color:#000

&#x20;   classDef obj fill:#f3f4f6,stroke:#6b7280,color:#000

```



| File | What it does |

|---|---|

| `01\_convert\_combined\_to\_scenic.R` | Exports `combined.rds` into `scenicplus\_input/` as matrices, barcodes, metadata, `peaks.bed`, fragment paths and the UMAP. |

| `02\_make\_genome\_annots.py` | Pulls the gene annotation from Ensembl BioMart and computes chromosome sizes from the GRCg6a FASTA. |

| `cistarget\_db/steps.txt` | How the cisTarget database was built: filtering motifs to chicken, filtering by expression, building the consensus FASTA, the database command, and downloading the genome annotation for Snakemake. |

| `cistarget\_db/run\_cistarget\_db\_osc.sh` | SLURM job that builds the cisTarget motif database on OSC. |



\---



\## Stage 3: SCENIC+ (`03\_scenicplus/`)



```mermaid

flowchart TD

&#x20;   inp\[("From Stages 1–2<br/>cisTopic model + RNA h5ad<br/>genome\_annotation.tsv<br/>da\_\*.csv + reconciled motif tables")]:::obj

&#x20;   v1\[("cisTarget DB v1<br/>from Stage 2")]:::obj



&#x20;   inp --> nb1(\["01\_scenicplus\_allcells.ipynb<br/>remaining cells + Snakemake"]):::notebook

&#x20;   v1 --> nb1

&#x20;   nb1 --> out1\[("scplusmdata.h5mu<br/>eRegulons\_extended.tsv<br/>MG\_TF\_crosslayer\_evidence.csv<br/>comparison tables, volcano plots")]:::obj



&#x20;   v2\[/"cisTarget DB v2<br/>motifs\_chicken\_for\_db\_v2\_final.lst"/]:::input

&#x20;   out1 --> nb2(\["02\_scenicplus\_MG\_v2db.ipynb<br/>MG-only cisTopic + Snakemake"]):::notebook

&#x20;   inp --> nb2

&#x20;   v2 --> nb2

&#x20;   nb2 --> final\[("results/<br/>Cytoscape network<br/>TF knockout summary<br/>top eRegulons + targets")]:::obj



&#x20;   classDef input fill:#fef3c7,stroke:#b45309,color:#000

&#x20;   classDef script fill:#dbeafe,stroke:#1d4ed8,color:#000

&#x20;   classDef notebook fill:#ede9fe,stroke:#6d28d9,color:#000

&#x20;   classDef obj fill:#f3f4f6,stroke:#6b7280,color:#000

```



Both notebooks run the SCENIC+ Snakemake pipeline from inside the notebook.



| File | What it does |

|---|---|

| `01\_scenicplus\_allcells.ipynb` | SCENIC+ run on all cell types. Cell 0 does the scanpy RNA preprocessing and cisTopic topic modelling (Stage 2). The remaining cells check and repair the feather database and region names, run the Snakemake pipeline, and compare activated vs resting MG across RNA, ATAC and eRegulon layers. |

| `02\_scenicplus\_MG\_v2db.ipynb` | \*\*Final analysis.\*\* Reruns cisTopic on MG cells only, builds per-cell-type region sets, and runs SCENIC+ against the rebuilt v2 database. Produces the Cytoscape network, simulated TF knockouts (MYC, RUNX1, THRB, ELF3) and the top eRegulon/target tables. |



\---



\## Side analysis: chicken–mouse SAMap (`04\_samap\_chicken\_mouse/`)



This maps chicken retina cells onto a mouse retina dataset with \[SAMap](https://github.com/atarashansky/SAMap), which aligns species using BLAST gene homology instead of shared gene names. It is related to the main project but doesn't feed into, or depend on, the stages above. There are two parallel runs: one with the full chicken dataset and one with a smaller chicken subset (`sub48`). The `b` files are the subset versions of the same steps.



```mermaid

flowchart TD

&#x20;   chrds\[/"Chicken Seurat object (.rds)"/]:::input

&#x20;   msrds\[/"mouse.rds"/]:::input

&#x20;   blast\[/"BLAST gene maps<br/>chicken\_and\_mouse/chms/maps/"/]:::input



&#x20;   a1\["01\_Seurat\_To\_Samap.R"]:::script

&#x20;   chrds --> a1

&#x20;   msrds --> a1

&#x20;   a1 --> h5in\[("samap\_adata\_sparse.h5ad (chicken)<br/>samap\_adata\_chicken\_sub48.h5ad (subset)<br/>samap\_adata\_sparse\_mouse.h5ad (mouse)")]:::obj



&#x20;   a2(\["02\_Chicken\_Mouse\_Samap\_Run.ipynb<br/>02b for subset"]):::notebook

&#x20;   h5in --> a2

&#x20;   blast --> a2

&#x20;   a2 --> res\[("samap\_object\_\*.pkl<br/>samap\_adata\_\*.pkl<br/>samap\_results\_h5ad\_\*.h5ad")]:::obj



&#x20;   a3(\["03\_Chicken\_Mouse\_Downstream.ipynb<br/>03b for subset"]):::notebook

&#x20;   res --> a3

&#x20;   a3 --> dout\[("Sankey + chord plots<br/>treatmentgenepairs.csv<br/>celltypegenepairs.csv")]:::obj



&#x20;   a4\["04\_Samap\_To\_Seurat.R"]:::script

&#x20;   res --> a4

&#x20;   a4 --> rds\[("SAMap results as Seurat .rds")]:::obj



&#x20;   classDef input fill:#fef3c7,stroke:#b45309,color:#000

&#x20;   classDef script fill:#dbeafe,stroke:#1d4ed8,color:#000

&#x20;   classDef notebook fill:#ede9fe,stroke:#6d28d9,color:#000

&#x20;   classDef obj fill:#f3f4f6,stroke:#6b7280,color:#000

```



| File | What it does |

|---|---|

| `01\_Seurat\_To\_Samap.R` | Converts a Seurat object to `.h5ad` with sceasy so SAMap can read it, storing the cluster identities in `active\_ident`. Set up for the mouse object, and rerun with the chicken object (or subset) to make the chicken input. |

| `02\_Chicken\_Mouse\_Samap\_Run.ipynb` | Loads the chicken and mouse `.h5ad` files and the BLAST homology maps, runs SAMap, shows the joint UMAP, and saves the SAMap object, the combined AnnData and a results `.h5ad`. |

| `02b\_Chicken\_Mouse\_Samap\_Run\_Sub48Chicken.ipynb` | Same run using the chicken subset. |

| `03\_Chicken\_Mouse\_Downstream.ipynb` | Loads the saved SAMap object and compares the species twice: chicken treatment against mouse cell types, then cell type against cell type. Each comparison gives mapping scores, Sankey and chord plots, and the gene pairs driving the alignment. |

| `03b\_Chicken\_Mouse\_Downstream\_Sub48Chicken.ipynb` | Same analysis on the subset run (groups chicken cells by `orig.ident` for the first comparison). |

| `04\_Samap\_To\_Seurat.R` | Loads the SAMap pickles through reticulate and converts the results `.h5ad` back into a Seurat object for further analysis in R. |



The BLAST maps come from SAMap's `map\_genes.sh` run on the chicken and mouse proteomes; that step isn't included here.



\---



\## Software



R: Seurat v5, Signac, SeuratWrappers, BSgenome.Ggallus.UCSC.galGal6, JASPAR2020, TFBSTools, monocle3, cicero, biomaRt, AnnotationHub.

Python: scanpy, pycisTopic, pycistarget, scenicplus (Snakemake pipeline), mudata, pyranges.

cisTarget database: `create\_cisTarget\_databases` with Cluster-Buster and the v10nr\_clust motif collection.

SAMap: samap, sam-algorithm (Python 3.10), sceasy and reticulate for converting between Seurat and AnnData.

