# RNA-seq analysis code — tomato leaf transcriptome

Scripts and commands behind the RNA-seq bioinformatic analyses and figures of
the Scientific Data descriptor *Multi-omics analyses of tomato rhizosphere and
leaves treated with a biowaste-derived biostimulant extract* (Bona *et al.*,
Fondazione Edmund Mach).

Everything from raw reads to the published figures is here: the HPC pipeline
that produced the count matrix, the technical-validation QC, the variance
decomposition that quantifies the contributions of experimental trial, time
point and treatment, and the assembly of Figure 2 itself.

The data are deposited separately and openly:

| Resource | Accession |
|---|---|
| Raw sequencing reads | NCBI SRA BioProject **PRJNA1425279** |
| Count matrices, physiology, soil and leaf chemistry | Zenodo **[10.5281/zenodo.19336336](https://doi.org/10.5281/zenodo.19336336)** |

---

## Layout

```
config/paths.R                 single source of truth for every file path
metadata/targets.txt           60-library sample sheet (Condition, Time, Experiment)
data/                          inputs; the large matrix is on Zenodo (data/README.md)
results/                       every figure and table the scripts write
env/                           R package versions for the archived run
  record_environment.R         regenerates r_packages.tsv and sessionInfo.txt

01_preprocessing_alignment/    HPC pipeline (SLURM + Apptainer)
  01_run_fastqc.slurm          read QC
  02_run_fastp.slurm           adapter and quality trimming
  03_run_star_samtools.slurm   STAR alignment to SL5.0 / ITAG5.0, sort, index
  04_run_featurecounts.slurm   featureCounts -s 2, gene-level counts
  legacy/                      an alternative HISAT2 route explored but NOT used

02_rnaseq_qc/                  technical validation; Figure 2a and 2b
  plot_read_retention.R            read fate through the pipeline
  distance_heatmap_functions.R     builders for the distance heatmap
  plot_sample_distance_heatmap.R   blind VST, Euclidean distances, clustering
  rnaseq_qc_table_and_plots.R      supplementary QC table and QC panels

03_batch_effect/               technical validation; Figure 2c and 2d
  variance_partition.R         variancePartition: contribution of trial, time
                               point and treatment, per gene and per principal
                               component, with a moderated F-test per term

04_figures/
  build_figure2.R              assembles Figure 2 from the four panel objects
```

## Running it

Every script resolves its paths through `config/paths.R`; no absolute path
appears anywhere in the R code. Run one step from any directory:

```bash
Rscript 03_batch_effect/variance_partition.R
```

or the whole downstream analysis in dependency order (about a minute and a
half, ending with the assembled Figure 2):

```bash
./run_all.sh
```

Override locations if the data live elsewhere:

```bash
export PERCIVAL_ROOT=/path/to/this/repo
export PERCIVAL_DATA=/path/to/large/matrices
export PERCIVAL_RESULTS=/path/to/output
```

The `.slurm` scripts ran on the FEM HPC cluster and still carry that cluster's
paths in their header variables. They are deposited as the executed record of
the pipeline; edit those variables to re-run elsewhere. See
`01_preprocessing_alignment/README.md`.

## Pipeline order

1. `01_run_fastqc.slurm` → `02_run_fastp.slurm` → `03_run_star_samtools.slurm`
   → `04_run_featurecounts.slurm`, producing the raw count matrix deposited on
   Zenodo as `RNAseq_rawcounts.tsv`.
2. `02_rnaseq_qc/*` and `03_batch_effect/*` each depend only on step 1 and are
   independent of one another.
3. `04_figures/build_figure2.R` depends on both and only composes: each panel
   is saved by the script that computes it, so the assembled figure cannot
   drift from the analysis behind it.

## Figure mapping

| Manuscript panel | Produced by | File under `results/` |
|---|---|---|
| Figure 2a | `02_rnaseq_qc/plot_read_retention.R` | `rnaseq_qc/read_retention.pdf` |
| Figure 2b | `02_rnaseq_qc/plot_sample_distance_heatmap.R` | `rnaseq_qc/sample_distance_heatmap.pdf` |
| Figure 2c, upper | `03_batch_effect/variance_partition.R` | `batch_effect/Fig_PCA.pdf` |
| Figure 2c, lower | `03_batch_effect/variance_partition.R` | `batch_effect/Fig_variance_per_axis.pdf` |
| Figure 2d | `03_batch_effect/variance_partition.R` | `batch_effect/Fig_variancePartition_violin.pdf` |
| Figure 2, assembled | `04_figures/build_figure2.R` | `figures/Figure2_complete.pdf` |
| Figure 2c–d alone | `04_figures/build_figure2.R` | `figures/Figure2_panels_c_d.pdf` |
| Supplementary | `03_batch_effect/variance_partition.R` | `batch_effect/Fig_canCorPairs.pdf`, `Fig_percentBars_top_genes.pdf` |
| Supplementary Table (QC) | `02_rnaseq_qc/rnaseq_qc_table_and_plots.R` | `rnaseq_qc/Supplementary_Table_RNAseq_QC.tsv` |

Two PCA plots exist and are not interchangeable. `rnaseq_qc/figures_QC/Fig_2D_PCA.pdf`
is an exploratory QC plot on the top 2,000 most-variable genes. The manuscript
PCA is `batch_effect/Fig_PCA.pdf`, computed on all 18,425 expressed genes with
no most-variable-gene selection; the percentages differ because the gene sets
differ. `rnaseq_qc/figures_QC/` also holds mapping-rate, marker-gene and
gene-detection-saturation panels, kept as the supporting QC evidence behind the
technical validation.

## Software

R 4.6.1 with Bioconductor: principally variancePartition 1.42.0, limma 3.68.5,
edgeR 4.10.3, DESeq2 1.52.0, ComplexHeatmap 2.28.0, cowplot 1.2.0 and
ggplot2 4.0.3. The full list is in `env/r_packages.tsv`, the complete session
including transitive dependencies in `env/sessionInfo.txt`, and the session for
the archived variance-decomposition run in
`results/batch_effect/sessionInfo.txt`. Regenerate the first two with
`Rscript env/record_environment.R`.

Command-line tools, as invoked in `01_preprocessing_alignment/`: FastQC 0.12.1,
fastp 0.23.4, STAR 2.7.11b, Subread/featureCounts 2.1.1, samtools 1.22.1, all
inside Apptainer images.

## Reference genome

*Solanum lycopersicum* v5.0 (assembly SL5.0, ITAG5.0 annotation), Phytozome
JGI: <https://phytozome-next.jgi.doe.gov/info/Slycopersicum_ITAG5_0>

## Licence

Code: MIT (`../LICENSE`). Data on Zenodo: CC-BY 4.0.
