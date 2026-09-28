# Analysis code for *Multi-omics analyses of tomato rhizosphere and leaves treated with a biowaste-derived biostimulant extract*

Bona *et al.*, Fondazione Edmund Mach — Scientific Data descriptor.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22232657.svg)](https://doi.org/10.5281/zenodo.22232657)

This repository holds the scripts and commands used for the bioinformatic
analyses and for generating the figures of the descriptor: the leaf
transcriptome (RNA-seq) and the rhizosphere prokaryotic and fungal communities
(16S and ITS metabarcoding). It is code and results; the data themselves are
deposited separately and openly:

| Resource | Accession |
|---|---|
| Raw sequencing reads | NCBI SRA BioProject **PRJNA1425279** |
| Count matrices, physiology, soil and leaf chemistry | Zenodo **[10.5281/zenodo.19336336](https://doi.org/10.5281/zenodo.19336336)** |
| This repository, archived | Zenodo **[10.5281/zenodo.22232657](https://doi.org/10.5281/zenodo.22232657)** |

The two Zenodo records are distinct: `19336336` is the data deposit,
`22232657` is the archived code. Both are concept DOIs and always resolve to
the most recent version of their record.

---

## What is here

```
RNA-seq/         everything behind the leaf transcriptome: raw reads to
                 Figure 2 and Supplementary Table 6
metabarcoding/   the rhizosphere 16S and ITS pipeline (MICCA) and the
                 validation plots behind Supplementary Figures S1 and S2
```

`RNA-seq/README.md` documents the layout, how to run it and which script
produces which panel. `metabarcoding/README.md` gives every
MICCA command with its parameters, for both markers.

## What the RNA-seq code does, end to end

**1. Reads to counts** — `RNA-seq/01_preprocessing_alignment/`
SLURM batch scripts as executed on the FEM HPC cluster, every tool inside an
Apptainer image: FastQC 0.12.1 for read QC, fastp 0.23.4 for adapter and
quality trimming, STAR 2.7.11b against *Solanum lycopersicum* SL5.0 with the
ITAG5.0 annotation, samtools 1.22.1, and featureCounts 2.1.1 with
`-p -s 2 --countReadPairs -t exon -g Parent` for reverse-stranded gene-level
counting. The scripts are deposited as executed, cluster paths and all, so the
deposited record matches what produced the data. Two points of provenance that
a reader could otherwise only guess at — which of two `featureCounts`
invocations produced the published matrix, and which version was used — are
documented in that directory's README.

**2. Technical validation** — `RNA-seq/02_rnaseq_qc/`
Read fate through the pipeline for all 60 libraries (Figure 2a), and a
sample-to-sample distance heatmap: blind variance-stabilising transformation of
the raw counts, Euclidean distances between libraries, rows and columns ordered
by complete-linkage clustering of the distance matrix (Figure 2b). Also the
supplementary QC table and the supporting QC panels — mapping rate, replicate
correlation, marker-gene expression, gene-detection saturation.

**3. Quantifying trial, time point and treatment** — `RNA-seq/03_batch_effect/`
Added at revision. The Bioconductor package **variancePartition** decomposes
the transcriptome variance under the model `~ trial * time point + condition`,
with sum-to-zero contrasts, on 18,425 genes retained by `filterByExpr` and with
no most-variable-gene selection. Four products: the canonical correlations
between the factors, the decomposition gene by gene, the same decomposition
applied axis by axis to the principal components, and a moderated *F*-test per
model term. `build_supplementary_table.R` assembles these into Supplementary
Table 6, one workbook with a legend sheet. The directory README explains the
method, why the contrast coding is load-bearing, and why significance is
reported per gene and not per axis.

**4. Figure assembly** — `RNA-seq/04_figures/`
Composes Figure 2 from the panels the analysis scripts save, in two layouts:
the single-page version at 180 × 202 mm submitted with the revision
(`build_figure2_pagefit.R`), and a tall alternative (`build_figure2.R`) with
different panel lettering. Both read the same panel objects, so neither can
drift from the analysis behind it or from the other.

## What the RNA-seq analysis found

The design is fully crossed and balanced: 3 conditions × 2 time points ×
2 trials × 5 biological replicates = 60 libraries, exactly five per cell. All
canonical correlations between trial, time point and condition are zero (the
largest off-diagonal value is 6 × 10⁻¹⁷), so the trial is crossed with the
design factors, not confounded with them, and its effect is separable however
large it is.

Variance explained, across all 18,425 genes:

| Source | Median gene | Variance-weighted | Genes significant at FDR < 0.05 |
|---|---:|---:|---:|
| Time point | 11.7% | 17.6% | 12,333 (66.9%) |
| Experimental trial | 9.9% | 13.1% | 11,859 (64.4%) |
| Trial × time point | 9.2% | 17.4% | 11,471 (62.3%) |
| Treatment (condition) | 1.6% | 2.3% | **0 (0.0%)** |
| Residual, between plants | 50.1% | 49.5% | — |

Axis by axis, the first three principal components carry 35.9%, 13.6% and 10.8%
of total transcriptome variance. Trial, time point and their interaction account
for 91% of PC1 and 97% of PC2; treatment accounts for less than 1% of either.

The trial effect is therefore substantial, is larger at the first sampling
point than at the second — which is why the interaction term is modelled — and
is fully separable from the design. Treatment accounts for a small and, gene by
gene, undetectable share of leaf transcriptome variance at either time point.

## What the metabarcoding code does

`metabarcoding/` holds the rhizosphere pipeline, built on MICCA 1.7.2 and run
separately for the prokaryotic (16S rRNA, V4–V5) and fungal (ITS1) libraries:
paired-end merging, primer trimming, length and expected-error filtering,
UNOISE denoising to amplicon sequence variants, RDP Classifier taxonomy (2.14
for 16S; 2.13 trained on UNITE+INSD 8.3 for ITS), rarefaction to 68,279 (16S) and
44,822 (ITS) reads, NAST or MUSCLE alignment and a rooted tree, and export to
BIOM. `plots_validation.R` draws the read-depth distributions and ASV
rarefaction curves from the exported BIOM, tree and sequence files.

## Reproducing the RNA-seq analysis

R 4.6.1 with Bioconductor. From `RNA-seq/`, after placing the count matrix from
Zenodo in `data/` (see `RNA-seq/data/README.md`):

```bash
./run_all.sh
```

runs the whole downstream analysis in dependency order in about a minute and a
half and ends by writing the assembled Figure 2. No absolute path appears
anywhere in the R code: every script resolves its paths through
`RNA-seq/config/paths.R`. Package versions for the archived run are recorded in
`RNA-seq/env/r_packages.tsv` and `RNA-seq/env/sessionInfo.txt`, both
regenerated by `RNA-seq/env/record_environment.R`.

The outputs are committed alongside the code, so the figures and tables the
descriptor reports can be checked against the scripts that made them without
running anything.

## Licence and citation

Code: MIT (`LICENSE`), covering both `RNA-seq/` and `metabarcoding/`. Data on
Zenodo: CC-BY 4.0. Citation metadata for
this repository is in `CITATION.cff`.

Each release of this repository is archived on Zenodo. Cite the concept DOI,
[10.5281/zenodo.22232657](https://doi.org/10.5281/zenodo.22232657), which
resolves to the latest release; the release pages carry their own version DOIs
if you need to pin an exact snapshot.
