# Input data

Large matrices are not versioned here. Download them from Zenodo
(<https://doi.org/10.5281/zenodo.19336336>) into this directory:

| File expected here | Zenodo file | Used by |
|---|---|---|
| `RNAseq_rawcounts.tsv` | `RNAseq_rawcounts.tsv` | 02, 03, 04 |
| `multiqc_data/multiqc_data.json` | not on Zenodo — regenerate with MultiQC over the FastQC/fastp/STAR/featureCounts logs | 02 |
| `read_retention_merged.tsv` | derived; small, versioned here | 02 |

Raw reads are at NCBI SRA, BioProject **PRJNA1425279**.

> `RNAseq_rawcounts.tsv` and `multiqc_data/` are present in this working copy
> for convenience, so the scripts run out of the box. They are **deliberately
> excluded from version control** (see `.gitignore`) because they are already
> archived on Zenodo: the deposit is code, not a second copy of the data. Their
> absence from a fresh `git clone` is expected, not a gap in the deposit.

`read_retention_merged.tsv` is a small join of the per-sample pipeline counters
(raw reads, reads surviving fastp, uniquely mapped reads from STAR, alignments
assigned by featureCounts) and is kept in the repository because it is the
direct input to Figure 2a.

Point `PERCIVAL_DATA` elsewhere if you keep the matrices outside the repository:

```bash
export PERCIVAL_DATA=/path/to/matrices
```
