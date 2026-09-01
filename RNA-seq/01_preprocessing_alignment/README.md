# HPC preprocessing and alignment

SLURM batch scripts as executed on the Fondazione Edmund Mach cluster. All
tools run inside Apptainer (formerly Singularity) images; there is no module
system on that cluster.

| Script | Tool | Purpose |
|---|---|---|
| `01_run_fastqc.slurm` | FastQC 0.12.1 | raw read QC |
| `02_run_fastp.slurm` | fastp 0.23.4 | adapter and quality trimming (default parameters) |
| `03_run_star_samtools.slurm` | STAR 2.7.11b, samtools 1.22.1 | splice-aware alignment to SL5.0 / ITAG5.0, coordinate sort, index |
| `04_run_featurecounts.slurm` | Subread 2.1.1 | `featureCounts -p -s 2 --countReadPairs -t exon -g Parent`, gene-level counts |
| `legacy/run_fastp_hisat2_samtools.slurm` | HISAT2 2.2.1 | **not used for the published data** |

## Re-running elsewhere

Each script declares its cluster paths as variables in the header block. To run
on another system, edit only these:

```
GENOME_DIR / SUBREAD_INDEX   built index for SL5.0
ANNOTATION                   Slycopersicum_796_ITAG5.0.gene_exons.gff3
FASTQ_DIR                    input reads
OUT_DIR                      output directory
CONTAINER                    path to the Apptainer image
```

They are kept as executed rather than rewritten, so that the deposited record
matches what actually produced the data.

## Two caveats, recorded honestly

1. **`04_run_featurecounts.slurm` is the array-job template**, and its
   `SUBREAD_INDEX` variable still points at a grapevine (PN40024) index left
   over from an earlier project, while `expected_samples` is set to 6. The
   published count matrix was **not** produced by this array path. It came from
   a single combined `featureCounts` invocation over the 60 STAR BAM files —
   the block at the end of the same script, writing
   `star_results/fastp/counts/counts.txt`. The `featureCounts` flags in both
   blocks are identical (`-p -s 2 --countReadPairs -t exon -g Parent`), which
   is what determines the counts; only the driver differs. Verified against
   `report_data_sources` in the MultiQC report.

2. **featureCounts version.** The counts were produced with **v2.1.1**. The
   count file itself records it on its first line
   (`# Program:featureCounts v2.1.1; Command:"featureCounts" "-T" "14" "-t"
   "exon" "-g" "Parent" "-p" "-s" "2" "--countReadPairs" ...`), and the
   Apptainer image used carries the same version in its name
   (`hisat2.2.1_samtools_1.22.1_subread_2.1.1.sif`).

## Reference

*Solanum lycopersicum* v5.0, assembly SL5.0, annotation ITAG5.0, from
Phytozome JGI:
<https://phytozome-next.jgi.doe.gov/info/Slycopersicum_ITAG5_0>
