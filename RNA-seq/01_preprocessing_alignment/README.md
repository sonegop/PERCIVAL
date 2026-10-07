# HPC preprocessing and alignment

SLURM batch scripts as executed on the Fondazione Edmund Mach cluster. All
tools run inside Apptainer (formerly Singularity) images; there is no module
system on that cluster.

| Script | Tool | Purpose |
|---|---|---|
| `01_run_fastqc.slurm` | FastQC 0.12.1 | raw read QC |
| `02_run_fastp.slurm` | FastQC 0.12.1 | despite its name, a second FastQC run on the sequencing files; **it does not run fastp** (see below) |
| `legacy/run_fastp_hisat2_samtools.slurm` | fastp 0.23.4, HISAT2 2.2.1 | **its fastp step produced the clean reads of the published data**; its HISAT2 alignment step was not used |
| `03_run_star_samtools.slurm` | STAR 2.7.11b, samtools 1.22.1 | two-pass alignment to SL5.0 / ITAG5.0 on the fastp-cleaned reads, coordinate sort, index |
| `04_run_featurecounts.slurm` | Subread 2.1.1 | `featureCounts -p -s 2 --countReadPairs -t exon -g Parent`, gene-level counts |

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

## Parameters of the published run

Neither fastp nor STAR was run with default parameters.

**fastp 0.23.4** (read pre-processing; step run by the fastp block of
`legacy/run_fastp_hisat2_samtools.slurm`, whose HISAT2 block was not used). The
command recorded in all 60 fastp reports of the MultiQC report is:

```
fastp --in1 <R1> --in2 <R2> --out1 <clean_1> --out2 <clean_2> --html <report> --json <report>
      --thread 16 --length_required 50 --qualified_quality_phred 20
      --cut_window_size 4 --cut_mean_quality 20 --trim_poly_g --trim_poly_x
      --correction --detect_adapter_for_pe
```

That is: paired-end adapter detection, sliding-window quality trimming (window
4, mean Phred 20), poly-G and poly-X tail trimming, overlap-based base
correction, qualified-base threshold Phred 20 and a minimum read length of
50 bp. The same block in the legacy script sets exactly these options. Its
outputs, `*_clean_1.fq.gz` and `*_clean_2.fq.gz`, are what `03_run_star_samtools.slurm`
reads. On average 97.6% of read pairs were retained (range 94.1-98.5%, 60
libraries).

**STAR 2.7.11b** (`03_run_star_samtools.slurm`), splice junctions from the
ITAG5.0 annotation:

```
STAR --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate
     --outSAMunmapped Within --outSAMattributes NH HI AS nM MD
     --outFilterType BySJout --outFilterMultimapNmax 20
     --outMultimapperOrder Random --outSAMmultNmax 1
     --alignSJoverhangMin 8 --alignSJDBoverhangMin 1
     --alignIntronMin 20 --alignIntronMax 20000 --alignMatesGapMax 1000000
     --sjdbGTFfile <ITAG5.0 gff3> --sjdbGTFtagExonParentTranscript Parent
     --sjdbGTFtagExonParentGene ID --sjdbOverhang 149
     --quantMode GeneCounts --twopassMode Basic
```

On average 93.7% of the clean read pairs mapped uniquely (range 85.0-96.6%).
The comment "for Vitis" next to `MAX_INTRON` in that script is a leftover from the
project it was adapted from; the value in force is 20,000 bp.

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
