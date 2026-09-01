# Sample metadata

`targets.txt` — tab-separated, 60 rows, one per RNA-seq library.

| Column | Meaning |
|---|---|
| `Sample_name` | human-readable label, e.g. `EXP1_CTR_T2_r1` |
| `Condition` | `Cntr` (water only), `E1_D3`, `E1_D4` (two doses of extract E1) |
| `Time` | `T1` = 3 weeks / 21 d.p.t.; `T2` = 6 weeks / 42 d.p.t. |
| `Experiment` | `1` (trial started November 2024) or `2` (February 2025) |
| `Sample` | library ID, matches the column names of `RNAseq_rawcounts.tsv` |

The design is fully crossed and balanced: 3 Condition x 2 Time x 2 Experiment
x 5 biological replicates (independent plants) = 60 libraries. Every cell holds
exactly 5 replicates; `results/batch_effect/design_balance.txt` is the
machine-generated proof.

**T1 and T2 are not a clean time course.** E1 was applied at 2, 9, 23 and
30 d.p.t. T1 plants (sampled 21 d.p.t.) had received only the first two
applications; T2 plants received all four. T1 and T2 therefore differ in both
developmental stage and cumulative exposure.

Mapping between these library IDs and the SRA BioSample / Run accessions is in
Supplementary Table 1 of the data descriptor.
