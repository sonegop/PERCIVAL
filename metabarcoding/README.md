# PERCIVAL - Rhizosphere Metabarcoding Pipeline (16S & ITS)

This section contains the bioinformatic pipelines used for processing the rhizosphere metabarcoding data for the PERCIVAL project, targeting both the Prokaryotic (16S rRNA) and Fungal (ITS) communities. The pipelines are built using **MICCA (v.1.7.2)** and integrate **VSEARCH**, **NAST**, **MUSCLE**, and **RDP Classifier**.

## Prerequisites & Software Versions
* **MICCA** v.1.7.2
* **RDP Classifier** v.2.14, updated to training set n.19 (used for 16S), and v.2.13 (custom-trained for ITS)
* **SILVA Database** v.132 (Prokaryotes reference)
* **UNITE Database** Full UNITE+INSD dataset for eukaryotes, version 10.05.2021 (Fungi reference)
* **Java** (for running the standalone RDP Classifier)

---

## 1. Quality Control & Sequence Processing

### Step 1: Paired-End Reads Merging (Common to 16S and ITS)
Merge forward (`R1`) and reverse (`R2`) Illumina MiSeq reads with a minimum overlap length of 100 bp and a maximum of 32 allowed mismatches:
```bash
micca mergepairs \
  -i *_R1*.fastq \
  -o merged.fastq \
  --notmerged-fwd notmerged_fwd.fastq \
  --notmerged-rev notmerged_rev.fastq \
  -l 100 -d 32
```

### Step 2: Primer Trimming
Trim specific universal primers from the merged datasets.

*   **For Prokaryotes (16S V4–V5 region):**
    ```bash
    micca trim \
      -i merged.fastq \
      -o trimmed.fastq \
      -w GTGCCAGCMGCCGCGGTAA \
      -r CCCCGYCAATTCMTTTRAGT \
      -W -R -c
    ```
*   **For Fungi (ITS1 region):**
    ```bash
    micca trim \
      -i merged.fastq \
      -o trimmed.fastq \
      -w CTTGGTCATTTAGAGGAAGTAA \
      -r GCTGCGTTCTTCATCGATGC \
      -W -R -c
    ```

### Step 3: Quality Filtering
Generate pre-filtering stats and remove low-quality sequences using custom parameters for length and expected error rate (`-e`).

*   **For Prokaryotes (16S):** Reads shorter than 350 bp or with an expected error rate higher than 0.25% are discarded.
    ```bash
    micca filterstats -i fastq/trimmed.fastq -o filterstats
    micca stats -i fastq/trimmed.fastq -o stats -n 10000
    micca filter -i fastq/trimmed.fastq -o fastq/filtered.fasta -e 0.25 -m 350
    ```
*   **For Fungi (ITS):** Reads shorter than 120 bp or with an expected error rate higher than 0.5% are discarded.
    ```bash
    micca filterstats -i fastq/trimmed.fastq -o filterstats
    micca stats -i fastq/trimmed.fastq -o stats -n 10000
    micca filter -i fastq/trimmed.fastq -o fastq/filtered.fasta -e 0.5 -m 120
    ```

---

## 2. Denoising & ASV Generation (Common to 16S and ITS)

Identification of Amplicon Sequence Variants (ASVs) from the filtered sequences using the **UNOISE** de novo algorithm available in MICCA (utilizing VSEARCH):
```bash
micca otu \
  -m denovo_unoise \
  -i fastq/filtered.fasta \
  -o denovo_greedy_unoise \
  -t 10
```

---

## 3. Taxonomic Assignment

### For Prokaryotes (16S)
Taxonomic classification executed using **RDP Classifier (v2.14)**, updated to training set n.19, against the SILVA reference database with a confidence cutoff of 0.5:
```bash
java -Xms20g -Xmx30g -jar /micca/db/rdp_classifier_2.14/dist/classifier.jar classify \
  -g 16srrna \
  -c 0.5 \
  -f fixrank \
  -o rdp.output \
  /micca/PERCIVAL/16S/denovo_greedy_unoise/otus.fasta
```

### For Fungi (ITS)
Taxonomic classification executed using **RDP Classifier (v2.13)** custom-trained on the Full UNITE+INSD v8.3 database for eukaryotes (Version 10.05.2021, Porter 2021):
```bash
java -Xms20g -Xmx30g -jar /micca/db/rdp_classifier_2.13/dist/classifier.jar \
  -t /micca/db/mydata_trained/rRNAClassifier.properties \
  -o rdp.output \
  /micca/PERCIVAL/ITS/denovo_greedy_unoise/otus.fasta
```

---

## 4. Downstream Table Processing & Rarefaction

Compute abundance table statistics and rarefy to standard baseline sampling depths based on dataset-specific technical thresholds.

*   **For Prokaryotes (16S):** Subsampled to a uniform depth of **68,279** reads.
    ```bash
    micca tablestats -i denovo_greedy_unoise/otutable.txt -o tablestats
    micca tablerare -i denovo_greedy_unoise/otutable.txt -o denovo_greedy_unoise/otutable_rare.txt -d 68279
    micca tabletotax -i denovo_greedy_unoise/otutable_rare.txt -t denovo_greedy_unoise/taxa.txt -o taxtables
    micca tablebar -i taxtables/taxtable5.txt -o taxtables/taxtable5.png
    ```
*   **For Fungi (ITS):** Subsampled to a uniform depth of **44,822** reads.
    ```bash
    micca tablestats -i denovo_greedy_unoise/otutable.txt -o tablestats
    micca tablerare -i denovo_greedy_unoise/otutable.txt -o denovo_greedy_unoise/otutable_rare.txt -d 44822
    micca tabletotax -i denovo_greedy_unoise/otutable_rare.txt -t denovo_greedy_unoise/taxa.txt -o taxtables
    micca tablebar -i taxtables/taxtable6.txt -o taxtables/taxtable6.png
    ```

---

## 5. Phylogenetic Inference

Build a rooted phylogenetic tree for downstream phylogenetic metrics (e.g., UniFrac distances).

*   **For Prokaryotes (16S):** Multiple Sequence Alignment executed via the **NAST** algorithm.
    ```bash
    micca msa -m nast -i denovo_greedy_unoise/otus.fasta -o denovo_greedy_unoise/msa.fasta --nast-template /micca/db/core_set_aligned.fasta.imputed --nast-threads 10
    micca tree -i denovo_greedy_unoise/msa.fasta -o denovo_greedy_unoise/tree.tree
    micca root -i denovo_greedy_unoise/tree.tree -o denovo_greedy_unoise/tree_rooted.tree
    ```
*   **For Fungi (ITS):** Multiple Sequence Alignment executed via the **MUSCLE** algorithm.
    ```bash
    micca msa -m muscle -i denovo_greedy_unoise/otus.fasta -o denovo_greedy_unoise/msa.fasta
    micca tree -i denovo_greedy_unoise/msa.fasta -o denovo_greedy_unoise/tree.tree
    micca root -i denovo_greedy_unoise/tree.tree -o denovo_greedy_unoise/tree_rooted.tree
    ```

---

## 6. Exporting to BIOM Format (Common to 16S and ITS)

Combine the processed ASV abundance tables, specific taxonomy classifications, and the unified sample metadata sheet (`sample_metadata.txt`) into a standardized **BIOM** format for downstream statistical computing in R (e.g., `phyloseq`):
```bash
micca tobiom \
  -i otutable.txt \
  -o output.biom \
  -t taxa.txt \
  -s sample_metadata.txt
```

---

## 7. Downstream Quality Visualization in R

To reproduce the read-depth distribution graphs and the ASV rarefaction curves presented in the technical validation section of the manuscript, you can run the provided R script available in this repository:
*   File: `plots_validation.R`
*   Required R packages: `phyloseq`, `vegan`, `ggplot2`

Simply open the script, configure your local working directory path, and execute it to generate high-resolution validation plots (`Supplementary_Figure_S1` and `S2`) for both 16S and ITS datasets.

