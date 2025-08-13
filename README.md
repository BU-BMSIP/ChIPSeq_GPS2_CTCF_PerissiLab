# ChIPSeq_GPS2_CTCF_PerissiLab

**BU-BMSIP / Perissi Lab (Boston University Chobanian & Avedisian School of Medicine)**  
Reproducible multi-omics analysis framework for ChIP-seq and RNA–chromatin interaction (iMARGI) data, focusing on GPS2-mediated mitochondrial retrograde signaling and its interaction with ATF4, CTCF, and related transcription factors.

---

## 1. Project Overview

Mitochondria depend on nuclear-encoded proteins, requiring tight coordination of gene expression via **anterograde** (nucleus→mitochondria) and **retrograde** (mitochondria→nucleus) signaling.  
GPS2 (G-protein Pathway Suppressor 2) is a key retrograde mediator that translocates to chromatin under mitochondrial stress, potentially cooperating with stress-response TFs such as **ATF4/ATF5** and the chromatin organizer **CTCF**.

This repository contains:

- **Snakemake-based workflows** for scalable, reproducible processing of:
  - GPS2, ATF4, CTCF ChIP-seq datasets (in-house + public)
  - iMARGI RNA–DNA interaction data (including mtRNA–DNA contacts)
- **Integrated analyses**:
  - Peak calling, intersection, and annotation
  - Motif discovery and promoter topology analysis
  - Signal profiling and heatmaps (deepTools)
  - Functional enrichment (GO/KEGG, GSEA)
  - Mapping of mtRNA–DNA interactions and overlap with TF peaks

---

## 2. Repository Layout

```
├── docs/                # Documentation, diagrams, protocol notes
├── envs/                # Conda environment YAMLs (bedtools, deeptools, macs3, homer, etc.)
├── notebooks/           # Jupyter/R notebooks (QC, annotation, enrichment, motif)
├── profile/             # Snakemake profiles (threads, paths, runtime configs)
├── scripts/             # Python/R/Bash utilities & Snakemake workflow files
│   ├── ChIP_CTCF_mouse_cleaned.smk
│   ├── iMargi.smk
│   └── ...
├── adapters_and_annotations/  # Reference FASTA/BED, genome annotations, blacklists
├── Peak_Calls_with_Control.csv
├── Peak_Calls_without_Control.csv
├── Sample_Data_for_CTCF.csv
├── Updated_Sample_Sheet.csv
└── README.md
```

---

## 3. Prerequisites

- Access to **BU Shared Computing Cluster (SCC)** or equivalent HPC
- Conda (≥4.10)
- Snakemake (≥7.0)
- Required tools: `bedtools`, `deeptools`, `macs3`, `homer`, `samtools`, `pairtools`, `bwa`

Activate environment before running:

```bash
conda activate chip_seq
```

---

## 4. Running the Pipeline

**Dry run first (recommended):**
```bash
snakemake -s scripts/ChIP_CTCF_mouse_cleaned.smk --profile profile -np
```
**Actual run:**
```bash
snakemake -s scripts/ChIP_CTCF_mouse_cleaned.smk --profile profile
```
For iMARGI workflow:
```bash
snakemake -s scripts/iMargi.smk --profile iMargi
```

---

## 5. Adding New Data

### ChIP-seq
1. **Place FASTQ files**:  
   - Raw (unmerged): `CTCF_3T3L1/raw_samples/`  
   - Pre-merged: `CTCF_3T3L1/samples/`
2. **Or** add FTP links + SRR IDs to `Sample_Data_for_CTCF.csv` for auto-download & merge
3. Register sample in `Updated_Sample_Sheet.csv`
4. Update peak calling configs:  
   - With control → `Peak_Calls_with_Control.csv`  
   - Without control → `Peak_Calls_without_Control.csv`

### iMARGI
- Edit `scripts/iMargi.smk` → update "Global config" with new SRA IDs  
- Follow iMARGI preprocessing steps (cleaning, mapping, parsing)

---

## 6. Key Outputs

- **Alignment & QC**  
  `*.bam`, `multiqc_report.html`
- **Coverage tracks**  
  `*.bw`
- **Peak calls**  
  `*_peaks.narrowPeak`
- **Annotation**  
  `results/annotation/*.bed` & gene lists
- **Motifs**  
  `motifs/**/knownResults.txt`
- **Signal profiling**  
  `matrix/*.gz`, `plots/*.png`
- **Enrichment analysis**  
  `results/Enrichment/*.tsv`
- **iMARGI outputs**  
  `output/final_*.pairs.gz` + promoter overlap stats

---

## 7. Example Analyses

- GPS2–ATF4 co-binding during adipocyte differentiation  
- CTCF motif enrichment flanking GPS2 peaks  
- Condition-specific GPS2 occupancy shifts (day 0 vs day 6)  
- mtRNA-binding sites overlapping GPS2/NCOR peaks (T263 endothelial cells)

---

## 8. References

1. Cardamone et al., *Mol Cell*, 2018 — GPS2 retrograde signaling mechanism  
2. Chen et al., *Cell Biol Toxicol*, 2022 — ATF4–CTCF cooperation in adipogenesis  
3. iMARGI Pipeline — [Yan Lab Documentation](https://sysbiocomp.ucsd.edu/public/frankyan/imargi_pipeline/index.html)
