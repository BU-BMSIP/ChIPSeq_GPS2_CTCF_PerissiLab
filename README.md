# ChIPSeq_GPS2_CTCF_PerissiLab

**BU-BMSIP / Perissi Lab (Boston University School of Medicine)**  
Repository for ChIP‑seq and RNA–chromatin interaction analysis focused on GPS2, ATF4, and CTCF under mitochondrial stress.

## Overview

Mitochondrial stress triggers retrograde signaling to the nucleus. GPS2 relocates to chromatin, but the recruitment mechanics—especially in conjunction with ATF4/CTCF—remain unclear. This repository implements:

- reproducible Snakemake workflows for ChIP‑seq and iMARGI
- integrated analysis of in-house GPS2 with public ATF4/CTCF datasets
- peak intersection, motif discovery, signal profiling, and RNA–DNA interaction mapping

## Repository Layout  
├── docs/ # Documentation, diagrams, protocol notes  
├── envs/ # Conda environments for workflows  
├── notebooks/ # Analysis Jupyter notebooks  
├── profile/ # Profiling/logs for workflow runs  
├── scripts/ # Python/R/Bash utilities (e.g., annotation, parsing)  
├── Peak_Calls_with_Control.csv  
├── Peak_Calls_without_Control.csv  
├── Sample_Data_for_CTCF.csv  
├── Updated_Sample_Sheet.csv  
└── README.md  


## Prerequisites
- pending

Key Outputs
*.bam, *.bw – aligned reads and coverage tracks

*_peaks.narrowPeak – GPS2 / ATF4 / CTCF peak calls

common_genes/…_overlap.txt – peak intersection gene lists

plots/…_profile_matrix.gz, .png – deepTools signal matrices and profile plots

motifs/…/knownResults.txt – HOMER motif enrichment results

iMargi/output/*.pairs.gz – parsed RNA–DNA interaction files

Findings(pending)
GPS2 and ATF4 co-bind stress and metabolism-related genes

CTCF motifs are enriched at GPS2 flanking regions, suggesting structural roles

Condition-specific differences observed during adipocyte differentiation (e.g., day 0 vs day 6)

iMARGI output .pairs.gz files generated and undergoing downstream filtering