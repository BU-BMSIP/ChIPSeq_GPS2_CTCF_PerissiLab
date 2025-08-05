# Bioinformatics Pipeline Scripts

The table below lists each script (using the updated file names) and a brief description of its purpose.

| File Name                                          | Purpose                                                                                                           |
|----------------------------------------------------|-------------------------------------------------------------------------------------------------------------------|
| **generate_common_peak_flanks.ipynb**              | Extend “common” ChIP-seq peaks by a series of flank sizes (e.g. 100–500 bp) using `bedtools slop`, then subtract the original peaks to retain only the flanking regions. |
| **calculate_FRiP_scores.ipynb**                    | Compute FRiP (Fraction of Reads in Peaks) for CTCF BAM replicates: count total vs. mapped reads (with `pysam`), count reads in MACS3 peaks (with `pybedtools`), and summarize in a Pandas table. |
| **CommonGeneFilter.R**                             | Filter an annotated peak table to keep only rows whose gene symbol appears in a “common genes” list, then write both a filtered annotation table and a minimal BED (chr/start/end) for visualization. |
| **convert_narrowPeak_mm10_to_mm39.ipynb**          | “Lift” a MACS3 narrowPeak file from mm10 to mm39 coordinates using CrossMap and a chain file, with automatic installation of CrossMap if needed. |
| **filter_bed_by_ensembl_ids.ipynb**                | Translate RefSeq transcript IDs in column 4 of a BED file to Ensembl gene IDs (using NCBI’s `mouse_gene2ensembl.tsv`), stripping version suffixes and saving a new BED. |
| **annotate_and_find_common_genes.R**               | Annotate promoter peaks for GPS2, CTCF, GPS2 day 6, and ATF4 with ChIPseeker (mm39), extract gene SYMBOLs, and compute pairwise & three-way gene set intersections. |
| **filter_promoter_peaks_by_common_genes.R**        | From an annotated peak table (e.g. GPS2 day 6), select only promoter‐annotated rows whose SYMBOL is in a “common genes” list, and output a BED of seqnames/start/end. |
| **extract_promoter_peaks.R**                       | Batch‐process any annotation file (CTCF t1–t4, GPS2 silencing, ATF4, etc.) to extract “Promoter” peaks and write promoter-only BEDs via a reusable `extract_promoter_peaks()` function. |
| **convert_refseq_to_ensembl_bed.ipynb**            | Replace RefSeq IDs in the fourth column of a genome-wide BED with Ensembl gene IDs by mapping via the `mouse_gene2ensembl.tsv` table and saving the updated BED. |
