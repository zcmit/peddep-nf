# ATAC-seq (`--assay atac`)

Paired-end ATAC-seq, from FASTQ to peaks, coverage tracks and a graded QC report. See the [README](../README.md) for installation, references and running at your site.

## Samplesheet

```csv
sample,fastq_1,fastq_2,outdir
SAMPLE_A,/data/SAMPLE_A_S1_L001_R1_001.fastq.gz,/data/SAMPLE_A_S1_L001_R2_001.fastq.gz,
SAMPLE_A,/data/SAMPLE_A_S1_L002_R1_001.fastq.gz,/data/SAMPLE_A_S1_L002_R2_001.fastq.gz,
SAMPLE_B,s3://bucket/SAMPLE_B_R1.fastq.gz,s3://bucket/SAMPLE_B_R2.fastq.gz,/projects/cellline_B/ATAC/SAMPLE_B
```

## What the workflow does

| Step | Tool |
|---|---|
| Raw read QC | FastQC 0.13 |
| Adapter trimming (Nextera by default), post-trim QC | Trim Galore 0.6.10 / cutadapt |
| Alignment, per lane | bowtie2 2.5.4 `--very-sensitive -X 2000 --no-mixed --no-discordant` |
| Lane merge, duplicate marking | samtools 1.24, Picard 3.5 MarkDuplicates |
| Library complexity (NRF, PBC1, PBC2) | ENCODE definitions, on fragments |
| Filtering: proper pairs, MAPQ >= 30, no duplicates, mito, blacklist or non-canonical contigs, no orphan mates | samtools + pysam |
| Tn5 shift (+4 / -5) | deepTools `alignmentSieve --ATACshift` |
| Peaks: narrow (q < 0.01, summits) and broad | MACS3 3.0.4, `BAMPE` (fragment-aware) |
| FRiP, peak count | samtools |
| Coverage: fragments (CPM, 10 bp) and Tn5 insertion sites (CPM, 1 bp) | deepTools `bamCoverage` |
| TSS enrichment, mito fraction, fragment periodicity | ataqv 1.3.1 (+ `mkarv` HTML viewer) |
| TSS profile/heatmap, fingerprint | deepTools |
| Fragment-length distribution, NFR fraction, nucleosome pattern | `bin/fragment_qc.py` |
| PEDDEP graded QC report and run summary | `bin/atac_qc_report.py`, `bin/qc_summary.py` |
| Everything in one report | MultiQC 1.35 |

### Outputs

```
<outdir>/
├── atac_qc_summary.tsv / .html     # one row per sample, PEDDEP grades
├── multiqc/multiqc_report.html
├── ataqv_report/index.html         # interactive ataqv viewer
├── pipeline_info/                  # execution report, timeline, trace, DAG, software versions
└── <sample>/                       # or the samplesheet's outdir
    ├── <sample>.atac_qc_report.txt / .html / .atac_qc_metrics.json
    ├── BAM_files/                  # <sample>.filtered.bam, <sample>.filtered.tn5.bam (+ .bai)
    ├── PEAKS/NARROW_peaks/         # _peaks.narrowPeak, _summits.bed, _peaks.xls, _peaks.bed.gz + .tbi (browser track)
    ├── PEAKS/BROAD_peaks/          # _peaks.broadPeak, _peaks.gappedPeak, _peaks.xls
    ├── COVERAGE_files/             # .fragments.CPM.bw, .tn5_insertions.CPM.bw
    └── QC/
        ├── FastQC/{raw,trimmed}/   trim_report/   alignment/   duplicates/
        ├── library_complexity/     peaks/ (FRiP)  Fragments/   TSS/
        └── ataqv/                  fingerprint/
```

### Genome-browser tracks

Each sample has two files ready for IGV or a ProteinPaint display: `COVERAGE_files/<sample>.fragments.CPM.bw` (bigwig track) and `PEAKS/NARROW_peaks/<sample>_peaks.bed.gz` with its `.tbi` (bedj track, sorted 3-column BED). They take the place of the WDL pipeline's `<sample>.sorted.bklist.rmdup.w50.FPM.bw` and `<sample>.sorted.bklist_p9_kd-auto_peaks.bed.gz`. For a new sample, copy the two files (plus the `.tbi`) and point the display's `file` entries at the new names. `<sample>.tn5_insertions.CPM.bw` is an extra 1-bp Tn5 insertion track for footprint-level views.

Optional: `--save_markdup_bam`, `--save_trimmed`, `--save_reference`; `--save_tn5_bam false` skips publishing the shifted BAM.

### PEDDEP graded QC report

This keeps the grading scheme of the PEDDEP/SEAseq reports. Each graded metric scores from -2 to +2, and `Overall_Quality` is the mean score: `>= 2` EXCELLENT, `>= 1` GOOD, `>= 0` AVERAGE, `>= -1` BELOW_AVERAGE, otherwise POOR. Sequencing-quality metrics come first.

| Metric | Source | Cutoffs for -1 / 0 / +1 / +2 |
|---|---|---|
| Raw_Read_Pairs | FastQC (raw) | 15M / 20M / 25M / 30M |
| Base_Quality | FastQC per-base quality, worst FASTQ | fail -2, warn 0, pass +2 |
| Sequence_Diversity | FastQC overrepresented sequences, worst FASTQ | fail -2, warn 0, pass +2 |
| Pct_Reads_Q30 | FastQC per-sequence quality, worst FASTQ | 60 / 70 / 75 / 80 % |
| Aligned_Percent | bowtie2 | 50 / 60 / 70 / 80 % |
| Mito_Percent (lower is better) | samtools idxstats | 40 / 20 / 10 / 5 % |
| NRF | library complexity | 0.5 / 0.6 / 0.7 / 0.8 |
| PBC1 | library complexity | 0.5 / 0.66 / 0.75 / 0.9 |
| PBC2 | library complexity | < 1 → -2, >= 1 → 0, >= 3 → +2 |
| Nucleosome_Pattern | fragment QC | fail -2, warn 0, pass +2 |
| TSS_Enrichment | ataqv | 3 / 5 / 7 / 10 |
| FRiP | narrow peaks | 0.10 / 0.15 / 0.20 / 0.30 |
| Total_Peaks | narrow peaks | 10k / 20k / 50k / 100k |

`Mean_Read_Length`, `Adapter_Trimmed_Pct`, `Duplication_Pct`, `Final_Read_Pairs`, `Median_Fragment` and `NFR_Fraction` are reported but not graded. The cutoffs live in [`assets/atac_qc_thresholds.json`](../assets/atac_qc_thresholds.json); pass your own file with `--qc_thresholds` (for example for small genomes or shallow pilot runs).

## Changes from the PEDDEP WDL ATAC pipeline

| | WDL (`atacpipeline.wdl`) | peddep-nf |
|---|---|---|
| Runs on | St. Jude LSF + Cromwell only | anywhere Nextflow runs |
| Input | one `python sj-atacdev.py` call per sample | one samplesheet for many samples and lanes |
| Trimming minimum length | max read length − adapter length (≈57 bp): drops short nucleosome-free fragments | 20 bp |
| Aligner | bowtie1 `-k 2 -m 2` (keeps reads with two hits) | bowtie2, then MAPQ >= 30 |
| chrM reads | kept | removed (reported as Mito_Percent) |
| Tn5 shift | not done | `alignmentSieve --ATACshift` |
| Peak caller | MACS 1.4.2 (reads, three variants) + SICER | MACS3 BAMPE narrow + broad |
| "single-end mode" re-analysis of R1 and R2 | yes | dropped |
| ROSE, BAM2GFF, peaksanno, TDF tracks | yes | dropped (deepTools TSS plots instead of BAM2GFF) |
| ATAC QC | NSC/RSC (ChIP metrics), FRiP | TSS enrichment, fragment periodicity, mito %, NRF/PBC, FRiP, MultiQC |
| Resume after failure | no (call caching off) | `-resume` |

## Comparison with the WDL pipeline

The same four PEDDEP ATAC-seq libraries (NextSeq 2000, 2×76 bp, 3–27 M read pairs) were run through both pipelines on GRCh38. Peak overlap uses the WDL's MACS 1.4 `p9_kd-auto` narrow peaks against MACS3 narrowPeak. Coverage correlation compares the WDL's 50-bp fragment FPM bigWig with the new fragment CPM bigWig, in 5-kb bins with chrM excluded.

| Library | Raw pairs (both) | WDL peaks recovered | New peaks in WDL set | Coverage Spearman ρ | New peaks |
|---|---|---|---|---|---|
| 1 | 3.0 M | 67.8 % | 58.7 % | 0.889 | 14,314 |
| 2 | 3.7 M | 74.2 % | 62.1 % | 0.891 | 18,517 |
| 3 | 12.7 M | 80.8 % | 83.7 % | 0.912 | 35,107 |
| 4 | 26.6 M | — (WDL run failed in MACS 1.4 peak model) | — | — | 41,698 |

- Raw read-pair counts agree exactly. The WDL report's own "Raw_Reads" column is the count *after* trimming: its minimum-length rule (≈57 bp) had already discarded ~22 % of pairs.
- WDL peaks the new pipeline does not call are the weaker ones: median MACS score 110–120, against 145–254 for recovered peaks. Agreement rises with sequencing depth.
- The new QC shows what the WDL report could not: 35–59 % of mapped reads in these libraries are mitochondrial. TSS enrichment is 5.4–6.0, FRiP 0.35–0.55, and all four show a clear nucleosome pattern.
