# HiChIP (`--assay hichip`)

From FASTQ to contact pairs, a `.hic` contact map, a coverage track, 1D peaks (loop anchors), FitHiChIP loops and a graded QC report, for paired-end HiChIP libraries. See the [README](../README.md) for installation, references and running at your site.

The design follows the Dovetail/4DN HiChIP route, `bwa-mem2 mem -5SP` then `pairtools`, in place of HiC-Pro. Every library, targets and their INPUT HiChIP controls alike, is aligned and deduplicated **once**. Targets are then paired with their INPUT for peak calling. Loops are called with FitHiChIP using the HiLOW settings, so loop calls stay comparable with earlier PEDDEP runs.

## Samplesheet

```csv
sample,fastq_1,fastq_2,mark,control,peaks,outdir
LINE1_H3K27ac,/data/L1_H3K27ac_L001_R1.fastq.gz,/data/L1_H3K27ac_L001_R2.fastq.gz,H3K27ac,LINE1_INPUT,,
LINE1_H3K27ac,/data/L1_H3K27ac_L002_R1.fastq.gz,/data/L1_H3K27ac_L002_R2.fastq.gz,H3K27ac,LINE1_INPUT,,
LINE1_CTCF,/data/L1_CTCF_R1.fastq.gz,/data/L1_CTCF_R2.fastq.gz,CTCF,LINE1_INPUT,/data/L1_CTCF_ChIP_peaks.bed,
LINE1_INPUT,/data/L1_INPUT_R1.fastq.gz,/data/L1_INPUT_R2.fastq.gz,INPUT,,,
```

- `mark`: the antibody target (`H3K27ac`, `CTCF`, `INPUT`, ...).
- `control`: the `sample` name of the target's INPUT HiChIP library. Leave it empty for the INPUT itself. A target with no control gets peaks called without one.
- `peaks` (optional): a BED of loop anchors for this target, for example ChIP-seq peaks of the same cell line (from `--assay chip`). FitHiChIP then uses these instead of the library's own MACS3 peaks; MACS3 peaks are still called for QC and the browser. A relative path is read from the samplesheet's folder.
- Both FASTQs are required. Lanes (repeated `sample` rows) and the optional `outdir` column work as for the other assays.

## What the workflow does

| Step | Tool |
|---|---|
| Raw read QC | FastQC 0.13 |
| Split each lane into chunks of `--split_fastq` read pairs (50M) that align side by side | seqkit 2.13 |
| Alignment: each read on its own, 5'-most piece of a ligation-junction read as primary | bwa-mem2 `mem -5SP -T0` |
| Contact pairs: MAPQ >= 30, junction reads resolved with `--walks-policy 5unique`, sorted | pairtools 1.1.3 `parse`, `sort` |
| Per library: merge chunks and lanes, remove duplicates, keep contacts between the kept contigs | pairtools `merge`, `dedup`, `select` |
| Contact map | Juicer Tools 1.22.01 `pre` (the build behind the existing PEDDEP `.hic` files) |
| Coverage: every read end of every unique contact, extended to 147 bp, CPM, 10-bp bins | deepTools 3.5.6 `bamCoverage` |
| 1D peaks (loop anchors) on the read ends, against the INPUT: q < 0.01, `--nomodel --extsize 147` | MACS3 3.0.4 |
| FRiP: fraction of read ends in the peaks | samtools |
| Loops: peak-to-all, 5-kb bins, 5 kb – 4 Mb, coverage bias correction, nearby loops merged; tracks at q < 0.01, 0.05 and 0.1 | FitHiChIP 10 (HiLOW settings) |
| PEDDEP graded QC report and run summary | `bin/hichip_qc_report.py`, `bin/qc_summary.py` |
| Everything in one report | MultiQC 1.35 |

**Contigs and blacklist.** Contacts are kept when both ends lie on contigs matching `--canonical_chroms` (minus the mitochondrial contig, unless `--keep_mito`). The `.hic` map and the pairs file keep blacklisted regions, as is usual for contact maps. Peaks, coverage and FRiP skip read ends in the blacklist. FitHiChIP skips contacts with an end in the blacklist and anchors that overlap it, as HiLOW did.

**INPUT libraries.** INPUT HiChIP libraries get the contact map, coverage track and QC report, and serve as the control for peak calling. They get no peaks or loops of their own.

**No anchors.** A target with no peaks outside the blacklist (for example a failed IP) gets empty loop files and 0 loops in its report, instead of stopping the run.

**Containers.** Juicer Tools and FitHiChIP run from the HiLOW image (`ghcr.io/stjude/abralab/hilow:v1.0`), as they did in HiLOW. These two steps need a container engine; they are not on conda.

## Outputs

```
<outdir>/
├── hichip_qc_summary.tsv / .html    # one row per library, PEDDEP grades
├── multiqc/multiqc_report.html
├── pipeline_info/                   # execution report, timeline, trace, DAG, software versions
└── <sample>/                        # or the samplesheet's outdir; targets and INPUTs alike
    ├── <sample>.hichip_qc_report.txt / .html / .hichip_qc_metrics.json
    ├── CONTACTS/                    # <sample>.hic; <sample>.pairs.gz (4DN pairs, --save_pairs)
    ├── COVERAGE_files/              # <sample>.CPM.bw
    ├── PEAKS/                       # targets: _peaks.narrowPeak, _summits.bed, _peaks.xls, _peaks.bed.gz + .tbi
    ├── LOOPS/                       # targets only
    │   ├── <sample>.loops.q{0.01,0.05,0.1}.bedpe        # chr1 s1 e1 chr2 s2 e2 contacts
    │   ├── <sample>.loops.q{0.01,0.05,0.1}.bed.gz + .tbi  # both directions, for ProteinPaint arcs
    │   ├── <sample>.FitHiChIP.interactions.bed.gz       # every tested bin pair with p and q
    │   ├── <sample>.FitHiChIP.Q0.01.bed / .Q0.01.merged.bed
    │   ├── <sample>.anchors.bed     # the anchors FitHiChIP used
    │   └── <sample>.loops.txt       # anchor and loop counts
    └── QC/
        ├── FastQC/raw/   pairtools/ (<sample>.dedup.stats)   peaks/ (FRiP)
```

For ProteinPaint, use `CONTACTS/<sample>.hic`, `COVERAGE_files/<sample>.CPM.bw`, `PEAKS/<sample>_peaks.bed.gz` and `LOOPS/<sample>.loops.q*.bed.gz` (each with its `.tbi`). They replace HiLOW's `*_fastq.allValidPairs.hic`, `*_fastqbedgraph.bw`, `*_fastqpeaks.bed.gz` and `*_Qvalue.<q>.3.0.4.5000.4000000.bed.gz`. The loop tracks have the same columns as HiLOW's.

## PEDDEP graded QC report

The grading works as for the other assays: each graded metric scores -2 to +2, and the overall grade is the mean.

| Metric | Source | Cutoffs for -1 / 0 / +1 / +2 |
|---|---|---|
| Raw_Read_Pairs | FastQC (raw) | 50M / 100M / 150M / 200M |
| Base_Quality, Sequence_Diversity | FastQC, worst FASTQ | fail -2, warn 0, pass +2 |
| Pct_Reads_Q30 | FastQC, worst FASTQ | 60 / 70 / 75 / 80 % |
| Mapped_Pairs_Pct: both ends mapped uniquely | pairtools | 50 / 60 / 70 / 80 % |
| Duplication_Pct: of mapped pairs (lower is better) | pairtools | <= 70 / 60 / 50 / 40 % |
| Unique_Contacts | pairtools | 25M / 50M / 75M / 100M |
| Cis_Pct: of unique contacts | pairtools | 40 / 50 / 60 / 70 % |
| FRiP: read ends in peaks | MACS3 peaks | 0.01 / 0.02 / 0.05 / 0.1 |
| Total_Peaks | MACS3 | 1k / 2k / 5k / 10k |
| Loops_Q0.01 | FitHiChIP | 500 / 1k / 5k / 10k |

`Mean_Read_Length`, `Cis_Under_1kb_Pct`, `Cis_Over_20kb_Pct`, `Trans_Pct` and `Merged_Loops_Q0.01` are reported but not graded. INPUT libraries are graded on sequencing and contact metrics only. Every contact fraction is a share of the deduplicated contacts. The HiLOW report divided deduplicated cis counts by all valid pairs, duplicates included, which made more-duplicated libraries look worse. The cutoffs live in [`assets/hichip_qc_thresholds.json`](../assets/hichip_qc_thresholds.json); pass your own with `--qc_thresholds`.

HICHIP_COMPARISON_PLACEHOLDER
