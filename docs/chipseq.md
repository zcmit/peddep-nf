# ChIP-seq and CUT&RUN (`--assay chip`, `--assay cutrun`)

From FASTQ to peaks, coverage tracks, ROSE stitched regions and a graded QC report, for paired-end or single-end libraries. See the [README](../README.md) for installation, references and running at your site.

The layout follows [pyflow-ChIPseq](https://github.com/crazyhottommy/pyflow-ChIPseq): every library, cases and their Input/IgG controls alike, is trimmed, aligned and filtered **once**. Cases are then paired with their control for peak calling and case-vs-control tracks. A control shared by several cases is processed one time, not once per case as in the SEAseq WDL.

`--assay chip` and `--assay cutrun` run the same steps. They differ in duplicate handling and in the read-depth cutoffs of the QC grade (see below).

## Samplesheet

```csv
sample,fastq_1,fastq_2,mark,control,outdir
LINE1_H3K27ac,/data/L1_H3K27ac_R1.fastq.gz,/data/L1_H3K27ac_R2.fastq.gz,H3K27ac,LINE1_Input,
LINE1_H3K27me3,/data/L1_H3K27me3_R1.fastq.gz,/data/L1_H3K27me3_R2.fastq.gz,H3K27me3,LINE1_Input,
LINE1_Input,/data/L1_Input_L001_R1.fastq.gz,/data/L1_Input_L001_R2.fastq.gz,Input,,
LINE1_Input,/data/L1_Input_L002_R1.fastq.gz,/data/L1_Input_L002_R2.fastq.gz,Input,,
LINE2_CTCF_SE,/data/L2_CTCF.fastq.gz,,CTCF,,/projects/LINE2/CTCF
```

- `mark`: the antibody target (`H3K27ac`, `CTCF`, `Input`, `IgG`, ...). It is matched without regard to case against `--broad_marks` and `--rose_marks`.
- `control`: the `sample` name of the case's Input or IgG. Leave it empty for the controls themselves. A case with no control gets peaks called without one.
- Leave `fastq_2` empty for single-end reads. All lanes of one sample must be the same type.
- Lanes (repeated `sample` rows) and the optional `outdir` column work as for ATAC-seq.

## What the workflow does

| Step | Tool |
|---|---|
| Raw read QC | FastQC 0.13 |
| Adapter trimming (adapter auto-detected unless `--adapter`), post-trim QC | Trim Galore 0.6.10 / cutadapt 4.9 |
| Alignment, per lane | bowtie2 2.5.4 `--very-sensitive` (paired-end: `-X 2000 --no-mixed --no-discordant`) |
| Lane merge, duplicate marking | samtools, Picard 3.5 MarkDuplicates |
| Library complexity (NRF, PBC1, PBC2) | ENCODE definitions |
| Filtering: MAPQ >= 30, proper pairs (paired-end); drop mito, blacklist and non-canonical contigs; drop duplicates (see below) | samtools |
| Strand cross-correlation: NSC, RSC, fragment length | phantompeakqualtools 1.2.2, on at most 15M reads (read 1 for paired-end, taken at the median insert size) |
| Peaks vs control: narrow (q < 0.05, summits) and broad (`--broad-cutoff 0.1`) | MACS3 3.0.4: `BAMPE` for paired-end; `--nomodel --extsize 200` for single-end |
| FRiP and peak count, on the mark's primary peak set: broad for `--broad_marks`, otherwise narrow | samtools |
| Coverage: fragments, CPM, 10 bp (every library); case minus control, CPM, 50 bp | deepTools 3.5.6 `bamCoverage`, `bamCompare` |
| Fingerprint (case + control), TSS profile and heatmap | deepTools |
| Stitched and super-enhancer-like regions for `--rose_marks` (default H3K27ac), from narrow peaks with p <= 1e-9 | ROSE 1.3.2 |
| PEDDEP graded QC report and run summary | `bin/chip_qc_report.py`, `bin/qc_summary.py` |
| Everything in one report | MultiQC 1.35 |

**Duplicates.** `--assay chip` removes duplicates from every library. `--assay cutrun` removes them from IgG controls only and keeps them in targets: in CUT&RUN, MNase cuts at the same positions on purpose, so identical fragments are often independent molecules. Add `--cutrun_dedup_targets` to remove them from targets as well. Duplicates are always counted (Duplication_Pct, NRF, PBC).

**Broad marks.** `--broad_marks` (default `H3K27me3,H3K36me3,H3K9me3,H3K9me2,H3K79me2,H3K4me1,H4K20me1`) lists the marks whose FRiP and peak count come from broad peaks. Both peak sets are produced for every case.

**ROSE.** ROSE runs only for the marks in `--rose_marks` (`''` turns it off). Its input is the narrow peaks with p <= `--rose_peak_pvalue` (1e-9, ROSE's usual MACS input and what SEAseq used), merged into `unionpeaks.gff`. Regions are stitched within `--rose_stitch` (12.5 kb), excluding `--rose_tss` (±2 kb) around TSSs, and ranked by case signal minus control. Output file names match the old `STITCHED_peaks` folder. Two slow steps of ROSE are replaced:

- **Read mapping.** `bin/rose_bamtogff.py` replaces `ROSE_bamToGFF.py`. It computes the same per-region densities from one pass over the BAM. Values match ROSE to the 4th decimal in all but a handful of regions.
- **Gene mapping.** `bin/rose_genemapper.py` runs ROSE's own `ROSE_geneMapper.py` with two lookups sped up: the transcript list and the 50 Mb TSS search. On the human data below its tables are byte-identical to ROSE's. It also fixes a ROSE crash on regions with no gene within 50 Mb, such as contigs without annotated genes; ROSE stopped there and wrote no gene tables.

On the 44M-pair H3K27ac library below, ROSE takes 13 minutes. The same input takes 1 h 49 min with ROSE's own gene mapping, and hours more with its own read mapping. ROSE needs a container engine; it is not on conda.

## Outputs

```
<outdir>/
├── chip_qc_summary.tsv / .html      # or cutrun_qc_summary.*; one row per library, PEDDEP grades
├── multiqc/multiqc_report.html
├── pipeline_info/                   # execution report, timeline, trace, DAG, software versions
└── <sample>/                        # or the samplesheet's outdir; cases and controls alike
    ├── <sample>.chip_qc_report.txt / .html / .chip_qc_metrics.json
    ├── BAM_files/                   # <sample>.filtered.bam (+ .bai)
    ├── COVERAGE_files/              # <sample>.fragments.CPM.bw; cases: <sample>.minus_control.CPM.bw
    ├── PEAKS/                       # cases only
    │   ├── NARROW_peaks/            # _peaks.narrowPeak, _summits.bed, _peaks.xls, _narrow_peaks.bed.gz + .tbi
    │   ├── BROAD_peaks/             # _peaks.broadPeak, _peaks.gappedPeak, _peaks.xls, _broad_peaks.bed.gz + .tbi
    │   └── STITCHED_peaks/          # --rose_marks only: unionpeaks_{All,Super}Stitched.table.txt, _withGENES, plots
    └── QC/
        ├── FastQC/{raw,trimmed}/    trim_report/   alignment/   duplicates/   library_complexity/
        └── phantompeakqualtools/    peaks/ (FRiP)  fingerprint/ TSS/
```

For a genome browser, use `COVERAGE_files/<sample>.fragments.CPM.bw` (or the control-subtracted `.minus_control.CPM.bw`) and `PEAKS/*/<sample>_{narrow,broad}_peaks.bed.gz` with its `.tbi`. They replace the WDL's `*.sorted.bklist.rmdup.w50.FPM.bw` and `*-p9_kd-auto_peaks.bed`.

## PEDDEP graded QC report

The grading works as for ATAC-seq: each graded metric scores -2 to +2, and the overall grade is the mean. The cutoffs come from the PEDDEP/SEAseq `evaluation-statistics.pl` report.

| Metric | Source | Cutoffs for -1 / 0 / +1 / +2 |
|---|---|---|
| Raw_Reads (read pairs for paired-end) | FastQC (raw) | ChIP 15M / 20M / 25M / 30M; CUT&RUN 3M / 5M / 8M / 10M |
| Base_Quality, Sequence_Diversity | FastQC, worst FASTQ | fail -2, warn 0, pass +2 |
| Pct_Reads_Q30 | FastQC, worst FASTQ | 60 / 70 / 75 / 80 % |
| Aligned_Percent | bowtie2 | 50 / 60 / 70 / 80 % |
| NRF | library complexity | 0.5 / 0.6 / 0.7 / 0.8 |
| PBC1 | library complexity | 0.5 / 0.66 / 0.75 / 0.9 |
| PBC2 | library complexity | < 1 → -2, >= 1 → 0, >= 3 → +2 |
| NSC | phantompeakqualtools | < 1.045 → -2, otherwise +2 |
| RSC | phantompeakqualtools | < 0.75 → -2, >= 0.75 → +1, >= 1 → +2 |
| FRiP | primary peak set | 0.0075 / 0.01 / 0.02 / 0.05 |
| Total_Peaks | primary peak set | 1k / 2k / 5k / 10k |

`Mean_Read_Length`, `Adapter_Trimmed_Pct`, `Duplication_Pct`, `Final_Reads`, `Fragment_Length`, `Peak_Type`, `Stitched_Regions` and `SE_Like_Regions` are reported but not graded. Input/IgG controls are graded on sequencing and library metrics only, since a good control shows no enrichment. Their NSC and RSC are still shown. The cutoffs live in [`assets/chip_qc_thresholds.json`](../assets/chip_qc_thresholds.json) and [`assets/cutrun_qc_thresholds.json`](../assets/cutrun_qc_thresholds.json); pass your own with `--qc_thresholds`.

## Changes from the PEDDEP SEAseq WDL

| | SEAseq (`sj-seaseq.py`, PEASEQ) | peddep-nf |
|---|---|---|
| Runs on | St. Jude LSF + Cromwell only | anywhere Nextflow runs |
| Input | one call per case, with its control | one samplesheet for a whole batch, any number of lanes |
| Controls | re-processed for every case that uses them | processed once and shared |
| Paired-end data | analysed three times: paired-end, then R1 and R2 again as "single-end mode" | analysed once; single-end only for single-end libraries |
| Trimming | none | Trim Galore |
| Aligner | bowtie1 `-k 2 -m 2 -X 600` | bowtie2 `-X 2000`, then MAPQ >= 30 |
| Duplicates | removed | removed for ChIP and for controls; kept in CUT&RUN targets |
| Peak calling | MACS 1.4.2, three runs per case (p 1e-9 with and without duplicates, fixed shift) + SICER | MACS3 narrow + broad, one run each, fragment-aware for paired-end |
| ROSE | every mark, hours per run | `--rose_marks` only, about 15 min per run |
| BAM2GFF metagene plots, peaksanno, TDF / wig tracks | yes | dropped (deepTools TSS plots instead of BAM2GFF) |
| Motifs (MEME / HOMER) | optional, off in PEDDEP runs | dropped |
| QC report "Raw_Reads" | BAM records after bowtie1 `-k 2`, so multi-mapped reads count twice | read pairs (or reads) in the FASTQ |
| Resume after failure | no | `-resume` |

## Comparison with the SEAseq WDL

Two ChIP-seq libraries (H3K27ac, H3K27me3; 44–57M pairs, sharing one Input) and ten CUT&RUN libraries (CTCF and H3K4me3 in five cell lines, 24–38M pairs, each line with its own IgG) were run through both pipelines on GRCh38. All are NextSeq 2000, 2×76 bp. Old narrow peaks are the WDL's paired-end MACS 1.4 `p9_kd-auto` set. Old broad peaks are SICER islands (FDR 0.01). Coverage correlation compares the WDL's 50-bp fragment FPM bigWig with the new fragment CPM bigWig, in 5-kb bins with chrM excluded.

| Library | Old narrow peaks recovered | New narrow peaks | SICER bp in new broad peaks | New broad bp in SICER | Coverage Spearman ρ |
|---|---|---|---|---|---|
| ChIP H3K27ac | 98.3 % of 28,544 | 62,682 | 64.4 % | 97.5 % | 0.937 |
| ChIP H3K27me3 | 98.1 % of 59,052 | 261,453 | 72.1 % | 94.9 % | 0.958 |
| CUT&RUN CTCF (5 lines) | 98.1–98.8 % | 59–87k | 39–62 % | 82–94 % | 0.80–0.88 |
| CUT&RUN H3K4me3 (5 lines) | 96.2–97.7 % | 17–23k | 55–67 % | 92–99 % | 0.78–0.86 |

- Raw read-pair counts agree exactly with the WDL's raw FastQC for all 12 libraries.
- Nearly every WDL peak is found again. The new pipeline calls more peaks because q < 0.05 is a looser cutoff than p < 1e-9. ROSE gets the p <= 1e-9 subset, which matches the old peak count (26,247 vs 28,544 for H3K27ac).
- New broad peaks lie almost entirely inside SICER islands. SICER's islands are wider.
- CUT&RUN coverage correlates less than ChIP for two reasons. Most 5-kb bins in CUT&RUN hold only background. The new tracks also keep duplicate fragments in CUT&RUN targets, where the WDL removed them.
- ROSE (H3K27ac) finds 558 super-enhancer-like regions against the WDL's 557. 94.8 % of the old ones overlap a new one (96.4 % the other way), and 99.6 % fall inside a new stitched region. There are 10,359 stitched regions against 12,577: 96.7 % of the new ones overlap an old one and 78.2 % of the old ones overlap a new one, because MACS3 and MACS 1.4 give somewhat different p <= 1e-9 peak sets.

### Compute cost

| | SEAseq WDL | peddep-nf |
|---|---|---|
| One CUT&RUN case + its IgG (SJ0242 CTCF) | 212 LSF jobs, 143 CPU-hours, 23 h from first to last job | 30 tasks, 23 CPU-hours |
| The batch of 10 cases + 5 IgG | 10 runs (each IgG processed twice): ≈ 1,430 CPU-hours (10 × the run above) | 251 tasks, 142 CPU-hours, 1 h 55 min wall time |

The largest costs in the WDL run were BAM2GFF metagene plots (95 CPU-hours), the second "single-end mode" pass, and ROSE (9.5 CPU-hours) for a CTCF run. In peddep-nf, bowtie2 alignment takes three quarters of the CPU time. The ChIP batch above (H3K27ac, H3K27me3 and their Input) took 2 h 51 min wall time and 49 CPU-hours, with 13 minutes for ROSE.
