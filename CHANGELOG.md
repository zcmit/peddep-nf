# Changelog

## v0.3.0 (unreleased)

HiChIP (`--assay hichip`), replacing the PEDDEP HiLOW WDL pipeline (see docs/hichip.md).

- Samplesheet with `mark` and `control` columns as for ChIP-seq; INPUT HiChIP libraries are
  processed once and serve as the peak-calling control; optional `peaks` column for loop anchors.
- bwa-mem2 `-5SP` and pairtools (parse, dedup) instead of HiC-Pro; lanes are split into chunks of
  `--split_fastq` read pairs that align in parallel. bwa-mem2 index built on first use.
- Juicer Tools 1.22.01 `.hic` maps and FitHiChIP 10 loops with the HiLOW settings, from the HiLOW
  image; loop tracks at q < 0.01 / 0.05 / 0.1 in the HiLOW layout.
- MACS3 peaks on the read ends against the INPUT replace HiChIP-Peaks; CPM coverage bigWig.
- PEDDEP graded HiChIP QC report: contact fractions use deduplicated contacts throughout, the
  duplicate rate is reported, and long-range cis contacts (>= 20 kb) are graded as in HiLOW.
- `test_hichip` profile (nf-core HiCAR reads on human chr22) and CI.
- ATAC-seq, ChIP-seq and CUT&RUN results are unchanged from v0.2.0.

## v0.2.0 (2026-10-20)

ChIP-seq and CUT&RUN (`--assay chip`, `--assay cutrun`), a rewrite of the PEDDEP SEAseq/PEASEQ
WDL pipeline in the style of pyflow-ChIPseq (see docs/chipseq.md).

- Samplesheet with `mark` and `control` columns; each library, including shared Input/IgG
  controls, is processed once; paired-end and single-end libraries.
- Trim Galore, bowtie2, Picard; duplicates removed for ChIP and controls, kept in CUT&RUN targets.
- MACS3 narrow and broad peaks against the control; FRiP on broad peaks for `--broad_marks`.
- Fragment CPM and control-subtracted bigWigs; bgzipped + tabix-indexed peaks for browsers.
- ROSE 1.3.2 for `--rose_marks` (default H3K27ac) on p <= 1e-9 peaks. Its read-mapping and
  gene-mapping steps are replaced by faster code with the same output (`bin/rose_bamtogff.py`,
  `bin/rose_genemapper.py`): 13 min instead of hours on a 44M-pair H3K27ac library. Regions
  with no gene within 50 Mb no longer stop ROSE's gene mapping.
- QC: phantompeakqualtools NSC/RSC (at the median insert size for paired-end), ENCODE library
  complexity, deepTools fingerprint and TSS plots, MultiQC, PEDDEP graded report with ChIP and
  CUT&RUN cutoffs.
- `test_chip` profile and CI; ATAC docs moved to docs/atac.md.
- ATAC results are unchanged from v0.1.0.

## v0.1.0 (2026-10-14)

First release: the ATAC-seq workflow (`--assay atac`), a portable rewrite of the PEDDEP
WDL/Cromwell ATAC pipeline (see the README's "Changes from the PEDDEP WDL ATAC pipeline").

- Samplesheet input for many samples; lanes merged after alignment; optional per-sample output folder.
- Trim Galore, bowtie2, Picard, MAPQ/blacklist/mito filtering, Tn5 shift, MACS3 narrow + broad peaks.
- ATAC QC: TSS enrichment (ataqv), fragment-length / nucleosome pattern, ENCODE library complexity,
  FRiP, deepTools TSS plots and fingerprint, MultiQC.
- PEDDEP graded QC report per sample and a run-level summary table.
- Browser-ready tracks: fragment CPM bigWig, Tn5 insertion bigWig, bgzipped + tabix-indexed peaks.
- Docker, Podman, Singularity/Apptainer, conda; nf-core institutional profiles; example SLURM/SGE/PBS
  and AWS/Google/Azure Batch configs; public GRCh38/GRCm38 references; yeast test profile and CI.
