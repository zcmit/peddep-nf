# Changelog

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
