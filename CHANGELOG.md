# Changelog

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
