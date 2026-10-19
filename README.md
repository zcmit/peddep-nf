# peddep-nf

Portable [Nextflow](https://www.nextflow.io) pipelines for PEDDEP epigenomic profiling.

| Assay | Option | Details |
|---|---|---|
| **ATAC-seq** (paired-end) | `--assay atac` | [docs/atac.md](docs/atac.md) |
| **ChIP-seq** (paired- or single-end) | `--assay chip` | [docs/chipseq.md](docs/chipseq.md) |
| **CUT&RUN** (paired- or single-end) | `--assay cutrun` | [docs/chipseq.md](docs/chipseq.md) |
| HiChIP / Micro-C | planned | |

The pipeline runs anywhere Nextflow runs: a laptop, any HPC scheduler (LSF, SLURM, SGE, PBS) or the cloud (AWS Batch, Google Batch, Azure Batch), with Docker, Podman, Singularity/Apptainer or conda. Nothing in this repository is tied to one site; local paths and queues go in a config file you keep yourself (see [Running at your site](#running-at-your-site)).

## Quick start

```bash
# Nextflow >= 24.10 and one container engine (or conda)
nextflow run zcmit/peddep-nf -profile test,docker --outdir results            # ATAC-seq
nextflow run zcmit/peddep-nf -profile test_chip,docker --outdir results_chip  # ChIP-seq
```

Both test profiles use small public nf-core yeast datasets (100k read pairs per FASTQ, with one sample split over two lanes) and finish in a few minutes on 2 CPUs. `test_chip` also covers a shared Input control, a single-end library and ROSE; it exercises the workflow, not the biology. Swap `docker` for `singularity`, `apptainer` or `podman` (or `conda`, except for ROSE, which needs a container).

## Running on your data

```bash
nextflow run zcmit/peddep-nf -r v0.2.0 \
    --assay chip \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    -profile singularity
```

### Samplesheet

| Assay | Columns |
|---|---|
| ATAC-seq | `sample,fastq_1,fastq_2[,outdir]` ([example](docs/atac.md#samplesheet)) |
| ChIP-seq, CUT&RUN | `sample,fastq_1,fastq_2,mark,control[,outdir]`: `control` names the case's Input/IgG row ([example](docs/chipseq.md#samplesheet)) |

- FASTQs may be gzipped or not, as local paths or `https://`, `s3://`, `gs://`, `az://` URLs.
- Rows that share a `sample` are lanes of that sample. Each lane is trimmed and aligned separately, then the lanes are merged.
- `outdir` (optional column) sends that sample's results to its own folder, for example a project tree organised by cell line. Without it, results go to `<--outdir>/<sample>/`. Run-level files (MultiQC, QC summary, `pipeline_info`) always go to `--outdir`.

### References

Pick one:

- `--genome GRCh38` / `hg38`, or `GRCm38` / `mm10`. These pull GENCODE primary-assembly FASTA and GTF files plus the Boyle-lab blacklist v2. The bowtie2 index is built on first use; add `--save_reference` to keep it under `<outdir>/reference/` and pass it as `--bowtie2_index` next time.
- Your own files: `--fasta`, `--gtf`, and optionally `--blacklist`, `--bowtie2_index` (directory or prefix), `--mito_name`, `--macs_gsize`, `--canonical_chroms`. They can be local or remote, and gzipped or not.

Reads are kept only on contigs that match `--canonical_chroms` (for the human and mouse genome keys this is `chr1-22, X, Y`). Reads on the mitochondrial contig and in blacklisted regions are also dropped.

## Running at your site

The **container engine** is chosen with `-profile docker | podman | singularity | apptainer | conda | mamba | wave`. Singularity/Apptainer images are downloaded as ready-made `.sif` files, so no image is ever built on a compute node. Set `NXF_SINGULARITY_CACHEDIR` (or `singularity.cacheDir`) to a shared folder so each image is downloaded only once.

The **scheduler** can come from three places:

1. **An nf-core institutional profile**, if your institution has one ([list](https://nf-co.re/configs)), e.g. `-profile stjude,singularity`. These profiles are loaded automatically from nf-core/configs. Set `NXF_OFFLINE=true` to skip that download.
2. **An example config**, edited for your site: [`conf/examples/`](conf/examples) has `slurm`, `sge`, `pbs`, `awsbatch`, `google_batch` and `azurebatch` configs. Pass one with `-c`.
3. **Your own config file.** Start from [`conf/examples/site_template.config`](conf/examples/site_template.config), which shows local reference paths, the queue, resource ceilings and the image cache. Keep it outside this repository.

```bash
# e.g. a SLURM cluster with a local copy of the references
nextflow run zcmit/peddep-nf -profile singularity -c my_site.config --input samples.csv --outdir results

# e.g. AWS Batch, everything on S3
nextflow run zcmit/peddep-nf -c conf/examples/awsbatch.config --genome GRCh38 \
    --input s3://bucket/samples.csv --outdir s3://bucket/results -work-dir s3://bucket/work
```

Per-task ceilings: `--max_cpus`, `--max_memory`, `--max_time` on the command line, or `process.resourceLimits` in a config. Re-running with `-resume` reuses every finished task.

## Testing

| Environment | What was run | Result |
|---|---|---|
| `nextflow lint`, `-stub` | full DAG for all three assays, incl. multi-lane merge, shared controls, single-end ChIP and per-sample `outdir`; Nextflow 24.10.2 and 25.10.2 | pass |
| Local executor + Singularity (one 8-core node) | `-profile test,singularity` and `-profile test_chip,singularity`; re-run with `-resume` | pass; all tasks cached on resume |
| LSF via the nf-core `stjude` institutional profile | `-profile test,stjude` | pass; results identical to the local run |
| LSF + Singularity, human data, GRCh38 | 4 PEDDEP ATAC libraries; 3 ChIP-seq libraries (2 marks + shared Input); 15 CUT&RUN libraries (10 cases + 5 IgG) | compared with the WDL pipelines: [ATAC](docs/atac.md#comparison-with-the-wdl-pipeline), [ChIP / CUT&RUN](docs/chipseq.md#comparison-with-the-seaseq-wdl) |
| GitHub Actions, Docker (`ubuntu-latest`) | lint; stub and Docker runs of `test` and `test_chip` (+ CUT&RUN stub) on Nextflow 24.10.2 and latest | pass |
| SLURM, SGE, PBS, AWS/Google/Azure Batch | example configs only | **not yet tested** |

## License

MIT
