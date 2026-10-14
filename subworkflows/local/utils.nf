/*
========================================================================================
    Shared helper functions (no processes)
========================================================================================
*/

// Value of a --genome entry in conf/genomes.config, or null
def getGenomeAttribute(String attribute) {
    if (params.genome && params.genomes && params.genomes.containsKey(params.genome)) {
        return params.genomes[params.genome][attribute]
    }
    return null
}

// Sample-level meta map: drops the lane fields so lanes of one sample join together
def toSampleMeta(Map meta) {
    return [id: meta.id, outdir: meta.outdir]
}

def validateParameters() {
    def supported = ['atac']
    if (!(params.assay in supported)) {
        error("--assay '${params.assay}' is not available yet. Supported: ${supported.join(', ')}")
    }
    if (!params.input) {
        error("Missing --input <samplesheet.csv>")
    }
    if (!params.outdir) {
        error("Missing --outdir <directory>")
    }
    if (params.genome && !params.genomes.containsKey(params.genome)) {
        error("--genome '${params.genome}' is not defined. Available: ${params.genomes.keySet().join(', ')}. Or give --fasta/--gtf directly.")
    }
    if (!params.fasta && !getGenomeAttribute('fasta')) {
        error("No reference: give --genome <${params.genomes.keySet().join('|')}> or --fasta (and --gtf)")
    }
    if (!params.gtf && !getGenomeAttribute('gtf')) {
        error("No annotation: --gtf is needed for TSS enrichment")
    }
}

/*
 * Samplesheet: sample,fastq_1,fastq_2[,outdir]
 * Rows sharing a sample name are lanes of that sample and are merged after alignment.
 * Returns [ meta, [fastq_1, fastq_2] ] per lane, meta = [id, outdir, lane, n_lanes].
 */
def parseSamplesheet(String path) {
    def rows = file(path, checkIfExists: true).splitCsv(header: true, strip: true)
    if (!rows) {
        error("Samplesheet '${path}' has no rows")
    }
    def missing = ['sample', 'fastq_1', 'fastq_2'].findAll { col -> !rows[0].containsKey(col) }
    if (missing) {
        error("Samplesheet '${path}' is missing column(s): ${missing.join(', ')}. Expected: sample,fastq_1,fastq_2[,outdir]")
    }

    def lane_counts = [:]
    def outdirs     = [:]
    def seen        = []
    rows.eachWithIndex { row, idx ->
        def line = idx + 2
        def sample = row.sample
        if (!sample || !(sample ==~ /^[A-Za-z0-9._+-]+$/)) {
            error("Samplesheet line ${line}: sample name '${sample}' must be non-empty and use only letters, digits, '.', '_', '+', '-'")
        }
        if (!row.fastq_1 || !row.fastq_2) {
            error("Samplesheet line ${line}: both fastq_1 and fastq_2 are required (paired-end ATAC-seq)")
        }
        if (row.fastq_1 == row.fastq_2) {
            error("Samplesheet line ${line}: fastq_1 and fastq_2 are the same file")
        }
        def key = "${sample}\t${row.fastq_1}"
        if (key in seen) {
            error("Samplesheet line ${line}: duplicate row for sample '${sample}' (${row.fastq_1})")
        }
        seen << key
        def outdir = row.outdir ?: null
        if (outdirs.containsKey(sample) && outdirs[sample] != outdir) {
            error("Samplesheet line ${line}: sample '${sample}' has different outdir values on different rows")
        }
        outdirs[sample] = outdir
        lane_counts[sample] = (lane_counts[sample] ?: 0) + 1
    }

    def lane_index = [:]
    return rows.collect { row ->
        def sample = row.sample
        lane_index[sample] = (lane_index[sample] ?: 0) + 1
        def meta = [
            id     : sample,
            outdir : outdirs[sample],
            lane   : lane_index[sample],
            n_lanes: lane_counts[sample]
        ]
        def reads = [
            file(row.fastq_1, checkIfExists: true),
            file(row.fastq_2, checkIfExists: true)
        ]
        [meta, reads]
    }
}

def helpMessage() {
    return """
    zcmit/peddep-nf  v${workflow.manifest.version}

    Usage:
      nextflow run zcmit/peddep-nf --assay atac --input samplesheet.csv --outdir results \\
          --genome GRCh38 -profile <docker|singularity|apptainer|podman|conda>[,<institution>]

    Required:
      --input             CSV with columns sample,fastq_1,fastq_2[,outdir]
                          (repeat a sample on several rows for multiple lanes)
      --outdir            Run-level output directory (MultiQC, QC summary, pipeline_info)
      --genome | --fasta  A key from conf/genomes.config (GRCh38, hg38, GRCm38, mm10)
                          or an explicit FASTA (local or remote, may be .gz)

    References (override or replace --genome):
      --gtf               Gene annotation, used for TSS enrichment
      --blacklist         BED of regions to exclude
      --bowtie2_index     Directory or prefix of a prebuilt index (built if absent)
      --mito_name         Mitochondrial contig name [chrM]
      --macs_gsize        Effective genome size for MACS3 (computed from FASTA if absent)
      --canonical_chroms  Regex of contigs to keep, e.g. '^chr([0-9]+|X|Y|M)\$'
      --save_reference    Publish the prepared reference files (incl. a built index)

    Processing:
      --adapter           Adapter for Trim Galore [${params.adapter}]
      --min_trim_length   Discard pairs with a read shorter than this after trimming [${params.min_trim_length}]
      --min_mapq          Minimum MAPQ [${params.min_mapq}]
      --keep_mito         Keep mitochondrial reads [${params.keep_mito}]
      --keep_dups         Keep duplicate reads [${params.keep_dups}]
      --macs_qvalue       MACS3 narrow-peak q-value [${params.macs_qvalue}]
      --skip_broad_peaks  Skip MACS3 broad peaks [${params.skip_broad_peaks}]
      --qc_thresholds     JSON with the grading thresholds [assets/atac_qc_thresholds.json]

    Resources:
      --max_cpus / --max_memory / --max_time   Per-task ceilings [${params.max_cpus} / ${params.max_memory} / ${params.max_time}]

    See README.md for outputs, site configuration and cloud/HPC examples.
    """.stripIndent()
}
