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
    return meta.findAll { key, _value -> !(key in ['lane', 'n_lanes']) }
}

// Collect per-lane files of each sample into one list: [ sample_meta, [ files ] ]
def groupLaneFiles(ch) {
    return ch
        .map { meta, files -> [ groupKey(toSampleMeta(meta), meta.n_lanes), files ] }
        .groupTuple()
        .map { key, files -> [ key.getGroupTarget(), files.flatten() ] }
}

def validateParameters() {
    def supported = ['atac', 'chip', 'cutrun']
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
      nextflow run zcmit/peddep-nf --assay <atac|chip|cutrun> --input samplesheet.csv --outdir results \\
          --genome GRCh38 -profile <docker|singularity|apptainer|podman|conda>[,<institution>]

    Required:
      --assay             atac (default), chip (ChIP-seq) or cutrun (CUT&RUN)
      --input             atac:         sample,fastq_1,fastq_2[,outdir]
                          chip/cutrun:  sample,fastq_1,fastq_2,mark,control[,outdir]
                          (repeat a sample on several rows for multiple lanes; fastq_2 empty =
                          single-end; control = sample name of the Input/IgG row)
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
      --adapter           Adapter for Trim Galore [ATAC: Nextera CTGTCTCTTATACACATCT; ChIP/CUT&RUN: auto-detect]
      --min_trim_length   Discard pairs with a read shorter than this after trimming [${params.min_trim_length}]
      --min_mapq          Minimum MAPQ [${params.min_mapq}]
      --keep_mito         Keep mitochondrial reads [${params.keep_mito}]
      --keep_dups         Keep duplicate reads [${params.keep_dups}]
      --macs_qvalue       MACS3 q-value [ATAC 0.01, ChIP/CUT&RUN 0.05]
      --skip_broad_peaks  ATAC only: skip MACS3 broad peaks [${params.skip_broad_peaks}]
      --qc_thresholds     JSON with the grading thresholds [assets/<assay>_qc_thresholds.json]

    ChIP-seq / CUT&RUN:
      --broad_marks       Marks graded on broad peaks [${params.broad_marks}]
      --rose_marks        Marks that get ROSE stitched regions ('' = none) [${params.rose_marks}]
      --se_fragment_length  Single-end fragment length for MACS3/bigWigs [${params.se_fragment_length}]
      --cutrun_dedup_targets  Also remove duplicates from CUT&RUN targets [${params.cutrun_dedup_targets}]

    Resources:
      --max_cpus / --max_memory / --max_time   Per-task ceilings [${params.max_cpus} / ${params.max_memory} / ${params.max_time}]

    See README.md for outputs, site configuration and cloud/HPC examples.
    """.stripIndent()
}

/*
 * ChIP-seq / CUT&RUN samplesheet: sample,fastq_1,fastq_2,mark,control[,outdir]
 *   - fastq_2 empty = single-end; rows sharing a sample name are lanes
 *   - control = the sample name of that library's Input/IgG row (empty for controls
 *     and for cases without a control); one control can serve many cases
 * Returns [ meta, [fastqs] ] per lane.
 */
def parseChipSamplesheet(String path, String assay) {
    def rows = file(path, checkIfExists: true).splitCsv(header: true, strip: true)
    if (!rows) {
        error("Samplesheet '${path}' has no rows")
    }
    def missing = ['sample', 'fastq_1', 'fastq_2', 'mark', 'control'].findAll { col -> !rows[0].containsKey(col) }
    if (missing) {
        error("Samplesheet '${path}' is missing column(s): ${missing.join(', ')}. Expected: sample,fastq_1,fastq_2,mark,control[,outdir]")
    }

    def broad_marks = (params.broad_marks ?: '').tokenize(',').collect { m -> m.trim().toLowerCase() }
    def info  = [:]
    def seen  = []
    rows.eachWithIndex { row, idx ->
        def line   = idx + 2
        def sample = row.sample
        if (!sample || !(sample ==~ /^[A-Za-z0-9._+-]+$/)) {
            error("Samplesheet line ${line}: sample name '${sample}' must be non-empty and use only letters, digits, '.', '_', '+', '-'")
        }
        if (!row.fastq_1) {
            error("Samplesheet line ${line}: fastq_1 is required")
        }
        if (!row.mark) {
            error("Samplesheet line ${line}: mark is required (e.g. H3K27ac, CTCF, Input, IgG)")
        }
        def key = "${sample}\t${row.fastq_1}"
        if (key in seen) {
            error("Samplesheet line ${line}: duplicate row for sample '${sample}' (${row.fastq_1})")
        }
        seen << key
        def current = [
            outdir    : row.outdir ?: null,
            mark      : row.mark,
            control   : row.control ?: '',
            single_end: !row.fastq_2
        ]
        if (info.containsKey(sample)) {
            def first = info[sample]
            ['outdir', 'mark', 'control', 'single_end'].each { field ->
                if (first[field] != current[field]) {
                    error("Samplesheet line ${line}: sample '${sample}' has a different ${field} than on its earlier row(s)")
                }
            }
            first.n_lanes += 1
        } else {
            current.n_lanes = 1
            info[sample] = current
        }
    }

    def controls = info.values().collect { v -> v.control }.findAll { c -> c } as Set
    controls.each { c ->
        if (!info.containsKey(c)) {
            error("Samplesheet: control '${c}' is not a sample in ${path}")
        }
        if (info[c].control) {
            error("Samplesheet: '${c}' is used as a control but has its own control '${info[c].control}'")
        }
    }

    def lane_index = [:]
    return rows.collect { row ->
        def sample = row.sample
        def s = info[sample]
        lane_index[sample] = (lane_index[sample] ?: 0) + 1
        def is_control = sample in controls
        def meta = [
            id        : sample,
            outdir    : s.outdir,
            mark      : s.mark,
            control   : s.control,
            is_control: is_control,
            single_end: s.single_end,
            peak_type : s.mark.toLowerCase() in broad_marks ? 'broad' : 'narrow',
            // ChIP: remove duplicates everywhere; CUT&RUN: only in controls unless --cutrun_dedup_targets
            dedup     : assay == 'chip' || is_control || params.cutrun_dedup_targets,
            lane      : lane_index[sample],
            n_lanes   : s.n_lanes
        ]
        def fastqs = [ file(row.fastq_1, checkIfExists: true) ]
        if (row.fastq_2) {
            fastqs << file(row.fastq_2, checkIfExists: true)
        }
        [meta, fastqs]
    }
}
