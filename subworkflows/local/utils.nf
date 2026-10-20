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
    def supported = ['atac', 'chip', 'cutrun', 'hichip']
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
    if (params.assay != 'hichip' && !params.gtf && !getGenomeAttribute('gtf')) {
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
      nextflow run zcmit/peddep-nf --assay <atac|chip|cutrun|hichip> --input samplesheet.csv --outdir results \\
          --genome GRCh38 -profile <docker|singularity|apptainer|podman|conda>[,<institution>]

    Required:
      --assay             atac (default), chip (ChIP-seq), cutrun (CUT&RUN) or hichip (HiChIP)
      --input             atac:         sample,fastq_1,fastq_2[,outdir]
                          chip/cutrun:  sample,fastq_1,fastq_2,mark,control[,outdir]
                          hichip:       sample,fastq_1,fastq_2,mark,control[,peaks][,outdir]
                          (repeat a sample on several rows for multiple lanes; fastq_2 empty =
                          single-end (not HiChIP); control = sample name of the Input/IgG/INPUT
                          row; peaks = optional BED of loop anchors, e.g. ChIP-seq peaks)
      --outdir            Run-level output directory (MultiQC, QC summary, pipeline_info)
      --genome | --fasta  A key from conf/genomes.config (GRCh38, hg38, GRCm38, mm10)
                          or an explicit FASTA (local or remote, may be .gz)

    References (override or replace --genome):
      --gtf               Gene annotation, used for TSS enrichment (not needed for HiChIP)
      --blacklist         BED of regions to exclude
      --bowtie2_index     Directory or prefix of a prebuilt index (built if absent)
      --bwamem2_index     HiChIP: directory or prefix of a prebuilt bwa-mem2 index (built if absent)
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
      --rose_peak_pvalue  Narrow peaks with p <= this go into ROSE [${params.rose_peak_pvalue}]
      --se_fragment_length  Single-end fragment length for MACS3/bigWigs [${params.se_fragment_length}]
      --cutrun_dedup_targets  Also remove duplicates from CUT&RUN targets [${params.cutrun_dedup_targets}]

    HiChIP:
      --split_fastq       Read pairs per alignment chunk, 0 = whole lanes [${params.split_fastq}]
      --read_length       Read length MACS3 assumes for the read-end records [${params.read_length}]
      --fithichip_bin_size / --fithichip_low_dist / --fithichip_upp_dist
                          FitHiChIP bin size and loop distance range [${params.fithichip_bin_size} / ${params.fithichip_low_dist} / ${params.fithichip_upp_dist}]
      --fithichip_int_type  3 = peak-to-all, 1 = peak-to-peak [${params.fithichip_int_type}]
      --fithichip_qvalues   Loop q-value cutoffs written as tracks [${params.fithichip_qvalues}]
      --save_pairs        Publish the deduplicated contact pairs [${params.save_pairs}]

    Resources:
      --max_cpus / --max_memory / --max_time   Per-task ceilings [${params.max_cpus} / ${params.max_memory} / ${params.max_time}]

    See README.md for outputs, site configuration and cloud/HPC examples.
    """.stripIndent()
}

/*
 * ChIP-seq / CUT&RUN / HiChIP samplesheet: sample,fastq_1,fastq_2,mark,control[,peaks][,outdir]
 *   - fastq_2 empty = single-end (ChIP-seq / CUT&RUN only); rows sharing a sample name are lanes
 *   - control = the sample name of that library's Input/IgG row (empty for controls
 *     and for cases without a control); one control can serve many cases
 *   - peaks (HiChIP only, optional) = BED of loop anchors to use instead of the
 *     library's own MACS3 peaks; a relative path is taken from the samplesheet's folder
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
        if (assay == 'hichip' && !row.fastq_2) {
            error("Samplesheet line ${line}: HiChIP needs both fastq_1 and fastq_2")
        }
        if (row.peaks && assay != 'hichip') {
            error("Samplesheet line ${line}: the peaks column is only used with --assay hichip")
        }
        if (!row.mark) {
            error("Samplesheet line ${line}: mark is required (e.g. H3K27ac, CTCF, Input, IgG, INPUT)")
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
            peaks     : row.peaks ?: '',
            single_end: !row.fastq_2
        ]
        if (info.containsKey(sample)) {
            def first = info[sample]
            ['outdir', 'mark', 'control', 'peaks', 'single_end'].each { field ->
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
        if (info[c].peaks) {
            error("Samplesheet: '${c}' is used as a control but has a peaks file; controls get no peaks or loops")
        }
    }

    // relative peaks paths are read from the samplesheet's folder
    def sheet_dir = file(path).parent
    info.each { _sample, v ->
        if (v.peaks && !(v.peaks ==~ /^[A-Za-z][A-Za-z0-9+.-]*:\/\/.*/) && !v.peaks.startsWith('/')) {
            v.peaks = sheet_dir.resolve(v.peaks).toUriString()
        }
        if (v.peaks) {
            file(v.peaks, checkIfExists: true)
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
        if (assay == 'hichip') {
            // kept out of ChIP/CUT&RUN metas so their task hashes (and -resume) are unchanged
            meta.peaks = s.peaks
        }
        def fastqs = [ file(row.fastq_1, checkIfExists: true) ]
        if (row.fastq_2) {
            fastqs << file(row.fastq_2, checkIfExists: true)
        }
        [meta, fastqs]
    }
}
