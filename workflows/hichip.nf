/*
========================================================================================
    HiChIP (--assay hichip)

    Every library -- targets and their INPUT HiChIP controls -- is split into chunks,
    aligned with bwa-mem2 (-5SP) and turned into contact pairs with pairtools; the
    chunks of a library are merged and deduplicated once. Each library gets a .hic
    contact map, a coverage track and a graded QC report. Targets also get MACS3
    peaks against their INPUT (the loop anchors, unless the samplesheet gives a
    `peaks` file) and FitHiChIP loops; INPUT libraries get neither.
========================================================================================
*/

include { PREPARE_GENOME                                } from '../subworkflows/local/prepare_genome'
include { FASTQC as FASTQC_RAW                          } from '../modules/local/fastqc'
include { SEQKIT_SPLIT2                                 } from '../modules/local/seqkit_split2'
include { BWAMEM2_MEM                                   } from '../modules/local/bwamem2_mem'
include { PAIRTOOLS_PARSE                               } from '../modules/local/pairtools_parse'
include { PAIRTOOLS_DEDUP                               } from '../modules/local/pairtools_dedup'
include { JUICER_PRE                                    } from '../modules/local/juicer_pre'
include { DEEPTOOLS_BAMCOVERAGE as BIGWIG_COVERAGE      } from '../modules/local/deeptools_bamcoverage'
include { MACS3_CALLPEAK as MACS3_NARROW                } from '../modules/local/macs3_callpeak'
include { BROWSER_PEAKS as BROWSER_PEAKS_NARROW         } from '../modules/local/browser_peaks'
include { FRIP                                          } from '../modules/local/frip'
include { FITHICHIP                                     } from '../modules/local/fithichip'
include { HICHIP_QC_REPORT                              } from '../modules/local/hichip_qc_report'
include { QC_SUMMARY as HICHIP_QC_SUMMARY               } from '../modules/local/qc_summary'
include { MULTIQC                                       } from '../modules/local/multiqc'
include { toSampleMeta; groupLaneFiles                  } from '../subworkflows/local/utils'

// Lane-level meta map: drops the chunk fields so the chunks of one lane join together
def toLaneMeta(Map meta) {
    return meta.findAll { key, _value -> !(key in ['chunk', 'n_chunks']) }
}

workflow HICHIP {
    take:
    ch_reads    // [ meta(id, outdir, lane, n_lanes, mark, control, is_control, peaks), [ fastq_1, fastq_2 ] ]

    main:
    ch_versions = channel.empty()

    PREPARE_GENOME()
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // Per lane: raw QC, then chunks of --split_fastq read pairs that align side by side
    //
    FASTQC_RAW(ch_reads)
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions)

    ch_chunks = channel.empty()
    if (params.split_fastq) {
        SEQKIT_SPLIT2(ch_reads, params.split_fastq)
        ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)
        ch_chunks = SEQKIT_SPLIT2.out.reads.flatMap { meta, files ->
            def r1 = files.findAll { f -> f.name.contains('_R1.part_') }.sort { f -> f.name }
            def r2 = files.findAll { f -> f.name.contains('_R2.part_') }.sort { f -> f.name }
            [r1, r2].transpose().withIndex().collect { pair, i ->
                [ meta + [chunk: i + 1, n_chunks: r1.size()], pair ]
            }
        }
    } else {
        ch_chunks = ch_reads.map { meta, reads -> [ meta + [chunk: 1, n_chunks: 1], reads ] }
    }

    BWAMEM2_MEM(ch_chunks, PREPARE_GENOME.out.index)
    PAIRTOOLS_PARSE(BWAMEM2_MEM.out.bam, PREPARE_GENOME.out.sizes)
    ch_versions = ch_versions
        .mix(BWAMEM2_MEM.out.versions)
        .mix(PAIRTOOLS_PARSE.out.versions)

    //
    // Per library (targets and INPUTs alike): merge chunks and lanes, deduplicate
    //
    ch_library_pairs = PAIRTOOLS_PARSE.out.pairs
        .map { meta, pairs -> [ groupKey(toLaneMeta(meta), meta.n_chunks), pairs ] }
        .groupTuple()
        .map { key, pairs -> [ groupKey(toSampleMeta(key.getGroupTarget()), key.getGroupTarget().n_lanes), pairs ] }
        .groupTuple()
        .map { key, pairs -> [ key.getGroupTarget(), pairs.flatten() ] }

    PAIRTOOLS_DEDUP(ch_library_pairs, PREPARE_GENOME.out.sizes, PREPARE_GENOME.out.include_bed)
    JUICER_PRE(PAIRTOOLS_DEDUP.out.pairs, PREPARE_GENOME.out.sizes)
    BIGWIG_COVERAGE(PAIRTOOLS_DEDUP.out.bam)
    ch_versions = ch_versions
        .mix(PAIRTOOLS_DEDUP.out.versions)
        .mix(JUICER_PRE.out.versions)
        .mix(BIGWIG_COVERAGE.out.versions)

    //
    // Targets: MACS3 peaks on the read ends against the INPUT, FRiP, FitHiChIP loops
    //
    ch_ends    = PAIRTOOLS_DEDUP.out.bam
    ch_by_id   = ch_ends.map { meta, bam, bai -> [ meta.id, bam, bai ] }
    ch_targets = ch_ends.filter { meta, _bam, _bai -> !meta.is_control }
    ch_with_control = ch_targets
        .filter { meta, _bam, _bai -> meta.control }
        .map { meta, bam, bai -> [ meta.control, meta, bam, bai ] }
        .combine(ch_by_id, by: 0)
        .map { _control_id, meta, bam, bai, control_bam, control_bai -> [ meta, bam, bai, control_bam, control_bai ] }
    ch_target_only = ch_targets
        .filter { meta, _bam, _bai -> !meta.control }
        .map { meta, bam, bai -> [ meta, bam, bai, [], [] ] }

    MACS3_NARROW(ch_with_control.mix(ch_target_only), PREPARE_GENOME.out.gsize)
    BROWSER_PEAKS_NARROW(MACS3_NARROW.out.peak)
    FRIP(ch_targets.join(MACS3_NARROW.out.peak))

    // Anchors: the samplesheet's peaks file when given, else the MACS3 peaks
    ch_anchors = MACS3_NARROW.out.peak.map { meta, peaks ->
        [ meta, meta.peaks ? file(meta.peaks, checkIfExists: true) : peaks ]
    }
    FITHICHIP(
        PAIRTOOLS_DEDUP.out.pairs.join(ch_anchors),
        PREPARE_GENOME.out.sizes,
        PREPARE_GENOME.out.blacklist
    )
    ch_versions = ch_versions
        .mix(MACS3_NARROW.out.versions)
        .mix(BROWSER_PEAKS_NARROW.out.versions)
        .mix(FRIP.out.versions)
        .mix(FITHICHIP.out.versions)

    //
    // PEDDEP graded report per library, then one summary table for the run
    //
    ch_report_inputs = groupLaneFiles(FASTQC_RAW.out.zip)
        .join(PAIRTOOLS_DEDUP.out.stats)
        .join(FRIP.out.txt, remainder: true)
        .join(FITHICHIP.out.counts, remainder: true)
        .map { meta, fastqc, stats, frip, loops -> [ meta, fastqc, stats, frip ?: [], loops ?: [] ] }

    HICHIP_QC_REPORT(ch_report_inputs, file(params.qc_thresholds ?: "${projectDir}/assets/hichip_qc_thresholds.json", checkIfExists: true))
    HICHIP_QC_SUMMARY(HICHIP_QC_REPORT.out.json.map { _meta, json -> json }.collect())
    ch_versions = ch_versions
        .mix(HICHIP_QC_REPORT.out.versions)
        .mix(HICHIP_QC_SUMMARY.out.versions)

    //
    // Software versions + MultiQC
    //
    ch_versions_yml = ch_versions
        .map { yml -> yml.text }
        .unique()
        .collectFile(name: 'peddep_software_mqc_versions.yml', storeDir: "${params.outdir}/pipeline_info", newLine: true, sort: true)

    ch_multiqc_files = channel.empty()
        .mix(FASTQC_RAW.out.zip.map { _meta, zip -> zip })
        .mix(PAIRTOOLS_DEDUP.out.stats.map { _meta, stats -> stats })
        .mix(MACS3_NARROW.out.xls.map { _meta, xls -> xls })
        .mix(HICHIP_QC_SUMMARY.out.multiqc)
        .mix(ch_versions_yml)
        .collect()

    MULTIQC(ch_multiqc_files, file(params.multiqc_config ?: "${projectDir}/assets/multiqc_config_hichip.yml", checkIfExists: true))

    emit:
    multiqc_report = MULTIQC.out.report
    qc_summary     = HICHIP_QC_SUMMARY.out.tsv
}
