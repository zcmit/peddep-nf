/*
========================================================================================
    ChIP-seq / CUT&RUN (--assay chip | cutrun)

    Every library -- cases and their Input/IgG controls -- is trimmed, aligned and
    filtered once; a control shared by several cases is reused, not re-aligned.
    Cases are then paired with their control for MACS3 narrow + broad peaks,
    case-vs-control tracks, fingerprints and (for --rose_marks) ROSE.
========================================================================================
*/

include { PREPARE_GENOME                                } from '../subworkflows/local/prepare_genome'
include { GTF2ROSE                                      } from '../modules/local/gtf2rose'
include { FASTQC as FASTQC_RAW                          } from '../modules/local/fastqc'
include { TRIMGALORE                                    } from '../modules/local/trimgalore'
include { BOWTIE2_ALIGN                                 } from '../modules/local/bowtie2_align'
include { SAMTOOLS_MERGE                                } from '../modules/local/samtools_merge'
include { PICARD_MARKDUPLICATES                         } from '../modules/local/picard_markduplicates'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_MARKDUP      } from '../modules/local/samtools_stats'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_FILTERED     } from '../modules/local/samtools_stats'
include { LIBRARY_COMPLEXITY                            } from '../modules/local/library_complexity'
include { BAM_FILTER                                    } from '../modules/local/bam_filter'
include { PHANTOMPEAKQUALTOOLS                          } from '../modules/local/phantompeakqualtools'
include { DEEPTOOLS_BAMCOVERAGE as BIGWIG_COVERAGE      } from '../modules/local/deeptools_bamcoverage'
include { DEEPTOOLS_BAMCOMPARE as BIGWIG_VS_CONTROL     } from '../modules/local/deeptools_bamcompare'
include { DEEPTOOLS_FINGERPRINT                         } from '../modules/local/deeptools_fingerprint'
include { DEEPTOOLS_TSS_PLOTS                           } from '../modules/local/deeptools_tss_plots'
include { MACS3_CALLPEAK as MACS3_NARROW                } from '../modules/local/macs3_callpeak'
include { MACS3_CALLPEAK as MACS3_BROAD                 } from '../modules/local/macs3_callpeak'
include { BROWSER_PEAKS as BROWSER_PEAKS_NARROW         } from '../modules/local/browser_peaks'
include { BROWSER_PEAKS as BROWSER_PEAKS_BROAD          } from '../modules/local/browser_peaks'
include { FRIP                                          } from '../modules/local/frip'
include { ROSE                                          } from '../modules/local/rose'
include { CHIP_QC_REPORT                                } from '../modules/local/chip_qc_report'
include { QC_SUMMARY as CHIP_QC_SUMMARY                 } from '../modules/local/qc_summary'
include { MULTIQC                                       } from '../modules/local/multiqc'
include { toSampleMeta; groupLaneFiles                  } from '../subworkflows/local/utils'

workflow CHIP {
    take:
    ch_reads    // [ meta(id, outdir, lane, n_lanes, single_end, mark, control, is_control, peak_type, dedup), fastqs ]

    main:
    ch_versions = channel.empty()

    PREPARE_GENOME()
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    //
    // Per lane: raw QC, adapter trimming, alignment
    //
    FASTQC_RAW(ch_reads)
    TRIMGALORE(ch_reads)
    BOWTIE2_ALIGN(TRIMGALORE.out.reads, PREPARE_GENOME.out.index)
    ch_versions = ch_versions
        .mix(FASTQC_RAW.out.versions)
        .mix(TRIMGALORE.out.versions)
        .mix(BOWTIE2_ALIGN.out.versions)

    //
    // Per library (cases and controls alike): merge lanes, mark duplicates, filter
    //
    ch_lane_bams = BOWTIE2_ALIGN.out.bam
        .map { meta, bam -> [ groupKey(toSampleMeta(meta), meta.n_lanes), bam ] }
        .groupTuple()
        .map { key, bams -> [ key.getGroupTarget(), bams ] }
        .branch { meta, bams ->
            single: bams.size() == 1
                return [ meta, bams[0] ]
            multiple: true
        }
    SAMTOOLS_MERGE(ch_lane_bams.multiple)
    ch_sample_bam = ch_lane_bams.single.mix(SAMTOOLS_MERGE.out.bam)

    PICARD_MARKDUPLICATES(ch_sample_bam)
    ch_markdup = PICARD_MARKDUPLICATES.out.bam
    SAMTOOLS_STATS_MARKDUP(ch_markdup)
    LIBRARY_COMPLEXITY(ch_markdup, PREPARE_GENOME.out.mito_name)

    BAM_FILTER(ch_markdup, PREPARE_GENOME.out.include_bed)
    ch_filtered = BAM_FILTER.out.bam
    SAMTOOLS_STATS_FILTERED(ch_filtered)
    PHANTOMPEAKQUALTOOLS(ch_filtered)
    BIGWIG_COVERAGE(ch_filtered)
    ch_versions = ch_versions
        .mix(SAMTOOLS_MERGE.out.versions)
        .mix(PICARD_MARKDUPLICATES.out.versions)
        .mix(SAMTOOLS_STATS_MARKDUP.out.versions)
        .mix(LIBRARY_COMPLEXITY.out.versions)
        .mix(BAM_FILTER.out.versions)
        .mix(SAMTOOLS_STATS_FILTERED.out.versions)
        .mix(PHANTOMPEAKQUALTOOLS.out.versions)
        .mix(BIGWIG_COVERAGE.out.versions)

    //
    // Pair each case with its control (one control can serve many cases)
    //
    ch_by_id = ch_filtered.map { meta, bam, bai -> [ meta.id, bam, bai ] }
    ch_cases = ch_filtered.filter { meta, _bam, _bai -> !meta.is_control }
    ch_with_control = ch_cases
        .filter { meta, _bam, _bai -> meta.control }
        .map { meta, bam, bai -> [ meta.control, meta, bam, bai ] }
        .combine(ch_by_id, by: 0)
        .map { _control_id, meta, bam, bai, control_bam, control_bai -> [ meta, bam, bai, control_bam, control_bai ] }
    ch_case_only = ch_cases
        .filter { meta, _bam, _bai -> !meta.control }
        .map { meta, bam, bai -> [ meta, bam, bai, [], [] ] }
    ch_peak_input = ch_with_control.mix(ch_case_only)

    //
    // Peaks, browser tracks, enrichment QC
    //
    MACS3_NARROW(ch_peak_input, PREPARE_GENOME.out.gsize)
    MACS3_BROAD(ch_peak_input, PREPARE_GENOME.out.gsize)
    BROWSER_PEAKS_NARROW(MACS3_NARROW.out.peak)
    BROWSER_PEAKS_BROAD(MACS3_BROAD.out.peak)
    BIGWIG_VS_CONTROL(ch_with_control)
    DEEPTOOLS_FINGERPRINT(
        ch_peak_input.map { meta, bam, bai, control_bam, control_bai ->
            control_bam ? [ meta, [ bam, control_bam ], [ bai, control_bai ] ] : [ meta, [ bam ], [ bai ] ]
        }
    )
    DEEPTOOLS_TSS_PLOTS(
        BIGWIG_COVERAGE.out.bigwig.filter { meta, _bw -> !meta.is_control },
        PREPARE_GENOME.out.tss
    )

    // FRiP on the mark's primary peak set (broad for --broad_marks, else narrow)
    ch_primary_peaks = MACS3_NARROW.out.peak
        .join(MACS3_BROAD.out.peak)
        .map { meta, narrow, broad -> [ meta, meta.peak_type == 'broad' ? broad : narrow ] }
    FRIP(ch_filtered.join(ch_primary_peaks))

    ch_versions = ch_versions
        .mix(MACS3_NARROW.out.versions)
        .mix(MACS3_BROAD.out.versions)
        .mix(BROWSER_PEAKS_NARROW.out.versions)
        .mix(BIGWIG_VS_CONTROL.out.versions)
        .mix(DEEPTOOLS_FINGERPRINT.out.versions)
        .mix(DEEPTOOLS_TSS_PLOTS.out.versions)
        .mix(FRIP.out.versions)

    //
    // ROSE stitched / super-enhancer-like regions for the chosen marks only
    //
    def rose_marks = (params.rose_marks ?: '').tokenize(',').collect { m -> m.trim().toLowerCase() }
    ch_rose_tables = channel.empty()
    if (rose_marks) {
        GTF2ROSE(PREPARE_GENOME.out.gtf)
        ROSE(
            ch_peak_input
                .filter { meta, _bam, _bai, _cbam, _cbai -> meta.mark.toLowerCase() in rose_marks }
                .join(MACS3_NARROW.out.peak),
            GTF2ROSE.out.refseq
        )
        ch_rose_tables = ROSE.out.tables
        ch_versions = ch_versions.mix(GTF2ROSE.out.versions).mix(ROSE.out.versions)
    }

    //
    // PEDDEP graded report per library, then one summary table for the run
    //
    ch_report_inputs = groupLaneFiles(FASTQC_RAW.out.zip)
        .join(groupLaneFiles(TRIMGALORE.out.log))
        .join(groupLaneFiles(BOWTIE2_ALIGN.out.log))
        .join(PICARD_MARKDUPLICATES.out.metrics)
        .join(LIBRARY_COMPLEXITY.out.txt)
        .join(SAMTOOLS_STATS_FILTERED.out.flagstat)
        .join(PHANTOMPEAKQUALTOOLS.out.spp)
        .join(FRIP.out.txt, remainder: true)
        .join(ch_rose_tables, remainder: true)
        .map { meta, fastqc, trim, bt2, dup, complexity, flagstat, spp, frip, rose ->
            [ meta, fastqc, trim, bt2, dup, complexity, flagstat, spp, frip ?: [], rose ?: [] ]
        }

    def default_thresholds = "${projectDir}/assets/${params.assay == 'cutrun' ? 'cutrun' : 'chip'}_qc_thresholds.json"
    CHIP_QC_REPORT(ch_report_inputs, file(params.qc_thresholds ?: default_thresholds, checkIfExists: true))
    CHIP_QC_SUMMARY(CHIP_QC_REPORT.out.json.map { _meta, json -> json }.collect())
    ch_versions = ch_versions
        .mix(CHIP_QC_REPORT.out.versions)
        .mix(CHIP_QC_SUMMARY.out.versions)

    //
    // Software versions + MultiQC
    //
    ch_versions_yml = ch_versions
        .map { yml -> yml.text }
        .unique()
        .collectFile(name: 'peddep_software_mqc_versions.yml', storeDir: "${params.outdir}/pipeline_info", newLine: true, sort: true)

    ch_multiqc_files = channel.empty()
        .mix(FASTQC_RAW.out.zip.map { _meta, zip -> zip })
        .mix(TRIMGALORE.out.log.map { _meta, log -> log })
        .mix(TRIMGALORE.out.zip.map { _meta, zip -> zip })
        .mix(BOWTIE2_ALIGN.out.log.map { _meta, log -> log })
        .mix(PICARD_MARKDUPLICATES.out.metrics.map { _meta, txt -> txt })
        .mix(SAMTOOLS_STATS_MARKDUP.out.flagstat.map { _meta, txt -> txt })
        .mix(SAMTOOLS_STATS_MARKDUP.out.stats.map { _meta, txt -> txt })
        .mix(SAMTOOLS_STATS_FILTERED.out.flagstat.map { _meta, txt -> txt })
        .mix(SAMTOOLS_STATS_FILTERED.out.stats.map { _meta, txt -> txt })
        .mix(PHANTOMPEAKQUALTOOLS.out.spp.map { _meta, txt -> txt })
        .mix(MACS3_NARROW.out.xls.map { _meta, xls -> xls })
        .mix(DEEPTOOLS_FINGERPRINT.out.raw.map { _meta, txt -> txt })
        .mix(DEEPTOOLS_FINGERPRINT.out.metrics.map { _meta, txt -> txt })
        .mix(DEEPTOOLS_TSS_PLOTS.out.profile_data.map { _meta, txt -> txt })
        .mix(CHIP_QC_SUMMARY.out.multiqc)
        .mix(ch_versions_yml)
        .collect()

    MULTIQC(ch_multiqc_files, file(params.multiqc_config ?: "${projectDir}/assets/multiqc_config_chip.yml", checkIfExists: true))

    emit:
    multiqc_report = MULTIQC.out.report
    qc_summary     = CHIP_QC_SUMMARY.out.tsv
}
