/*
========================================================================================
    ATAC-seq: per lane -> per sample -> per run
========================================================================================
*/

include { PREPARE_GENOME                                } from '../subworkflows/local/prepare_genome'
include { FASTQC as FASTQC_RAW                          } from '../modules/local/fastqc'
include { TRIMGALORE                                    } from '../modules/local/trimgalore'
include { BOWTIE2_ALIGN                                 } from '../modules/local/bowtie2_align'
include { SAMTOOLS_MERGE                                } from '../modules/local/samtools_merge'
include { PICARD_MARKDUPLICATES                         } from '../modules/local/picard_markduplicates'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_MARKDUP      } from '../modules/local/samtools_stats'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_FILTERED     } from '../modules/local/samtools_stats'
include { LIBRARY_COMPLEXITY                            } from '../modules/local/library_complexity'
include { BAM_FILTER                                    } from '../modules/local/bam_filter'
include { TN5_SHIFT                                     } from '../modules/local/tn5_shift'
include { MACS3_CALLPEAK as MACS3_NARROW                } from '../modules/local/macs3_callpeak'
include { MACS3_CALLPEAK as MACS3_BROAD                 } from '../modules/local/macs3_callpeak'
include { FRIP                                          } from '../modules/local/frip'
include { BROWSER_PEAKS                                 } from '../modules/local/browser_peaks'
include { DEEPTOOLS_BAMCOVERAGE as BIGWIG_FRAGMENTS     } from '../modules/local/deeptools_bamcoverage'
include { DEEPTOOLS_BAMCOVERAGE as BIGWIG_INSERTIONS    } from '../modules/local/deeptools_bamcoverage'
include { DEEPTOOLS_TSS_PLOTS                           } from '../modules/local/deeptools_tss_plots'
include { DEEPTOOLS_FINGERPRINT                         } from '../modules/local/deeptools_fingerprint'
include { FRAGMENT_QC                                   } from '../modules/local/fragment_qc'
include { ATAQV                                         } from '../modules/local/ataqv'
include { MKARV                                         } from '../modules/local/mkarv'
include { ATAC_QC_REPORT                                } from '../modules/local/atac_qc_report'
include { ATAC_QC_SUMMARY                               } from '../modules/local/atac_qc_summary'
include { MULTIQC                                       } from '../modules/local/multiqc'
include { toSampleMeta                                  } from '../subworkflows/local/utils'

workflow ATAC {
    take:
    ch_reads    // [ meta(id, outdir, lane, n_lanes), [ fastq_1, fastq_2 ] ] per lane

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
    // Per sample: merge lanes (a sample's group is released as soon as all its lanes finish)
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
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)
    ch_sample_bam = ch_lane_bams.single.mix(SAMTOOLS_MERGE.out.bam)

    PICARD_MARKDUPLICATES(ch_sample_bam)
    ch_markdup = PICARD_MARKDUPLICATES.out.bam
    SAMTOOLS_STATS_MARKDUP(ch_markdup)
    LIBRARY_COMPLEXITY(ch_markdup, PREPARE_GENOME.out.mito_name)
    ch_versions = ch_versions
        .mix(PICARD_MARKDUPLICATES.out.versions)
        .mix(SAMTOOLS_STATS_MARKDUP.out.versions)
        .mix(LIBRARY_COMPLEXITY.out.versions)

    //
    // Analysis BAM, Tn5-shifted BAM, peaks, coverage
    //
    BAM_FILTER(ch_markdup, PREPARE_GENOME.out.include_bed)
    ch_filtered = BAM_FILTER.out.bam
    SAMTOOLS_STATS_FILTERED(ch_filtered)
    TN5_SHIFT(ch_filtered)

    MACS3_NARROW(ch_filtered, PREPARE_GENOME.out.gsize)
    if (!params.skip_broad_peaks) {
        MACS3_BROAD(ch_filtered, PREPARE_GENOME.out.gsize)
        ch_versions = ch_versions.mix(MACS3_BROAD.out.versions)
    }

    FRIP(ch_filtered.join(MACS3_NARROW.out.peak))
    BROWSER_PEAKS(MACS3_NARROW.out.peak)
    BIGWIG_FRAGMENTS(ch_filtered)
    BIGWIG_INSERTIONS(TN5_SHIFT.out.bam)
    ch_versions = ch_versions
        .mix(BAM_FILTER.out.versions)
        .mix(SAMTOOLS_STATS_FILTERED.out.versions)
        .mix(TN5_SHIFT.out.versions)
        .mix(MACS3_NARROW.out.versions)
        .mix(FRIP.out.versions)
        .mix(BROWSER_PEAKS.out.versions)
        .mix(BIGWIG_FRAGMENTS.out.versions)
        .mix(BIGWIG_INSERTIONS.out.versions)

    //
    // ATAC-specific QC
    //
    DEEPTOOLS_TSS_PLOTS(BIGWIG_FRAGMENTS.out.bigwig, PREPARE_GENOME.out.tss)
    DEEPTOOLS_FINGERPRINT(ch_filtered)
    FRAGMENT_QC(ch_filtered)
    ATAQV(
        ch_markdup.join(MACS3_NARROW.out.peak),
        PREPARE_GENOME.out.tss,
        PREPARE_GENOME.out.excluded,
        PREPARE_GENOME.out.autosomes,
        PREPARE_GENOME.out.mito_name
    )
    MKARV(ATAQV.out.json.map { _meta, json -> json }.collect())
    ch_versions = ch_versions
        .mix(DEEPTOOLS_TSS_PLOTS.out.versions)
        .mix(DEEPTOOLS_FINGERPRINT.out.versions)
        .mix(FRAGMENT_QC.out.versions)
        .mix(ATAQV.out.versions)
        .mix(MKARV.out.versions)

    //
    // PEDDEP graded report per sample, then one summary table for the run
    //
    ch_raw_fastqc = groupLaneFiles(FASTQC_RAW.out.zip)
    ch_trim_logs  = groupLaneFiles(TRIMGALORE.out.log)
    ch_bt2_logs   = groupLaneFiles(BOWTIE2_ALIGN.out.log)

    ch_report_inputs = ch_raw_fastqc
        .join(ch_trim_logs)
        .join(ch_bt2_logs)
        .join(SAMTOOLS_STATS_MARKDUP.out.idxstats)
        .join(PICARD_MARKDUPLICATES.out.metrics)
        .join(LIBRARY_COMPLEXITY.out.txt)
        .join(SAMTOOLS_STATS_FILTERED.out.flagstat)
        .join(FRIP.out.txt)
        .join(ATAQV.out.json)
        .join(FRAGMENT_QC.out.metrics)

    ATAC_QC_REPORT(
        ch_report_inputs,
        file(params.qc_thresholds, checkIfExists: true),
        PREPARE_GENOME.out.mito_name
    )
    ATAC_QC_SUMMARY(ATAC_QC_REPORT.out.json.map { _meta, json -> json }.collect())
    ch_versions = ch_versions
        .mix(ATAC_QC_REPORT.out.versions)
        .mix(ATAC_QC_SUMMARY.out.versions)

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
        .mix(SAMTOOLS_STATS_MARKDUP.out.idxstats.map { _meta, txt -> txt })
        .mix(SAMTOOLS_STATS_MARKDUP.out.stats.map { _meta, txt -> txt })
        .mix(SAMTOOLS_STATS_FILTERED.out.flagstat.map { _meta, txt -> txt })
        .mix(SAMTOOLS_STATS_FILTERED.out.stats.map { _meta, txt -> txt })
        .mix(DEEPTOOLS_TSS_PLOTS.out.profile_data.map { _meta, txt -> txt })
        .mix(DEEPTOOLS_FINGERPRINT.out.raw.map { _meta, txt -> txt })
        .mix(DEEPTOOLS_FINGERPRINT.out.metrics.map { _meta, txt -> txt })
        .mix(ATAQV.out.json.map { _meta, json -> json })
        .mix(ATAC_QC_SUMMARY.out.multiqc)
        .mix(ch_versions_yml)
        .collect()

    MULTIQC(ch_multiqc_files, file(params.multiqc_config, checkIfExists: true))

    emit:
    multiqc_report = MULTIQC.out.report
    qc_summary     = ATAC_QC_SUMMARY.out.tsv
}

// Collect per-lane files of each sample into one list: [ sample_meta, [ files ] ]
def groupLaneFiles(ch) {
    return ch
        .map { meta, files -> [ groupKey(toSampleMeta(meta), meta.n_lanes), files ] }
        .groupTuple()
        .map { key, files -> [ key.getGroupTarget(), files.flatten() ] }
}
