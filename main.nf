#!/usr/bin/env nextflow
/*
========================================================================================
    zcmit/peddep-nf
    Portable Nextflow pipelines for PEDDEP epigenomic assays
    https://github.com/zcmit/peddep-nf
========================================================================================
*/

include { ATAC                                        } from './workflows/atac'
include { CHIP                                        } from './workflows/chip'
include { HICHIP                                      } from './workflows/hichip'
include { parseSamplesheet; parseChipSamplesheet; validateParameters; helpMessage } from './subworkflows/local/utils'

workflow {
    main:
    if (params.help) {
        log.info(helpMessage())
    } else {
        validateParameters()

        log.info """
        zcmit/peddep-nf v${workflow.manifest.version}  (assay: ${params.assay})
          input   : ${params.input}
          outdir  : ${params.outdir}
          genome  : ${params.genome ?: params.fasta}
          profile : ${workflow.profile}
          workDir : ${workflow.workDir}
        """.stripIndent()

        if (params.assay == 'atac') {
            ATAC(channel.fromList(parseSamplesheet(params.input)))
        } else if (params.assay == 'hichip') {
            HICHIP(channel.fromList(parseChipSamplesheet(params.input, params.assay)))
        } else {
            CHIP(channel.fromList(parseChipSamplesheet(params.input, params.assay)))
        }
    }
}
