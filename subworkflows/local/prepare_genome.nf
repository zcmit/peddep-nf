/*
========================================================================================
    Reference preparation: decompress, index, filter regions, TSS, genome size
========================================================================================
*/

include { GUNZIP as GUNZIP_FASTA     } from '../../modules/local/gunzip'
include { GUNZIP as GUNZIP_GTF       } from '../../modules/local/gunzip'
include { GUNZIP as GUNZIP_BLACKLIST } from '../../modules/local/gunzip'
include { SAMTOOLS_FAIDX             } from '../../modules/local/samtools_faidx'
include { BOWTIE2_BUILD              } from '../../modules/local/bowtie2_build'
include { GENOME_FILTER              } from '../../modules/local/genome_filter'
include { GTF2TSS                    } from '../../modules/local/gtf2tss'
include { GENOME_SIZE                } from '../../modules/local/genome_size'
include { getGenomeAttribute         } from './utils'

workflow PREPARE_GENOME {
    main:
    ch_versions = channel.empty()

    def fasta_path     = params.fasta            ?: getGenomeAttribute('fasta')
    def gtf_path       = params.gtf              ?: getGenomeAttribute('gtf')
    def blacklist_path = params.blacklist        ?: getGenomeAttribute('blacklist')
    def index_path     = params.bowtie2_index    ?: getGenomeAttribute('bowtie2_index')
    def gsize_value    = params.macs_gsize       ?: getGenomeAttribute('macs_gsize')
    def mito           = params.mito_name        ?: (getGenomeAttribute('mito_name') ?: 'chrM')
    def canonical      = params.canonical_chroms ?: (getGenomeAttribute('canonical_chroms') ?: '')

    // FASTA
    ch_fasta = channel.empty()
    if (fasta_path.toString().endsWith('.gz')) {
        GUNZIP_FASTA(file(fasta_path, checkIfExists: true))
        ch_fasta    = GUNZIP_FASTA.out.file
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = channel.value(file(fasta_path, checkIfExists: true))
    }

    // GTF
    ch_gtf = channel.empty()
    if (gtf_path.toString().endsWith('.gz')) {
        GUNZIP_GTF(file(gtf_path, checkIfExists: true))
        ch_gtf      = GUNZIP_GTF.out.file
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    } else {
        ch_gtf = channel.value(file(gtf_path, checkIfExists: true))
    }

    // Blacklist (optional)
    ch_blacklist = channel.value([])
    if (blacklist_path) {
        if (blacklist_path.toString().endsWith('.gz')) {
            GUNZIP_BLACKLIST(file(blacklist_path, checkIfExists: true))
            ch_blacklist = GUNZIP_BLACKLIST.out.file
            ch_versions  = ch_versions.mix(GUNZIP_BLACKLIST.out.versions)
        } else {
            ch_blacklist = channel.value(file(blacklist_path, checkIfExists: true))
        }
    }

    SAMTOOLS_FAIDX(ch_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // bowtie2 index: a directory, a prefix (files <prefix>*.bt2*), or built here
    ch_index = channel.empty()
    if (index_path) {
        def index_file = file(index_path)
        if (index_file.isDirectory()) {
            ch_index = channel.value(index_file)
        } else {
            ch_index = channel.fromPath("${index_path}*.bt2*", checkIfExists: true).collect()
        }
    } else {
        BOWTIE2_BUILD(ch_fasta)
        ch_index    = BOWTIE2_BUILD.out.index
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    }

    GENOME_FILTER(
        SAMTOOLS_FAIDX.out.fai,
        ch_blacklist,
        mito,
        params.keep_mito,
        canonical,
        params.sex_chroms
    )
    ch_versions = ch_versions.mix(GENOME_FILTER.out.versions)

    GTF2TSS(ch_gtf)
    ch_versions = ch_versions.mix(GTF2TSS.out.versions)

    ch_gsize = channel.empty()
    if (gsize_value) {
        ch_gsize = channel.value(gsize_value.toString())
    } else {
        GENOME_SIZE(ch_fasta)
        ch_gsize    = GENOME_SIZE.out.size.map { size_file -> size_file.text.trim() }
        ch_versions = ch_versions.mix(GENOME_SIZE.out.versions)
    }

    // ataqv takes the blacklist only when there is one
    ch_excluded = blacklist_path ? GENOME_FILTER.out.blacklist : channel.value([])

    emit:
    fasta       = ch_fasta
    fai         = SAMTOOLS_FAIDX.out.fai
    index       = ch_index
    include_bed = GENOME_FILTER.out.include_bed
    autosomes   = GENOME_FILTER.out.autosomes
    excluded    = ch_excluded
    tss         = GTF2TSS.out.tss
    gtf         = ch_gtf
    gsize       = ch_gsize
    mito_name   = channel.value(mito)
    versions    = ch_versions
}
