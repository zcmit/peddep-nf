process MACS3_CALLPEAK {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::macs3=3.0.4"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2f/2fb492856efb63a7f824f0801b1386d08468cd4b7819ddc4c21e7f10e09b4fda/data' :
        'community.wave.seqera.io/library/macs3:3.0.4--e0346d811b8b428e' }"

    input:
    tuple val(meta), path(bam), path(bai), path(control), path(control_bai)
    val gsize

    output:
    tuple val(meta), path("*_peaks.{narrowPeak,broadPeak}"), emit: peak
    tuple val(meta), path("*_peaks.xls")                   , emit: xls
    tuple val(meta), path("*_summits.bed")                 , emit: summits, optional: true
    tuple val(meta), path("*_peaks.gappedPeak")            , emit: gapped , optional: true
    path "versions.yml"                                    , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ctrl   = control ? "--control ${control}" : ''
    def format = meta.single_end ? 'BAM' : 'BAMPE'
    """
    mkdir -p macs_tmp
    macs3 callpeak \\
        --treatment ${bam} \\
        ${ctrl} \\
        --format ${format} \\
        --gsize ${gsize} \\
        --name ${prefix} \\
        --keep-dup all \\
        --tempdir macs_tmp \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        macs3: \$(macs3 --version | sed 's/^macs3 //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext    = (task.ext.args ?: '').contains('--broad') ? 'broadPeak' : 'narrowPeak'
    """
    printf 'chr1\\t100\\t500\\tp1\\t50\\t.\\t5\\t10\\t8\\t200\\n' > ${prefix}_peaks.${ext}
    touch ${prefix}_peaks.xls versions.yml
    """
}
