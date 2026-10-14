process PICARD_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::picard=3.5.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b474c6f12c0502377b95062b75ef4b7778864b1353c7fc1d5a3b1f3b3017fd2e/data' :
        'community.wave.seqera.io/library/picard:3.5.0--842d4c70c98af9b4' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.markdup.bam"), path("*.markdup.bam.bai"), emit: bam
    tuple val(meta), path("*.markdup.metrics.txt")                   , emit: metrics
    path "versions.yml"                                              , emit: versions

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def java_mem = task.memory ? (task.memory.toMega() * 0.8).intValue() : 3072
    """
    mkdir -p tmp
    picard -Xmx${java_mem}M MarkDuplicates \\
        ${args} \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}.markdup.bam \\
        --METRICS_FILE ${prefix}.markdup.metrics.txt \\
        --ASSUME_SORT_ORDER coordinate \\
        --REMOVE_DUPLICATES false \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY LENIENT \\
        --TMP_DIR tmp
    mv ${prefix}.markdup.bai ${prefix}.markdup.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard MarkDuplicates --version 2>&1 | grep -Eo '[0-9]+\\.[0-9]+\\.[0-9]+' | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.markdup.bam ${prefix}.markdup.bam.bai ${prefix}.markdup.metrics.txt versions.yml
    """
}
