process FASTQC {
    tag "${meta.id}_L${meta.lane}"
    label 'process_low'

    conda "bioconda::fastqc=0.13.0"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/816bda6dd014e810753f642611232443ddeabb2b47f9681ad43f9044816f00ff/data' :
        'community.wave.seqera.io/library/fastqc:0.13.0--b342db7b694c9af4' }"

    input:
    tuple val(meta), path(reads, stageAs: 'input*/*')

    output:
    tuple val(meta), path("*_fastqc.html"), emit: html
    tuple val(meta), path("*_fastqc.zip") , emit: zip
    path "versions.yml"                   , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext    = reads[0].name.endsWith('.gz') ? 'fastq.gz' : 'fastq'
    """
    ln -s ${reads[0]} ${prefix}_R1.${ext}
    ln -s ${reads[1]} ${prefix}_R2.${ext}
    fastqc ${args} --threads ${task.cpus} ${prefix}_R1.${ext} ${prefix}_R2.${ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/^FastQC v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for r in R1 R2; do touch ${prefix}_\${r}_fastqc.html ${prefix}_\${r}_fastqc.zip; done
    touch versions.yml
    """
}
