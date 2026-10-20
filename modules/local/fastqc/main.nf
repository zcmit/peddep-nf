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
    def fastqs = [reads].flatten()
    def ext    = fastqs[0].name.endsWith('.gz') ? 'fastq.gz' : 'fastq'
    def names  = meta.single_end ? ["${prefix}.${ext}"] : ["${prefix}_R1.${ext}", "${prefix}_R2.${ext}"]
    def links  = [fastqs, names].transpose().collect { pair -> "ln -s ${pair[0]} ${pair[1]}" }.join('\n    ')
    """
    ${links}
    fastqc ${args} --threads ${task.cpus} ${names.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed 's/^FastQC v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def names  = meta.single_end ? ["${prefix}"] : ["${prefix}_R1", "${prefix}_R2"]
    """
    for r in ${names.join(' ')}; do touch \${r}_fastqc.html \${r}_fastqc.zip; done
    touch versions.yml
    """
}
