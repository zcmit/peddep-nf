process DEEPTOOLS_FINGERPRINT {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::deeptools=3.5.6 bioconda::samtools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' :
        'biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.plotFingerprint.pdf")          , emit: pdf
    tuple val(meta), path("*.plotFingerprint.raw.txt")      , emit: raw
    tuple val(meta), path("*.plotFingerprint.qcmetrics.txt"), emit: metrics
    path "versions.yml"                                     , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib

    plotFingerprint \\
        --bamfiles ${bam} \\
        --labels ${meta.id} \\
        --extendReads \\
        --skipZeros \\
        --numberOfProcessors ${task.cpus} \\
        ${args} \\
        --plotFile ${prefix}.plotFingerprint.pdf \\
        --outRawCounts ${prefix}.plotFingerprint.raw.txt \\
        --outQualityMetrics ${prefix}.plotFingerprint.qcmetrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(deeptools --version | sed "s/^deeptools //")
        samtools: \$(samtools --version | head -1 | sed "s/^samtools //")
        python: \$(python3 --version | sed "s/^Python //")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.plotFingerprint.pdf ${prefix}.plotFingerprint.raw.txt ${prefix}.plotFingerprint.qcmetrics.txt versions.yml
    """
}
