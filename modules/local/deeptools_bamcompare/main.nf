/*
 * Case-vs-control signal track (e.g. input-subtracted CPM, as in pyflow-ChIPseq).
 */
process DEEPTOOLS_BAMCOMPARE {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::deeptools=3.5.6 bioconda::samtools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' :
        'biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(control), path(control_bai)

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    path "versions.yml"          , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bamCompare \\
        --bamfile1 ${bam} \\
        --bamfile2 ${control} \\
        --outFileName ${prefix}.bw \\
        --outFileFormat bigwig \\
        --numberOfProcessors ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(deeptools --version | sed "s/^deeptools //")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bw versions.yml
    """
}
