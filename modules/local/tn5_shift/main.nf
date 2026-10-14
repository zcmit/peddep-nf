/*
 * Tn5 shift (+4 bp on the forward strand, -5 bp on the reverse strand) so read
 * 5' ends mark the centre of the transposase insertion.
 */
process TN5_SHIFT {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::deeptools=3.5.6 bioconda::samtools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' :
        'biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.tn5.bam"), path("*.tn5.bam.bai"), emit: bam
    path "versions.yml"                                      , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    alignmentSieve --ATACshift --numberOfProcessors ${task.cpus} --bam ${bam} --outFile ${prefix}.unsorted.bam
    samtools sort -@ ${task.cpus} -T ${prefix}.sort -o ${prefix}.filtered.tn5.bam ${prefix}.unsorted.bam
    samtools index ${prefix}.filtered.tn5.bam
    rm ${prefix}.unsorted.bam

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
    touch ${prefix}.filtered.tn5.bam ${prefix}.filtered.tn5.bam.bai versions.yml
    """
}
