process SAMTOOLS_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::samtools=1.24 bioconda::htslib=1.24"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd' }"

    input:
    tuple val(meta), path(bams, stageAs: 'lanes/*')

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), emit: bam
    path "versions.yml"                           , emit: versions

    script:
    """
    samtools merge -@ ${task.cpus} -o ${meta.id}.sorted.bam ${bams}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam versions.yml
    """
}
