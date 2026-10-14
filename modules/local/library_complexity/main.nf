/*
 * ENCODE library-complexity metrics (NRF, PBC1, PBC2) from the
 * duplicate-marked (not removed) BAM.
 */
process LIBRARY_COMPLEXITY {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::deeptools=3.5.6 bioconda::samtools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' :
        'biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val mito_name

    output:
    tuple val(meta), path("*.library_complexity.txt"), emit: txt
    path "versions.yml"                              , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    library_complexity.py --bam ${bam} --sample ${meta.id} --exclude-contig '${mito_name}' ${args} --out ${prefix}.library_complexity.txt

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
    printf 'Sample\tTotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF\tPBC1\tPBC2\n${meta.id}\t100\t90\t80\t10\t0.9\t0.89\t8\n' > ${prefix}.library_complexity.txt
    touch versions.yml
    """
}
