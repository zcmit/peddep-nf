/*
 * Exact fragment-length histogram from the filtered BAM, nucleosome-free /
 * mono-nucleosome fractions and a check for the nucleosome banding pattern.
 */
process FRAGMENT_QC {
    tag "${meta.id}"
    label 'process_single'

    conda "bioconda::deeptools=3.5.6 bioconda::samtools=1.20"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' :
        'biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.fragment_qc.txt")     , emit: metrics
    tuple val(meta), path("*.fragment_lengths.txt"), emit: lengths
    tuple val(meta), path("*.fragment_sizes.png")  , emit: png
    path "versions.yml"                            , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLCONFIGDIR=\$PWD/.matplotlib
    fragment_qc.py --bam ${bam} --sample ${meta.id} --prefix ${prefix} ${args}

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
    printf 'Sample\tFragments\tMedian_Fragment\tNFR_Fraction\tMono_Fraction\tNucleosome_Pattern\n${meta.id}\t100\t120\t0.5\t0.2\tpass\n' > ${prefix}.fragment_qc.txt
    touch ${prefix}.fragment_lengths.txt ${prefix}.fragment_sizes.png versions.yml
    """
}
